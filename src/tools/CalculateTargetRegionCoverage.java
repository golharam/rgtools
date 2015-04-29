package tools;

/* Useful approach:
 * 1.  Download refGene.gtf from UCSC Genome Browser
 * 2.  Construct a list of accession numbers and put in text file.  
 *     Make sure accession numbers aren't versioned (ie .1, .2, etc)
 * 3.  for accession in `cat accession_numbers.txt`
 * 			do echo $accession
 * 			grep $accession ../../pcgc/reference/b37/refGene.hg19.gtf > $accession.refGene.gtf
 * 	   done 
 * 4.  Concatenate all GTF files into one file.
 * 5.  Run this program on BAM and GTF file.
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import reader.BEDReader;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.CoordMath;

public class CalculateTargetRegionCoverage extends CommandLineProgram {
	private static String VERSION = "0.1";
	
	private static final Log log = Log.getInstance(CalculateTargetRegionCoverage.class);
	
    @Usage public final String USAGE =
            "Calculates a set of coverage metrics from an aligned SAM or BAM file of all features in a BED file";

    @Option(shortName="B", doc="Input BAM file")
    public File BAM;

    @Option(shortName="G", doc="Input BED file")
    public File BED;

    @Option(shortName="O", doc="Output TXT file")
    public File OUTPUT;
    
    public static void main(final String[] argv) {
    	System.exit(new CalculateTargetRegionCoverage().instanceMain(argv));
    }
    
    /* Return true if the read should be filtered out
     * 
     */
    private boolean filterRead(SAMRecord rec) {
        // Just plain avoid records that are marked as not-primary
        if (rec.getNotPrimaryAlignmentFlag()) return true;
        // Check for PF reads
        if (rec.getReadFailsVendorQualityCheckFlag()) return true;
        // Check for unmapped reads
        if (rec.getReadUnmappedFlag()) return true;
        // Check for reads that are marked as duplicates
        if (rec.getDuplicateReadFlag()) return true;
        // Don't bother with reads that didn't align uniquely
        if (rec.getMappingQuality() == 0) return true;

        return false;
    }
    
	@SuppressWarnings("resource")
	@Override
	protected int doWork() {
		log.info("Version: " + VERSION);
		
        IOUtil.assertFileIsReadable(BAM);
        IOUtil.assertFileIsReadable(BED);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        AsciiWriter out = null;
        StringBuffer buffer = new StringBuffer();
        
        long totalReads = 0;
        long totalTargetBases = 0;
        long totalPFMappedReads = 0;
        long totalPFMappedOnTargetReads = 0;
        long totalPFMappedOnTargetBases = 0;
        long totalPFMappedOnTargetBases0x = 0;
        long totalPFMappedOnTargetBases10x = 0;
        long totalPFMappedOnTargetBases20x = 0;
        long totalPFMappedOnTargetBases50x = 0;
        long totalPFMappedOnTargetBases100x = 0;

        long bases_at_normalizedCoverage_greaterThan0_0 = 0;
        long bases_at_normalizedCoverage_greaterThan0_2 = 0;
        long bases_at_normalizedCoverage_greaterThan0_5 = 0;
        long bases_at_normalizedCoverage_greaterThan1_0 = 0;
        long bases_at_normalizedCoverage_greaterThan1_5 = 0;
        long bases_at_normalizedCoverage_greaterThan2_0 = 0;
        long bases_at_normalizedCoverage_greaterThan2_5 = 0;
        long bases_at_normalizedCoverage_greaterThan3_0 = 0;
        
        try {
			out = new AsciiWriter(new FileOutputStream(OUTPUT));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
			return -1;
		}
        
        final SAMFileReader samReader = new SAMFileReader(BAM);
        
		// output header, part 1
    	try {
            out.write("# CalculateTargetRegionCoverage: B=" + BAM.getCanonicalPath() + " BED=" + BED.getCanonicalPath() + " O=" + OUTPUT.getCanonicalPath() + "\n\n");
    	
	    	// count total # of PF reads mapping to reference
	        SAMRecordIterator records = samReader.iterator();
	        while (records.hasNext()) {
	            final SAMRecord rec = records.next();
	            // TBD: This value does not match whats in the fastq file.  Why?
	            totalReads++;	
	            
	            if (filterRead(rec)) continue;
	
	            totalPFMappedReads++;
	            
	            if (++totalReads % 1000000 == 0) {
	                log.info("Counting " + totalReads + " records so far.");
	            }
	        }
	        records.close();
	    	out.write("Total Reads\t" + totalReads + "\n");
	    	out.write("Total PF Mapped Reads\t" + totalPFMappedReads + "\n");
	    	log.info("Total Reads: " + totalReads + ", Total PF Mapped Reads:" + totalPFMappedReads);
    	} catch (IOException e) {
    		e.printStackTrace();
    		return -1;
    	}
    	
    	// collect feature information
        FeatureList listFeatures = null;
		try {
			listFeatures = BEDReader.read(BED.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
				
        buffer.append("Target\tChromosome\tStart\tEnd\tStrand\tLength\tTotal Reads\tMean Coverage\tAvgPerBaseCoverage\tBases10X\tPctBases10X\n");

        // iterate over each feature
    	int target = 1;
        for (FeatureI feature : listFeatures) {
        	// Make sure the feature chromosome is in header of the sam file
        	if (samReader.getFileHeader().getSequenceDictionary().getSequence(feature.seqname()) == null) {
        		log.warn("Feature " + feature.seqname() + " does not exist in the SAM reference. Skipping...");
        		continue;
        	}
        	
        	HashMap<String, String> userData = feature.userData();
        	int[] perBaseCoverage = new int[feature.location().length()];
        	double[] normalizedBaseCoverage = new double[feature.location().length()];
        	
        	// iterate over each sam record that overlaps this feature
        	SAMRecordIterator it = samReader.queryOverlapping(feature.seqname(), feature.location().bioStart(), feature.location().bioEnd());
        	int featureReadCount = 0;
        	int onTargetBases = 0;
        	while (it.hasNext()) {
        		SAMRecord rec = it.next();
        		
        		if (filterRead(rec)) continue;

        		totalPFMappedOnTargetReads++;
                featureReadCount++;
                
                // walk over each base of the feature and add 1 to each base this read covers.
                for (final AlignmentBlock block : rec.getAlignmentBlocks()) {
                    final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                    for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                        if (pos >= feature.location().bioStart() && pos <= feature.location().bioEnd()) {
                        	totalPFMappedOnTargetBases++;
                        	onTargetBases++;
                        	perBaseCoverage[pos - feature.location().bioStart()]++;
                        }
                    }
                }
        	}
        	it.close();
        	
        	// calculate some information regarding this feature
        	totalTargetBases += feature.location().length();
        	int bases10x = 0;
        	int bases20x = 0;
        	int bases50x = 0;
        	int bases100x = 0;
        	double coverage = featureReadCount / (double)feature.location().length();
        	userData.put("ReadCount", String.valueOf(featureReadCount));
        	userData.put("MeanCoverage", String.valueOf(coverage));
        	
        	double perBaseAverageCoverage = onTargetBases / (double) feature.location().length();
        	userData.put("perBaseCoverageArray", String.valueOf(perBaseCoverage));
        	userData.put("avgPerBaseCoverage", String.valueOf(perBaseAverageCoverage));
        	for (int i = 0; i < perBaseCoverage.length; i++) {
        		if (perBaseCoverage[i] == 0) {
        			totalPFMappedOnTargetBases0x++;
        		}
        		if (perBaseCoverage[i] >= 10) {
        			bases10x++;
        			totalPFMappedOnTargetBases10x++;
        		}
        		if (perBaseCoverage[i] >= 20) {
        			bases20x++;
        			totalPFMappedOnTargetBases20x++;
        		}
        		if (perBaseCoverage[i] >= 50) {
        			bases50x++;
        			totalPFMappedOnTargetBases50x++;
        		}
        		if (perBaseCoverage[i] >= 100) {
        			bases100x++;
        			totalPFMappedOnTargetBases100x++;
        		}
        		
        		normalizedBaseCoverage[i] = perBaseCoverage[i] / perBaseAverageCoverage;
        		
        		if (normalizedBaseCoverage[i] >= 0)
        			bases_at_normalizedCoverage_greaterThan0_0++;
        		if (normalizedBaseCoverage[i] >= 0.2)
        			bases_at_normalizedCoverage_greaterThan0_2++;
        		if (normalizedBaseCoverage[i] >= 0.5)
        			bases_at_normalizedCoverage_greaterThan0_5++;
        		if (normalizedBaseCoverage[i] >= 1.0)
        			bases_at_normalizedCoverage_greaterThan1_0++;
        		if (normalizedBaseCoverage[i] >= 1.5)
        			bases_at_normalizedCoverage_greaterThan1_5++;
        		if (normalizedBaseCoverage[i] >= 2.0)
        			bases_at_normalizedCoverage_greaterThan2_0++;
        		if (normalizedBaseCoverage[i] >= 2.5)
        			bases_at_normalizedCoverage_greaterThan2_5++;
        		if (normalizedBaseCoverage[i] >= 3.0)
        			bases_at_normalizedCoverage_greaterThan3_0++;	
        	}
        	userData.put("normalizedBaseCoverageArray", String.valueOf(normalizedBaseCoverage));
        	
        	buffer.append(feature.getAttribute("name") + "\t" + feature.seqname() + "\t" + feature.location().bioStart() + "\t" + feature.location().bioEnd() + 
        					"\t" + feature.location().bioStrand() + "\t" + feature.location().length() + 
        					"\t" + featureReadCount + "\t" + coverage + "\t" + perBaseAverageCoverage + "\t" + 
        					bases10x + "\t" + bases10x/(double)feature.location().length() + 
        					bases20x + "\t" + bases20x/(double)feature.location().length() + 
        					bases50x + "\t" + bases50x/(double)feature.location().length() + 
        					bases100x + "\t" + bases100x/(double)feature.location().length() + 
        					"\n");
        }

        // output header, part 2 
        // then the feature information
        try {
	        out.write("Total Reads Mapped to Targets\t" + totalPFMappedOnTargetReads + "\n");
	        out.write("Enrichment Efficiency\t" + totalPFMappedOnTargetReads/(double)totalPFMappedReads + "\n");
	        out.write("Total Target Bases\t" + totalTargetBases + "\n");
	        out.write("Total Bases Mapped on Target\t" + totalPFMappedOnTargetBases + "\n");
	        out.write("Total Bases Mapped on Target at 0x\t" + totalPFMappedOnTargetBases0x + "\n");
	        out.write("Total Bases Mapped on Target at 10x\t" + totalPFMappedOnTargetBases10x + "\n");
	        out.write("Total Bases Mapped on Target at 20x\t" + totalPFMappedOnTargetBases20x + "\n");
	        out.write("Total Bases Mapped on Target at 50x\t" + totalPFMappedOnTargetBases50x + "\n");
	        out.write("Total Bases Mapped on Target at 100x\t" + totalPFMappedOnTargetBases100x + "\n");
	        out.write("\n");
	        
	        // Write out summary coverage information
	        out.write("Normalized Coverage\t# Bases_at_Given_NormalizedCoverage\tFraction_Bases_at_Given_NormalizedCoverage\n");
	        out.write(">= 0\t" + bases_at_normalizedCoverage_greaterThan0_0 + "\t" + bases_at_normalizedCoverage_greaterThan0_0 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 0.2\t" + bases_at_normalizedCoverage_greaterThan0_2 + "\t" + bases_at_normalizedCoverage_greaterThan0_2 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 0.5\t" + bases_at_normalizedCoverage_greaterThan0_5 + "\t" + bases_at_normalizedCoverage_greaterThan0_5 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 1.0\t" + bases_at_normalizedCoverage_greaterThan1_0 + "\t" + bases_at_normalizedCoverage_greaterThan1_0 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 1.5\t" + bases_at_normalizedCoverage_greaterThan1_5 + "\t" + bases_at_normalizedCoverage_greaterThan1_5 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 2.0\t" + bases_at_normalizedCoverage_greaterThan2_0 + "\t" + bases_at_normalizedCoverage_greaterThan2_0 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 2.5\t" + bases_at_normalizedCoverage_greaterThan2_5 + "\t" + bases_at_normalizedCoverage_greaterThan2_5 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write(">= 3.0\t" + bases_at_normalizedCoverage_greaterThan3_0 + "\t" + bases_at_normalizedCoverage_greaterThan3_0 / (double)totalPFMappedOnTargetBases + "\n");
	        out.write("\n");
	        
	        // Write out per exon information
	        out.write(buffer.toString());

        } catch (IOException e) {
        	e.printStackTrace();
        	return -1;
        }
    
    	try {
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	return 0;
	}

}
