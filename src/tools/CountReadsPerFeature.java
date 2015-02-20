package tools;

import genome.features.GFFFeature;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.FormatUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.StringLineReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.Location;

import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import reader.GFF3Reader;

/*
 * Similar to htseq-count
 * Version 0.01: Initial implementation
 * TODO: 1.  Store seen read names in case read is encountered spanning over an intron
 *       2.  Handle read pairs and not double count       
 */
public class CountReadsPerFeature extends CommandLineProgram {
	public String VERSION = "0.01";
	
    @Option(shortName="G", doc="Input GTF file with Picard/SAM-style header (Picard Interval Header)")
    public File GTF;

    @Option(shortName="I", doc="Input BAM file")
    public File BAM;

    @Option(shortName="O", doc="Output Text file")
    public File OUTPUT;
    
    private final Log log = Log.getInstance(CountReadsPerFeature.class);
    
    public static void main(final String[] args) {
        new CountReadsPerFeature().instanceMainWithExit(args);
    }

    /**
     * Parses an GTF file from a reader in a stream based fashion.
     * @param GTF File
     * @return an IntervalList object that contains the headers and intervals from the file
     */
	private IntervalList GTFToPicardInterval(File GTF) {
		BufferedReader in = null;
        try {
        	in = new BufferedReader(new InputStreamReader(new FileInputStream(GTF)));

        	// Setup a reader and parse the header
            final StringBuilder builder = new StringBuilder(4096);
            String line = null;

            while ((line = in.readLine()) != null) {
                if (line.startsWith("@")) {
                    builder.append(line).append('\n');
                }
                else {
                    break;
                }
            }

            if (builder.length() == 0) {
                throw new IllegalStateException("GTF Interval list file must contain header. ");
            }

            final StringLineReader headerReader = new StringLineReader(builder.toString());
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            final IntervalList list = new IntervalList(codec.decode(headerReader, "BufferedReader"));
            final SAMSequenceDictionary dict = list.getHeader().getSequenceDictionary();

            // Then read in the intervals
            final FormatUtil format = new FormatUtil();
            do {
                if (line.trim().length() == 0) continue; // skip over blank lines

                // Make sure we have the right number of fields
                final String[] fields = line.split("\t");
                if (fields.length != 9) {
                    throw new PicardException("Invalid GTF record contains " +
                                              fields.length + " fields: " + line);
                }

                // Then parse them out
                final String seq = fields[0];
                final int start = format.parseInt(fields[3]);
                final int end   = format.parseInt(fields[4]);

                final boolean negative;
                if (fields[6].equals("-")) negative = true;
                else if (fields[6].equals("+")) negative = false;
                else throw new IllegalArgumentException("Invalid strand field: " + fields[6]);

                // TODO: Pull out gene_id from field.  Field looks like
                // gene_id "NM_183416"; transcript_id "NM_183416"; 
                final String name = fields[8].substring(fields[8].indexOf("gene_id")+9, fields[8].indexOf("\"; transcript_id"));

                final Interval interval = new Interval(seq, start, end, negative, name);
                if (dict.getSequence(seq) == null) {
                    log.warn("Ignoring interval for unknown reference: " + interval);
                }
                else {
                	list.add(interval);                    
                }
            }
            while ((line = in.readLine()) != null);

            return list;
        }
        catch (IOException ioe) {
            throw new PicardException("Error parsing interval list.", ioe);
        }
        finally {
            try { in.close(); } catch (Exception e) { /* do nothing */ }
        }
	}
    
    /**
     * Count the number of reads that overlap a feature in a sam file
     * @param feat Feature to compare
     * @param inputSam SAM file
     * @return number of overlapping reads
     */
	protected int countOverlappingReads(GFFFeature feature, SAMFileReader inputSam) {
		String chr = feature.getChromosome();
		int start = feature.getStart();
		int end = feature.getEnd();
	
		@SuppressWarnings("deprecation")
		SAMRecordIterator iterator = inputSam.queryOverlapping(chr, start, end);
		
		int count = 0;
		while (iterator.hasNext()) {
			SAMRecord sam = iterator.next();
			
            // Just plain avoid records that are marked as not-primary
            if (sam.getNotPrimaryAlignmentFlag())
            	continue;
            
            // Check for PF reads
            if (sam.getReadFailsVendorQualityCheckFlag()) {
                continue;
            }
            if (sam.getReadUnmappedFlag()){
            	continue;
            }
            // Check for reads that are marked as duplicates
            if (sam.getDuplicateReadFlag()) {
                continue;
            }
            // Don't bother with reads that didn't align uniquely
            if (sam.getReadUnmappedFlag() || sam.getMappingQuality() == 0) {
                continue;
            }

			count++;
		}
		iterator.close();

		return count;
	}

	@Override
	protected int doWork() {
		try {
			log.info("Version " + VERSION);
			IOUtil.assertFileIsReadable(BAM);
			IOUtil.assertFileIsReadable(GTF);
	        IOUtil.assertFileIsWritable(OUTPUT);
        
	        final SAMFileReader inputSam = new SAMFileReader(BAM);
	        GFF3Reader gffReader = new GFF3Reader(GTF);
	        GFFFeature feature;
	        //List<GFFFeature> featureList = GFF3Reader.read(GTF);
	        HashMap<String, Integer> readCountPerGene = new HashMap<String, Integer>();
	        

			// Iterating over each read to see which feature it overlaps will take a long time
	        // since the reads in teh file can be anywhere.  Instead its faster to query for a 
	        // bunch of reads that overlap a given feature.  So, iterate over each feature
	        // and count the reads that intersect with that feature.
	        while ((feature = gffReader.next()) != null) {	     
	        	if (feature.getType().equals("exon")) {
		        	int overlappingReads = countOverlappingReads(feature, inputSam);
		        	String gene = (String)feature.getAttribute("gene_id");
		        	
		        	int readCount = readCountPerGene.containsKey(gene) ? readCountPerGene.get(gene) : 0;
		        	readCountPerGene.put(gene, readCount + overlappingReads);	        	
	        	}
	        }
	        
	        AsciiWriter out = new AsciiWriter(new FileOutputStream(OUTPUT));
	        out.write("gene_id\traw_count\n");
	        for (String gene : readCountPerGene.keySet()) {
	        	int readCount = readCountPerGene.get(gene);
	        	out.write(gene + "\t" + readCount + "\n");
	        }
	        
	        gffReader.close();
	        out.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			
		}
		return 0;
	}
/*       	            
	        exoncounts.write("chr\tgene_id\ttranscript_id\texon_num\tstart\tend\tstrand\texon_len\tread_count\tbase_coverage\tavg_cov\n");
	        
	        long genes = 0;
	        String gene = null;
	        int exon = 1;
	        GFFFeature feature;
	        
	        HashMap<String, Integer> geneCounts = new HashMap<String, Integer>();
	        HashMap<String, Integer> geneLengths = new HashMap<String, Integer>();
	        Integer i;
	        
			while ((feature = gffReader.next()) != null) {
				
		            if (feature.getType().equals(FEATURE_TYPE)) {
		            	if (feature.getAttribute("gene_id").equals(gene)) {
		            		exon++;
		            	} else {
		            		gene = (String) feature.getAttribute("gene_id");
		            		exon = 1;
		            		genes++;
		                    if (++genes % 1000 == 0)
		                    	log.info("Processed " + genes + " genes.");
		            	}
		            	String gene_id = (String) feature.getAttribute("gene_id");
		            	String transcript_id = (String) feature.getAttribute("transcript_id");
		            	int exon_len = feature.getLength();
						int reads = countOverlappingReads(feature, inputSam);
						char strand = feature.getStrand();
						
						int basecoverage = 0;
						double avg_cov = 0;
//						if (reads > 0) {
//			            	int[] coverage = countReadsPerPosition(feature, inputSam);
//			            	for (int data : coverage) {
//			            		basecoverage += data;
//			            	}
//			            	avg_cov = basecoverage / (double)coverage.length;
//						}
					
						i = geneCounts.get(gene_id);
						if (i == null)
							i = reads;
						else
							i += reads;
						geneCounts.put(gene_id, i);
						
						i = geneLengths.get(gene_id);
						if (i == null)
							i = exon_len;
						else
							i += exon_len;
						geneLengths.put(gene_id, i);
						
						exoncounts.write(feature.seqname() + "\t" + gene_id + "\t" + transcript_id + "\t" + exon + "\t" + 
								  Math.abs(feature.location().start()) + "\t" + Math.abs(feature.location().end()) + "\t" + 
								  strand + "\t" + exon_len + "\t" + reads + "\t" + basecoverage + 
								  "\t" + avg_cov + "\n");
		            }
		          
	        }    
			exoncounts.close();

			log.info("Writing counts file...");
			genecounts.write("gene_id\tgene_len\tread_count\n");
			Iterator<String> it = geneCounts.keySet().iterator();
			while (it.hasNext()) {
				String gene_id = it.next();
				Integer readCount = geneCounts.get(gene_id);
				Integer geneLength = geneLengths.get(gene_id);
				
				genecounts.write(gene_id + "\t" + geneLength + "\t" + readCount + "\n");
			}
			genecounts.close();
			
	        log.info("Finished! Total of " + genes + " genes processed.");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
*/


	
	// unused but could come in handy.
	/*
	private int[] countReadsPerPosition(FeatureI feature, SAMFileReader inputSam) {
		//TODO: Check that this works properly for - strand
		int[] coverage = new int[feature.location().length()];
		String chr = feature.seqname();
		int start = Math.abs(feature.location().start()) + 1;
		
		for (int i = 0; i < coverage.length; i++) {
			SAMRecordIterator iterator = inputSam.queryOverlapping(chr, start + i, start + i);
			int count = 0;
			while (iterator.hasNext()) {
				count++;
				SAMRecord rec = iterator.next();
			}
			iterator.close();
			coverage[i] = count;
		}
		return coverage;
	}

	private in[]countReadsPerPositionNew(FeatureI feature, SAMFileReader inputSam) {
        // Prefetch the list of target overlaps here as they're needed multiple times.
        final Collection<Interval> targets;
        final Interval read = new Interval(feature.getReferenceName(), sam.getAlignmentStart(), sam.getAlignmentEnd());
        targets = targetDetector.getOverlaps(read);
        
        // Find the target overlaps
        if (targets != null && !targets.isEmpty()) {
            for (final Interval target : targets) {
                final Coverage coverage = coverageByTarget.get(target);

                for (final AlignmentBlock block : sam.getAlignmentBlocks()) {
                    final int end = CoordMath.getEnd(block.getReferenceStart(), block.getLength());
                    for (int pos=block.getReferenceStart(); pos<=end; ++ pos) {
                        if (pos >= target.getStart() && pos <= target.getEnd()) {
                            coverage.addBase(pos - target.getStart());
                        }
                    }
                }
            }
        }	
	}
*/

}
