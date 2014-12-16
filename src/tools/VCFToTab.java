package tools;

import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;
import ca.mcgill.mcb.pcingola.vcf.VcfHeader;
import ca.mcgill.mcb.pcingola.vcf.VcfInfo;

/*
 * Version 0.01: Initial version
 * 				 BUG: Per-sample Genotype data not being written out
 * Version 0.02a: Output standard INFO fields separately from per-sample fields.
 * Version 0.02b: Output each per-sample genotype field
 */

public class VCFToTab extends CommandLineProgram {
	public String VERSION = "0.02b";
	
    @Usage public final String USAGE = "Convert a VCF file to a tab-delimited file";

    @Option(shortName="V", doc="Input VCF file")
    public File VCF;

    @Option(shortName="O", doc="Output tab-delimited text file")
    public File OUTPUT;

    private final Log log = Log.getInstance(VCFToTab.class);

    public static void main(final String[] args) {
        new VCFToTab().instanceMainWithExit(args);
    }

	protected int doWork() {
		try { 
			log.info("Version " + VERSION);
			
			IOUtil.assertFileIsReadable(VCF);
	        IOUtil.assertFileIsWritable(OUTPUT);
	
	        // TODO: Redo this to make use to VCFFile
	        VcfFileIterator vcfFile = new VcfFileIterator(VCF.getAbsolutePath());
	        final AsciiWriter out = new AsciiWriter(new FileOutputStream(OUTPUT));
	        long total = 0;

	        // Read the header and construct a table of all the possible info fields
	        VcfHeader header = vcfFile.readHeader();
	    	HashMap<String, VcfInfo> vcfInfoById = new HashMap<String, VcfInfo>();
	    	HashMap<String, VcfInfo> vcfFormatById = new HashMap<String, VcfInfo>();
	    	ArrayList<String> keyList = new ArrayList<String>();
	    	ArrayList<String> perSampleKeyList = new ArrayList<String>();
			ArrayList<String> sampleNames = new ArrayList<String>();
	    	
	    	String headerLines[] = header.toString().split("\n");
			for (String line : headerLines) {
				if (line.startsWith("##INFO=")) {
					VcfInfo vcfInfo = new VcfInfo(line);
					vcfInfoById.put(vcfInfo.getId(), vcfInfo);
					keyList.add(vcfInfo.getId());
				} else
				if (line.startsWith("##FORMAT=")) {
					VcfInfo vcfInfo = new VcfInfo(line);
					vcfFormatById.put(vcfInfo.getId(), vcfInfo);
					perSampleKeyList.add(vcfInfo.getId());
				} else
				if (line.startsWith("#CHROM")) {
					String fields[] = line.split("\t");
					for (int i = 9; i < fields.length; i++)
						sampleNames.add(fields[i]);
				}
			}
	
			// TODO: Print out key
			
			// Print out header lines then data lines
			out.write("CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER");
    		for (String key : keyList) {
    			out.write("\t" + key);
    		}
    		// Now print out the genotype fields for each sample
    		for (String sample : sampleNames) {
    			for (String perSampleKey : perSampleKeyList) {
        			out.write("\t" + sample + "-" + perSampleKey);    				
    			}
    		}
    		out.write("\n");
			
	        // print out the variants with the info field values in the same order
	        for (VcfEntry ve : vcfFile) {
	        		total++;
	        		int pos = ve.getStart()+1;
	        		out.write(ve.getChromosomeNameOri() + "\t" + pos + "\t" + ve.getId() + "\t" + ve.getRef() + "\t" + 
	        				ve.getAltsStr() + "\t" + ve.getQuality() + "\t" + ve.getFilterPass());
	        		for (String key : keyList) {
	        			out.write("\t" + ve.getInfo(key));
	        		}
	        		List<VcfGenotype> genotypes = ve.getVcfGenotypes();
	        		for (VcfGenotype genotype : genotypes) {
	        			for (String perSampleKey : perSampleKeyList) {
	        				if (perSampleKey.equals("GT")) 
		        				out.write("\t" + genotype.getGenotypeStr());
	        				else
	        					out.write("\t" + genotype.get(perSampleKey));
	        			}
	        		}
	        		out.write("\n");
	        }
	        out.close();
	        log.info("Finished! Total of " + total + " lines processed.");
		} catch (Exception e) {
			log.error(e.getMessage());
			return -1;
		}

        return 0;
	}

}
