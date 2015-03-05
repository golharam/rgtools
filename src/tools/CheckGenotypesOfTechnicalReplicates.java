package tools;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;

/*
 * Version 0.01: Initial implementation - Scan VCF file and make sure the 
 * 	             genotypes of technical replicates agree with each other.
 */
public class CheckGenotypesOfTechnicalReplicates extends CommandLineProgram {
	public String VERSION = "0.01";

    @Usage public final String USAGE = "Check the genotypes of technical replicates within a VCF file agree with each other";
    
    @Option(shortName="V", doc="Input VCF file")
    public File VCF;

    @Option(shortName="T", doc="Input replicate file")
    public File REPLICATES;

    private final Log log = Log.getInstance(CheckGenotypesOfTechnicalReplicates.class);
    
    public static void main(final String[] args) {
    	new CheckGenotypesOfTechnicalReplicates().instanceMainWithExit(args);
    }
    
    private HashMap<String,String> readSamples(File sampleReplicateInfo) {
    	HashMap<String,String> repToSampleMap = new HashMap<String,String>();
		try {
	    	String line;
	    	BufferedReader in;
			in = new BufferedReader(new InputStreamReader(new FileInputStream(sampleReplicateInfo)));
		    	
	        while ((line = in.readLine()) != null) {
	            final String[] fields = line.split("\t");
	            log.debug(fields[1] + " -> " + fields[0]);
	            repToSampleMap.put(fields[1], fields[0]);
	        }
	        
	        in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

    	return repToSampleMap;
    }
    
	@Override
	protected int doWork() {
		log.info("Version " + VERSION);

		IOUtil.assertFileIsReadable(VCF);
		IOUtil.assertFileIsReadable(REPLICATES);

		HashMap samples = readSamples(REPLICATES);
		
		//VCFFileReader vcfFileReader = new VCFFileReader(VCF, false);
		return 0;
	}

}
