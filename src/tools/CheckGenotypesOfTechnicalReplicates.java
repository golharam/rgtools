package tools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

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
    
    private HashMap<String,String[]> readSampleReplicatesMap(File sampleReplicateInfo) {
    	HashMap<String,String[]> sampleToReplicates = new HashMap<String,String[]>();
		try {
	    	String line;
	    	BufferedReader in;
			in = new BufferedReader(new InputStreamReader(new FileInputStream(sampleReplicateInfo)));
		    	
	        while ((line = in.readLine()) != null) {
	            String[] fields = line.split("\t");
	            String sample = fields[0];   
	            String[] replicates = fields[1].split(",");
	            
	            sampleToReplicates.put(sample, replicates);
	        }
	        
	        in.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

    	return sampleToReplicates;
    }
    
	@Override
	protected int doWork() {
		log.info("Version " + VERSION);

		IOUtil.assertFileIsReadable(VCF);
		IOUtil.assertFileIsReadable(REPLICATES);

		HashMap<String,String[]> sampleToReplicates = readSampleReplicatesMap(REPLICATES);
		
		VCFFileReader vcfFileReader = new VCFFileReader(VCF, false);
		CloseableIterator<VariantContext> variantIterator = vcfFileReader.iterator();
		while (variantIterator.hasNext()) {
			VariantContext vc = variantIterator.next();
			
			// for each sample...
			for (String sample : sampleToReplicates.keySet()) {
				// get the genotypes for each of the replicates for this sample
				String[] replicates = sampleToReplicates.get(sample);
				Genotype[] genotypes = new Genotype[replicates.length];
				int i = 0;
				for (String replicate : replicates) {
					genotypes[i] = vc.getGenotype(replicate);
					// if this isn't the first replicate, make sure all the genotypes for the technical replicates 
					// agree with each other
					if (i > 0) {
						if (!genotypes[i].getGenotypeString().equals(genotypes[0].getGenotypeString())) {
							String msg = String.format("Genotype for replicate %s (%s) does not match first replicate %s (%s) for variant at %s:%d", replicate, genotypes[i].getGenotypeString(), replicates[0], genotypes[0].getGenotypeString(), vc.getChr(), vc.getStart());
							log.warn(msg);
						}
					}
					i++;
				}
				
			}
		}
		vcfFileReader.close();
		return 0;
	}

}
