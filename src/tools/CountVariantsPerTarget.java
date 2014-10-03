package tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;

import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import reader.BEDReader;
import tools.vcfbrowser.VarTableModel;

public class CountVariantsPerTarget extends CommandLineProgram {
	private static String VERSION = "0.1";

	private static final Log log = Log.getInstance(CountVariantsPerTarget.class);
	
    @Usage public final String USAGE =
            "Calculates a set of variant metrics from an aligned BAM file of all features in a BED file";

    @Option(shortName="V", doc="Input VCF file")
    public File VCF;

    @Option(shortName="B", doc="Input BED file")
    public File BED;

    @Option(shortName="O", doc="Output TXT file")
    public File OUTPUT;
    
	public static void main(String[] args) {
		System.exit(new CountVariantsPerTarget().instanceMain(args));
	}

	@Override
	protected int doWork() {
		VarTableModel varTableModel = null;
		String[] sampleNames;
        FeatureList bedIntervals = null;

        log.info("Version: " + VERSION);
		
        IOUtil.assertFileIsReadable(VCF);
        IOUtil.assertFileIsReadable(BED);
        IOUtil.assertFileIsWritable(OUTPUT);

        try {
			bedIntervals = BEDReader.read(BED.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
        
        // 1.  Determine how many samples are in the VCF file
        // 2.  For each target in bed file:
		// 3.  		For each sample
        //				Count the number of variants that exist in this sample in this target

        // Determine how many samples are in the VCF file
        varTableModel = new VarTableModel(VCF, -1);
        sampleNames = varTableModel.getSampleNames();
        if (sampleNames.length > 1) {
            log.error("Found more than 1 sample in VCF file.  This tool currently only supports 1 sample per VCF file.");
        	return -1;
        }

        // For each target in bed file
        for (FeatureI target : bedIntervals) {
			// Get a list of variants for this target for this sample
			ArrayList<VcfEntry> variants = varTableModel.getVariantsInRange(target.seqname(), target.location().bioStart(), target.location().bioEnd());
			System.out.println(target.toString() + "\t" + variants.size());
        }
		return 0;
	}

}
