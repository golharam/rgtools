package tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import org.biojava3.genome.parsers.gff.FeatureI;
import org.biojava3.genome.parsers.gff.FeatureList;

import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import htsjdk.samtools.util.AsciiWriter;
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
        AsciiWriter out = null;

        log.info("Version: " + VERSION);
		
        IOUtil.assertFileIsReadable(VCF);
        IOUtil.assertFileIsReadable(BED);
        IOUtil.assertFileIsWritable(OUTPUT);

        // Open the BED file
        try {
			bedIntervals = BEDReader.read(BED.getAbsolutePath());
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}

        // Open the output text file
        try {
			out = new AsciiWriter(new FileOutputStream(OUTPUT));
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
			return -1;
		}

        // Right now, handle 1 sample per VCF file
        varTableModel = new VarTableModel(VCF, -1);
        sampleNames = varTableModel.getSampleNames();
        if (sampleNames.length > 1) {
            log.error("Found more than 1 sample in VCF file.  This tool currently only supports 1 sample per VCF file.");
        	return -1;
        }

        // For each target in bed file
        for (FeatureI target : bedIntervals) {
			// Get a list of variants for this target
			ArrayList<VcfEntry> variants = varTableModel.getVariantsInRange(target.seqname(), target.location().bioStart(), target.location().bioEnd());
			
			try {
				out.write(target.toString() + "\t" + variants.size() + "\n");
			} catch (IOException e) {
				continue;
			}
        }

        // close the output file
    	try {
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return 0;
	}

}
