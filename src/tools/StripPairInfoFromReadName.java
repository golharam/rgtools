package tools;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.IOException;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

/**
 * Strip read names of /1, /2, /3
 */
public class StripPairInfoFromReadName extends CommandLineProgram {
	private static String VERSION = "0.1";

	private static final Log log = Log.getInstance(StripPairInfoFromReadName.class);

    @Usage public final String USAGE =
            "Strip read names of suffix";

    @Option(shortName="I", doc="Input BAM file")
    public File BAM;

    @Option(shortName="S", doc="Suffix")
    public String SUFFIX;

    @Option(shortName="O", doc="Output BAM file")
    public File OUT;

    public static void main(final String[] argv) {
    	System.exit(new StripPairInfoFromReadName().instanceMain(argv));
    }

	@Override
	protected int doWork() {
		log.info("Version: " + VERSION);
		
        IOUtil.assertFileIsReadable(BAM);
        IOUtil.assertFileIsWritable(OUT);
        
        SAMFileReader inputSam = new SAMFileReader(BAM);
        inputSam.setValidationStringency(ValidationStringency.SILENT);

        SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(), true, OUT);
        
        int count = 0;
        
        for (SAMRecord read : inputSam) {
        	String readName = read.getReadName();

        	int suffixLength = SUFFIX.length();
        	
        	if (readName.endsWith(SUFFIX)) {
        		readName = readName.substring(0,  readName.length() - suffixLength);
        		
        		read.setReadName(readName);
        	}
        	
        	writer.addAlignment(read);
        	
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
        }
        
        writer.close();
        
		return 0;
	}	
}
