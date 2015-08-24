package tools;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.AsciiWriter;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class CountMappedReadPairs extends CommandLineProgram {
	private static String VERSION = "0.2";

	private static final Log log = Log.getInstance(CountMappedReadPairs.class);

	@Usage public final String USAGE =
		"Count mapped/unmapped read pairs in a read name sorted BAM file";

	@Option(shortName="I", doc="Input BAM file")
	public File INPUT;

	@Option(shortName="O", doc="Output TXT file")
	public File OUT;

	public static void main(final String[] argv) {
		System.exit(new CountMappedReadPairs().instanceMain(argv));
	}

	@Override
	protected int doWork() {
		log.info("Version: " + VERSION);
		
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUT);

        // Make sure sort order of input BAM file is by read name
        SAMFileReader inputSam = new SAMFileReader(INPUT);
        final net.sf.samtools.SAMFileHeader fileHeader = inputSam.getFileHeader();
        final SortOrder inputSortOrder = fileHeader.getSortOrder();

        if (inputSortOrder != SortOrder.queryname) {
            log.error("Input BAM file must be sorted by query name");
            inputSam.close();
            return -1;
        }
        
        // Count the read pairs
        inputSam.setValidationStringency(ValidationStringency.LENIENT);
        int count = 0;
        int readPairs = 0, unmappedReadPairs = 0, mappedReadPairs = 0, unPairedReads = 0;
        SAMRecord first = null, second = null;
        boolean isFirstMapped, isSecondMapped;

        for (SAMRecord read : inputSam) {
            // only want 1 count per read so skip non primary alignments
        	if (read.isSecondaryOrSupplementary())
                continue;
        	
            if (read.getReadPairedFlag()) {
                if (read.getFirstOfPairFlag()) {
                	if (first != null)
                		log.warn("Read already encountered: " + read.getReadName());

               		first = read;
                } else {
                	if (second != null)
                		log.warn("Read already encountered: " + read.getReadName());

                	second = read;
                }
                
                if ((first != null) && (second != null)) {
                	if (first.getReadName().equals(second.getReadName())) {
                		readPairs++;
                		
                		isFirstMapped = !first.getReadUnmappedFlag();
                		isSecondMapped = !second.getReadUnmappedFlag();
                		
                		if (isFirstMapped && isSecondMapped) {
                			if ((first.getMateUnmappedFlag() == true) || (second.getMateUnmappedFlag() == true)) {
                				log.warn("Mate information not correct for " + first.getReadName());
                			}

                			mappedReadPairs++;
                		} else {
                			unmappedReadPairs++;
                		}
                		
                		first = null;
                		second = null;
                	} else {
                		log.error("Read names don't match: " + first.getReadName() + " and " + second.getReadName());
                		inputSam.close();
                		return -1;
                	}
                }
                	
            }
            else {
            	unPairedReads++;
            }
            
            if ((count++ % 1000000) == 0) {
                System.out.println("Processed " + count + " reads.");
            }
            
        }
        inputSam.close();
        
        try {
			AsciiWriter out = new AsciiWriter(new FileOutputStream(OUT));
			out.write("MAPPED_READPAIRS\tUNMAPPED_READPAIRS\tREADPAIRS\n");
			out.write(mappedReadPairs + "\t" + unmappedReadPairs + "\t" + readPairs + "\n");
			out.close();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
			return -1;
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
        
		return 0;
	}	
}
