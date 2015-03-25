package tools;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;

public class FixContigsAndVCFHeader extends CommandLineProgram {
	public String VERSION = "0.01";
	
    @Option(shortName="V", doc="Input VCF file")
    public File VCF;

    @Option(shortName="O", doc="Output VCF file")
    public File OUTPUT;
    
    @Option(shortName="S", doc="SAM Header of Reference to use")
    public File SAM;

    private final Log log = Log.getInstance(FixContigsAndVCFHeader.class);

	@Override
	protected int doWork() {
		log.info("FixContigsAndVCFHeader v" + VERSION);
		
		VCFFileReader vcfFileReader = new VCFFileReader(VCF, false);
		
        SAMFileReader reader = new SAMFileReader(SAM);
        SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();

        VariantContextWriterBuilder vcfBuilder = new VariantContextWriterBuilder();
        vcfBuilder.setOutputFile(OUTPUT);
        vcfBuilder.setOutputFileType(VariantContextWriterBuilder.OutputType.VCF);
        vcfBuilder.setReferenceDictionary(dict);        
        VariantContextWriter vcfWriter = vcfBuilder.build();
        
        VCFHeader vcfHeader = vcfFileReader.getFileHeader();
        vcfHeader.setSequenceDictionary(dict);
        vcfWriter.writeHeader(vcfHeader);
        
        log.info("Writing variants.");
        
        int i = 0;
        for (VariantContext vc : vcfFileReader) {
        	i++;

        	VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);
        	vcBuilder.chr("chr"+vc.getChr());
        	vcfWriter.add(vcBuilder.make());
        	
        	if ((i % 100000) == 0) {
        		log.info("Written " + i + " variants.");
        	}
        }
        
        vcfWriter.close();
        log.info("Wrote " + i + " variants.");
		return 0;
	}

    public static void main(final String[] argv) {
    	System.exit(new FixContigsAndVCFHeader().instanceMain(argv));
    }

}
