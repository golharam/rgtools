package workflow.examples;

import net.sf.picard.util.Log;

import workflow.Job;
import workflow.Workflow;
import workflow.WorkflowFactory;

public class KittleExomeWorkflow {
	public final static String VERSION = "0.01";
    private final Log log = Log.getInstance(KittleExomeWorkflow.class);
    
    static String java_gatk = "java -jar /mnt/isilon/cag/ngs/hiseq/gatk/GenomeAnalysisTK-2.1-8-g5efb575/GenomeAnalysisTK.jar";
    
    static String project_dir = "/mnt/isilon/cag/ngs/hiseq/golharr/Kittle";
    
    static String tlr_pathway_gtffile = project_dir + "/Annotation/TLR_Pathway.gtf";
    
    static String reference = "/mnt/isilon/cag/ngs/hiseq/gatk/hg19/hg19.fa";
    
	static String[] samples = {"SID5768975236", "SID7232589898", "SID7744627578"};
	static String[] sampleNames = {"dad", "mom", "child"};
	
	public static void calculateCoverage(Workflow workflow, String sample) {
		String bamfile = String.format("%s/BAM/%s.dedup.indelrealigner.recal.bam", project_dir, sample);
		String outputfile = String.format("%s/coverage/%s.coverageBed.txt", project_dir, sample);
		String cmd = String.format("coverageBed -abam %s -b %s > %s", bamfile, tlr_pathway_gtffile, outputfile);
		//workflow.createJob(sample + "_bedcoverage", 1, cmd, bamfile, outputfile);
		
		String baitsFile = "/mnt/isilon/cag/ngs/hiseq/golharr/Kittle/Annotation/TLR_Pathway.picard_interval";
		outputfile = project_dir + "/coverage/" + sample + ".calcHsMetrics.txt";
		String pertarget_file = project_dir + "/coverage/" + sample + ".per_target.calcHsMetrics.txt";
		// Sample information is not coming across in the SAMPLE field.  Is this because of PER_TARGET?
		cmd = String.format("java -jar ~/workspace/picard-tools/CalculateHsMetrics.jar BI=%s TI=%s I=%s O=%s R=%s PER_TARGET_COVERAGE=%s " +
							  "LEVEL=ALL_READS VALIDATION_STRINGENCY=LENIENT",
							  baitsFile, baitsFile, bamfile, outputfile, reference, pertarget_file);
		Job perTargetJob = workflow.createJob(sample + "_pertarget", 1, cmd, bamfile, outputfile);
		
		cmd = String.format("java -jar ~/workspace/picard-tools/CalculateHsMetrics.jar BI=%s TI=%s I=%s O=%s R=%s " +
							"LEVEL=null LEVEL=SAMPLE VALIDATION_STRINGENCY=LENIENT",
							baitsFile, baitsFile, bamfile, outputfile, reference);
		Job calcHsMetrics = workflow.createJob(sample + "_calcHsMetrics", 1, cmd, bamfile, outputfile);
		calcHsMetrics.addParent(perTargetJob);
	}

	public static void mergeVariants(Workflow workflow, String type) {
		String[] inputFiles = { project_dir + "/AnnotatedVCF/SID5768975236." + type + ".vcf",
								project_dir + "/AnnotatedVCF/SID7232589898." + type + ".vcf",
								project_dir + "/AnnotatedVCF/SID7744627578." + type + ".vcf" };
		String outputFile = project_dir + "/AnnotatedVCF/family." + type + ".vcf";
		String cmd = String.format("%s -T CombineVariants -R %s -V:%s %s " +
															   "-V:%s %s " +
															   "-V:%s %s " +
															   "-o %s", 
															   java_gatk, reference, 
															   sampleNames[0], inputFiles[0], 
															   sampleNames[1], inputFiles[1], 
															   sampleNames[2], inputFiles[2], 
															   outputFile);
		workflow.createJob("merge"+type, 1, cmd, inputFiles, outputFile);
		
		// convert the vcf to a text file
		/*
		String txtFile = project_dir + "/AnnotatedVCF/family." + type + ".txt";
		cmd = String.format("java -jar ~/workspace/rgtools/VCFToTab.jar V=%s O=%s", outputFile, txtFile);
		workflow.createJob("vcf2txt"+type, 1, cmd, outputFile, txtFile);
		*/
	}

	public static void selectSharedVariants(Workflow workflow, String type) {
		String inputVCF = project_dir + "/AnnotatedVCF/family." + type + ".vcf";
		String outputVCF = project_dir + "/AnnotatedVCF/family.momchild." + type + ".vcf";
		String cmd = "/home/golharr/workspace/rgtools/scripts/filterMomChildVCF.pl " + inputVCF + " > " + outputVCF;
		workflow.createJob("momchild"+type, 1, cmd, inputVCF, outputVCF);

		// convert the vcf to a text file
		String txtFile = project_dir + "/AnnotatedVCF/family.momchild." + type + ".txt";
		cmd = String.format("java -jar ~/workspace/rgtools/VCFToTab.jar V=%s O=%s", outputVCF, txtFile);
		workflow.createJob("vcf2txt"+type, 1, cmd, outputVCF, txtFile);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length != 0) {
            System.err.println("Usage: java -jar KittleExomeWorkflow.jar");
            System.exit(1);
        }

		Workflow workflow = WorkflowFactory.createNewWorkflow("Kittle Workflow");

		// 1.  Determine coverage of TLR Pathway genes
		for (String sample : samples) {
			calculateCoverage(workflow, sample);
		}
		
		// 2.  Look at SNPs 
		mergeVariants(workflow, "snp");
		selectSharedVariants(workflow, "snp");
		
		// 3.  Look at Indels
		mergeVariants(workflow, "indel");
		selectSharedVariants(workflow, "indel");
				
		workflow.writeToSTDOUT();

	}
}
