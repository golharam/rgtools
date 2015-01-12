package workflow.examples;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.picard.util.Log;

import org.apache.commons.lang3.StringUtils;

import workflow.Job;
import workflow.Workflow;
import workflow.WorkflowFactory;

public class CreateExomeWorkflow {
	final static String VERSION = "0.01";
    private final Log log = Log.getInstance(CreateExomeWorkflow.class);

	final static String chromosomes[] = { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
										  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
										  "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"};
	
	final static String java = "/usr/java/latest/bin/java";
	final static String gatk = "/share/apps/GenomeAnalysisTK-1.4-14-g2e47336";
	final static String gatk_old = "/share/apps/GenomeAnalysisTK-1.0.5974";
	
	final static String picard = "/share/apps/picard/current/dist";
	final static String java_gatk = java + " -jar " + gatk + "/GenomeAnalysisTK.jar";
	final static String java_gatk_old = java + " -jar " + gatk_old + "/GenomeAnalysisTK.jar";
	final static String java_picard = java + " -jar " + picard;
	
	final static String dbSNP132vcf = gatk + "/dbsnp_132.b37.vcf";
	final static String project_dir = "/share/ngs/ryang/exome_pipeline";
	final static String bam_dir = project_dir + "/BAM";
	final static String exome_dir = project_dir + "/exome";
	
	final static String bwa_reference = "/share/ngs/references/hg19/bwa/hg19.fa";
	
	/* Annovar settings */
	final static String annovar_dir = "/share/apps/annovar";
	final static String annovar_db = annovar_dir + "/humandb";
	final static String annovar_opts = "--genetype knowngene --ver1000g 1000g2011may --buildver hg19";
	
	/* Instance specific variables */
	Workflow workflow;
	Job lastJob;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 3) {
            System.err.println("Usage: java -jar CreateExomeWorkflow.jar <sample name> <forward.fastq.gz> <reverse.fastq.gz> [output file]");
            for (int i = 0; i < args.length; i++) {
            	System.err.println("Arg " + i + ": " + args[i]);
            }
            System.exit(1);
        }
		
		String sampleName = args[0];
		File forwardReads = new File(args[1]);
		File reverseReads = new File(args[2]);		

		CreateExomeWorkflow workflowCreator = new CreateExomeWorkflow();
		Workflow workflow = workflowCreator.createWorkflow_v1(sampleName, forwardReads, reverseReads);

		if (args.length == 4) {
			try {
				workflow.writeToWriter(new FileWriter(new File(args[3])), true);
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else
			workflow.writeToSTDOUT();

	}
	
	public Workflow createWorkflow_v1(String sampleName, File forwardReads, File reverseReads) {
		/* 
		 * Exome Pipeline Steps:
		 * 1. Count Reads
		 * 2. Mapping
		 *    a.  Run 'bwa aln' on the forward and reverse reads
		 *    b.  Run 'bwa sampe' on the sai files from (1a)
		 *    c. Convert/Sort the sam (1b) to bam
		 *    d. Mark Duplicates
		 * 3. Local Realignment
		 *    a. Create Realign Targets
		 *    b. Realign
		 * 3. Fix Mate Information
		 * 4. Recalibrate Quality Scores
		 * 	  a.  Count Covariates
		 * 	  b.  Write new QVs
		 * 5. Depth of Coverage
		 * 6. Reduce BAMs
		 */	
		log.info("VERSION" + VERSION);
		log.info("Creating Workflow Exome_Pipeline");
		workflow = WorkflowFactory.createNewWorkflow("Exome_Pipeline_v1");
		lastJob = null;
		
		String fastqFiles[] = new String[2];
		fastqFiles[0] = forwardReads.getAbsolutePath();
		fastqFiles[1] = reverseReads.getAbsolutePath();
		
		// 1. Count Reads - This can be done via Map/Reduce, but how will that work in the XML format?
		// For now, just submit as its own job
		countReads(sampleName, fastqFiles);
		
/*		// 1. Mapping
		String bamFile = map(sampleName, fastqFiles);
				
		// 2. Realign around indels
		bamFile = realign(bamFile);
				
		// 3. Mark Duplicates
		bamFile = markDuplicates(bamFile);
		
		// 4. Recalibrate QVs
		bamFile = recalibrateQualityScores(bamFile);
	
		// 5.  Call SNPs/Indels
		String vcfFile = callSNPs(bamFile);
		
		// 6.  Annotate Variants
		annotateVCF(workflow, vcfFile, lastJob);
*/		
		return workflow;
	}

	private String countReads(String sampleName, String[] fastqFiles) {
		String outputINI = sampleName + "/" + sampleName + ".ini";
		String cmd = String.format("gunzip -c %s | wc -l | awk '{print $1/4}' > %s; gunzip -c %s | wc -l | awk '{print $1/4}' >> %s", fastqFiles[0], outputINI, fastqFiles[1], outputINI);
		workflow.createJob(sampleName+".countReads", 1, cmd, fastqFiles, outputINI);
		return outputINI;
	}
	 
    /*
     * Maps fragment/paired-end fastq files using BWA.  Steps:
     * a.  Run 'bwa aln' on the forward and reverse reads
	 * b.  Run 'bwa samse/sampe' on the sai file(s) from (1a)
	 * c. Convert the sam (b) to bam
	 * d. Sort the bam (and index)
     * Also sets the reference to last job in the task
	 *  
	 * @param sampleName Name of Sample
     * @param fastqFiles Absolute path to forward (and reverse) fastq file(s)
     * @return outputBAM absolute path of resulting BAM file
     */
	private String map(String sampleName, String[] fastqFiles) {
		String outputBAM = bam_dir + "/" + sampleName + ".bam";

		String saiFiles[] = new String[2];
		saiFiles[0] = bam_dir + "/" + sampleName + "_1.sai";
		saiFiles[1] = bam_dir + "/" + sampleName + "_2.sai";

		String readGroup = String.format("@RG\\tID:1\\tLB:%s\\tPL:Illumina\\tPU:Hiseq\\tSM:%s", sampleName, sampleName);

		// a.  Run 'bwa aln' on the forward (and reverse) reads
		Job align_left, align_right = null;
		
		String cmd = String.format("bwa aln -t 8 %s %s > %s", bwa_reference, fastqFiles[0], saiFiles[0]);
		align_left = workflow.createJob("align_left", 8, cmd, fastqFiles[0], saiFiles[0]);	// (name, cores, cmd, inputfile, outputfile)
		
		if (fastqFiles.length == 2) {
			cmd = String.format("bwa aln -t 8 %s %s > %s", bwa_reference, fastqFiles[1], saiFiles[1]);
			align_right = workflow.createJob("align_right", 8, cmd, fastqFiles[1], saiFiles[1]);
		}
		
		// b.  Run 'bwa samse/sampe' on the sai file(s) from the previous step
		// c.  Convert sam to bam
		// d.  Sort the bam
		Job sai2bam;
		
		if (fastqFiles.length == 1) {
			cmd = String.format("bwa samse -r \"%s\" %s %s %s | samtools view -S -b - | samtools sort - %s", 
					readGroup, bwa_reference, saiFiles[0], fastqFiles[0], bam_dir + "/" + sampleName);
			sai2bam = workflow.createJob("samse", 1, cmd, saiFiles, outputBAM);
			sai2bam.addParent(align_left);
		} else {
			cmd = String.format("bwa sampe -r \"%s\" %s %s %s %s %s | samtools view -S -b - | samtools sort - %s", 
								readGroup, bwa_reference, saiFiles[0], saiFiles[1], fastqFiles[0], 
								fastqFiles[1], bam_dir + "/" + sampleName);
			sai2bam = workflow.createJob("sampe", 1, cmd, saiFiles, outputBAM);
			sai2bam.addParent(align_left);
			sai2bam.addParent(align_right);
		}
		
		// e.  Index the bam file (this can be added in the prior step, but because we already have the BAM
		//     file created, we'll just create a new job here.
		cmd = String.format("samtools index %s", outputBAM);
		Job bamindex = workflow.createJob("bamindex", 1, cmd, outputBAM, outputBAM+".bai");
		bamindex.addParent(sai2bam);
		
		lastJob = bamindex;
		return outputBAM;
	}

	/* 
	 * Perform GATK Indel Realigner
	 * a.  Create local realignment target intervals
	 * b.  Realign target
     * Also sets the reference to last job in the task
     *
	 * @param inputBAM absolute path to BAM file
     * @return outputBAM absolute path of resulting BAM file
	 */
	private String realign(String inputBAM) {
		String cmd;
		
		// 4a.  Create local realignment target intervals
		//	TODO: Split by chromosome, then realign by chromosome, then merge resulting bam file
		/*
		Job[] chr_realignertargetcreator = new Job[chromosomes.length];
		String[] chr_realignertargets = new String[chromosomes.length];
		for (int i = 0; i < chromosomes.length; i++) {
			chr_realignertargets[i] = bam_dir + "/" + sampleName + "." + chromosomes[i] + ".realign.intervals";
			cmd = String.format("%s -T RealignerTargetCreator -nt 8 -R %s -I %s -o %s -L %s", java_gatk, bwa_reference, bamFile, chr_realignertargets[i], chromosomes[i]);
			chr_realignertargetcreator[i] = workflow.createJob("realignertargetcreator_"+chromosomes[i], 8, cmd, bamFile, chr_realignertargets[i]);
			chr_realignertargetcreator[i].addParent(bamindex);
		}
		*/
		
		String intervalsFile = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + ".intervals";
		cmd = String.format("%s -T RealignerTargetCreator -nt 8 -R %s -I %s -o %s", java_gatk, bwa_reference, inputBAM, intervalsFile);
		Job realignertargetcreator = workflow.createJob("realignertargetcreator", 8, cmd, inputBAM, intervalsFile);
		realignertargetcreator.addParent(lastJob);

		// 4b.  Realign Indel targets
		// Can this be divided up by chromosome and distributed?
		String outputBAM = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + ".realign.bam";
		cmd = String.format("%s -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s", java_gatk, bwa_reference, inputBAM, intervalsFile, outputBAM);
		Job realigner = workflow.createJob("realigner", 1, cmd, intervalsFile, outputBAM);
		realigner.addParent(realignertargetcreator);
		
		lastJob = realigner;
		return outputBAM;
	}

	/* 
	 * Run Picard MarkDuplicates
     * Also sets the reference to last job in the task
	 *
	 * @param inputBAM absolute path to BAM file
     * @return absolute path of resulting BAM file
	 */
	private String markDuplicates(String inputBAM) {
		String prefix = inputBAM.substring(0, inputBAM.lastIndexOf(".bam"));
		String outputBAM =  prefix + ".markdups.bam";
		String markDupsMetricsFile = prefix + ".markdups.metrics";
		String cmd = String.format("%s/MarkDuplicates.jar INPUT=%s OUTPUT=%s M=%s VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true ASSUME_SORTED=true", 
							java_picard, inputBAM, outputBAM, markDupsMetricsFile);
		Job markdups = workflow.createJob("markduplicates", 1, cmd, inputBAM, outputBAM);
		markdups.addParent(lastJob);
		
		lastJob = markdups;
		return outputBAM;
	}

	/* 
	 * Perform GATK Quality Score Recalibration
	 * a.  CountCovariates
	 * b.  TableRecalibration
     * Also sets the reference to last job in the task
	 *
	 * @param inputBAM absolute path to BAM file
     * @return absolute path of resulting BAM file
	 */
	private String recalibrateQualityScores(String inputBAM) {
		String recalFile = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + ".recal.csv";
		String cmd = String.format("%s -T CountCovariates -nt 8 -l INFO -R %s -knownSites %s -I %s --standard_covs " +
						    "-cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile %s", 
							 java_gatk, bwa_reference, dbSNP132vcf, inputBAM, recalFile);
		Job count_covariates = workflow.createJob("count_covariates", 8, cmd, inputBAM, recalFile);
		count_covariates.addParent(lastJob);
		
		// Recalibrate
		String recalBamFile = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + ".recal.bam";
		cmd = String.format("%s -T TableRecalibration -I %s -l INFO -R %s --out %s -recalFile %s", 
							java_gatk, inputBAM, bwa_reference, recalBamFile, recalFile);
		Job recalibrate = workflow.createJob("recalibrate", 1, cmd, recalFile, recalBamFile);
		recalibrate.addParent(count_covariates);
		
		lastJob = recalibrate;
		return recalBamFile;
	}

	/* 
	 * Run GATK Unified Genotyper
     * Also sets the reference to last job in the task
	 *
	 * @param inputBAM absolute path to BAM file
     * @return absolute path of resulting VCF file
	 */
	private String callSNPs(String inputBAM) {
		// 7. Call SNPs via UnifiedGenotyper
		// This can be split by chromosome to speed it up, then combined later
		/*
		String samplevcf = exome_dir + "/" + sampleName + ".vcf";
		cmd = String.format("%s -T UnifiedGenotyper -nt 8 -glm BOTH -dcov 1000 -A DepthOfCoverage -A AlleleBalance -A FisherStrand -A ReadPosRankSumTest " +
							"-B:dbsnp,vcf %s -R %s -I %s -o %s",
							java_gatk, dbSNP132vcf, bwa_reference, recalBamFile, samplevcf);
		*/
		String cmd;
		Job[] chrsnp = new Job[chromosomes.length];
		String[] chrvcf = new String[chromosomes.length];
		for (int i = 0; i < chromosomes.length; i++) {
			chrvcf[i] = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + "." + chromosomes[i] + ".vcf";
			String metrics = chrvcf[i].substring(0, chrvcf[i].lastIndexOf(".vcf")) + ".metrics";
			cmd = String.format("%s -T UnifiedGenotyper -nt 8 -glm BOTH -dcov 1000 -A DepthOfCoverage -A AlleleBalance -A FisherStrand -A ReadPosRankSumTest " +
					"-R %s -I %s -o %s -L %s -metrics %s",
					java_gatk, bwa_reference, inputBAM, chrvcf[i], chromosomes[i], metrics);
			chrsnp[i] = workflow.createJob("unified_genotyper_" + chromosomes[i], 8, cmd, inputBAM, chrvcf[i]);
			chrsnp[i].addParent(lastJob);
		}

		// now join the chromosomes
		// this cannot be parallelized
		String samplevcf = inputBAM.substring(0, inputBAM.lastIndexOf(".bam")) + ".vcf";
		String chrparts = new String();
		for (int i = 0; i < chromosomes.length; i++) {
			chrparts = String.format("%s -V:%s %s ", chrparts, chromosomes[i], chrvcf[i]);
		}
		String priority = StringUtils.join(chromosomes, ",");
		cmd = String.format("%s -T CombineVariants %s -priority %s -R %s -o %s", java_gatk, chrparts, priority, bwa_reference, samplevcf);
		Job mergechrvcf = workflow.createJob("mergechrvcf", 1, cmd, chrvcf, samplevcf);
		for (int i = 0; i < chromosomes.length; i++) {
			mergechrvcf.addParent(chrsnp[i]);
		}

		lastJob = mergechrvcf;
		return samplevcf;
	}

	/* 
	 * Run Annovar to annotate VCF file
	 *
	 * @param workflow to add jobs to
	 * @param inputVCF absolute path to VCF file
     * @return nothing
	 */
	private static void annotateVCF(Workflow workflow, String inputVCF, Job lastDependentJob) {
		String cmd;
		
		// 1. Convert SNP VCF to Annovar
		String avFile = inputVCF.substring(0, inputVCF.lastIndexOf(".vcf")) + ".snp.av";
		cmd = String.format("%s/convert2annovar.pl -format vcf4 %s > %s", annovar_dir, inputVCF, avFile);
		Job convertToAnnovar = workflow.createJob("convertToAnnovar", 1, cmd, inputVCF, avFile);
		convertToAnnovar.addParent(lastDependentJob);
		
		// 2. Run auto_annovar.pl
		String geneFile = avFile + ".genelist";
		cmd = String.format("%s/auto_annovar.pl %s -step 1,4,7-9 -model recessive --verdbsnp 132 %s %s",
							 annovar_dir, annovar_opts, avFile, annovar_db);
		Job autoAnnovar = workflow.createJob("autoAnnovar", 1, cmd, avFile, geneFile);
		autoAnnovar.addParent(convertToAnnovar);
	}
}
