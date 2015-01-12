package workflow.examples;

import workflow.Job;
import workflow.Workflow;
import workflow.WorkflowFactory;
import workflow.comet.Comet;

public class RNASeq_AWS_Workflow {
	private Workflow workflow = null;
	private Comet comet = new Comet("sge"); // should be specified in config file
	private int cores = 8; // should be specified in config file
	
	private String bwa_index = null;
	private int combinedReadLength = 100; // should be estimated at runtime
	
	private String picard_tools = "java -jar /share/apps/picard-tools/current";
	
	private String sampleName = null;
	
	private String forwardReads = null;
	private String reverseReads = null;
	
	private long forwardReadCount = -1;
	private long reverseReadCount = -1;
	
	public RNASeq_AWS_Workflow(String name, String fqForward, String fqReverse) {
		sampleName = name;
		
		forwardReads = fqForward;
		reverseReads = fqReverse;

		workflow = WorkflowFactory.createNewWorkflow(name);
}

	private void countReads() {
		if ((forwardReadCount != -1) && (reverseReadCount != -1)) {
			return;
		}
		
		String[] cmds = new String[3];
		Job[] jobs = new Job[3];
		
		String[] countFiles = {"forward.readcount", "reverse.readcount"};
		cmds[0] = String.format("gunzip -c %s | wc -l | awk '{print $1/4}' > %s", forwardReads, countFiles[0]);
		cmds[1] = String.format("gunzip -c %s | wc -l | awk '{print $1/4}' > %s", reverseReads, countFiles[1]);
	
		jobs[0] = workflow.createJob("countForwardReads_"+sampleName, 1, cmds[0], forwardReads, countFiles[0]);
		jobs[1] = workflow.createJob("countReverseReads_"+sampleName, 1, cmds[1], reverseReads, countFiles[1]);
		
		cmds[2] = String.format("echo foward.readcount >> %s.ini ; echo reverse.readcount >> %s.ini", sampleName, sampleName);
		jobs[2] = workflow.createJob("collectReadCounts_"+sampleName, 1, cmds[2], countFiles, sampleName+".ini");
		jobs[2].addParent(jobs[0]);
		jobs[2].addParent(jobs[1]);
	}
	
	/* Estimate the insert size of paired-end reads
	 * Input: fastq files
	 * Output: mean inner distance and std dev
	 */
	private void estimateInsertSize() {
		String[] cmds = new String[7];
		Job[] jobs = new Job[7];
		
		// 1.  Create subset of reads (cut first 5 million reads)
		String[] subsetFastq = {sampleName+"_1_bwa.fastq", sampleName+"_2_bwa.fastq"};
		
		cmds[0] = String.format("gunzip -c %s | head -n 20000000 > %s", forwardReads, subsetFastq[0]);
		cmds[1] = String.format("gunzip -c %s | head -n 20000000 > %s", reverseReads, subsetFastq[1]);
		
		jobs[0] = workflow.createJob("subsetForward_"+sampleName, 1, cmds[0], forwardReads, subsetFastq[0]);
		jobs[1] = workflow.createJob("subsetReverse_"+sampleName, 1, cmds[1], reverseReads, subsetFastq[1]);
		
		// 2.  Align the reads
		String[] saiFiles = {sampleName+"_1_bwa.sai", sampleName+"_2_bwa.sai"};
		
		cmds[2] = String.format("bwa aln -t %d %s %s > %s", cores, bwa_index, subsetFastq[0], saiFiles[0]);
		cmds[3] = String.format("bwa aln -t %d %s %s > %s", cores, bwa_index, subsetFastq[1], saiFiles[1]);
		
		jobs[2] = workflow.createJob("bwa_aln_1_"+sampleName, cores, cmds[2], subsetFastq[0], saiFiles[0]);
		jobs[3] = workflow.createJob("bwa_aln_2_"+sampleName, cores, cmds[3], subsetFastq[1], saiFiles[1]);
		
		jobs[2].addParent(jobs[0]);
		jobs[1].addParent(jobs[1]);
		
		// 3.  sai 2 sam, sam 2 bam, sort bam
		String bamFile = sampleName + "_bwa.bam";
		cmds[4] = String.format("bwa sampe %s %s %s %s %s | samtools view -bS - | samtools sort - %s", bwa_index, saiFiles[0], saiFiles[1], subsetFastq[0], subsetFastq[1], bamFile);
		jobs[4] = workflow.createJob("bwa_sampe_"+sampleName, 3, cmds[4], saiFiles, bamFile);
		jobs[4].addParent(jobs[2]);
		jobs[4].addParent(jobs[3]);
				
		// 4. Determine insert metrics using PICARD
		String insertMetricsTxtFile = sampleName + "_Picard_Transcriptome_InsertMetrics.txt";
		cmds[5] = String.format("java -Xmx2g -jar %s/CollectInsertSizeMetrics.jar INPUT=%s OUTPUT=%s HISTOGRAM_FILE=%s_Picard_Transcriptome_InsertMetrics.pdf VALIDATION_STRINGENCY=SILENT", picard_tools, bamFile, insertMetricsTxtFile, sampleName);
		jobs[5] = workflow.createJob("CollectInsertSizeMetrics_"+sampleName, 1, cmds[5], bamFile, insertMetricsTxtFile);
		jobs[5].addParent(jobs[4]);
		
		// 5. Get the inner distance
		cmds[6] = String.format("head -n 8 %s | tail -n 1 | cut -f5 | awk -v var1=\"%d\" '{printf \"%.0f\\n\", $1-var1}' > %s.innerdistance", insertMetricsTxtFile, combinedReadLength, sampleName);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length != 2) {
			System.err.println("Usage: RNASeq_AWS_Workflow <sample same> <forward reads file> <reverse reads file>");
			System.exit(-1);
		}

		RNASeq_AWS_Workflow aws_workflow = new RNASeq_AWS_Workflow(args[0], args[1], args[2]);
				
		aws_workflow.countReads();
		
		aws_workflow.estimateInsertSize();
	}

}
