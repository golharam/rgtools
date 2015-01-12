package workflow.examples;

import java.io.File;
import java.util.ArrayList;

import workflow.Job;
import workflow.Workflow;
import workflow.WorkflowFactory;

public class SHWTSWorkflow {
	final static String picard_dir = "~/workspace/picard-tools";
	final static String picard = "java -jar " + picard_dir;
	
	final static String project_dir = "/mnt/isilon/cag/ngs/hiseq/golharr/SHWTS";
	final static String bam_dir = project_dir + "/BAM";
	
	private static Job mergeLanes(Workflow workflow, String sampleName) {
		String sampleBamFile = bam_dir + "/" + sampleName + ".bam";
        String sampleBamFilesDir = bam_dir + "/SOLiD_BAM_Files/" + sampleName;

        ArrayList<String> bamFiles = new ArrayList<String>();
        String inputLaneFiles = "";
        File dir = new File(sampleBamFilesDir);
        for (File laneBamFile : dir.listFiles()) {
        	inputLaneFiles += "I=" + laneBamFile.getAbsolutePath() + " ";
        	bamFiles.add(laneBamFile.getAbsolutePath());
        }

        String cmd = String.format("%s/MergeSamFiles.jar %s O=%s SO=coordinate MSD=true USE_THREADING=true", picard, inputLaneFiles, sampleBamFile);
        return workflow.createJob(sampleName + "_MergeSamFiles", 2, cmd, bamFiles.toArray(new String[0]), sampleBamFile);        
	}
	
	private static void fixReadGroup(Workflow workflow, String sampleName, Job lastJob) {
		String sampleBamFile = bam_dir + "/" + sampleName + ".bam";
		String outputFile = bam_dir + "/" + sampleName + ".rg.bam";
		String cmd = String.format("%s/AddOrReplaceReadGroups.jar I=%s O=%s SO=coordinate LB=%s PL=illumina PU=101pe SM=%s", 
									picard, sampleBamFile, outputFile, sampleName, sampleName);
		Job fixRGJob = workflow.createJob(sampleName + "_FixReadGroup", 1, cmd, sampleBamFile, outputFile);
		fixRGJob.addParent(lastJob);
	}
		
	/* This needs to be run on respublica */
	public static void main(String[] args) {
		Workflow workflow = WorkflowFactory.createNewWorkflow("SHWTS_Pipeline");
		
		for (int i = 1; i <= 8; i++) {
			String sampleName = String.format("SHWTS0%s", i);
			
            // 1.  Merge the BAM files per lane into 1 sample BAM file
    		Job mergeLanesJob = mergeLanes(workflow, sampleName);

    		// 2.  Fix Read Group information
    		fixReadGroup(workflow, sampleName, mergeLanesJob);
		}
	            	    
	    workflow.writeToSTDOUT();
	}
}
 