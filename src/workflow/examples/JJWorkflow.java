package workflow.examples;

import java.io.File;
import java.util.ArrayList;

import workflow.Job;
import workflow.Workflow;
import workflow.WorkflowFactory;

public class JJWorkflow {
	final static String java = "/usr/java/latest/bin/java";
	final static String picard = "/share/apps/picard/current/dist";
	final static String java_picard = java + " -jar " + picard;

	final static String project_dir = "/share/ngs/ryang/JJ";
	final static String bam_dir = "/share/ngs/ryang/JJ/bam";

	final static String refFlat = "/share/ngs/references/mm9/refFlat.txt";
	final static String reference = "/share/ngs/references/mm9/mm9.fa";
	
	public static void compareV2_V3(Workflow workflow) {
		File[] directories = new File(bam_dir).listFiles();
	    for (File dir : directories) {
	        if (dir.isDirectory()) {
	            
	            for (File bamFile : dir.listFiles()) {
            		if (bamFile.getName().endsWith(".bam")) {

	    	            if (bamFile.getName().contains("_v2")) {
		    	            System.err.println("Sample: " + dir.getName());
		    	            System.err.println("\t" + bamFile.getName());

		    	            // Collect some mapping stats
		    	            String sampleMetrics = bam_dir + "/" + dir.getName() + ".v2.alignment.metrics";
		    	            String cmd = String.format("%s/CollectAlignmentSummaryMetrics.jar LEVEL=ALL_READS BS=false I=%s O=%s R=%s", java_picard, bamFile.getAbsoluteFile(), sampleMetrics, reference);
		    	            Job collectAlignmentMetrics = workflow.createJob("alnmetrics_V2_"+dir.getName(), 1, cmd, bamFile.getAbsolutePath(), sampleMetrics);
		    	            
		    	            // Collect some RNA-Seq metrics (this is not dependent on indexing)
		    	            sampleMetrics = bam_dir + "/" + dir.getName() + ".v2.rg.metrics";
		    	    		String sampleMetricsPdf = bam_dir + "/" + dir.getName() + ".v2.rg.metrics.pdf";
		    	    		cmd = String.format("%s/CollectRnaSeqMetrics.jar REF_FLAT=%s CHART=%s LEVEL=SAMPLE I=%s O=%s R=%s STRAND=NONE",
		    	    							java_picard, refFlat, sampleMetricsPdf, bamFile, sampleMetrics, reference);
		    	    		Job collectRnaSeqMetrics = workflow.createJob("rnaseqmetrics_V2_"+dir.getName(), 1, cmd, bamFile.getAbsolutePath(), sampleMetrics);
	    	            }
		            	if (bamFile.getName().contains("_v3")) {
		    	            System.err.println("Sample: " + dir.getName());
		    	            System.err.println("\t" + bamFile.getName());

		    	            // Collect some mapping stats
		    	            String sampleMetrics = bam_dir + "/" + dir.getName() + ".v3.alignment.metrics";
		    	            String cmd = String.format("%s/CollectAlignmentSummaryMetrics.jar LEVEL=ALL_READS BS=false I=%s O=%s R=%s", java_picard, bamFile.getAbsoluteFile(), sampleMetrics, reference);
		    	            Job collectAlignmentMetrics = workflow.createJob("alnmetrics_V3_"+dir.getName(), 1, cmd, bamFile.getAbsolutePath(), sampleMetrics);

		    	            // Collect some RNA-Seq metrics (this is not dependent on indexing)
		    	    		sampleMetrics = bam_dir + "/" + dir.getName() + ".v3.rg.metrics";
		    	    		String sampleMetricsPdf = bam_dir + "/" + dir.getName() + ".v3.rg.metrics.pdf";
		    	    		cmd = String.format("%s/CollectRnaSeqMetrics.jar REF_FLAT=%s CHART=%s LEVEL=SAMPLE I=%s O=%s R=%s STRAND=NONE",
		    	    							java_picard, refFlat, sampleMetricsPdf, bamFile, sampleMetrics, reference);
		    	    		Job collectRnaSeqMetrics = workflow.createJob("rnaseqmetrics_V3_"+dir.getName(), 1, cmd, bamFile.getAbsolutePath(), sampleMetrics);
		    	    		
		            	}
	    	            
            		}
	            }
	        }
	    }
	    workflow.writeToSTDOUT();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Workflow workflow = WorkflowFactory.createNewWorkflow("J&J_Transcriptome_Coverage_Pipeline");

		//compareV2_V3(workflow);
		
		// 1.  For each sample
		File[] files = new File(bam_dir).listFiles();
	    for (File dir : files) {
	        if (dir.isDirectory()) {
	            System.err.println("Sample: " + dir.getName());
	            
	            // 2. Merge lane into a single bam file
	            String sampleBamFile = bam_dir + "/" + dir.getName() + ".bam";
	            
	            ArrayList<String> laneBams = new ArrayList<String>();
	            String inputFiles = new String();
	            for (File bamFile : dir.listFiles()) {
	            	if (bamFile.getName().endsWith(".bam")) {
	            		System.err.println("\tBAM File: " + bamFile.getName());
	            		inputFiles += "I="+bamFile.getAbsolutePath() + " ";
	            		laneBams.add(bamFile.getAbsolutePath());
	            	}
	            }  
	    		String cmd = String.format("%s/MergeSamFiles.jar %s O=%s MSD=true USE_THREADING=true", java_picard, inputFiles, sampleBamFile);
	    		Job mergeLanesJob = workflow.createJob("mergeLane_"+dir.getName(), 2, cmd, laneBams.toArray(new String[0]), sampleBamFile);

	    		// 3. Fix the Read Group information
	    		String sampleRGBamFile = bam_dir + "/" + dir.getName() + ".rg.bam";
	    		cmd = String.format("%s/AddOrReplaceReadGroups.jar I=%s O=%s SORT_ORDER=coordinate RGLB=1 RGPL=solid RGPU=1 RGSM=%s CREATE_INDEX=true", 
	    							java_picard, sampleBamFile, sampleRGBamFile, dir.getName());
	    		Job fixReadGroups = workflow.createJob("fixReadGroups_"+dir.getName(), 1, cmd, sampleBamFile, sampleRGBamFile);
	    		fixReadGroups.addParent(mergeLanesJob);
	    		
	    		// 4. Index the BAM file
//	    		cmd = String.format("samtools index %s", sampleRGBamFile, sampleRGBamFile+".bai");
//	    		Job index = workflow.createJob("index_"+dir.getName(), 1, cmd, sampleRGBamFile, sampleRGBamFile+".bai");
//	    		index.addParent(fixReadGroups);
	    		
	    		// 5. Collect alignment and RNA-Seq metrics
	    		String sampleMetrics = bam_dir + "/" + dir.getName() + ".alignment.metrics";
	    		cmd = String.format("%s/CollectAlignmentSummaryMetrics.jar LEVEL=SAMPLE I=%s O=%s R=%s", java_picard, sampleRGBamFile, sampleMetrics, reference);
	    		Job alignmentMetrics = workflow.createJob("alignmentmetrics_"+dir.getName(), 1, cmd, sampleRGBamFile, sampleMetrics);
	    		alignmentMetrics.addParent(fixReadGroups);
	    		
	    		// 6. Collect some RNA-Seq metrics (this is not dependent on indexing)
	    		sampleMetrics = bam_dir + "/" + dir.getName() + ".rnaseq.metrics";
	    		String sampleMetricsPdf = bam_dir + "/" + dir.getName() + ".rnaseq.metrics.pdf";
	    		
	    		cmd = String.format("%s/CollectRnaSeqMetrics.jar REF_FLAT=%s CHART=%s LEVEL=SAMPLE I=%s O=%s R=%s STRAND=NONE",
	    							java_picard, refFlat, sampleMetricsPdf, sampleRGBamFile, sampleMetrics, reference);
	    		Job collectRnaSeqMetrics = workflow.createJob("rnaseqmetrics_"+dir.getName(), 1, cmd, sampleRGBamFile, sampleMetrics);
	    		collectRnaSeqMetrics.addParent(fixReadGroups);
	    		
	        }
	    }
	    
	    workflow.writeToSTDOUT();
	}
}
