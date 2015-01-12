package workflow.comet;

import java.io.File;
import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;

import com.lifetechnologies.solid.wt.cluster.ClusterInterface;
import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;

import workflow.Job;
import workflow.Parser;
import workflow.Workflow;

/* Version 0.01: Initial version
 * 				 BUG: Multi-core requests only work on PBS, not SGE
 * 				 BUG: Job dpeendencies not working with SGE
 * Version 0.02a: Only specify multi-core requests for PBS, not SGE.
 * Version 0.02b: Fixed job dependency with SGE
 */
public class Comet {
	public static String VERSION = "0.02b";

	Job[] jobArray;
	String clusterType;

	public Comet(String clusterType) {
		this.clusterType = clusterType;
		this.jobArray = null;		
	}

	public Comet(String clusterType, Job[] jobArray) {
		this.clusterType = clusterType;
		this.jobArray = jobArray;
		
		// Get the cluster interface
//		JobSubmissionParameters parms = new JobSubmissionParameters();
//		parms.setEnvironment(clusterType);
//		clusterInterface = ClusterInterface.getClusterInterface(parms);
	}
	
	private void printExecutionTree() {
		int executionlevel = 1;
		int currentlevel = jobArray[0].getOrder();
		
		System.out.print("Execution Level " + executionlevel + ": ");
		for (int i = 0; i < jobArray.length; i++) {
			if (jobArray[i].getOrder() == currentlevel) {
				System.out.print(jobArray[i].getName() + " ");
			} else {
				System.out.println();
				executionlevel++;
				System.out.print("Execution Level " + executionlevel + ": ");
				currentlevel = jobArray[i].getOrder();
				System.out.print(jobArray[i].getName() + " ");
			}			
		}
		System.out.println();
	}
	
	/* Submit a command immediately to the cluster and return the job id */
	public String submitJob(String jobName, int cores, String cmd) {
		JobSubmissionParameters params = new JobSubmissionParameters();
		params.setEnvironment(clusterType);

		// Now set the cpu/core requirement
		if (clusterType.equals("pbs"))
			params.setResourceString("nodes=1:ppn="+cores);
		if (clusterType.equals("sge")) {
			params.setResourceString("pe_smp="+cores);
		}

		ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);

		ArrayList<String> jobCmds = new ArrayList<String>();
		jobCmds.add(cmd);
		
		String jobID = null;
		try {
			jobID = clusterInterface.executeJob(new File(jobName + ".sh"), jobCmds);
			//if (clusterType.equals("pbs"))
			//	return jobID;
			if (clusterType.equals("sge")) {
				jobID = jobID.substring(9, jobID.indexOf('('));
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return jobID;
	}
	
	/* 
	 * Algorithm:
	 * 1.  Determine if the jobs needs to be run based on timestamp of input/output files
	 *     as well as lack of existence of output files
	 * 2a. If the job does not have to be run, mark it as such
	 * 2b. If the job has to be run, submit job to cluster.
	 *     In doing so, make sure parent jobs have also been submitted, else we have a problem.
	 * 
	 */
	private boolean submitJob(Job job) {
		if (job.isRunnable() == false) {
			// 2a.  Job does not need to be run as the output files are newer than the input files
			//      or the output files don't exist.
			return true;
		}
		
		// 2b. Submit job to cluster 
		JobSubmissionParameters params = new JobSubmissionParameters();
		ArrayList<String> depends = new ArrayList<String>();

		params.setEnvironment(clusterType);
		
		// First, check if jobs parents have been submitted and if so, record their jobIDs
		for (Job parentJob: job.getParentJobs()) {
			// make sure the parent was submitted if it was runnable.  If it wasn't submit it and get its job id.
			if (parentJob.isRunnable()) {
				if (parentJob.getJobID() == null) {
					System.err.println("Error: Job " + parentJob.getName() + " should have been submitted before " + job.getName());
					return false;
					//submitJob(parentJob);
				}
				depends.add(parentJob.getJobID());
			}
		}

		if (depends.size() > 0) {
			// I shouldn't have to specify the parameter here.  This should be in the ClusterInface specific implementation
			if (clusterType.equals("pbs"))
				params.setAdditionalOptions("-W depend=afterok:" + StringUtils.join(depends, ":"));
			if (clusterType.equals("sge")) {
				params.setAdditionalOptions("-hold_jid " + StringUtils.join(depends, ":"));
			}
		}
		
		// Now set the cpu/core requirement
		if (clusterType.equals("pbs"))
			params.setResourceString("nodes=1:ppn="+job.getCores());
		if (clusterType.equals("sge")) {
			// TODO: Implement multiple core specs here 
		}

			
		try {
			ClusterInterface clusterInterface = ClusterInterface.getClusterInterface(params);
			// NOTE: If the parent job(s) finish before we submit this job, this job will never
			// execute because PBS does not know who the parent jobid is.  This can happen if a 
			// parent job fails on submission.
			String jobID = clusterInterface.executeJob(new File(job.getName() + ".sh"), job.getCommandStrings());
			if (clusterType.equals("pbs"))
				job.setJobID(jobID);
			if (clusterType.equals("sge")) {
				job.setJobID(jobID.substring(9, jobID.indexOf('(')));
			}
		} catch (Exception e) {
			e.printStackTrace();
			return false;
		}
		return true;
	}

	private int runExecutionTree() {		
		// Iterate through the jobs and execute them in order
		int executionlevel = 1;
		int currentlevel = jobArray[0].getOrder();

		System.out.print("Execution Level " + executionlevel + ": ");
		for (int i = 0; i < jobArray.length; i++) {
			Job job = jobArray[i];

			// this is really just for printing purposes
			if (job.getOrder() != currentlevel) {
				System.out.println();
				executionlevel++;
				System.out.print("Execution Level " + executionlevel + ": ");
				currentlevel = job.getOrder();
			}
			
			// run this job and record its job id
			if (submitJob(job))
				System.out.print(job.getName() + "(" + job.getJobID() + ") ");
			else {
				System.err.println("Faied to submit job: " + job.getName());
				return -1;
			}
		}
		System.out.println();
		
		return 0;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 4) {
           System.err.println("Usage: java Comet -c <pbs|sge|lsf> -w <workflow.xml> [-p]");
            System.exit(1);
        }
		System.out.println("Version " + VERSION);
		
		String clusterType = args[1];
		String workflowFile = args[3];
		System.out.println("Cluster Type: " + clusterType);
		System.out.println("Workflow File: " + workflowFile);
		
		Parser parser = new Parser();
		Workflow workflow = parser.parse(workflowFile);
		Job[] jobArray = workflow.getOrderedJobExecutionArray();

		// To check for cycles, use
		// http://en.wikipedia.org/wiki/Topological_sorting
		// For job scheduling, use
		// http://en.wikipedia.org/wiki/Coffman–Graham_algorithm
		// In the end, use org.apache.avalon.fortress.util.dag to create a DAG
		Comet comet = new Comet(clusterType, jobArray);
		
		if ((args.length == 5) && (args[4].equalsIgnoreCase("-p")))
			comet.printExecutionTree();
		else {
			int retVal = comet.runExecutionTree();
			if (retVal == 0)
				System.out.println("Workflow successful");
			else
				System.err.println("Workflow failed: " + retVal);
		}
	}
}
