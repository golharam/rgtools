package workflow.comet;

import java.io.File;
import java.util.ArrayList;

import com.lifetechnologies.solid.wt.cluster.JobSubmissionParameters;
import com.lifetechnologies.solid.wt.cluster.PBSClusterInterface;

public class PBSTest {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		JobSubmissionParameters jobParams = new JobSubmissionParameters();
		PBSClusterInterface pbs = new PBSClusterInterface(jobParams);
		
		ArrayList<String> commands = new ArrayList<String>();
		commands.add("uname -n");
		commands.add("sleep 20");
		commands.add("echo hello");
		String jobId = pbs.executeJob(new File("hello.sh"), commands);
		System.out.println("Job ID: " + jobId);
		
		Thread.sleep(5000);
		
		String jobStatus =  pbs.getJobStatus(jobId);
		System.out.println("Job Status: " + jobStatus);

		while (!jobStatus.equals("C")) {
			Thread.sleep(5000);
			jobStatus =  pbs.getJobStatus(jobId);
			System.out.println("Job Status: " + jobStatus);
		}
	}

}
