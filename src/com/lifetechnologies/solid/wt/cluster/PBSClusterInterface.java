package com.lifetechnologies.solid.wt.cluster;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * User: tuchbb
 * Date: Aug 14, 2008
 * Time: 9:56:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class PBSClusterInterface extends ClusterInterface {

	private static final Logger logger = Logger.getLogger(PBSClusterInterface.class.getName());
	
    public PBSClusterInterface(JobSubmissionParameters params) {
		super(params, "qsub", "qdel");
	}

    @Override
    protected String[] getCommand(File scriptFile) {
        return new String[] { this.pathToClusterSubmissionExe, scriptFile.getPath() };
    }

    protected File generateScript(File fileScript, ArrayList<String> commandStrings) throws Exception {

        if (fileScript == null)
            throw new Exception("You must specify a script file.");

        
        File fileScriptOutput = new File(fileScript.getPath() + SUFFIX_SCRIPT_OUTPUT_FILES);

        if (fileScriptOutput.exists())
            fileScriptOutput.delete();

        BufferedWriter writerScriptFile = new BufferedWriter(new FileWriter(fileScript));

        writeScriptHeader(writerScriptFile);        

        writerScriptFile.write("## define job name");
        writerScriptFile.newLine();
        writerScriptFile.write("#PBS -N " + fileScript.getName());
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        writerScriptFile.write("## export environment variables");
        writerScriptFile.newLine();
        writerScriptFile.write("#PBS -V ");
        writerScriptFile.newLine();
        writerScriptFile.newLine();
    
        if (this.jobSubmissionParameters.isRerunnable()) {
            writerScriptFile.write("## Declare job re-runnable");
            writerScriptFile.newLine();
            writerScriptFile.write("##PBS -r y");
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        } else {
            writerScriptFile.write("## Declare job not re-runnable");
            writerScriptFile.newLine();
            writerScriptFile.write("##PBS -r n");
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }

        writerScriptFile.write("## redirect output to file");
        writerScriptFile.newLine();
        writerScriptFile.write("#PBS -o '" + fileScriptOutput.getAbsolutePath() +  "'");
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        writerScriptFile.write("## join stdout and stderr");
        writerScriptFile.newLine();
        writerScriptFile.write("#PBS -j oe");
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        if (this.jobSubmissionParameters.getQueueName() != null) {
            writerScriptFile.write("## run on a particular queue");
            writerScriptFile.newLine();
            writerScriptFile.write("#PBS -q " + this.jobSubmissionParameters.getQueueName());
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }
        String resourceString = jobSubmissionParameters.getResourceString();
        if (resourceString != null) {
        	Long memoryRequirement = jobSubmissionParameters.getMemoryRequirement();
        	if (resourceString.contains("${pbs.mem}") || resourceString.contains("${pbs.vmem}")) {
        		if (memoryRequirement == null || memoryRequirement < 1 ) throw new Exception("No memory requirement specified.");
        		Long memoryRequirementMB = memoryRequirement / Constants.BYTES_PER_MEGABYTE;
        		if (memoryRequirementMB % Constants.BYTES_PER_MEGABYTE > 0) memoryRequirementMB++;
        		Long vmemRequirementMB = memoryRequirementMB * 2;
        		resourceString = resourceString.replaceAll("\\$\\{pbs.mem\\}", memoryRequirementMB + "mb");
        		resourceString = resourceString.replaceAll("\\$\\{pbs.vmem\\}", vmemRequirementMB + "mb");
        	}
        	writerScriptFile.write("## specify the resource string");
        	writerScriptFile.newLine();
        	writerScriptFile.write("#PBS -l " + resourceString);
        	writerScriptFile.newLine();
        	writerScriptFile.newLine();
        } 	

        String additionalOptions = jobSubmissionParameters.getAdditionalOptions();
        if (additionalOptions != null) {
        	additionalOptions = additionalOptions.trim();
        	if(additionalOptions.isEmpty() == false ) {
        		writerScriptFile.write("## Additional Options");
        		writerScriptFile.newLine();
        		writerScriptFile.write("#PBS "+ additionalOptions);
        		writerScriptFile.newLine();
        		writerScriptFile.newLine();
        	}
        }
        
        writerScriptFile.write("echo Directory is `pwd`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Running on host $PBS_O_HOST / `hostname`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Job $PBS_O_JOBID - $PBS_JOBNAME in Queue $PBS_QUEUE");
        writerScriptFile.newLine();

        writeScriptFooter(writerScriptFile, commandStrings);

        writerScriptFile.close();


        return fileScriptOutput;
    }

    protected String submitScript(File fileScriptFile) throws Exception {

        String idOfJob = null;

        logger.fine("Submitting script " + fileScriptFile.getPath() + " to cluster.");

        String[] cmd = getCommand(fileScriptFile);

        BufferedReader inputReader = null;
        BufferedReader errorReader = null;
        try {
            Process p = Runtime.getRuntime().exec(cmd);

            inputReader = new BufferedReader(new InputStreamReader(p.getInputStream()));

            String readLine;
            StringBuffer message = new StringBuffer();
            while ((readLine = inputReader.readLine()) != null) {
                message.append(readLine.concat("\n"));
                if (idOfJob == null)
                    idOfJob = readLine;
            }
            logger.fine(message.toString());
            
            message = new StringBuffer();
            errorReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((readLine = errorReader.readLine()) != null) {
            	message.append(readLine.concat("\n"));
            }
            errorReader.close();
        	if (message.length() > 0) logger.warning(message.toString());
        	
            inputReader.close();
            p.waitFor();
            if (p.exitValue() != 0) throw new Exception("Script submission failed.");
        } catch (Exception e) {
            if (inputReader != null)    inputReader.close();
            if (errorReader != null)    errorReader.close();
            throw e;
        }

        return idOfJob;
    }

    public String getJobStatus(String jobId) throws Exception {
    	String jobStatus = null;
    	String [] cmd = { "qstat", "-f", jobId };
    	
        BufferedReader inputReader = null;
        BufferedReader errorReader = null;
        try {
            Process p = Runtime.getRuntime().exec(cmd);

            inputReader = new BufferedReader(new InputStreamReader(p.getInputStream()));

            String readLine;
            while ((readLine = inputReader.readLine()) != null) {
            	if (readLine.contains("job_state")) {
            		jobStatus = readLine.substring(readLine.length()-1, readLine.length());
            		break;
            	}
            }

            StringBuffer message = new StringBuffer();
            errorReader = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((readLine = errorReader.readLine()) != null) {
            	message.append(readLine.concat("\n"));
            }
            errorReader.close();
        	if (message.length() > 0) logger.warning(message.toString());
        	
            inputReader.close();
            p.waitFor();
            if (p.exitValue() != 0) throw new Exception("Unable to check job state.");
        } catch (Exception e) {
            if (inputReader != null)    inputReader.close();
            if (errorReader != null)    errorReader.close();
            throw e;
        }

        return jobStatus;
    	
    }
}
