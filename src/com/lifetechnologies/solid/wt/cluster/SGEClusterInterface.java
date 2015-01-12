package com.lifetechnologies.solid.wt.cluster;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.logging.Logger;

/**
 * User: tuchbb
 * Date: Aug 14, 2008
 * Time: 9:56:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class SGEClusterInterface extends ClusterInterface {

	private Logger logger = Logger.getLogger(SGEClusterInterface.class.getName());
	
    public SGEClusterInterface(JobSubmissionParameters params) {
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
        writerScriptFile.write("## Force bash shell.");
        writerScriptFile.newLine();
        writerScriptFile.write("#$ -S /bin/bash");
        writerScriptFile.newLine();
        writerScriptFile.write("## define job name");
        writerScriptFile.newLine();
        writerScriptFile.write("#$ -N " + fileScript.getName());
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        writerScriptFile.write("## export environment variables");
        writerScriptFile.newLine();
        writerScriptFile.write("#$ -V ");
        writerScriptFile.newLine();
        writerScriptFile.newLine();
    
        if (this.jobSubmissionParameters.isRerunnable()) {
            writerScriptFile.write("## Declare job re-runnable");
            writerScriptFile.newLine();
            writerScriptFile.write("#$ -r y");
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        } else {
            writerScriptFile.write("## Declare job not re-runnable");
            writerScriptFile.newLine();
            writerScriptFile.write("#$ -r n");
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }

        writerScriptFile.write("## redirect output to file");
        writerScriptFile.newLine();
        writerScriptFile.write("#$ -o '" + fileScriptOutput.getAbsolutePath() +  "'");
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        writerScriptFile.write("## join stdout and stderr");
        writerScriptFile.newLine();
        writerScriptFile.write("#$ -j yes");
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        if (this.jobSubmissionParameters.getQueueName() != null) {
            writerScriptFile.write("## run on a particular queue");
            writerScriptFile.newLine();
            writerScriptFile.write("#$ -q " + this.jobSubmissionParameters.getQueueName());
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }
        String resourceString = jobSubmissionParameters.getResourceString();
        if (resourceString != null) {
        	Long memoryRequirement = jobSubmissionParameters.getMemoryRequirement();
        	if (resourceString.contains("${sge.mem}") || resourceString.contains("${sge.vmem}")) {
        		if (memoryRequirement == null || memoryRequirement < 1 ) throw new Exception("No memory requirement specified.");
        		Long memoryRequirementMB = memoryRequirement / Constants.BYTES_PER_MEGABYTE;
        		if (memoryRequirement % Constants.BYTES_PER_MEGABYTE > 0) memoryRequirementMB++;
        		Long vmemRequirementMB = memoryRequirementMB * 2;
        		resourceString = resourceString.replaceAll("\\$\\{sge.mem\\}", memoryRequirementMB + "M");
        		resourceString = resourceString.replaceAll("\\$\\{sge.vmem\\}", vmemRequirementMB + "M");
        	}
        	writerScriptFile.write("## specify the resource string");
        	writerScriptFile.newLine();
        	writerScriptFile.write("#$ -l " + resourceString);
        	writerScriptFile.newLine();
        	writerScriptFile.newLine();
        } 	

        String additionalOptions = jobSubmissionParameters.getAdditionalOptions();
        if (additionalOptions != null) {
        	additionalOptions = additionalOptions.trim();
        	if(additionalOptions.isEmpty() == false ) {
        		writerScriptFile.write("## Additional Options");
        		writerScriptFile.newLine();
        		writerScriptFile.write("#$ "+ additionalOptions);
        		writerScriptFile.newLine();
        		writerScriptFile.newLine();
        	}
        }
        
        writerScriptFile.write("echo Directory is `pwd`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Running on host $SGE_O_HOST / `hostname`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Job $JOB_ID in Queue $QUEUE");
        writerScriptFile.newLine();

        writeScriptFooter(writerScriptFile, commandStrings);

        writerScriptFile.close();


        return fileScriptOutput;
    }

    protected String submitScript(File fileScriptFile) throws Exception {

        String idOfJob = null;

        logger.info("Submitting script " + fileScriptFile.getPath() + " to cluster.");

        String[] cmd = getCommand(fileScriptFile);

        //if (verbose) System.out.println(cmd);
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
            logger.info(message.toString());
            
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
    
    /**
     * A job is complete if the script output file (copy of stdout) has been written.
     * and contains the "+ RETVAL= string
     *
     * @return
     */
    @Override
    public boolean checkIfLoggedJobsComplete()  {
        boolean allJobsComplete = true;
        Iterator<File> iteratorOverLoggedScriptOutputFiles = logOfScriptOutputFiles.iterator();
        while (allJobsComplete && iteratorOverLoggedScriptOutputFiles.hasNext()) {
            File fileScriptOutput = iteratorOverLoggedScriptOutputFiles.next();
            if (!fileScriptOutput.exists()) return false;
            BufferedReader reader = null;
            boolean foundRetval = false;
            try {
	            try {
	            	reader = new BufferedReader(new InputStreamReader(new FileInputStream(fileScriptOutput)));
	            	for (String line = reader.readLine(); line != null; line = reader.readLine()) {
	            		if (line.startsWith("+ RETVAL")) foundRetval = true;
	            	}
	            } finally {
	            	if (reader != null) reader.close();
	            }
            } catch (IOException e) {
            	e.printStackTrace();
            	return false;
            }
            if (foundRetval == false) return false;
        }
        return allJobsComplete;
    }
}
