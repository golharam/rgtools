package com.lifetechnologies.solid.wt.cluster;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * User: tuchbb
 * Date: Aug 14, 2008
 * Time: 9:56:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public class LSFClusterInterface extends ClusterInterface {

	private static Logger logger = Logger.getLogger(LSFClusterInterface.class.getName());
	
    public LSFClusterInterface(JobSubmissionParameters jobSubmissionParameters) {
    	//TODO implement and test job deletion.
        super(jobSubmissionParameters, "bsub", "bkill");
    }

    @Override
    protected String[] getCommand(File scriptFile) throws Exception {
    	File f = new File(scriptFile.getParent(), scriptFile.getName()+".runner.sh");
    	try {
    		PrintWriter writer = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new FileOutputStream(f))));
    		writer.println("#!/bin/bash");
    		writer.println(this.pathToClusterSubmissionExe + " < " + scriptFile.getPath());
    		writer.close();
    		f.setExecutable(true);
    	} catch (FileNotFoundException e) {
    		throw new Exception("Failed to create runner script.");
    	}
    	return new String[] { f.getPath() };
    }


    protected File generateScript(File fileScript,
                               ArrayList<String> commandStrings) throws Exception {

        if (fileScript == null)
            throw new Exception("You must specify a script file.");

        
        File fileScriptOutput = new File(fileScript.getPath() + SUFFIX_SCRIPT_OUTPUT_FILES);

        if (fileScriptOutput.exists())
            fileScriptOutput.delete();

        BufferedWriter writerScriptFile = new BufferedWriter(new FileWriter(fileScript));

        writeScriptHeader(writerScriptFile);        

        writerScriptFile.write("## define job name");
        writerScriptFile.newLine();
        writerScriptFile.write(String.format("#BSUB -J '%s'", fileScript.getName()));
        writerScriptFile.newLine();
        writerScriptFile.newLine();
    
        if (this.jobSubmissionParameters.isRerunnable()) {
            writerScriptFile.write("## Declare job re-runnable");
            writerScriptFile.newLine();
            writerScriptFile.write("##BSUB -r");
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }

        writerScriptFile.write("## redirect output to file");
        writerScriptFile.newLine();
        writerScriptFile.write(String.format("#BSUB -o '%s'", fileScriptOutput.getAbsolutePath()));
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        String queueName = jobSubmissionParameters.getQueueName();
        if (queueName != null) {
            writerScriptFile.write("## run on a particular queue");
            writerScriptFile.newLine();
            writerScriptFile.write(String.format("#BSUB -q '%s'", queueName));
            writerScriptFile.newLine();
            writerScriptFile.newLine();
        }
               
        if (jobSubmissionParameters.getResourceString() != null) {
        	String resourceString = jobSubmissionParameters.getResourceString();
        	if (resourceString.contains("${lsf.select.mem}") ||
        		resourceString.contains("${lsf.select.swp}") ||
        		resourceString.contains("${lsf.rusage.mem}") ||
        		resourceString.contains("${lsf.rusage.swp}") ){
        		Long memoryRequirement = jobSubmissionParameters.getMemoryRequirement();
        		if (memoryRequirement == null && memoryRequirement < 1) throw new Exception("No Memory requirement specified");
        		Long memoryRequirementMB = memoryRequirement / Constants.BYTES_PER_MEGABYTE;
        		if (memoryRequirement % Constants.BYTES_PER_MEGABYTE > 0) memoryRequirementMB++;
        		Long swapRequirementMB = memoryRequirementMB * 2;
        		resourceString = resourceString.replaceAll("\\$\\{lsf.select.mem\\}", memoryRequirementMB.toString());
        		resourceString = resourceString.replaceAll("\\$\\{lsf.select.swp\\}", swapRequirementMB.toString());
        		resourceString = resourceString.replaceAll("\\$\\{lsf.rusage.mem\\}", memoryRequirementMB.toString());
        		resourceString = resourceString.replaceAll("\\$\\{lsf.rusage.swp\\}", swapRequirementMB.toString());
        	}
        	writerScriptFile.write("## specify the resource string");
        	writerScriptFile.newLine();
        	writerScriptFile.write("#BSUB -R \"" + resourceString + "\"");
        	writerScriptFile.newLine();
        	writerScriptFile.newLine();

        }
        
        String additionalOptions = jobSubmissionParameters.getAdditionalOptions();
        if (additionalOptions != null) {
        	additionalOptions = additionalOptions.trim();
        	if(additionalOptions.isEmpty() == false ) {
        		writerScriptFile.write("## Additional Options");
        		writerScriptFile.newLine();
        		writerScriptFile.write("#BSUB "+ additionalOptions);
        		writerScriptFile.newLine();
        		writerScriptFile.newLine();
        	}
        }

        writerScriptFile.write("echo Directory is `pwd`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Running on host $LSB_HOSTS / `hostname`");
        writerScriptFile.newLine();
        writerScriptFile.write("echo Job $LSB_JOBID - $LSB_JOBNAME in Queue $LSB_QUEUE");
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
              //Grab first string of digits as the job id.
                if (idOfJob == null) {
                	Matcher matcher = Pattern.compile("\\d+").matcher(readLine);
                	if (matcher.find())           		
                		idOfJob = matcher.group();
                }
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

}
