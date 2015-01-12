package com.lifetechnologies.solid.wt.cluster;

import java.util.ArrayList;
import java.util.Iterator;
import java.io.*;


/**
 * User: tuchbb
 * Date: Aug 14, 2008
 * Time: 9:56:39 AM
 * Revision: $Rev$
 * This code was originally developed as part of the SOLiD Whole Transcriptome package.
 */
public abstract class ClusterInterface {

	protected JobSubmissionParameters jobSubmissionParameters;
	
    protected String pathToClusterSubmissionExe;
    protected String pathToClusterRemovalExe;

    public static final String SUFFIX_SCRIPT_OUTPUT_FILES = ".out";

    protected ArrayList<File> logOfScriptFilesGenerated = new ArrayList<File>();
    protected ArrayList<File> logOfScriptOutputFiles = new ArrayList<File>();
    protected ArrayList<String> logOfJobSubmissionIds = new ArrayList<String>();

    private String pathToTemporaryLocalFolderWhereCommandsWillBeExecuted = null;

    public ClusterInterface(JobSubmissionParameters params, String pathToClusterSubmissionExe, String pathToClusterRemovalExe) {
    	this.jobSubmissionParameters = params;
    	this.pathToClusterSubmissionExe = pathToClusterSubmissionExe;
    	this.pathToClusterRemovalExe = pathToClusterRemovalExe;
    }
    
    public static ClusterInterface getClusterInterface(JobSubmissionParameters jobSubmissionParameters) {
    	String environment = jobSubmissionParameters.getEnvironment();
    	if (environment.equalsIgnoreCase("pbs")) {
    		return new PBSClusterInterface(jobSubmissionParameters);
    	} else if (environment.equalsIgnoreCase("lsf")) {
    		return new LSFClusterInterface(jobSubmissionParameters);
    	} else if (environment.equalsIgnoreCase("sge")) {
    		return new SGEClusterInterface(jobSubmissionParameters);
    	} else {
    		throw new IllegalArgumentException("Unrecognized Cluster environment: "+ environment);
    	}
    }
    
    public void setLocalTemporaryWorkingPath(String pathToFolderWhereCommandsWillBeExecuted) {
        this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted = pathToFolderWhereCommandsWillBeExecuted;
    }


    /**
     * @param fileScriptFile
     * @param commandStrings
     * @return the File where script output will be written when script is executed (this file should allow one to check for job completion)
     * @throws Exception
     */
    protected abstract File generateScript(File fileScriptFile, ArrayList<String> commandStrings) throws Exception;

    protected abstract String[] getCommand(File scriptFile) throws Exception;
    
    protected abstract String submitScript(File fileScriptFile) throws Exception ;

    protected void writeScriptHeader(BufferedWriter writerScriptFile) throws IOException {
        writerScriptFile.write("#!/bin/bash");
        writerScriptFile.newLine();
        writerScriptFile.newLine();
    }

    protected void writeScriptFooter(BufferedWriter writerScriptFile, ArrayList<String> commandStrings) throws IOException {
        writerScriptFile.write("date '+TS[JOB_START]: %Y-%m-%d %k:%M:%S.%N'");
        writerScriptFile.newLine();
        writerScriptFile.write("echo StartTime is `date`");
        writerScriptFile.newLine();
        writerScriptFile.write("set +e -x");
        writerScriptFile.newLine();
        if (this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted != null) {
            writerScriptFile.write("mkdir " + this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted);
            writerScriptFile.newLine();
            writerScriptFile.write("cd " + this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted);
            writerScriptFile.newLine();
        }
        writerScriptFile.newLine();

        writerScriptFile.write("## =============== COMMANDS ARE RUN HERE ============================");
        writerScriptFile.newLine();
        Iterator<String> commandStringsIterator = commandStrings.iterator();
        while (commandStringsIterator.hasNext()) {
            writerScriptFile.write(commandStringsIterator.next());
            writerScriptFile.newLine();
        }
        writerScriptFile.write("## ==================================================================");
        writerScriptFile.newLine();
        writerScriptFile.newLine();

        writerScriptFile.write("RETVAL=$?");
        writerScriptFile.newLine();
        if (this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted != null) {
            writerScriptFile.write("rm -rf " + this.pathToTemporaryLocalFolderWhereCommandsWillBeExecuted);
            writerScriptFile.newLine();
        }
        writerScriptFile.write("set -e +x");
        writerScriptFile.newLine();
        writerScriptFile.write("date '+TS[JOB_END]: %Y-%m-%d %k:%M:%S.%N'");
        writerScriptFile.newLine();
        writerScriptFile.write("echo EndTime is `date`");
        writerScriptFile.newLine();
        writerScriptFile.write("exit $RETVAL");
        writerScriptFile.newLine();
    }

    public void writeMasterJobSubmissionFileFromLog(File fileMasterJobSubmissionScript) throws IOException {
        //fileMasterJobSubmissionScript.setExecutable(true, true);
        BufferedWriter writerMasterScriptFile = new BufferedWriter(new FileWriter(fileMasterJobSubmissionScript));
        Iterator<File> iteratorOverLoggedScriptFiles = this.logOfScriptFilesGenerated.iterator();
        while (iteratorOverLoggedScriptFiles.hasNext()) {
            File fileScript = iteratorOverLoggedScriptFiles.next();
            writerMasterScriptFile.write(pathToClusterSubmissionExe + " " + fileScript.getPath());
            writerMasterScriptFile.newLine();
        }
        writerMasterScriptFile.close();
    }

    public void writeMasterJobRemovalFileFromLog(File fileMasterJobRemovalScript) throws IOException {
        //fileMasterJobRemovalScript.setExecutable(true, true);
        BufferedWriter writerMasterScriptFile = new BufferedWriter(new FileWriter(fileMasterJobRemovalScript));
        Iterator<String> iteratorOverLoggedScriptFiles = this.logOfJobSubmissionIds.iterator();
        while (iteratorOverLoggedScriptFiles.hasNext()) {
            String idOfJob = iteratorOverLoggedScriptFiles.next();
            writerMasterScriptFile.write(pathToClusterRemovalExe + " " + idOfJob);
            writerMasterScriptFile.newLine();
        }
        writerMasterScriptFile.close();
    }

    public void writeMasterListOfJobOutputFilesFromLog(File fileMasterListOfJobOutputFiles) throws IOException {
        BufferedWriter writerMasterListOfScriptOutputFiles = new BufferedWriter(new FileWriter(fileMasterListOfJobOutputFiles));
        Iterator<File> iteratorOverLoggedScriptOutputFiles = this.logOfScriptOutputFiles.iterator();
        while (iteratorOverLoggedScriptOutputFiles.hasNext()) {
            File fileScriptOutput = iteratorOverLoggedScriptOutputFiles.next();
            writerMasterListOfScriptOutputFiles.write(fileScriptOutput.getPath());
            writerMasterListOfScriptOutputFiles.newLine();
        }
        writerMasterListOfScriptOutputFiles.close();

    }

    public void clearLoggedJobs() {
        this.logOfScriptFilesGenerated.clear();
        this.logOfScriptOutputFiles.clear();
    }   

    public ArrayList<File> getLogOfScriptFilesGenerated() {
        return logOfScriptFilesGenerated;
    }

    public String executeJob(File fileScript, ArrayList<String> commandStrings) throws Exception {
        File fileJobOutput = generateScript(fileScript, commandStrings);
        String idOfSubmittedJob = submitScript(fileScript);

        this.logOfScriptFilesGenerated.add(fileScript);
        this.logOfScriptOutputFiles.add(fileJobOutput);
        this.logOfJobSubmissionIds.add(idOfSubmittedJob);

        //return fileJobOutput;
        return idOfSubmittedJob;
    }

    /**
     * A job is complete if the script output file (copy of stdout) has been written.
     *
     * @return
     */
    public boolean checkIfLoggedJobsComplete()  {
        boolean allJobsComplete = true;
        Iterator<File> iteratorOverLoggedScriptOutputFiles = logOfScriptOutputFiles.iterator();
        while (allJobsComplete && iteratorOverLoggedScriptOutputFiles.hasNext()) {
            File fileScriptOutput = iteratorOverLoggedScriptOutputFiles.next();
            if (!fileScriptOutput.exists())
                allJobsComplete = false;
        }
        return allJobsComplete;
    }

    /**
     * A job has completed successfully if the script output file (copy of stdout) has been written and indicates
     * success.  If one or more script output files do not exist or do not indicate success (which is the
     * case if the job submitted fails or if it is deleted) then false is returned.
     *
     * @return
     * @throws Exception
     */
    public boolean checkIfLoggedJobsCompletedSuccessfully() throws Exception {
        boolean allJobsCompletedSuccessfully = true;
        Iterator<File> iteratorOverLoggedScriptOutputFiles = logOfScriptOutputFiles.iterator();
        while (allJobsCompletedSuccessfully && iteratorOverLoggedScriptOutputFiles.hasNext()) {
            File fileScriptOutput = iteratorOverLoggedScriptOutputFiles.next();
            if (fileScriptOutput.exists()) {
                BufferedReader readerScriptOutput = new BufferedReader(new FileReader(fileScriptOutput));
                boolean currentJobComplete = false;
                String line = null;
                while (!currentJobComplete && (line = readerScriptOutput.readLine()) != null)
                    if (line.startsWith("+ RETVAL=0"))  // RETVAL=0 is our inidcator that the submitted job succeeded
                        currentJobComplete = true;

                readerScriptOutput.close();
                allJobsCompletedSuccessfully &= currentJobComplete;
            } else
                allJobsCompletedSuccessfully = false;
        }

        return allJobsCompletedSuccessfully;
    }

    /**
     * A job has failed if the script output file (copy of stdout) has been written and it does not cotain
     * an indicator of success.
     *
     * @return A list of logged script output file that exist, but which did not contain the indicator
     * of success.
     *
     */
    public ArrayList<File> getListOfScriptOutputFilesForLoggedJobsThatIndicateFailure() throws Exception {

        return getListOfScriptOutputFilesIndicatingFailure(this.logOfScriptOutputFiles);

    }


    /**
     * A job has failed if the script output file (copy of stdout) has been written and it does not cotain
     * an indicator of success.  Each script output file that exists, but which does not contain the indicator
     * of success is returned.
     *
     * @param listOfScriptOutputFiles
     * @return A list, which is a sublist of the list passed to this method, of script output files that exist,
     * but which did not contain the indicator of success.
     * @throws IOException
     */
    public static ArrayList<File> getListOfScriptOutputFilesIndicatingFailure(ArrayList<File> listOfScriptOutputFiles) throws Exception {

        ArrayList<File> listOfJobsThatFailed = new ArrayList<File>();
        Iterator<File> iteratorOverLoggedScriptOutputFiles = listOfScriptOutputFiles.iterator();
        while (iteratorOverLoggedScriptOutputFiles.hasNext()) {
            File fileScriptOutput = iteratorOverLoggedScriptOutputFiles.next();
            if (fileScriptOutput.exists()) {

                String statementRETVAL = extractRetvalStatementFromFile(fileScriptOutput);

                if (statementRETVAL == null)
                    Thread.sleep(10000);

                statementRETVAL = extractRetvalStatementFromFile(fileScriptOutput);

                // RETVAL=0 is our inidcator that the submitted job succeeded
                if (statementRETVAL == null || !statementRETVAL.startsWith("+ RETVAL=0"))
                    listOfJobsThatFailed.add(fileScriptOutput);
            }
        }


        return listOfJobsThatFailed;
    }

    private static String extractRetvalStatementFromFile(File file) throws IOException {
        String statementRETVAL = null;
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = null;
        while (statementRETVAL == null && (line = reader.readLine()) != null)
            if (line.startsWith("+ RETVAL="))
                statementRETVAL = line.trim();
        reader.close();

        return statementRETVAL;
    }

    
    public JobSubmissionParameters getJobSubmissionParameters() {
    	return jobSubmissionParameters;
    }
}
