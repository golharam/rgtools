#!/usr/bin/env python
import drmaa
import os
import re
import shutil
import sys
import subprocess

fastQValidator_app = "/ngs/ngs15/golharr/apps/fastQValidator/bin/fastQValidator"
fastQC_app = "/ngs/ngs15/golharr/apps/FastQC-0.11.2/fastqc"
tophat2_app = "/ngs/ngs15/golharr/apps/tophat-2.0.13.Linux_x86_64/tophat2"

# Pipeline functions
def removeFiles(files):
    for file in files:
        subprocess.call(["/bin/rm", "-rf", file])

def runCommand(command, outputFiles):
    retCode = 0
    try:
        outputFile = outputFiles[0]
        print ' '.join(command) + " 2>&1 > " + outputFile
        with open(outputFile, 'w') as outfile:
            retCode = subprocess.call(command, stdout=outfile, stderr=subprocess.STDOUT)
    except KeyboardInterrupt:
        print "\nYou pressed Ctrl-C!!! Cleaning up."
        removeFiles(outputFiles)
        sys.exit(-1)

    print "Return Code: %d" % retCode
    if retCode != 0:
        removeFiles(outputFiles[1:])
 
    return retCode

def validateFastQ(inputFile, outputFile):
    retCode = 0
    if (os.access(outputFile, os.F_OK) == False):
        retCode = runCommand([fastQValidator_app, "--file", inputFile], [outputFile])
    return retCode

def fastQC(fq1, fq2, outputFile):
    # strip off .gz and and add _fastqc.html and _fastqc.zip
    html_outputfile = fq1.strip('.gz') + "_fastqc.html"
    zip_outpufile = fq1.strip('.gz') + "_fastqc.zip"

    outputFiles = [outputFile, fq1.strip('.gz') + "_fastqc.html", fq1.strip('.gz') + "_fastqc.zip", fq2.strip('.gz') + "_fastqc.html", fq2.strip('.gz') + "_fastqc.zip"]
    retCode = 0
    if (os.access(outputFile, os.F_OK) == False):
        retCode = runCommand([fastQC_app, "-t", "2", "-o", "analysis/", fq1, fq2], outputFiles)
    return retCode

def mapReads(fq1, fq2, sample):
    retCode = 0
    outdir = "analysis/tophat2_"+sample
    tmpdir = "/scratch2/"+sample
    outtxt = "analysis/" + sample + ".tophat2.txt"
    if (os.access(outdir, os.F_OK) == False):
        retCode = runCommand([tophat2_app, "-o", outdir, "-p 4", "--tmp-dir", tmpdir, "/ngs/ngs15/golharr/NGS/reference/hg19/reference/bowtie_index/hg19", fq1, fq2], [outtxt, outdir])

    if retCode == 0:
        os.symlink("analysis/tophat2_"+sample+"/accepted_hits.bam", "analysis/"+sample+".bam")

    return retCode

def main():
    print "RNASeq Pipeline v0.01\n"

    # 1.  Get the cohort of samples to process.  The name of the file is either 
    #     passed on the command line as the first argument, or is assumed to be
    #     in the current directory and named samples.txt
    samplesFile = "samples.txt"
    if (os.access(samplesFile, os.R_OK) == False):
        print "ERROR: Unable to read %s" % samplesFile
        print "Usage: RNASeqPipeline.py [-cohort <samples.txt>]\n"
        sys.exit(-1)
    print "    -cohort: %s\n" % samplesFile

    # 2.  Read in the list of samples including sample name, fastq1, fastq2 and
    #     build a sample list
    samples = {}
    with open(samplesFile) as f:
        for line in f:
            fields = re.split('\t', line.strip('\n'))
	    samples[fields[0]] = {'fq1' : fields[1], 'fq2' : fields[2]}

    # Sanity Check - Print out samples
    #for s in samples:
    #    print "%s\t%s\t%s" % (s, samples[s]['fq1'], samples[s]['fq2'])

    # Run the following pipeline on each sample
    for sample in samples:
        # 1.  Validate FastQ file
        output = "analysis/%s_1.fastQValidator.txt" % sample
        samples[sample]['fq1_validate'] = validateFastQ(samples[sample]['fq1'], output)
        output = "analysis/%s_2.fastQValidator.txt" % sample
        samples[sample]['fq2_validate'] = validateFastQ(samples[sample]['fq2'], output)

        # 2. Run FastQC
        output = "analysis/%s.fastqc.txt" % sample
        samples[sample]['fastqc'] = fastQC(samples[sample]['fq1'], samples[sample]['fq2'], output)

        # 3. If fastq validation is successful, map reads
        if samples[sample]['fq1_validate'] == 0 and samples[sample]['fq2_validate'] == 0:
            samples[sample]['mapReads'] = mapReads(samples[sample]['fq1'], samples[sample]['fq2'], sample)

        return 0

if __name__ == '__main__':
    with drmaa.Session() as SGE:
        main()

