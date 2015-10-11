#!/bin/sh
# vim:et:ft=sh:sts=2:sw=2
#
# Author: ryan.golhar@bms.com
#
# Unit test for rnaseqqc_pipeline.sh
#

oneTimeSetUp()
{
  echo "SHUNIT_TMPDIR: ${SHUNIT_TMPDIR}"
  echo

  script="/home/golharr/workspace/rgtools/scripts/NGS_pipelines/RNASeq/QCPipeline/rnaseqqc_pipeline.sh"
  testDataDir="/home/golharr/workspace/rgtools/test_data/scripts/NGS_pipelines/RNASeq/QCPipeline"

  outputDir="${SHUNIT_TMPDIR}/output"
  mkdir "${outputDir}"
  stdoutF="${outputDir}/stdout"
  stderrF="${outputDir}/stderr"

  testDir="${SHUNIT_TMPDIR}/testRun_HumanSample_hg19ERCC_refSeq"
  outDir="${SHUNIT_TMPDIR}/analysis"
  mkdir ${outDir}

  # If you want to submit a previous analysis as a test,
  # set sampleName, copyResults=1, and resultsDir
  # else this script will re-run the analysis which can
  # take longer.
  sampleName="test1.hg19ERCC.refSeq"
  copyResults=0
  if [ ${copyResults} -eq 1 ]; then
    resultsDir="/home/golharr/ngsprojects/test/analysis/${sampleName}"
    echo "Copying ${resultsDir} to ${outDir}/"
    cp -r ${resultsDir} ${outDir}/
    rtrn_sampleRun=$?
  else 
    echo "Running ${sampleName}:"
    echo "${script} -s ${sampleName} \
            --fastq1 ${testDataDir}/test1_1.fastq.gz \
            --fastq2 ${testDataDir}/test1_2.fastq.gz \
            --sample-dir ${testDir} --outdir ${outDir} \
            > ${stdoutF} 2> ${stderrF}"
    ${script} -s ${sampleName} \
            --fastq1 ${testDataDir}/test1_1.fastq.gz \
            --fastq2 ${testDataDir}/test1_2.fastq.gz \
            --sample-dir ${testDir} --outdir ${outDir} \
            > ${stdoutF} 2> ${stderrF}
    rtrn_sampleRun=$?
  fi
}

oneTimeTearDown()
{
  echo "oneTimeTearDown() called"
  #rm -fr "${testDir}"
}

#-----------------------------------------------------------------------------
# suite tests
# All tests should check:
#  return value
#  stdout, stderr
#  and any files created

testRunWithNoParameters()
{
  ${script} > ${stdoutF} 2> ${stderrF}
  rtrn=$?

  # check return value
  assertTrue 'Help should be printed with no parameters' "[ ${rtrn} -eq 2 ]"

  # check stdout
  md1=`md5sum ${stdoutF} | cut -f 1 -d ' '`
  md2=`md5sum ${testDataDir}/helpscreen.stdout | cut -f 1 -d ' '`
  assertTrue 'Help screen don''t match' "[ '$md1' == '$md2' ]"
  unset md1 md2

  # check stderr
  assertFalse 'unexpected output to STDERR' "[ -s '${stderrF}' ]"

  # check output files
  ## no files created
}

test_exitCode()
{
  # check return value
  assertTrue 'Script exited with non-zero exit code.' "[ ${rtrn_sampleRun} -eq 0 ]"
}

test_fastQValidator()
{
  th_assertFileExists ${outDir}/${sampleName}/${sampleName}.fq1Validator.txt
  readCount=`head -n 1 ${outDir}/${sampleName}/${sampleName}.fq1Validator.txt | cut -f 8 -d ' '`
  testReadCount=`head -n 1 ${testDataDir}/${sampleName}.fq1Validator.txt | cut -f 8 -d ' '`
  assertTrue 'FQ1 Read count incorrect' "[ $readCount -eq $testReadCount ]"
  unset readCount testReadCount

  th_assertFileExists ${outDir}/${sampleName}/${sampleName}.fq2Validator.txt
  readCount=`head -n 1 ${outDir}/${sampleName}/${sampleName}.fq2Validator.txt | cut -f 8 -d ' '`
  testReadCount=`head -n 1 ${testDataDir}/${sampleName}.fq2Validator.txt | cut -f 8 -d ' '`
  assertTrue 'FQ1 Read count incorrect' "[ $readCount -eq $testReadCount ]"
  unset readCount testReadCount
}

test_FastQC()
{
  th_assertFileExists ${outDir}/${sampleName}/${sampleName}_1.subsampled_fastqc.zip
  th_assertFileExists ${outDir}/${sampleName}/${sampleName}_2.subsampled_fastqc.zip
}

test_ContaminationAlignmentMetrics()
{
  lines=`tail -n 7 ${outDir}/${sampleName}/contamination/${sampleName}.contaminated.alnMetrics.txt`
  testLines=`tail -n 7 ${testDataDir}/${sampleName}.contaminated.alnMetrics.txt`
  assertTrue 'Contamination Metrics do not match' "[ '${lines}' == '${testLines}' ]"
  unset lines testLines
}

test_AlignmentMetrics()
{
  lines=`tail -n 7 ${outDir}/${sampleName}/${sampleName}.alnMetrics.txt`
  testLines=`tail -n 7 ${testDataDir}/${sampleName}.alnMetrics.txt`
  assertTrue 'Alignment Metrics do not match' "[ '${lines}' == '${testLines}' ]"
  unset lines testLines
}

test_InsertSizeMetrics()
{
  th_assertFileIsNonZeroSize ${outDir}/${sampleName}/${sampleName}.insertSizeHistogram.pdf
  lines=`head -n 8 ${outDir}/${sampleName}/${sampleName}.insertSizeMetrics.txt | tail -n 1`
  testLines=`head -n 8 ${testDataDir}/${sampleName}.insertSizeMetrics.txt | tail -n 1`
  assertTrue 'InsertSize Metrics do not match' "[ '${lines}' == '${testLines}' ]"
  unset lines testLines
}

test_RNASeqMetrics()
{
  th_assertFileIsNonZeroSize ${outDir}/${sampleName}/${sampleName}.rnaseq.pdf
  lines=`head -n 8 ${outDir}/${sampleName}/${sampleName}.rnaseqMetrics.txt | tail -n 1`
  testLines=`head -n 8 ${testDataDir}/${sampleName}.rnaseqMetrics.txt | tail -n 1`
  assertTrue 'RNASeq Metrics do not match' "[ '${lines}' == '${testLines}' ]"
  unset lines testLines
}

test_HemoglobinMetrics()
{
  lines=`tail -n 90 ${outDir}/${sampleName}/${sampleName}.hemoglobinMetrics.txt`
  testLines=`tail -n 90 ${testDataDir}/${sampleName}.hemoglobinMetrics.txt`
  assertTrue 'Hemoglobin Metrics do not match' "[ '${lines}' == '${testLines}' ]"
  unset lines testLines
}

test_IndexStats()
{
  md1=`md5sum ${outDir}/${sampleName}/${sampleName}.idxStats.txt | cut -f 1 -d ' '`
  md2=`md5sum ${testDataDir}/${sampleName}.idxStats.txt | cut -f 1 -d ' '`
  assertTrue 'IndexStats don''t match' "[ '$md1' == '$md2' ]"
  unset md1 md2
}

test_BAMFileExists()
{
  assertTrue 'BAM File does not exist' "[ -e ${outDir}/${sampleName}/${sampleName}.bam ]" 
}

#testRun_HumanSample_hg19ERCC_ensembl()
#{
#  testDir="${SHUNIT_TMPDIR}/testRun_HumanSample_hg19ERCC_ensembl"
#  ${script} --refmodel hg19ERCC.ensembl -s test1 --fastq1 ${testDataDir}/test1_1.fastqgz --fastq2 ${testDataDir}/test1_2.fastq.gz --sample-dir ${testDir} > ${stdoutF} 2> ${stderrF}
#  rtrn=$?
  
  # check return value
#  assertTrue 'Script exited with non-zero exit code.' "[ ${rtrn} -eq 0 ]"
  
  # check output files
#}

#testRun_MouseSample_mm10ERCC_refSeq()
#{
  #${script} --refmodel mm10ERCC.refSeq -s 
#}

#testRun_MouseSample_mm10ERCC_ensembl()
#{
  #${script} --refmodel mm10ERCC.ensembl
#}

#testRun_RatSample_rn6ERCC_refSeq()
#{
  #${script} --refmodel rn6ERCC.refSeq
#}

#testRun_RatSample_rn6ERCC_ensembl()
#{
  #${script} --refmodel rn6ERCC.ensembl
#}

#-----------------------------------------------------------------------------
# suite helper functions
#

th_assertTrueWithNoOutput()
{
  th_return_=$1
  th_stdout_=$2
  th_stderr_=$3

  assertFalse 'unexpected output to STDOUT' "[ -s '${th_stdout_}' ]"
  assertFalse 'unexpected output to STDERR' "[ -s '${th_stderr_}' ]"

  unset th_return_ th_stdout_ th_stderr_
}

th_assertFileExists()
{
  th_testFile=$1
  
  assertTrue "${th_testFile} does not exist" "[ -e '${th_testFile}' ]"

  unset th_testFile
}

th_assertFileIsNonZeroSize()
{
  th_testFile=$1

  assertTrue '${th_testFile} does not exist' "[ -e '${th_testFile}' ]"
  assertTrue '${th_testfile} is zero size' "[ -s '${th_testFile}' ]" 

  unset th_testFile
}

# load and run shUnit2
[ -n "${ZSH_VERSION:-}" ] && SHUNIT_PARENT=$0
. /home/golharr/workspace/shunit2/source/2.1/src/shunit2
