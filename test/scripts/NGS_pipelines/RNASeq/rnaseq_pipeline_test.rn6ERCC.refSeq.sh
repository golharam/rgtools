#!/bin/sh
# vim:et:ft=sh:sts=2:sw=2
#
# Author: ryan.golhar@bms.com
#
# Unit test for rnaseq_pipeline.sh
#

oneTimeSetUp()
{
  echo "SHUNIT_TMPDIR: ${SHUNIT_TMPDIR}"
  echo

  script="/home/golharr/workspace/rgtools/scripts/NGS_pipelines/RNASeq/rnaseq_pipeline.sh"
  testDataDir="/home/golharr/workspace/rgtools/test_data/scripts/NGS_pipelines/RNASeq"

  outputDir="${SHUNIT_TMPDIR}/output"
  mkdir "${outputDir}"
  stdoutF="${outputDir}/stdout"
  stderrF="${outputDir}/stderr"

  testDir="${SHUNIT_TMPDIR}/testRun_RatSample_rn6ERCC_refSeq"
  outDir="${SHUNIT_TMPDIR}/analysis"
  mkdir ${outDir}

  # If you want to submit a previous analysis as a test,
  # set sampleName, copyResults=1, and resultsDir
  # else this script will re-run the analysis which can
  # take longer.
  sampleName="rat1.rn6ERCC.refSeq"
  copyResults=0
  if [ ${copyResults} -eq 1 ]; then
    resultsDir="/home/golharr/ngsprojects/test/analysis/${sampleName}"
    echo "Copying ${resultsDir} to ${outDir}/"
    cp -r ${resultsDir} ${outDir}/
    rtrn_sampleRun=$?
  else 
    echo "Running ${sampleName}:"
    echo "${script} -s ${sampleName} \
            --fastq1 ${testDataDir}/QCPipeline/SRR1169893.fastq.gz \
            --sample-dir ${testDir} --outdir ${outDir} \
            --refmodel rn6ERCC.refSeq \
            > ${stdoutF} 2> ${stderrF}"
    ${script} -s ${sampleName} \
            --fastq1 ${testDataDir}/QCPipeline/SRR1169893.fastq.gz \
            --sample-dir ${testDir} --outdir ${outDir} \
            --refmodel rn6ERCC.refSeq \
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

test_genes()
{
  th_assertFileExists ${outDir}/${sampleName}/${sampleName}.genes.results

  lineCount=`wc -l ${outDir}/${sampleName}/${sampleName}.genes.results | cut -f1 -d ' '`
  testLineCount=`wc -l ${testDataDir}/${sampleName}/${sampleName}.genes.results | cut -f1 -d ' '`
  assertTrue 'Line count incorrect' "[ $lineCount -eq $testLineCount ]"
  unset lineCount testLineCount

  header=`head -n 1 ${outDir}/${sampleName}/${sampleName}.genes.results`
  testHeader=`head -n 1 ${testDataDir}/${sampleName}/${sampleName}.genes.results`
  assertTrue 'Header lines do not match' "[ '${header}' == '${testHeader}' ]"
  unset header testHeader

  md5=`md5sum ${outDir}/${sampleName}/${sampleName}.genes.results | cut -f1 -d ' '`
  testmd5=`md5sum ${testDataDir}/${sampleName}/${sampleName}.genes.results | cut -f1 -d ' '`
  assertTrue 'File contents does not match' "[ '${md5}' == '${testmd5}' ]"
  unset md5 testmd5
}

test_isoforms()
{
  th_assertFileExists ${outDir}/${sampleName}/${sampleName}.isoforms.results

  lineCount=`wc -l ${outDir}/${sampleName}/${sampleName}.isoforms.results | cut -f1 -d ' '`
  testLineCount=`wc -l ${testDataDir}/${sampleName}/${sampleName}.isoforms.results | cut -f1 -d ' '`
  assertTrue 'Line count incorrect' "[ $lineCount -eq $testLineCount ]"
  unset lineCount testLineCount

  header=`head -n 1 ${outDir}/${sampleName}/${sampleName}.isoforms.results`
  testHeader=`head -n 1 ${testDataDir}/${sampleName}/${sampleName}.isoforms.results`
  assertTrue 'Header lines do not match' "[ '${header}' == '${testHeader}' ]"
  unset header testHeader

  md5=`md5sum ${outDir}/${sampleName}/${sampleName}.isoforms.results | cut -f1 -d ' '`
  testmd5=`md5sum ${testDataDir}/${sampleName}/${sampleName}.isoforms.results | cut -f1 -d ' '`
  assertTrue 'File contents does not match' "[ '${md5}' == '${testmd5}' ]"
  unset md5 testmd5
}

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
