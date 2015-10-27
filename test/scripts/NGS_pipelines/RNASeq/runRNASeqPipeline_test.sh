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

  script="/home/golharr/workspace/rgtools/scripts/NGS_pipelines/RNASeq/runRNASeqPipeline.pl"
  testDataDir="/home/golharr/workspace/rgtools/test_data/scripts/NGS_pipelines/RNASeq"

  outputDir="${SHUNIT_TMPDIR}/output"
  mkdir "${outputDir}"
  stdoutF="${outputDir}/stdout"
  stderrF="${outputDir}/stderr"

  # Sample specific stuff.  change as needed: 
#  sampleName="rat1.rn6ERCC.refSeq"
  # Where the test analysis should run:
#  testDir="${SHUNIT_TMPDIR}/testRun_${sampleName}"
  # Where the test analysis should write final test data
#  outDir="${SHUNIT_TMPDIR}/analysis"
#  mkdir ${outDir}

  # If you want to submit a previous analysis as a test,
  # set sampleName, copyResults=1, and resultsDir
  # else this script will re-run the analysis which can
  # take longer.
#  copyResults=1

#  if [ ${copyResults} -eq 1 ]; then
#    resultsDir="/home/golharr/ngsprojects/test/analysis/${sampleName}"
#    echo "Copying ${resultsDir} to ${outDir}/"
#    cp -r ${resultsDir} ${outDir}/
#    rtrn_sampleRun=$?
#  else 
#    echo "Running ${sampleName}:"
#    echo "${script} -s ${sampleName} \
#            --fastq1 ${testDataDir}/SRR1169893.fastq.gz \
#            --sample-dir ${testDir} --outdir ${outDir} \
#            --refmodel rn6ERCC \
#            > ${stdoutF} 2> ${stderrF}"
#    ${script} -s ${sampleName} \
#            --fastq1 ${testDataDir}/SRR1169893.fastq.gz \
#            --sample-dir ${testDir} --outdir ${outDir} \
#            --refmodel rn6ERCC \
#            > ${stdoutF} 2> ${stderrF}
#    rtrn_sampleRun=$?
#  fi
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
  assertFalse 'unexpected output to STDOUT' "[ -s '${stdoutF}' ]"

  # check stderr
  md1=`md5sum ${stderrF} | cut -f 1 -d ' '`
  md2=`md5sum ${testDataDir}/runRNASeqPipeline.helpscreen.err | cut -f 1 -d ' '`
  assertTrue 'Help screen don''t match' "[ '$md1' == '$md2' ]"
  unset md1 md2

  # check output files
  ## no files created
}

testDryRunWith2Samples()
{
  ${script} --dryRun --outdir /scratch/analysis --refmodel hg19ERCC.refSeq --tmpdir /scratch --samples ${testDataDir}/QCPipeline/samples.txt > ${stdoutF} 2> ${stderrF}
  rtrn=$?

  # check return value
  assertTrue 'Return value non-zero' "[ ${rtrn} -eq 0 ]"

  # check stdout
  md1=`md5sum ${stdoutF} | cut -f 1 -d ' '`
  md2=`md5sum ${testDataDir}/runRNASeqPipeline.testDryRunWith2Samples.stdout | cut -f 1 -d ' '`
  assertTrue 'Unexpected output to STDOUT' "[ '$md1' == '$md2' ]"
  unset md1 md2

  # check stderr
  md1=`md5sum ${stderrF} | cut -f 1 -d ' '`
  md2=`md5sum ${testDataDir}/runRNASeqPipeline.testDryRunWith2Samples.stderr | cut -f 1 -d ' '`
  assertTrue 'Unexpected output to STDERR' "[ '$md1' == '$md2' ]"
  unset md1 md2
}

testRunWith2Samples()
{
  oldWD=`pwd`
  cd ${SHUNIT_TMPDIR}
 
  analysisDir="${SHUNIT_TMPDIR}/analysis"
  mkdir ${analysisDir}

  # We can't actually run the samples here becasuse the test directory
  # is not shared.  So, let's assume that tests execute correctly and
  # we are just testing the matrix building portion. The test to run
  # is from testDryRunWith2Samples() without --dryRun
  tests=( test1 test2 )
  for test in ${tests[@]}
  do
   # Comparison results are in ${testDataDir}/hg19ERCC.refSeq.testresults
    resultsDir="${testDataDir}/twoSampleAnalysisOutput/${test}"
    echo "cp -r ${resultsDir} ${analysisDir}/"
    cp -r ${resultsDir} ${analysisDir}/
  done
  unset test tests

  ${script} --outdir ${analysisDir} --refmodel hg19ERCC.refSeq --tmpdir /scratch --samples ${testDataDir}/QCPipeline/samples.txt > ${stdoutF} 2> ${stderrF}
  rtrn=$?

  # check return value
  assertTrue 'Return value non-zero' "[ ${rtrn} -eq 0 ]"

  # check for matrix files 
  genesIsoforms=( genes isoforms )
  matrixType=( counts FPKM TPM UQ )
  for gi in ${genesIsoforms[@]}; do
    for matrix in ${matrixType[@]}; do
      md1=`sort ${analysisDir}/${gi}.${matrix}.matrix | md5sum | cut -f 1 -d ' '`
      md2=`sort $testDataDir/twoSampleAnalysisOutput/${gi}.${matrix}.matrix | md5sum | cut -f 1 -d ' '`
      assertTrue "${gi}.${matrix}.matrix do not match" "[ '$md1' == '$md2' ]"
    done
  done

  cd $oldWD
  unset oldWD
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
