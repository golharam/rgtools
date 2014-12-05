#!/bin/bash

# SGE Grid Engine Parameters
#$ -S /bin/sh
#$ -pe pe 4
#$ -j y
#$ -cwd

# Example to run:
# qsub -v PARAMETERS_FILE=`pwd`/exome_common_parameters.ini `pwd`/ExomePipeline_Step2.sh

VERSION="v0.1"

echo "Running Exome Analysis Pipeline Step 2 ($VERSION) on `uname -n`"
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
##  Exome Analysis Pipeline 		                    ##
##############################################################

# We can define this as an environment variable.
if [ -z ${PARAMETERS_FILE} ]
then
#	Because each project may have different parameters, its better load 
#	a project specific file from the current directory, if its not
#       specified.
	PARAMETERS_FILE=./exome_common_parameters.ini
fi

if [ ! -e ${PARAMETERS_FILE} ] 
then
	echo "Unable to locate ${PARAMETERS_FILE}"
	exit -1
fi

echo "Using parameters file $PARAMETERS_FILE"
date '+%m/%d/%y %H:%M:%S' 
echo
. ${PARAMETERS_FILE}

##############################################################
##  Go to analysis directory                                ##
##############################################################

if [ -z "$BASE_DIR" ]
then
	echo "BASE_DIR is not set"
	exit 1
fi

if [ ! -d "${BASE_DIR}" ] 
then
	echo "$BASE_DIR does not exist"
	exit -1
fi

cd ${BASE_DIR}

echo Using Analysis Directory: `pwd` 
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
## Print out version of tools used in this analysis         ##
##############################################################

echo The following applications will be used:
echo fastqValidator: ${FASTQ_VALIDATOR}
echo
echo FastQC: ${FASTQC}
echo
echo BWA: ${BWA}
 ${BWA} 2>&1 | head -n 3 | tail -n 1
echo BWA Reference: ${BWA_REFERENCE}
echo
echo Java: ${JAVA}
 ${JAVA} -version 2>&1 | head -n 1
echo
echo PicardTools: ${PICARD_TOOLS}
echo
echo GATK: ${GATK_JAR}
echo GATK Version `${JAVA} -jar $GATK_JAR -version`
echo
echo SAMTools: ${SAM_TOOLS}
echo SAMTools `${SAM_TOOLS} 2>&1 | head -n 3 | tail -n 1`
echo

# Alignment summary metrics cannot be specified on the whole lot of BAM files at once.
# So instead, we have to collect them all at the end :/
if [ ! -e alignment_summary_metrics.txt ]
then
        echo Merging alignment summary metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

	for f in `ls *.alignment_summary_metrics.txt`
	do
		cat $f >> alignment_summary_metrics.txt
	done

        if [ $? -ne 0 ]
        then
                echo "Encountered error merging alignment summary metrics"
                rm alignment_summary_metrics.txt 2>/dev/null
                exit -1
        fi

        echo "Finished merging alignment summary metrics"
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

if [ ! -e dedup_metrics.txt ]
then
        echo Merging dedup metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

        for f in `ls *.dedup_metrics.txt`
        do
                head -8 $f >> dedup_metrics.txt
        done

        if [ $? -ne 0 ]
        then
                echo "Encountered error merging dedup metrics"
                rm dedup_metrics.txt 2>/dev/null
                exit -1
        fi

        echo "Finished merging dedup metrics"
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

if [ ! -e indelrealigner.list ]
then
	echo Creating indel realigner targets
	date '+%m/%d/%y %H:%M:%S'
	echo

	# This next line is project dependent.  Some samples may want to be 
	# excluded, such as no-template controls, etc.
	#find . -name '*.dedup.bam' > bamfiles.list
	${JAVA} -jar ${GATK_JAR} \
		-T RealignerTargetCreator \
		-nt ${THREADS} \
		-R ${REFERENCE} \
		-I bamfiles.list \
		-o indelrealigner.list
	
	if [ $? -ne 0 ]
	then
		echo "Encountered error running indel realigner target creator."
		rm indelrealigner.list 2>/dev/null
		exit -1
	fi

	echo Finished creating indel realigner targets
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e indelrealigner.done ]
then
	echo Realigning around indels
	date '+%m/%d/%y %H:%M:%S'
	echo

	# 1:  bamfiles.list contains the full path to the bam files
	# 2:  bamfiles.realigned.amp contains two columns:
	#	first col is just the filename of the bam file
	#	second col is the full path to the output bam file
	# This line only works in conjunction with find . -name '*.dedup.bam' > bamfiles.list above
	paste bamfiles.list bamfiles.list > bamfiles.realigned.map
	${JAVA} -jar ${GATK_JAR} \
		-T IndelRealigner \
		-R ${REFERENCE} \
		-I bamfiles.list \
		-targetIntervals indelrealigner.list \
		-nWayOut bamfiles.realigned.map

	if [ $? -ne 0 ]
	then
		echo "Encountered error running indel realigner."
		rm indelrealigner.done 2>/dev/null
		exit -1
	fi

	touch indelrealigner.done
	echo Finished realigning around indels
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e recal.table ]
then
	echo Counting covariates using BaseCalibrator
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -jar ${GATK_JAR} \
		-T BaseRecalibrator \
		-I bamfiles.realigned.list \
		-R ${REFERENCE} \
		-o recal.table \
		-knownSites $GATK_KNOWN_SNPS

	if [ $? -ne 0 ]
	then
		echo Error counting covariates
		rm recal.table 2>/dev/null
		exit -1
	fi

	echo Finished counting covariates
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e post_recal.table ]
then
	echo 2nd pass to look for covariation using BaseCalibrator
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -jar ${GATK_JAR} \
		-T BaseRecalibrator \
		-I bamfiles.realigned.list \
		 -knownSites $GATK_KNOWN_SNPS \
		-R ${REFERENCE} \
		-BQSR recal.table \
		-o post_recal.table \

	if [ $? -ne 0 ]
	then
		echo Error counting covariates 2nd pass
		post_recal.table 2>/dev/null
		exit -1
	fi

	echo Finished 2nd pass to look for covariation using BaseCalibrator
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e recalibration_plots.pdf ]
then
	echo Creating recalibration plots
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -jar ${GATK_JAR} \
		-T AnalyzeCovariates \
		-R ${REFERENCE} \
		-before recal.table \
		-after post_recal.table \
		-plots recalibration_plots.pdf

	if [ $? -ne 0 ]
	then
		echo Error creating plots
		recalibration_plots.pdf 2>/dev/null
		exit -1
	fi

	echo Finished creating recalibration plots
	date '+%m/%d/%y %H:%M:%S'
	echo	
fi

echo
echo Finished step 2 analysis
date '+%m/%d/%y %H:%M:%S' 
echo
