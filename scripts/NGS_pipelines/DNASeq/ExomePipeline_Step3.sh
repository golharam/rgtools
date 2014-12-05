#!/bin/bash

# SGE Grid Engine Parameters
#$ -S /bin/bash
#$ -cwd
#$ -j y

# Example to run:
# FASTQ files are assumed to be in $BASEDIR as specified in exome_common_parameters.ini
# qsub -v PARAMETERS_FILE=`pwd`/exome_common_parameters.ini,SAMPLE_NAME=...,FASTQ_1=... `pwd`/ExomePipeline.sh

VERSION="v0.1"

echo "Running Exome Analysis Pipeline ($VERSION) on `uname -n`"
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
##  Check for SAMPLE_NAME env var                           ##
##############################################################

if [ -z "$SAMPLE_NAME" ]
then
	echo "SAMPLE_NAME is not set.  Set as environment variable."
	exit 1
fi

echo Analyzing Sample ${SAMPLE_NAME}
date '+%m/%d/%y %H:%M:%S' 
echo

##############################################################
##  Make analysis directory                                 ##
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

#if [ ! -d "${SAMPLE_NAME}" ] 
#then
#	mkdir ${SAMPLE_NAME}
#fi

#cd ${SAMPLE_NAME}

echo Using Analysis Directory: `pwd` 
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
## Print out version of tools used in this analysis         ##
##############################################################

echo The following applications will be used:
echo
echo Java: ${JAVA}
 ${JAVA} -version 2>&1 | head -n 1
echo
echo GATK: ${GATK_JAR}
echo GATK Version `${JAVA} -jar $GATK_JAR -version`
echo

if [ ! -e ${SAMPLE_NAME}.dedup.realigned.recal.bam ]
then
	echo Recalibrating quality scores
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -jar ${GATK_JAR} \
		-T PrintReads \
		-R ${REFERENCE} \
		-I ${SAMPLE_NAME}.dedup.realigned.bam \
		-BQSR recal.table \
		-o ${SAMPLE_NAME}.dedup.realigned.recal.bam

	if [ $? -ne 0 ]
	then
		echo Error running PrintReads
		rm ${SAMPLE_NAME}.dedup.realigned.recal.ba? 2>/dev/null
		exit -1
	fi

	echo Finished recalibrating quality scores
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

echo
echo Finished analyzing ${SAMPLE_NAME} 
date '+%m/%d/%y %H:%M:%S' 
echo
