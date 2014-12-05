#!/bin/bash

# SGE Grid Engine Parameters
#$ -S /bin/sh
#$ -pe pe 4
#$ -j y
#$ -cwd

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

##############################################################
##  Create/Read ini file (if it doesn't already exist)      ##
##############################################################

if [ -e ${SAMPLE_NAME}.ini ]
then
        echo "Reading previous sample settings: "
	echo `pwd`/${SAMPLE_NAME}.ini
        date '+%m/%d/%y %H:%M:%S'
        echo ---
        cat ${SAMPLE_NAME}.ini
        echo ---
	echo 

        source ./${SAMPLE_NAME}.ini
fi

##############################################################
##  Check for FASTQ files and read group information        ##
##############################################################

if [ -z "${FASTQ_1}" ]
then
	echo "FASTQ_1 not defined."
	date '+%m/%d/%y %H:%M:%S'
	exit -1
fi

if [ ! -e ${FASTQ_1} ] 
then
	echo Unable to locate FASTQ_1: $FASTQ_1 
	echo in relation to analysis directory. Specify full path to FASTQ file.
	date '+%m/%d/%y %H:%M:%S'
	echo
	exit -1	
fi

##############################################################
## Count # of Reads (assuming fastq is compressed with gz)  ##
##############################################################

if [ -z "$COUNT_READ1" ]
then
	echo "Validating forward reads (${FASTQ_1})"
	date '+%m/%d/%y %H:%M:%S'
	
	lines=`$FASTQ_VALIDATOR --file ${FASTQ_1}`
	if [ $? -ne 0 ]
	then
		echo $lines
		echo "$FASTQ_1 is not valid"
		exit -1
	fi
	
	COUNT_READ1=`echo "$lines" | cut -d ' ' -f 8`
	
	echo FASTQ_1=${FASTQ_1} >> ${SAMPLE_NAME}.ini
	echo COUNT_READ1=${COUNT_READ1} >> ${SAMPLE_NAME}.ini

	echo The forward or read1 fastq contains ${COUNT_READ1} reads
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ -z "$COUNT_READ2" ]
then
        echo "Validating reverse reads (${FASTQ_2})"
        date '+%m/%d/%y %H:%M:%S'

        lines=`$FASTQ_VALIDATOR --file ${FASTQ_2}`
        if [ $? -ne 0 ]
        then
                echo $lines
                echo "$FASTQ_2 is not valid"
                exit -1
        fi

        COUNT_READ2=`echo "$lines" | cut -d ' ' -f 8`

        echo FASTQ_2=${FASTQ_2} >> ${SAMPLE_NAME}.ini
        echo COUNT_READ2=${COUNT_READ2} >> ${SAMPLE_NAME}.ini

        echo The reverse or read2 fastq contains ${COUNT_READ2} reads
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

##############################################################
## Run FastQC                                               ##
##############################################################

FASTQC_OUTFILE=`basename $FASTQ_1 $EXT`
FASTQC_OUTFILE="${FASTQC_OUTFILE}_fastqc.html"
if [ ! -e $FASTQC_OUTFILE ]
then
	echo "Running FastQC"
	date '+%m/%d/%y %H:%M:%S'

	# TODO: Add --outfile
	$FASTQC -t ${THREADS} --outdir=. $FASTQ_1 $FASTQ_2

	echo
	echo "Finished running FastQC"
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

##############################################################
## Map Reads						    ##
##############################################################

if [ ! -e ${SAMPLE_NAME}_1.sai ]
then
	echo Aligning forward reads 
	date '+%m/%d/%y %H:%M:%S'
	echo

	echo ${BWA} aln -t ${THREADS} -f ${SAMPLE_NAME}_1.sai ${BWA_REFERENCE} ${FASTQ_1}
	${BWA} aln -t ${THREADS} -f ${SAMPLE_NAME}_1.sai ${BWA_REFERENCE} ${FASTQ_1}

	if [ $? -ne 0 ]
	then
		echo "Encountered error aligning forward reads"
		rm ${SAMPLE_NAME}_1.sai
		exit -1
	fi

	echo Finished aligning forward reads
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e ${SAMPLE_NAME}_2.sai ]
then
        echo Aligning reverse reads
        date '+%m/%d/%y %H:%M:%S'
        echo

        echo ${BWA} aln -t ${THREADS} -f ${SAMPLE_NAME}_2.sai ${BWA_REFERENCE} ${FASTQ_2}
        ${BWA} aln -t ${THREADS} -f ${SAMPLE_NAME}_2.sai ${BWA_REFERENCE} ${FASTQ_2}

        if [ $? -ne 0 ]
        then
                echo "Encountered error aligning reverse reads"
                rm ${SAMPLE_NAME}_2.sai
                exit -1
        fi

        echo Finished aligning reverse reads
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

if [ ! -e ${SAMPLE_NAME}.sam ]
then
	echo Converting sai to sam
	date '+%m/%d/%y %H:%M:%S'
	echo

	${BWA} sampe -f ${SAMPLE_NAME}.sam -r "@RG\tID:${SAMPLE_NAME}\tLB:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:illumina" ${BWA_REFERENCE} ${SAMPLE_NAME}_1.sai ${SAMPLE_NAME}_2.sai ${FASTQ_1} ${FASTQ_2}

	if [ $? -ne 0 ]
	then
		echo "Encountered error converting sai to sam"
		rm ${SAMPLE_NAME}.sam
		exit -1
	fi

	echo Finished converting to sam
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e ${SAMPLE_NAME}.bam ]
then
	echo Converting sam to bam
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -Xmx4G -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/SortSam.jar \
		SO=coordinate \
		INPUT=${SAMPLE_NAME}.sam \
		OUTPUT=${SAMPLE_NAME}.bam \
		VALIDATION_STRINGENCY=LENIENT \
		CREATE_INDEX=true

	if [ $? -ne 0 ]
	then
		echo "Encountered error converting sam to BAM"
		rm ${SAMPLE_NAME}.ba? 2>/dev/null
		exit -1
	fi

	echo Finished converting sam to bam
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

if [ ! -e ${SAMPLE_NAME}.alignment_summary_metrics.txt ]
then
	echo Collection alignment summary metrics
	date '+%m/%d/%y %H:%M:%S'
        echo

	$JAVA -jar $PICARD_TOOLS/CollectAlignmentSummaryMetrics.jar \
		I=${SAMPLE_NAME}.bam \
		O=${SAMPLE_NAME}.alignment_summary_metrics.txt \
		LEVEL=SAMPLE \
		R=$REFERENCE

        if [ $? -ne 0 ]
        then
                echo "Encountered error collecting alignment summary metrics"
                rm ${SAMPLE_NAME}.alignment_summary_metrics.txt 2>/dev/null
                exit -1
        fi

        echo "Finished collecting alignment summary metrics"
        date '+%m/%d/%y %H:%M:%S'
        echo

fi

if [ ! -e ${SAMPLE_NAME}.dedup.bam ]
then
	echo Removing duplicates
	date '+%m/%d/%y %H:%M:%S'
	echo

	${JAVA} -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/MarkDuplicates.jar \
		INPUT=${SAMPLE_NAME}.bam \
		OUTPUT=${SAMPLE_NAME}.dedup.bam \
		METRICS_FILE=${SAMPLE_NAME}.dedup_metrics.txt \
		VALIDATION_STRINGENCY=LENIENT \
		CREATE_INDEX=true

	if [ $? -ne 0 ]
	then
		echo "Encountered removing duplicates"
		rm ${SAMPLE_NAME}.dedup.ba? 2>/dev/null
		exit -1
	fi

	echo Finished removing duplicates
	date '+%m/%d/%y %H:%M:%S'
	echo
fi

echo
echo Finished analyzing ${SAMPLE_NAME} 
date '+%m/%d/%y %H:%M:%S' 
echo
