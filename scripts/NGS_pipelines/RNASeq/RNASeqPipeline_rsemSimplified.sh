#!/bin/bash

# SGE Grid Engine Parameters
#$ -S /bin/bash
#$ -pe pe 2
#$ -j y
#$ -cwd

# To run:
# qsub -v PARAMETERS_FILE=./rnaseq_rsem_parameters.ini,SAMPLE_NAME=...,FASTQ_1=...,FASTQ_2=... ./RNASeqPipeline.sh

VERSION="v0.1"

##############################################################
##  Check for SAMPLE_NAME env var                           ##
##############################################################

if [ -z "$SAMPLE_NAME" ]
then
        echo "SAMPLE_NAME is not set.  Set as environment variable."
        exit 1
fi

echo "Running RNA-Seq Analysis Pipeline ($VERSION) on Sample ${SAMPLE_NAME} on `uname -n`"
echo "FASTQ_1: $FASTQ_1"
echo "FASTQ_2: $FASTQ_2"
date '+%m/%d/%y %H:%M:%S'
analysis_date_started=$(date +"%s")
echo

##############################################################
##  Check for parameters file 		                    ##
##############################################################

# We can define this as an environment variable.
if [ -z ${PARAMETERS_FILE} ]
then
#	Because each project may have different parameters, its better load 
#	a project or sample specific file from the current directory, if its not
#       specified.
	PARAMETERS_FILE=./rnaseq_rsem_parameters.ini
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
echo "THREADS=$THREADS"
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

#cd ${BASE_DIR}
cd /scratch2

if [ ! -d "${SAMPLE_NAME}" ] 
then
	mkdir ${SAMPLE_NAME}
fi

cd ${SAMPLE_NAME}

echo Using Analysis Directory: `pwd` 
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
## Print out version of tools used in this analysis         ##
##############################################################

echo The following applications will be used:
echo bowtie: `which bowtie`
bowtie --version
echo fastqValidator: ${FASTQ_VALIDATOR}
echo FastQC: ${FASTQC}
echo FastQ_Validator: ${FASTQ_VALIDATOR}
echo RSEM: ${RSEM_DIR}
echo RSEM Reference: ${RSEM_REFERENCE}
echo

##############################################################
##  Create/Read ini file (if it doesn't already exist)      ##
##############################################################

if [ -e ${SAMPLE_NAME}.ini ]
then
        echo "Reading previous sample settings"
        date '+%m/%d/%y %H:%M:%S'
        echo ---
        cat ${SAMPLE_NAME}.ini
        echo ---
	echo 

        source ${SAMPLE_NAME}.ini
fi

##############################################################
##  Check for FASTQ files                                   ##
##############################################################

# FASTQ_1 MUST be defined and exist.  If FASTQ_2 is defined, 
# make sure it exists.
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

if [ -n "${FASTQ_2}" ]
then
	# FASTQ_2 is defined, make sure it exists
	if [ ! -e ${FASTQ_2} ]
	then
	        echo Unable to locate FASTQ_2: $FASTQ_2 
	        echo in relation to analysis directory. Specify full path to FASTQ file.
	        date '+%m/%d/%y %H:%M:%S'
	        echo
	        exit -1
	fi
fi

##############################################################
## Count # of Reads (assuming fastq is compressed with gz)  ##
##############################################################

if [ -z "$COUNT_READ1" ]
then
	echo "Validating forward reads (${FASTQ_1})"
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

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
	date2=$(date +"%s")
	diff=$(($date2-$date1))
	echo "Validating forward reads took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

if [ -n "${FASTQ_2}" ] && [ -z "$COUNT_READ2" ]
then
        echo "Validating reverse reads (${FASTQ_2})"
        date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

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
	date2=$(date +"%s")
	diff=$(($date2-$date1))
        echo "Validating reverse reads took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
        echo
fi

##############################################################
## Map Reads						    ##
##############################################################

if [ ! -e ${SAMPLE_NAME}.genes.results ]
then
	echo "Running RSEM Calculate Expression"
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	if [ -n "${FASTQ_2}" ] 
	then
		$RSEM_DIR/rsem-calculate-expression -p $THREADS --time --paired-end ${FASTQ_1} ${FASTQ_2} $RSEM_REFERENCE $SAMPLE_NAME
	else
		$RSEM_DIR/rsem-calculate-expression -p $THREADS --time ${FASTQ_1} $RSEM_REFERENCE $SAMPLE_NAME
	fi

	if [ $? -ne 0 ]
	then
		echo "Error running RSEM"
		rm -f ${SAMPLE_NAME}.genes.results 2>/dev/null
		exit -1
	fi

	echo
	echo "Finished running RSEM"
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
	diff=$(($date2-$date1))
	echo "Running RSEM took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

##############################################################
## Clean Up                                                 ##
##############################################################

cd ..
mv ${SAMPLE_NAME} ${BASE_DIR}/

echo
echo Finished analyzing ${SAMPLE_NAME} 
date '+%m/%d/%y %H:%M:%S' 
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
echo
