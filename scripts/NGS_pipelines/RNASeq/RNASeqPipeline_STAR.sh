#!/bin/bash

# SGE Grid Engine Parameters
#$ -pe orte 4
#$ -j y

# PBS Parameters
#PBS -l nodes=1:ppn=4
#PBS -j oe

# Example to run:
# FASTQ files are assumed to be in FASTQ_1 and FASTQ_2 as env var or specified in rnaseq_common_parameters.ini
# qsub -v PARAMETERS_FILE=`pwd`/rnaseq_common_parameters.ini,SAMPLE_NAME=...,FASTQ_1=...,FASTQ_2=... `pwd`/RNASeqPipeline.sh

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
	PARAMETERS_FILE=./rnaseq_common_parameters.ini
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

cd ${BASE_DIR}

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
echo fastqValidator: ${FASTQ_VALIDATOR}
echo FastQC: ${FASTQC}
echo
echo Tophat: `which tophat`
tophat -v
echo Tophat Reference: ${TOPHAT_REFERENCE}
echo
echo Java: ${JAVA}
 ${JAVA} -version 2>&1 | head -n 1
echo
echo PicardTools: ${PICARD_TOOLS}
echo
echo GATK: ${GATK_JAR}
#echo GATK Version `${JAVA} -jar $GATK_JAR -version`
echo
echo SAMTools: ${SAM_TOOLS}
echo SAMTools `${SAM_TOOLS} 2>&1 | head -n 3 | tail -n 1`
echo
echo STAR: ${STAR}
echo STAR `$STAR --version`
echo

#############################################################
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
## Run FastQC                                               ##
##############################################################

if [ ! -d fastqc ]
then
	echo "Running FastQC "
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	mkdir fastqc
	$FASTQC --outdir=fastqc $FASTQ_1 $FASTQ_2

        if [ $? -ne 0 ]
        then
                echo "Error running FastQC"
		rm -rf fastqc 2>/dev/null
                exit -1
        fi

	echo
	echo "Finished running FastQC"
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
	diff=$(($date2-$date1))
        echo "Runng FastQC took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

##############################################################
## Map Reads						    ##
##############################################################

#if [ ! -d tophat_out ]
#then
#	echo Running tophat2
#	date '+%m/%d/%y %H:%M:%S'
#	date1=$(date +"%s")

	# TBD: Read group information isn't taking.  Why?
#	tophat2 -p $THREADS --rg-id $SAMPLE_NAME --rg-sample $SAMPLE_NAME --rg-library $SAMPLE_NAME --rg-platform illumina --rg-platform-unit MiSeq320 $TOPHAT_REFERENCE $FASTQ_1 $FASTQ_2 

#	if [ $? -ne 0 ]
#	then
#		echo "Encountered error running tophat2"
#		rm -rf tophat_out 2>/dev/null
#		exit -1
#	fi

#	ln -s tophat_out/accepted_hits.bam ./${SAMPLE_NAME}.bam

#	echo Finished running tophat2
#	date '+%m/%d/%y %H:%M:%S'
#	date2=$(date +"%s")
#	diff=$(($date2-$date1))
#       echo "Runng tophat2 took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
#	echo
#fi

if [ ! -e Aligned.out.sam ]
then
	echo Running STAR
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	$STAR --runThreadN $THREADS --runMode alignReads --genomeDir $REF_STAR \
		--readFilesIn $FASTQ_1 $FASTQ_2 --readFilesCommand gzip -dc

	if [ $? -ne 0 ]
        then
        	echo "Encountered error running STAR"
                rm -f Aligned.out.sam 2>/dev/null
                exit -1
        fi

	ln -s Aligned.out.sam $SAMPLE_NAME.sam

	echo Finished running STAR
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
	diff=$(($date2-$date1))
	echo Running STAR took $(($diff / 60)) minutes and $(($diff % 60)) seconds.
	echo
fi

if [ ! -e ${SAMPLE_NAME}.bam ] 
then
	echo "Converting SAM to BAM"
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	#${JAVA} -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/SortSam.jar \
	#	INPUT=$SAMPLE_NAME.sam OUTPUT=${SAMPLE_NAME}.bam SORT_ORDER=coordinate \
	#	VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

	$JAVA -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/AddOrReplaceReadGroups.jar \
		INPUT=$SAMPLE_NAME.sam OUTPUT=${SAMPLE_NAME}.bam \
		RGID=$SAMPLE_NAME RGSM=$SAMPLE_NAME RGLB=$SAMPLE_NAME \
		RGPL=illumina RGPU=MiSeq320 \
		SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true

        if [ $? -ne 0 ]
        then
                echo "Encountered error converting SAM to BAM"
                rm ${SAMPLE_NAME}.ba? 2>/dev/null
                exit -1
        fi

        echo Finished converting SAM to BAM
        date '+%m/%d/%y %H:%M:%S'
        date2=$(date +"%s")
        diff=$(($date2-$date1))
        echo "Runng SAM->BAM took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
        echo
fi

if [ ! -e ${SAMPLE_NAME}.dedup.bam ]
then
	echo Mark duplicates
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	${JAVA} -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/MarkDuplicates.jar \
		INPUT=$SAMPLE_NAME.bam \
		OUTPUT=${SAMPLE_NAME}.dedup.bam \
		METRICS_FILE=${SAMPLE_NAME}.dedup_metrics.txt \
		VALIDATION_STRINGENCY=LENIENT \
		CREATE_INDEX=true

	if [ $? -ne 0 ]
	then
		echo "Encountered error marking duplicates"
		rm ${SAMPLE_NAME}.dedup.ba? 2>/dev/null
		exit -1
	fi

	echo Finished marking duplicates
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
        diff=$(($date2-$date1))
        echo "Runng MarkDuplicates took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

if [ ! -e ${SAMPLE_NAME}.dedup.karyotype.bam ]
then
	ln -s ${SAMPLE_NAME}.dedup.bam ${SAMPLE_NAME}.dedup.karyotype.bam
	ln -s ${SAMPLE_NAME}.dedup.bai ${SAMPLE_NAME}.dedup.karyotype.bai
#	echo Reordering BAM file
#        date '+%m/%d/%y %H:%M:%S'
#	date1=$(date +"%s")

#	${JAVA} -Djava.io.tmpdir=$TMPDIR -jar ${PICARD_TOOLS}/ReorderSam.jar \
#		I=${SAMPLE_NAME}.dedup.bam \
#		O=${SAMPLE_NAME}.dedup.karyotype.bam \
#		R=$REFERENCE \
#		CREATE_INDEX=true

#        if [ $? -ne 0 ]
#        then
#                echo "Encountered error reordering BAM file"
#                rm ${SAMPLE_NAME}.dedup.karyotype.* 2>/dev/null
#                exit -1
#        fi

#        echo Finished reordering BAM file
#        date '+%m/%d/%y %H:%M:%S'
#        date2=$(date +"%s")
#        diff=$(($date2-$date1))
#        echo "Runng ReorderSam took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
#        echo
fi

if [ ! -e ${SAMPLE_NAME}.dedup.karyotype.trimmed.bam ]
then
	echo Running GATK SplitNTrim
        date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")

	$JAVA -jar $GATK_JAR \
		-T SplitNCigarReads \
		-I ${SAMPLE_NAME}.dedup.karyotype.bam \
		-o ${SAMPLE_NAME}.dedup.karyotype.trimmed.bam \
		-R $REFERENCE \
		-U ALLOW_N_CIGAR_READS

        if [ $? -ne 0 ]
        then
                echo "Encountered error running SplitNTrim"
                rm ${SAMPLE_NAME}.dedup.karyotype.trimmed.* 2>/dev/null
                exit -1
        fi

        echo Finished running SplitNTrim
        date '+%m/%d/%y %H:%M:%S'
        date2=$(date +"%s")
        diff=$(($date2-$date1))
        echo "Runng GATK::SplitNCigarReads took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
        echo
fi

##############################################################
## Count Reads                                              ##
##############################################################

if [ ! -e $SAMPLE_NAME.rawcounts.txt ]
then
        echo Count Reads
        date '+%m/%d/%y %H:%M:%S'
        date1=$(date +"%s")

	# See http://www-huber.embl.de/users/anders/HTSeq/doc/count.html for parameters
	$HTSEQ_COUNT -f bam -r pos --stranded=$HTSEQ_STRANDED ${SAMPLE_NAME}.dedup.karyotype.trimmed.bam $REF_GFF > $SAMPLE_NAME.raw_counts.txt

        if [ $? -ne 0 ]
        then
                echo "Encountered error counting reads"
                rm ${SAMPLE_NAME}.raw_counts.txt 2>/dev/null
                exit -1
        fi

        echo Finished counting reads
        date '+%m/%d/%y %H:%M:%S'
        date2=$(date +"%s")
        diff=$(($date2-$date1))
        echo "Running htseq-count took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
        echo
fi

echo
echo Finished analyzing ${SAMPLE_NAME} 
date '+%m/%d/%y %H:%M:%S' 
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
echo
