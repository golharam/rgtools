#!/bin/bash

#$ -S /bin/sh
#$ -pe orte 8
#$ -j y
#$ -cwd

##############################################################################
# Set default values.  These can be set via command-line options only.  
# Setting them here overrides environment variables
###############################################################################
VERSION="0.1"
HELP=0
THREADS=8

if [ $# -eq 0 ]
then
	if [ -z "$SAMPLE" ]
	then
		HELP=1
	fi
fi

###############################################################################
# Get command line options.
# a.  Use > 1 to consume two arguments per pass in the loop (e.g. each argument
#     has a corresponding value to go with it).
# b.  Use > 0 to consume one or more arguments per pass in the loop (e.g. some
#     arguments don't have a corresponding value to go with it such as --help).
###############################################################################
while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		-1|--fastq1)
		FASTQ1="$2"
		shift # past argument
		;;

		-2|--fastq2)
		FASTQ2="$2"
		shift # past argument
		;;

		-a|--aws)
		AWS=1
		;;

		-h|--help)
		HELP=1
		;;

		-s|--sample)
		SAMPLE="$2"
		shift # past argument
		;;

		--tmpdir)
		TMP_DIR="$2"
		shift # past argument
		;;

		*)
		echo "Unknown option: $key"
		exit -1
		;;
	esac
	shift	# past argument or value
done

if [ $HELP == 1 ]; then
	echo "Usage 1: $0 <options> <usage>"
	echo "where usage:"
	echo "Usage 1: <options> [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>]"
	echo "Options:"
	echo "	-a|--aws"
	echo "	-h|--help"
	echo "  --tmpdir <scratch space>"
	exit 0
fi

# If variables are not set as an environment variable or as a parameter, set them to default here.
if [ -z "$AWS" ]; then
	AWS=0
fi
if [ -z "$TMP_DIR" ]; then
	TMP_DIR=/scratch
fi

if [ $AWS -eq 1 ]; then
	EXT_PKGS_DIR=/ngs/apps
	REFERENCE_DIR=/ngs/reference
else
	EXT_PKGS_DIR=/apps/sys/galaxy/external_packages
	REFERENCE_DIR=/ng18/galaxy/reference_genomes	
fi

# Applications / Programs
BWA=$EXT_PKGS_DIR/bwa-0.7.12/bwa
FASTQC=$EXT_PKGS_DIR/FastQC-0.11.3-RG/fastqc
FASTQVALIDATOR=$EXT_PKGS_DIR/fastQValidator-0.1.1a/bin/fastQValidator

# Reference data
SPECIES=hg19
REFERENCE="$REFERENCE_DIR/$SPECIES/$SPECIES.fa"
BWA_REFERENCE="$REFERENCE_DIR/$SPECIES/bwa_index/$SPECIES"

# Start Analysis
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "Pipeline: $VERSION"
echo "TmpDir: $TMP_DIR"
echo "AWS: $AWS"
echo

# Make tmp working directory aka SAMPLE_DIR
SAMPLE_DIR=/scratch/$SAMPLE
#if [ -z "$SAMPLE_DIR" ]
#then
#	cd $TMP_DIR
#	SAMPLE_DIR=`mktemp -d --tmpdir=${TMP_DIR} ${SAMPLE}_XXXXXX`
	if [ ! -d $SAMPLE_DIR ]
	then
	        mkdir $SAMPLE_DIR
	fi
#fi
cd $SAMPLE_DIR
echo "Working in `uname -n`:`pwd`"
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################################
# Step 1: Download the file to the local scratch space
# If AWS, copy from S3, FASTQ1, FASTQ2 
# Else assume file is local and unknown if compressed or not.  
# Regardless, $FASTQ1 and $FASTQ2 point to the files.
##############################################################################

if [ $AWS -eq 1 ] 
then
	if [ ! -e ${SAMPLE}_1.fastq.gz ]
	then
		echo "Downloading from S3: $FASTQ1"
		date '+%m/%d/%y %H:%M:%S'
		echo

		aws s3 cp $FASTQ1 ${SAMPLE}_1.fastq.gz
	
		if [ $? -ne 0 ]
		then
		        echo "Error downloading $FASTQ1"
		        rm -f ${SAMPLE}_1.fastq.gz 2>/dev/null
		        exit -1
		fi
	fi
	FASTQ1=`pwd`/${SAMPLE}_1.fastq.gz
	
	if [ ! -e ${SAMPLE}_2.fastq.gz ]
	then
		echo "Downloading from S3: $FASTQ2"
		date '+%m/%d/%y %H:%M:%S'
		echo
	
		aws s3 cp $FASTQ2 ${SAMPLE}_2.fastq.gz
	
		if [ $? -ne 0 ]
		then
		        echo "Error downloading $FASTQ2"
		        rm -f ${SAMPLE}_2.fastq.gz 2>/dev/null
		        exit -1
		fi
	fi
	FASTQ2=`pwd`/${SAMPLE}_2.fastq.gz
fi

# If we didn't do any of the above, then FASTQ1 and FASTQ2 are local files.
# They are either links to the original files or are the original files.
# In either case, make a link to them so we don't accidentially delete them.


echo "Using Fastq file(s):"
if [ ! -e $FASTQ1 ]
then
	echo "Unable to locate $FASTQ1.  Make sure full path is provided"
	exit -1
fi
echo FASTQ1: $FASTQ1
if [ -n "$FASTQ2" ]
then
	if [ ! -e $FASTQ2 ]
	then
		echo "Unable to locate $FASTQ2.  Make sure full path is provided"
		exit -1
	fi
	echo FASTQ2: $FASTQ2
fi

date '+%m/%d/%y %H:%M:%S'
analysis_date_started=$(date +"%s")
echo

##############################################################################
# Step 2: Run FastQ Validator (make sure the fastq [zip] file is valid)
##############################################################################
if [ ! -e $SAMPLE.fq1Validator.txt ]
then
	echo "Validating forward reads ($FASTQ1)"
	date '+%m/%d/%y %H:%M:%S'
	echo

	$FASTQVALIDATOR --disableSeqIDCheck --file $FASTQ1 > $SAMPLE.fq1Validator.txt

	if [ $? -ne 0 ]
	then
		echo "$FASTQ1 is not valid"
		rm $SAMPLE.fq1Validator.txt
		exit -1
	fi
fi

if [ -n "$FASTQ2" ] && [ ! -e $SAMPLE.fq2Validator.txt ]
then
        echo "Validating reverse reads ($FASTQ2)"
        date '+%m/%d/%y %H:%M:%S'
        echo

	$FASTQVALIDATOR --disableSeqIDCheck --file $FASTQ2 > $SAMPLE.fq2Validator.txt

        if [ $? -ne 0 ]
        then
                echo "$FASTQ2 is not valid"
		rm $SAMPLE.fq2Validator.txt
                exit -1
        fi

fi

##############################################################################
# Step 2: Run FastQC
##############################################################################
if [ ! -e "${SAMPLE}_1_fastqc.zip" ]
then
	echo "Running FastQC"
	date '+%m/%d/%y %H:%M:%S'
	echo

	$FASTQC --quiet --outdir=. --extract -t 2 $FASTQ1 $FASTQ2

	if [ $? -ne 0 ]
	then
		echo "Error running FastQC"
		rm -rf *_fastqc *_fastqc.zip *_fastqc.html
	        exit -1
	fi
fi

##############################################################
## Step 3: Map the sample				    ##
##############################################################

if [ ! -e ${SAMPLE}.sam ]
then
	echo "Running BWA"
	date '+%m/%d/%y %H:%M:%S'
	echo

	$BWA mem -t 8 -M -R '@RG\tID:$SAMPLE\tSM:$SAMPLE\tLB:$SAMPLE\tPL:illumina' $BWA_REFERENCE $FASTQ1 $FASTQ2 > ${SAMPLE}.sam

	if [ $? -ne 0 ]
	then
		echo "Error running BWA"
		rm $SAMPLE.sam
	        exit -1
	fi
fi

##############################################################
## Step 4.  Sort the BAM file	
##############################################################

if [ ! -e ${SAMPLE}.bam ]
then
	$JAVA -Xmx4G -Djava.io.tmpdir=. -jar $PICARD SortSam SO=coordinate INPUT=${SAMPLE}.sam OUTPUT=${SAMPLE}.bam VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
fi

##############################################################
## Step 5.  Collect alignment metrics		
##############################################################
if [ ! -e ${SAMPLE}.alignment_summary_metrics.txt ]
then
	echo Collection alignment summary metrics
	date '+%m/%d/%y %H:%M:%S'
        echo

	$JAVA -Xmx4G -Djava.io.tmpdir=. -jar $PICARD CollectAlignmentSummaryMetrics \
		I=${SAMPLE}.bam \
		O=${SAMPLE}.alignment_summary_metrics.txt \
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

##############################################################
## Step 5.  Mark Duplicates
##############################################################
if [ ! -e ${SAMPLE}.dedup.bam ]
then
	echo Marking duplicates
	date '+%m/%d/%y %H:%M:%S'
	echo

	$JAVA -Djava.io.tmpdir=. -jar $PICARD MarkDuplicates \
		INPUT=${SAMPLE}.bam \
		OUTPUT=${SAMPLE}.dedup.bam \
		METRICS_FILE=${SAMPLE}.dedup_metrics.txt \
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

##############################################################
## Step 6.  Collect target region coverage
##############################################################
if [ ! -e $SAMPLE.targetRegionCoverage.txt ]
then
	$JAVA -Xmx4G -Djava.io.tmpdir=. -jar $RGTOOLS/CalculateTargetRegionCoverage.jar B=${SAMPLE}.dedup.bam G=$TARGETREGIONS_BED O=$SAMPLE.targetRegionCoverage.txt
fi

##############################################################
# Step 7. Collect hsMetrics and per target metrics
##############################################################
if [ ! -e $SAMPLE.hsmetrics.txt ]
then
	$JAVA -Xmx4G -jar $PICARD CalculateHsMetrics BI=$BAITREGIONS_INTERVAL TI=$TARGETREGIONS_INTERVAL REFERENCE_SEQUENCE=$REFERENCE I=${SAMPLE}.dedup.bam O=$SAMPLE.hsmetrics.txt PER_TARGET_COVERAGE=$SAMPLE.per_target_metrics.txt
fi

##############################################################################
# 13.  Transfer data to final resting place
##############################################################################
cd ..
if [ $AWS -eq 0 ]
then
	echo "Moving ${SAMPLE_DIR} to ${OUTDIR}/${SAMPLE}"
	mv ${SAMPLE_DIR} ${OUTDIR}/${SAMPLE}
else
	aws s3 cp --recursive $SAMPLE_DIR s3://bmsrd-ngs-projects/${PROJECT}/${SAMPLE}
	
	if [ $? -ne 0 ]
	then
		echo Error copying $SAMPLE_DIR
		exit -1
	fi

	#rm -rf $SAMPLE_DIR 2>/dev/null
fi

echo
echo Finished analyzing ${SAMPLE} 
date '+%m/%d/%y %H:%M:%S' 
echo
