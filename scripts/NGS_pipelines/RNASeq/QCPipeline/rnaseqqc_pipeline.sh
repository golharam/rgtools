#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

###############################################################################
# Set default values.  These can be set via command-line options only.  
# Setting them here overrides environment variables
###############################################################################
VERSION=0.1
HELP=0
AWS=0
OUTDIR=`pwd`
DELETE_INTERMEDIATE=0
THREADS=8
REFERENCE=hg19
SUBSAMPLE=0

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
		--delete-intermediate)
		DELETE_INTERMEDIATE=1
		;;

		-s|--sample)
		SAMPLE="$2"
		shift # past argument
		;;

		-1|--fastq1)
		FASTQ1="$2"
		shift # past argument
		;;

		-2|--fastq2)
		FASTQ2="$2"
		shift # past argument
		;;

		-t|--threads)
		THREADS="$2"
		shift # past argument
		;;

		-a|--aws)
		AWS=1
		;;

		--sraftp)
		SRAFTP="$2"
		shift # past argument
		;;

		--subsample)
		SUBSAMPLE="$2"
		shift # past argument
		;;

		-h|--help)
		HELP=1
		;;

                -o|--outdir)
                OUTDIR="$2"
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
	echo "Usage 2: $0 <options> <usage>"
	echo "where usage:"
	echo "Usage 1: [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>]"
	echo "Usage 2: [-s|--sample <SRA ID>] [--sraftp <ftp location>]" 
	echo "Options:"
	echo "	-a|--aws"
	echo "  --delete-intermediate (delete intermediate files, not including source fq.gz) (not yet implemented)"
	echo "	-h|--help"
	echo "  --subsample [# of reads]"
	echo "  -t|--threads"
	echo "  -o|--outdir <output directory> [default=current working directory]"
	exit 0
fi

if [ $AWS == 0 ]; then
	EXT_PKGS_DIR=/apps/sys/galaxy/external_packages
	REFERENCE_DIR=/ng18/galaxy/reference_genomes
	TMP_DIR=/scratch
	
else
	EXT_PKGS_DIR=/ngs/apps
	REFERENCE_DIR=/ngs/reference/
	TMP_DIR=/scratch
fi

# Applications / Programs
BEDTOOLS=$EXT_PKGS_DIR/bedtools-2.23.0/bin
BOWTIE2=$EXT_PKGS_DIR/bowtie2-2.2.6/bowtie2
FASTQC=$EXT_PKGS_DIR/FastQC-0.11.2/fastqc
FASTQVALIDATOR=$EXT_PKGS_DIR/fastQValidator-0.1.1a/bin/fastQValidator
MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.137/picard.jar
RGTOOLS_DIR=$EXT_PKGS_DIR/rgtools
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
RSEM_DIR=$EXT_PKGS_DIR/RSEM-1.2.20
SAMTOOLS=$EXT_PKGS_DIR/samtools-0.1.19/samtools
SRA_DIR=$EXT_PKGS_DIR/sratoolkit-2.5.2/bin
TOPHAT2=$EXT_PKGS_DIR/tophat-2.1.0/tophat2
UBU_DIR=$EXT_PKGS_DIR/ubu
UBU_JAR=$EXT_PKGS_DIR/ubu-1.2-jar-with-dependencies.jar

# Reference data
CONTAMINATION_REFERENCE=$REFERENCE_DIR/contamination/bowtie2_index/contamination

# Start Analysis
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "Subsample: $SUBSAMPLE"
echo

# Make tmp working directory aka SAMPLE_DIR
if [ -z "$SAMPLE_DIR" ]
then
	cd $TMP_DIR
	SAMPLE_DIR=`mktemp -d --tmpdir=${TMP_DIR} ${SAMPLE}_XXXXXX`
	if [ ! -d $SAMPLE_DIR ]
	then
	        mkdir $SAMPLE_DIR
	fi
fi
cd $SAMPLE_DIR
echo "Working in `uname -n`:`pwd`"
echo "Output dir is $OUTDIR"
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################################
# Step 1: Download the file to the local scratch space
# If SRA, download from SRA
# If FTP, download from FTP
##############################################################################
if [ -n "$SRAFTP" ]
then
	if [ ! -e ${SAMPLE}.sra ]
	then
		echo "Downloading SRA Sample $SAMPLE"
		date '+%m/%d/%y %H:%M:%S'
		echo

		wget -nv $SRAFTP/${SAMPLE}.sra 

		if [ $? -ne 0 ]; then
			echo "Error downloading $SRAFTP/${SAMPLE}.sra"
			rm -f ${SAMPLE}.sra
			exit -1
		fi

		echo "Finished download."
		date '+%m/%d/%y %H:%M:%S'
		echo
	fi

	if [ ! -e ${SAMPLE}_1.fastq.gz ]
	then	
		echo "Extracting FastQ file(s) from ${SAMPLE}.sra"
		date '+%m/%d/%y %H:%M:%S'
		echo

		$SRA_DIR/fastq-dump --split-3 --gzip ${SAMPLE}.sra

		if [ $? -ne 0 ]; then
			rm -f *.gz
			echo "Error extracting FASTQ file(s)"
			exit -1
		fi

		if [ ! -e ${SAMPLE}_1.fastq.gz ]
		then
			mv ${SAMPLE}.fastq.gz ${SAMPLE}_1.fastq.gz
		fi

		rm ${SAMPLE}.sra
		touch ${SAMPLE}.sra

		echo "Finished extracting fastq files"
	        date '+%m/%d/%y %H:%M:%S'
        	echo
	fi

	FASTQ1=${SAMPLE}_1.fastq.gz
	if [ -e ${SAMPLE}_2.fastq.gz ]
	then
		FASTQ2=${SAMPLE}_2.fastq.gz
	fi
fi
# If FTP location, then download from FTP
if [[ $FASTQ1 == ftp* ]]
then
	echo "Downloading $FASTQ1"
	date '+%m/%d/%y %H:%M:%S'
	echo

	FILENAME=`basename $FASTQ1`
	wget -nv $FASTQ1
	
	if [ $? -ne 0 ]
	then
	        echo "Error downloading $FASTQ1"
	        rm -f $FASTQ1 2>/dev/null
        	exit -1
	fi
	FASTQ1=$FILENAME

	# Repeat for FastQ2
	echo "Downloading $FASTQ2"
        date '+%m/%d/%y %H:%M:%S'
        echo

        FILENAME=`basename $FASTQ2`
        wget -nv $FASTQ2

        if [ $? -ne 0 ]
        then
                echo "Error downloading $FASTQ2"
                rm -f $FASTQ2 2>/dev/null
                exit -1
        fi
        FASTQ2=$FILENAME
fi

if [ ! -e $FASTQ1 ]
then
	echo "Unable to locate $FASTQ1.  Make sure full path is provided"
	exit -1
fi
echo "Using Fastq file(s):"
echo FASTQ1: $FASTQ1
if [ -n "$FASTQ2" ]
then
	echo FASTQ2: $FASTQ2
fi
date '+%m/%d/%y %H:%M:%S'
analysis_date_started=$(date +"%s")
echo

##############################################################################
# Step 2: Run FastQ Validator
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
# Step 3: Run FastQC
##############################################################################
FASTQC1=`basename $FASTQ1 .fastq.gz`
if [ ! -e "${FASTQC1}_fastqc.zip" ]
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

##############################################################################
# Step 4: Perform contamination detection
##############################################################################
if [ ! -d contamination ]
then
	echo "Checking for bacterial/viral contamination"
	date '+%m/%d/%y %H:%M:%S'
	echo

	mkdir contamination
	cd contamination
	$BOWTIE2 --no-mixed --un-conc-gz $SAMPLE.uncontaminated.fastq.gz \
		 --al-conc-gz $SAMPLE.contaminated.fastq.gz \
		 -p $THREADS -1 $FASTQ1 -2 $FASTQ2 \
		 --no-unal --rg-id $SAMPLE \
		 --rg "SM:$SAMPLE\tLB:$SAMPLE\tPL:illumina" \
		 -S $SAMPLE.contaminated.sam -x $CONTAMINATION_REFERENCE 2>$SAMPLE.contamination.log

	if [ $? -ne 0 ]; then
		echo "Error running bowtie2 for contamination"
		cd ..
		rm -rf contamination
		exit -1
	fi
	
	mv ${SAMPLE}.uncontaminated.fastq.1.gz ../${SAMPLE}_1.fastq.gz
	mv ${SAMPLE}.uncontaminated.fastq.2.gz ../${SAMPLE}_2.fastq.gz
	
	cd ..
fi

##############################################################################
# Step 5: Subsample
##############################################################################
if [ "$SUBSAMPLE" -ne "0" ]; then
	if [ ! -e ${SAMPLE}_1.fastq ]; then
		echo "Subsampling for $SUBSAMPLE reads..."
		date '+%m/%d/%y %H:%M:%S'
		echo

		if [ -n "$FASTQ2" ]; then
			RandomSubFq -w $SUBSAMPLE -i $FASTQ1 -i $FASTQ2 -o ${SAMPLE}_1.fastq -o ${SAMPLE}_2.fastq
		else
			RandomSubFq -w $SUBSAMPLE -i $FASTQ1 -o ${SAMPLE}_1.fastq
		fi

	        if [ $? -ne 0 ]; then
	                echo "Error subsampling"
	                rm ${SAMPLE}_?.fastq
	                exit -1
	        fi

		gzip ${SAMPLE}_1.fastq
		if [ -n "$FASTQ2" ]; then
			gzip ${SAMPLE}_2.fastq
		fi
	fi
else
	echo "Not subsampling...using entire set of reads..."
	echo
	# Nothing to do here since we've got the uncontaminated set of reads in
	# ${SAMPLE}_1.fastq.gz and ${SAMPLE}_2.fastq.gz
fi

##############################################################################
# Step 6: Run tophat2
##############################################################################
if [ ! -d tophat_out ]
then
	echo Running tophat2
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")
	echo

	# TBD: Output unaligned reads as well, or else CollectAlignmentMetrics thinks all the reads aligned.
	$TOPHAT2 -p $THREADS \
		 --rg-id 1 --rg-sample $SAMPLE --rg-library $SAMPLE --rg-platform illumina \
		$REFERENCE_DIR/$REFERENCE/bowtie2_index/$REFERENCE ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz

	if [ $? -ne 0 ] && [ ! -e tophat_out/accepted_hits.bam ]
	then
		echo "Encountered error running tophat2"
		rm -rf tophat_out 2>/dev/null
		exit -1
	fi

	echo Finished running tophat2
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
	diff=$(($date2-$date1))
        echo "Runng tophat2 took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

# 5. Re-sort sequenences
if [ ! -e ${SAMPLE}.bam ]
then
	echo Resorting BAM file
	date '+%m/%d/%y %H:%M:%S'
	echo

	java -Xmx4G -jar $PICARD_JAR ReorderSam \
		R=$REFERENCE_DIR/$REFERENCE/$REFERENCE.fa \
		INPUT=tophat_out/accepted_hits.bam \
		OUTPUT=${SAMPLE}.bam \
		CREATE_INDEX=true \
		TMP_DIR=.

        if [ $? -ne 0 ]
        then
                echo "Error resorting bam file"
                rm ${SAMPLE}.bam
                exit -1
        fi
fi

# 6.  Collect Alignment Summary Metrics
if [ ! -e ${SAMPLE}.alnMetrics.txt ]
then
	echo Collect Alignment Summary Metrics
	date '+%m/%d/%y %H:%M:%S'
	echo

	java -Xmx4G -jar $PICARD_JAR CollectAlignmentSummaryMetrics \
		R=$REFERENCE_DIR/$REFERENCE/$REFERENCE.fa \
		INPUT=${SAMPLE}.bam \
		OUTPUT=${SAMPLE}.alnMetrics.txt

	if [ $? -ne 0 ]
	then
		echo "Error collecting alignment summary metrics"
		rm ${SAMPLE}.alnMetrics.txt
		exit -1
	fi
fi

# 7.  Collect Insert Size Metrics
if [ ! -e ${SAMPLE}.insertSizeMetrics.txt ]
then
        echo Collect InsertSize Metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

        java -Xmx4G -jar $PICARD_JAR CollectInsertSizeMetrics \
                INPUT=${SAMPLE}.bam \
                OUTPUT=${SAMPLE}.insertSizeMetrics.txt \
		HISTOGRAM_FILE=$SAMPLE.insertSizeHistogram.pdf \
		TMP_DIR=.

        if [ $? -ne 0 ]
        then
                echo "Error collecting insert size metrics"
                rm ${SAMPLE}.insertSize*
                exit -1
        fi
fi

# 7.  Collect RNA Seq Metrics
if [ ! -e ${SAMPLE}.rnaseqMetrics.txt ]
then
        echo Collect RNASeq Metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

        java -Xmx4G -jar $PICARD_JAR CollectRnaSeqMetrics \
                INPUT=${SAMPLE}.bam \
		OUTPUT=$SAMPLE.rnaseqMetrics.txt \
		REF_FLAT=$REFERENCE_DIR/$REFERENCE/annotation/refSeq/refFlat.txt \
		RIBOSOMAL_INTERVALS=$REFERENCE_DIR/$REFERENCE/annotation/refSeq/rRNA.list \
		STRAND_SPECIFICITY=NONE \
		CHART_OUTPUT=$SAMPLE.rnaseq.pdf

        if [ $? -ne 0 ]
        then
                echo "Error collecting rnaseq metrics"
                rm ${SAMPLE}.rnaseq*
                exit -1
        fi

fi

# 19. Cleanup large intermediate output
if [ $AWS == 1 ]
then
	rm $FASTQ1 $FASTQ2
fi
rm -rf *.fastq tophat_out/ $SAMPLE.ba? 

cd ..
if [ $AWS == 0 ]
then
	echo "Moving ${SAMPLE_DIR} to ${OUTDIR}/${SAMPLE}"
	mv ${SAMPLE_DIR} ${OUTDIR}/${SAMPLE}
else
	# TBD: This is project specific and needs to be parameterized
	aws s3 cp --recursive $SAMPLE_DIR s3://bmsrd-ngs-M2GEN/

	if [ $? -ne 0 ]
	then
		echo Error copying $SAMPLE_DIR to s3
		exit -1
	fi

	rm -rf $SAMPLE_DIR 2>/dev/null
fi

echo
echo $SAMPLE Complete
date '+%m/%d/%y %H:%M:%S'
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."

