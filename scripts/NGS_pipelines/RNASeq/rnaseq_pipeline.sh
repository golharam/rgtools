#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

# runRNASeqPipeline will pass the following env vars:
# SAMPLE, FASTQ1, FASTQ2, TMP_DIR, REFMODEL

###############################################################################
# Set default values.  These can be set via command-line options only.  
# Setting them here overrides environment variables
###############################################################################
VERSION=0.5.4
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

		-h|--help)
		HELP=1
		;;

                -o|--outdir)
                OUTDIR="$2"
                shift # past argument
                ;;
	
		-s|--sample)
		SAMPLE="$2"
		shift # past argument
		;;

                --sample-dir)
                SAMPLE_DIR="$2"
                shift # past argument
                ;;

                --refmodel)
                REFMODEL="$2"
                shift # past argument
                ;;

                --tmpdir)
                TMP_DIR="$2"
                shift # past argument
                ;;

		-t|--threads)
		THREADS="$2"
		shift # past argument
		;;

                --use-star)
                USE_STAR=1
                ;;

		*)
		echo "Unknown option: $key"
		exit -1
		;;
	esac
	shift	# past argument or value
done

if [ $HELP == 1 ]; then
        echo "Usage: rnaseq_pipeline.sh [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>]"
        echo "Options:"
        echo "  -h|--help"
        echo "  -o|--outdir <output directory> [default=current working directory]"
        echo "  --refmodel <reference model> (hg19ERCC/mm10ERCC/rn6ERCC) (default=hg19ERCC)"
        echo "  --sample-dir (analysis directory to use. Use same tmp dir, for debugging"
        exit 2
fi

# If variables are not set as an environment variable or as a parameter, set them to default here.
if [ -z "$AWS" ]; then
	AWS=0
fi
if [ -z "$OUTDIR" ]; then
	OUTDIR=`pwd`
fi
if [ ! -e $OUTDIR ]; then
	echo "$OUTDIR does not exist"
	exit -1
fi
if [ -z "$REFMODEL" ]; then
	REFMODEL="hg19ERCC.refSeq"
fi
if [ -z "$TMP_DIR" ]; then
	TMP_DIR=/scratch
fi
if [ -z "$USE_STAR" ]; then
	USE_STAR=0
fi

if [ -e /ngs/apps ]; then
	EXT_PKGS_DIR=/ngs/apps
	REFERENCE_DIR=/ngs/reference
else
	EXT_PKGS_DIR=/apps/sys/galaxy/external_packages
	REFERENCE_DIR=/ng18/galaxy/reference_genomes	
fi

# Applications / Programs
BEDTOOLS=$EXT_PKGS_DIR/bedtools-2.24.0/bin/bedtools
BOWTIE_DIR=$EXT_PKGS_DIR/bowtie-1.1.1
BOWTIE2=$EXT_PKGS_DIR/bowtie2-2.2.6/
FASTQC=$EXT_PKGS_DIR/FastQC-0.11.3-RG/fastqc
FASTQVALIDATOR=$EXT_PKGS_DIR/fastQValidator-0.1.1a/bin/fastQValidator
MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.137/picard.jar
R_DIR=$EXT_PKGS_DIR/R-3.2.2/bin
RGTOOLS_DIR=$EXT_PKGS_DIR/rgtools
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
#RSEM_DIR=$EXT_PKGS_DIR/RSEM-1.2.20
RSEM_DIR=$EXT_PKGS_DIR/rsem-1.2.22
SAMTOOLS=$EXT_PKGS_DIR/samtools-0.1.19/samtools
SRA_DIR=$EXT_PKGS_DIR/sratoolkit-2.5.2/bin
STAR=$EXT_PKGS_DIR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR
TOPHAT2=$EXT_PKGS_DIR/tophat-2.1.0/
UBU_DIR=$EXT_PKGS_DIR/ubu
UBU_JAR=$EXT_PKGS_DIR/ubu-1.2-jar-with-dependencies.jar

export PATH=$BOWTIE2:$R_DIR:$TOPHAT2:$PATH

# Reference data
CONTAMINATION_REFERENCE=$REFERENCE_DIR/contamination/contamination.fa
CONTAMINATION_REFERENCE_BOWTIE2=$REFERENCE_DIR/contamination/bowtie2_index/contamination
CONTAMINATION_REFERENCE_STAR=$REFERENCE_DIR/contamination/STAR_index
REFERENCE=`echo $REFMODEL | cut -f 1 -d '.'`
GENEMODEL=`echo $REFMODEL | cut -f 2 -d '.'`

###############################################################################
# Start analysis
###############################################################################
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "RNA-Seq Pipeline v$VERSION"
echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"
echo "THREADS: $THREADS"
echo "EXT_PKGS_DIR: $EXT_PKGS_DIR"
echo "REFERENCE: $REFERENCE"
echo "GENEMODEL: $GENEMODEL"
echo "SAMPLE_DIR: $SAMPLE_DIR"
echo "OUTDIR: $OUTDIR"
echo
analysis_date_started=$(date +"%s")

# Make tmp working directory aka SAMPLE_DIR
if [ -z "$SAMPLE_DIR" ]
then
	cd $TMP_DIR
	SAMPLE_DIR=`mktemp -d --tmpdir=${TMP_DIR} ${SAMPLE}_XXXXXX`
fi
if [ ! -d $SAMPLE_DIR ]
then
        mkdir $SAMPLE_DIR
fi
cd $SAMPLE_DIR

echo "Working in `uname -n`:`pwd`"
if [ $AWS -eq 0 ]
then
	echo "Output dir is $OUTDIR/$SAMPLE"
else
	echo "Output dir is ${S3BUCKET}/${PROJECT}/${SAMPLE}"
fi
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################################
# RSEM
##############################################################################
echo
echo "RSEM"
date '+%m/%d/%y %H:%M:%S'
echo

gunzip -c $FASTQ1 > ${SAMPLE}_1.fq
if [ ! -n "$FASTQ2" ]; then
	$RSEM_DIR/rsem-calculate-expression \
		-p $THREADS --no-bam-output --bowtie-path $BOWTIE_DIR \
		${SAMPLE}_1.fq $REFERENCE_DIR/$REFERENCE/rsem_${GENEMODEL,,}_index/$REFERENCE $SAMPLE
else
	gunzip -c $FASTQ2 > ${SAMPLE}_2.fq
	$RSEM_DIR/rsem-calculate-expression \
		-p $THREADS --no-bam-output --bowtie-path $BOWTIE_DIR \
		--paired-end ${SAMPLE}_1.fq ${SAMPLE}_2.fq $REFERENCE_DIR/$REFERENCE/rsem_${GENEMODEL,,}_index/$REFERENCE $SAMPLE 
	rm ${SAMPLE}_2.fq
fi
rm ${SAMPLE}_1.fq

if [ $? -ne 0 ]
then
	echo "Error running RSEM"
	exit -1
fi

echo Completed RSEM
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################################
# 2. Cleanup large intermediate output
##############################################################################

##############################################################################
# 13.  Transfer data to final resting place
##############################################################################
cd ..
#if [ $AWS -eq 0 ]; then
echo "Moving ${SAMPLE_DIR} to ${OUTDIR}/${SAMPLE}"
mv ${SAMPLE_DIR} ${OUTDIR}/${SAMPLE}
#else
#	aws s3 cp --recursive $SAMPLE_DIR ${BUCKET}/${SAMPLE}

#	if [ $? -ne 0 ]
#	then
#		echo Error copying $SAMPLE_DIR to s3
#		exit -1
#	fi

#	rm -rf $SAMPLE_DIR 2>/dev/null
#fi

##############################################################################
# Done.
##############################################################################
echo
echo $SAMPLE Complete
date '+%m/%d/%y %H:%M:%S'
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
