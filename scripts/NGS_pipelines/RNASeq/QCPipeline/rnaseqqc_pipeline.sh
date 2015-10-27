#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

# runRNASeqQCPipeline will pass the following env vars:
# SAMPLE, FASTQ1, FASTQ2, SUBSAMPLE, TMP_DIR, REFMODEL

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
		
		--sraftp)
		SRAFTP="$2"
		shift # past argument
		;;

		--subsample)
		SUBSAMPLE="$2"
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
	echo "Usage: rnaseqqc_pipeline.sh [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>]"
	echo "Options:"
	echo "  -h|--help"
	echo "  -o|--outdir <output directory> [default=current working directory]"
	echo "  --refmodel <reference model> (hg19ERCC/mm10ERCC/rn6ERCC) (default=hg19ERCC)"
	echo "  --sample-dir (analysis directory to use. Use same tmp dir, for debugging"
	echo "  --subsample [# of reads]"
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
if [ -z "$SUBSAMPLE" ]; then
	SUBSAMPLE=0
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
BOWTIE2=$EXT_PKGS_DIR/bowtie2-2.2.6/
FASTQC=$EXT_PKGS_DIR/FastQC-0.11.3-RG/fastqc
FASTQVALIDATOR=$EXT_PKGS_DIR/fastQValidator-0.1.1a/bin/fastQValidator
MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.137/picard.jar
R_DIR=$EXT_PKGS_DIR/R-3.2.2/bin
RGTOOLS_DIR=$EXT_PKGS_DIR/rgtools
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
RSEM_DIR=$EXT_PKGS_DIR/RSEM-1.2.20
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

# Start Analysis
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "Pipeline: $VERSION"
echo "Reference Model: $REFMODEL"
echo "Subsample: $SUBSAMPLE"
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
# Step 1: Download the file to the local scratch space
# If FTP, download from FTP -> [${SAMPLE}_1.fastq.gz, ${SAMPLE}_2.fastq.gz]
# If AWS, copy  from kraken -> [${SAMPLE}_1.fastq.gz, ${SAMPLE}_2.fastq.gz]
# Else assume file is local and unknown if compressed or not.  
# Regardless, $FASTQ1 and $FASTQ2 point to the files.
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

		# Single-end reads are extracted to ${SAMPLE}.fastq.gz
		# Paired-end reads are in ${SAMPLE}_1/2.fastq.gz 
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

	FASTQ1=`pwd`/${SAMPLE}_1.fastq.gz
	if [ -e ${SAMPLE}_2.fastq.gz ]
	then
		FASTQ2=`pwd`/${SAMPLE}_2.fastq.gz
	fi
elif [[ $FASTQ1 == ftp* ]] # If FTP location, then download from FTP
then
	echo "Downloading $FASTQ1"
	date '+%m/%d/%y %H:%M:%S'
	echo

	FILENAME=`basename $FASTQ1`
	wget -nv $FASTQ1
	
	if [ $? -ne 0 ]
	then
	        echo "Error downloading $FASTQ1"
	        rm -f $FILENAME 2>/dev/null
        	exit -1
	fi
	
	if [ -e ${SAMPLE}_1.fastq.gz ]
	then
		mv $FILENAME ${SAMPLE}_1.fastq.gz
	fi
	FASTQ1=`pwd`/${SAMPLE}_1.fastq.gz

	# Repeat for FastQ2
	echo "Downloading $FASTQ2"
        date '+%m/%d/%y %H:%M:%S'
        echo

        FILENAME=`basename $FASTQ2`
        wget -nv $FASTQ2

        if [ $? -ne 0 ]
        then
                echo "Error downloading $FASTQ2"
                rm -f $FILENAME 2>/dev/null
                exit -1
        fi

	if [ -e ${SAMPLE}_2.fastq.gz ]
	then
		mv $FILENAME ${SAMPLE}_2.fastq.gz
	fi
	FASTQ2=`pwd`/${SAMPLE}_2.fastq.gz
elif [ $AWS -eq 1 ] # Copy the data from kraken
then
	echo "Downloading from kraken: $FASTQ1"
	date '+%m/%d/%y %H:%M:%S'
	echo
	
	scp golharr@kraken.pri.bms.com:$FASTQ1 ./${SAMPLE}_1.fastq.gz
	
        if [ $? -ne 0 ]
        then
                echo "Error downloading $FASTQ1"
                rm -f ${SAMPLE}_1.fastq.gz 2>/dev/null
                exit -1
        fi
	FASTQ1=`pwd`/${SAMPLE}_1.fastq.gz
	
	echo "Downloading from kraken: $FASTQ2"
	date '+%m/%d/%y %H:%M:%S'
	echo
	
	scp golharr@kraken.pri.bms.com:$FASTQ2 ./${SAMPLE}_2.fastq.gz
	
        if [ $? -ne 0 ]
        then
                echo "Error downloading $FASTQ2"
                rm -f ${SAMPLE}_2.fastq.gz 2>/dev/null
                exit -1
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
# Step 3: Unzip if not already unzipped
# Output [${SAMPLE}_1.fastq, ${SAMPLE}_2.fastq]
##############################################################################
if [ ! -e ${SAMPLE}_1.fastq ]; then
	if [ $(file $FASTQ1 | cut -d' ' -f2) == "gzip" ]; then
		echo "Unzipping $FASTQ1 > ${SAMPLE}_1.fastq"
		date '+%m/%d/%y %H:%M:%S'
		echo

		gunzip -dc $FASTQ1 > ${SAMPLE}_1.fastq

		if [ $? -ne 0 ]; then
			echo "Error Unzipping $FASTQ1 > ${SAMPLE}_1.fastq"
			rm ${SAMPLE}_1.fastq
			exit -1
		fi
	else
		echo
		echo "$FASTQ1 is not gzip compressed:"
		echo $(file $FASTQ1)
		echo

		ln -s $FASTQ1 ${SAMPLE}_1.fastq
	fi
fi

if [ -n "$FASTQ2" ] && [ ! -e ${SAMPLE}_2.fastq ]; then
	if [ $(file $FASTQ2 | cut -d' ' -f2) == "gzip" ]; then
		echo "Unzipping $FASTQ2 > ${SAMPLE}_2.fastq"
		date '+%m/%d/%y %H:%M:%S'
		echo

		gunzip -dc $FASTQ2 > ${SAMPLE}_2.fastq

		if [ $? -ne 0 ]; then
			echo "Error Unzipping $FASTQ2 > ${SAMPLE}_2.fastq"
			rm ${SAMPLE}_2.fastq
			exit -1
		fi
	else
		echo
		echo "$FASTQ2 is not gzip compressed"
		echo $(file $FASTQ2)
		echo

		ln -s $FASTQ2 ${SAMPLE}_2.fastq
	fi
fi

##############################################################################
# Step 4: Subsample (if not subsampling, link to original fastq files)
# Output: [${SAMPLE}_1.subsampled.fastq, ${SAMPLE}_2.subsampled.fastq]
##############################################################################
if [ "$SUBSAMPLE" -ne "0" ]; then
	echo
	echo "Subsampling for $SUBSAMPLE reads..."
	date '+%m/%d/%y %H:%M:%S'
	echo

	NUMSEQS=`head -n 1 $SAMPLE.fq1Validator.txt | cut -f 8 -d ' '`
	echo "Sample has $NUMSEQS sequences"
	if [ -n "$FASTQ2" ]; then
		RandomSubFq -t $NUMSEQS -w $SUBSAMPLE -i $FASTQ1 -i $FASTQ2 -o ${SAMPLE}_1.subsampled.fastq -o ${SAMPLE}_2.subsampled.fastq
	else
		RandomSubFq -t $NUMSEQS -w $SUBSAMPLE -i $FASTQ1 -o ${SAMPLE}_1.subsampled.fastq
	fi

	if [ $? -ne 0 ]; then
		echo "Error subsampling"
		rm ${SAMPLE}_?.subsampled.fastq
		exit -1
	fi

else
	echo "Not subsampling...using entire set of reads..."
	echo
	
	if [ ! -e ${SAMPLE}_1.subsampled.fastq ]; then
		ln -s ${SAMPLE}_1.fastq ${SAMPLE}_1.subsampled.fastq
	fi
	if [ -e ${SAMPLE}_2.fastq ] && [ ! -e ${SAMPLE}_2.subsampled.fastq ]; then
		ln -s ${SAMPLE}_2.fastq ${SAMPLE}_2.subsampled.fastq
	fi
fi

##############################################################################
# Step 5: Run FastQC
##############################################################################
if [ ! -e "${SAMPLE}_1.subsampled_fastqc.zip" ]
then
	echo "Running FastQC"
	date '+%m/%d/%y %H:%M:%S'
	echo

	if [ ! -e "${SAMPLE}_2.subsampled.fastq" ]; then
		$FASTQC --quiet --outdir=. --extract -t 2 ${SAMPLE}_1.subsampled.fastq
	else 
		$FASTQC --quiet --outdir=. --extract -t 2 ${SAMPLE}_1.subsampled.fastq ${SAMPLE}_2.subsampled.fastq 
	fi

	if [ $? -ne 0 ]
	then
		echo "Error running FastQC"
		rm -rf *_fastqc *_fastqc.zip *_fastqc.html
	        exit -1
	fi
fi

##############################################################################
# Step 6: Perform contamination detection
# Output -> [${SAMPLE}.uncontaminated.fastq.1.gz, ${SAMPLE}.uncontaminated.fastq.2.gz]
##############################################################################
if [ ! -d contamination ]
then
	echo "Checking for bacterial/viral contamination (USE_STAR=$USE_STAR)"
	date '+%m/%d/%y %H:%M:%S'
	echo

	mkdir contamination
	cd contamination
	
	if [ $USE_STAR -eq 1 ]
	then
		$STAR --runThreadN $THREADS --genomeDir $CONTAMINATION_REFERENCE_STAR \
			--readFilesIn ../${SAMPLE}_1.subsampled.fastq ../${SAMPLE}_2.subsampled.fastq \
			--outSAMtype BAM SortedByCoordinate \
			--outReadsUnmapped Fastx \
			--genomeLoad LoadAndKeep 

		if [ $? -ne 0 ]; then
			echo "Error running STAR for contamination:"
			cat Log.out
			cd ..
			rm -rf contamination
			exit -1
		fi
		
		mv Aligned.sortedByCoord.out.bam $SAMPLE.contaminated.bam
		$SAMTOOLS index $SAMPLE.contaminated.bam
		mv Unmapped.out.mate1 ../${SAMPLE}_1.uncontaminated.fastq
		mv Unmapped.out.mate2 ../${SAMPLE}_2.uncontaminated.fastq
	else
		if [ ! -e ../${SAMPLE}_2.subsampled.fastq ]; then
			$BOWTIE2/bowtie2 \
				--un-gz $SAMPLE.uncontaminated.fastq.gz \
				--al-gz $SAMPLE.contaminated.fastq.gz \
				-p $THREADS \
				-x $CONTAMINATION_REFERENCE_BOWTIE2 \
				-U ../${SAMPLE}_1.subsampled.fastq \
				-S $SAMPLE.contaminated.sam 2>$SAMPLE.contamination.log
		else
			$BOWTIE2/bowtie2 --no-mixed --un-conc-gz $SAMPLE.uncontaminated.fastq.gz \
				--al-conc-gz $SAMPLE.contaminated.fastq.gz \
				-p $THREADS -1 ../${SAMPLE}_1.subsampled.fastq -2 ../${SAMPLE}_2.subsampled.fastq \
				--no-unal --rg-id $SAMPLE \
				--rg "SM:$SAMPLE\tLB:$SAMPLE\tPL:illumina" \
				-S $SAMPLE.contaminated.sam -x $CONTAMINATION_REFERENCE_BOWTIE2 2>$SAMPLE.contamination.log
		fi
	
		if [ $? -ne 0 ]; then
			echo "Error running bowtie2 for contamination:"
			cat $SAMPLE.contamination.log
			cd ..
			rm -rf contamination
			exit -1
		fi

		echo Sorting contamination BAM file
		date '+%m/%d/%y %H:%M:%S'
		echo

		java -Xmx4G -jar $PICARD_JAR SortSam \
			INPUT=$SAMPLE.contaminated.sam \
			OUTPUT=$SAMPLE.contaminated.bam \
			CREATE_INDEX=true \
			SO=coordinate \
			TMP_DIR=.

		if [ $? -ne 0 ]
		then
			echo "Error sorting contamination bam file"
			rm $SAMPLE.contaminated.bam
			exit -1
		fi

		echo Collect Alignment Summary Metrics on Contamination Data
		date '+%m/%d/%y %H:%M:%S'
		echo

		java -Xmx4G -jar $PICARD_JAR CollectAlignmentSummaryMetrics \
			R=$CONTAMINATION_REFERENCE \
			INPUT=$SAMPLE.contaminated.bam \
			OUTPUT=${SAMPLE}.contaminated.alnMetrics.txt

		if [ $? -ne 0 ]
		then
			echo "Error collecting contamination alignment summary metrics"
			rm ${SAMPLE}.contaminated.alnMetrics.txt
			exit -1
		fi		

		rm $SAMPLE.contaminated.sam
		if [ ! -e ../${SAMPLE}_2.subsampled.fastq ]; then
			mv ${SAMPLE}.uncontaminated.fastq.gz ../${SAMPLE}.uncontaminated.fastq.1.gz
		else
			mv ${SAMPLE}.uncontaminated.fastq.1.gz ../
			mv ${SAMPLE}.uncontaminated.fastq.2.gz ../
		fi
	fi
	
	cd ..
fi

##############################################################################
# Step 7: Run tophat2 or STAR
##############################################################################
if [ $USE_STAR -eq 1 ]
then
	echo
	echo Running STAR
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")
	echo
	
	mkdir star
	cd star
	
	$STAR --runThreadN $THREADS --genomeDir $REFERENCE_STAR \
		--readFilesIn ../${SAMPLE}_1.uncontaminated.fastq ../${SAMPLE}_2.uncontaminated.fastq \
		--outSAMtype BAM SortedByCoordinate --genomeLoad LoadAndKeep
		
	if [ $? -ne 0 ]; then
		echo "Error running STAR:"
		cat Log.out
		cd ..
		rm -rf star
		exit -1
	fi
	
	cd ..
	
	echo Finished running STAR
	date '+%m/%d/%y %H:%M:%S'
	date2=$(date +"%s")
	diff=$(($date2-$date1))
        echo "Running STAR took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
elif [ ! -d tophat_out ]
then
	echo
	echo Running tophat2
	date '+%m/%d/%y %H:%M:%S'
	date1=$(date +"%s")
	echo

	# TBD: Output unaligned reads as well, or else CollectAlignmentMetrics thinks all the reads aligned.
	if [ ! -e $SAMPLE.uncontaminated.fastq.2.gz ]
	then
		$TOPHAT2/tophat2 -p $THREADS \
			--rg-id 1 --rg-sample $SAMPLE --rg-library $SAMPLE --rg-platform illumina \
			$REFERENCE_DIR/$REFERENCE/bowtie2_index/$REFERENCE \
			$SAMPLE.uncontaminated.fastq.1.gz 
	else
                $TOPHAT2/tophat2 -p $THREADS \
                        --rg-id 1 --rg-sample $SAMPLE --rg-library $SAMPLE --rg-platform illumina \
                        $REFERENCE_DIR/$REFERENCE/bowtie2_index/$REFERENCE \
                        $SAMPLE.uncontaminated.fastq.1.gz $SAMPLE.uncontaminated.fastq.2.gz
	fi

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
        echo "Running tophat2 took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
	echo
fi

##############################################################################
# Step 7: Re-sort sequenences
##############################################################################
if [ ! -e ${SAMPLE}.bam ]
then
	echo Resorting BAM file
	date '+%m/%d/%y %H:%M:%S'
	echo

	if [ $USE_STAR -eq 1 ]
	then
		mv star/Aligned.sortedByCoord.out.bam ${SAMPLE}.bam
		$SAMTOOLS index ${SAMPLE}.bam
	else
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
		rm tophat_out/accepted_hits.bam
	fi	
fi

##############################################################################
# Step 8: Collect Alignment Summary Metrics
##############################################################################
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

##############################################################################
# 9.  Collect Insert Size Metrics
#     TBD: Is this really needed?  This is mainly for DNA-Seq, I believe.
##############################################################################
if [ ! -e ${SAMPLE}.insertSizeMetrics.txt ]
then
        echo Collect InsertSize Metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

	# For single-end data, this will product no output, so don't bother
	# running it.
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

##############################################################################
# 10.  Collect RNA Seq Metrics
##############################################################################
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
		CHART_OUTPUT=$SAMPLE.rnaseq.pdf \
		TMP_DIR=.

        if [ $? -ne 0 ]
        then
                echo "Error collecting rnaseq metrics"
                rm ${SAMPLE}.rnaseq*
                exit -1
        fi

fi

##############################################################################
# 10b.  Collect RNA Seq Metrics for hemoglobin
# Good article about hemoglobin genes: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3543078/
##############################################################################
if [ ! -e ${SAMPLE}.hemoglobinMetrics.txt ]
then
        echo "Collecting Hemoglobin Metrics"
        date '+%m/%d/%y %H:%M:%S'
        echo

	java -Xmx4G -jar $RGTOOLS_DIR/CalculateTargetRegionCoverage.jar \
                B=${SAMPLE}.bam \
		BED=$REFERENCE_DIR/$REFERENCE/annotation/refSeq/hemoglobin.bed \
		OUTPUT=$SAMPLE.hemoglobinMetrics.txt

        if [ $? -ne 0 ]
        then
                echo "Error collecting hemoglobin metrics"
                rm ${SAMPLE}.hemoglobinMetrics*
                exit -1
        fi

fi

##############################################################################
# 11.  Collect ERCC Metrics
##############################################################################
if [ ! -e ${SAMPLE}.idxStats.txt ]
then
	echo
        echo Collecting Index Metrics
        date '+%m/%d/%y %H:%M:%S'
        echo

	$SAMTOOLS idxstats ${SAMPLE}.bam > ${SAMPLE}.idxStats.txt
	
        if [ $? -ne 0 ]
        then
                echo "Error collecting index metrics:"
                cat ${SAMPLE}.idxStats.txt
                rm ${SAMPLE}.idxStats.txt
                exit -1
        fi

fi

##############################################################################
# 12. Cleanup intermediate output
##############################################################################
if [ $AWS -eq 1 ]
then
	rm *.gz
fi

# Removing *.gz files on local files 
rm *.fastq

##############################################################################
# 13.  Transfer data to final resting place
##############################################################################
cd ..
if [ $AWS -eq 0 ]
then
	echo "Moving ${SAMPLE_DIR} to ${OUTDIR}/${SAMPLE}"
	mv ${SAMPLE_DIR} ${OUTDIR}/${SAMPLE}
else
	aws s3 cp --recursive $SAMPLE_DIR ${S3BUCKET}/${PROJECT}/${SAMPLE}

	if [ $? -ne 0 ]
	then
		echo Error copying $SAMPLE_DIR
		exit -1
	fi

	rm -rf $SAMPLE_DIR 2>/dev/null
fi

##############################################################################
# Done.
##############################################################################
echo
echo $SAMPLE Complete
date '+%m/%d/%y %H:%M:%S'
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
