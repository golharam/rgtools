#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

# Pipeline is from: https://raw.githubusercontent.com/cc2qe/sandbox/master/unc_rnaseqV2_pipeline.sh
###############################################################################
# Set default values.  These can be set via command-line options only.  
# Setting them here overrides environment variables
###############################################################################
VERSION=0.04
HELP=0
DELETE_INTERMEDIATE=0
THREADS=8

if [ $# -eq 0 ]
then
	if [ -z "$SAMPLE" ]
	then
		HELP=1
	fi
fi

# Command line arguments and where to download is Ryan Golhar.
# Use > 1 to consume two arguments per pass in the loop (e.g. each argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g. some arguments don't have a corresponding value to go with it such as --help).
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

		-b|--bucket)
		BUCKET="$2"
		shift # past argument
		;;
		
		--delete-intermediate)
		DELETE_INTERMEDIATE=1
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

		--sraftp)
		SRAFTP="$2"
		shift # past argument
		;;

		--reference)
		REFERENCE="$2"
		shift # past argument
		;;

		-t|--threads)
		THREADS="$2"
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
	echo "Usage 1: $0 <options> [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>]"
	echo "Usage 2: $0 <options> [-s|--sample <SRA ID>] [--sraftp <ftp location>]" 
	echo "Options:"
	echo "	-a|--aws"
	echo "  -b|--bucket <bucket to transfer results to if running in AWS>"
	echo "  --delete-intermediate (delete intermediate files, not including source fq.gz) (not yet implemented)"
	echo "	-h|--help"
	echo "  -o|--outdir <output directory> [default=current working directory]/analysis"
	echo "  --reference <path to reference directory>"
	echo "  -t|--threads <threads to use> [default=8]"
	exit 0
fi

# If variables are not set as an environment variable or as a parameter, set them to default here.
if [ -z "$AWS" ]; then
	AWS=0
fi
if [ -z "$OUTDIR" ]; then
	OUTDIR=`pwd`/analysis
fi
if [ -z "$REFERENCE" ]; then
	REFERENCE=hg19_M_rCRS
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
BEDTOOLS=$EXT_PKGS_DIR/bedtools-2.23.0/bin
MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.129/picard.jar
RGTOOLS_DIR=$EXT_PKGS_DIR/rgtools
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
RSEM_DIR=$EXT_PKGS_DIR/rsem-1.2.22
SAMTOOLS=$EXT_PKGS_DIR/samtools-0.1.19/samtools
SRA_DIR=$EXT_PKGS_DIR/sratoolkit-2.5.2/bin
UBU_DIR=$EXT_PKGS_DIR/ubu
UBU_JAR=$EXT_PKGS_DIR/ubu-1.2-jar-with-dependencies.jar

# Reference data
#REFERENCE=$REFERENCE_DIR/$REFERENCE


# Start Analysis
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "Pipeline Version: $VERSION"
echo "AWS: $AWS"
echo "TmpDir: $TMP_DIR"
echo "Reference: $REFERENCE"
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
# If SRA, download from SRA -> [${SAMPLE}_1.fastq.gz, ${SAMPLE}_2.fastq.gz]
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
analysis_date_started=$(date +"%s")
echo

##############################################################################
# Step 0: Unzip if not already unzipped
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
		echo "$FASTQ1 is not gzip compressed:"
		echo $(file $FASTQ1)
		echo "Assuming uncompresed."

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
                echo "$FASTQ2 is not gzip compressed:"
                echo $(file $FASTQ2)
                echo "Assuming uncompresed."

                ln -s $FASTQ2 ${SAMPLE}_2.fastq
        fi
fi

##############################################################################
# 1. Format fastq 1 for Mapsplice
##############################################################################
if [ ! -e prep_1.fastq ]
then
	echo "Prepping prep_1.fastq"
        date '+%m/%d/%y %H:%M:%S'
        echo

	java -Xmx512M -jar $UBU_JAR fastq-format --phred33to64 --strip --suffix /1 --in ${SAMPLE}_1.fastq --out prep_1.fastq

	if [ $? -ne 0 ]; then
		echo "Error Prepping prep_1.fastq"
       		rm prep_1.fastq
                exit -1
        fi

	rm ${SAMPLE}_1.fastq
	touch ${SAMPLE}_1.fastq
fi

##############################################################################
# 2. Format fastq 2 for Mapsplice
##############################################################################
if [ -n "$FASTQ2" ] && [ ! -e prep_2.fastq ] 
then
	echo
	echo "Prepping prep_2.fastq"
        date '+%m/%d/%y %H:%M:%S'
        echo

	java -Xmx512M -jar $UBU_JAR fastq-format --phred33to64 --strip --suffix /2 --in ${SAMPLE}_2.fastq --out prep_2.fastq

        if [ $? -ne 0 ]; then
                echo "Error Prepping prep_2.fastq"
                rm prep_2.fastq
                exit -1
        fi

	rm ${SAMPLE}_2.fastq
	touch ${SAMPLE}_2.fastq
fi

##############################################################################
# 3. Mapsplice
##############################################################################
if [ ! -e $SAMPLE.mapsplice/alignments.bam ]
then

	## example command. this is for Mapsplice 12_07 though, i'm modifying it for Mapsplice 2.0.1.9, which they used in more recent samples.
	#python mapsplice_multi_thread.py --fusion --all-chromosomes-files hg19_M_rCRS/hg19_M_rCRS.fa --pairend -X 8 -Q fq --chromosome-filesdir hg19_M_rCRS/chromosomes --Bowtieidx hg19_M_rCRS/ebwt/humanchridx_M_rCRS -1 working/prep_1.fastq -2 working/prep_2.fastq -o SAMPLE_BARCODE 2> working/mapsplice.log
	## actual command
	#time python mapsplice.py --fusion --bam -p $THREADS -c ~/refdata/genomes/unc_tcga_hg19/chromosomes --qual-scale phred33 -x ~/refdata/genomes/unc_tcga_hg19/ebwt/humanchridx_M_rCRS -1 ${SAMPLE}_1.fastq -2 ${SAMPLE}_2.fastq -o $SAMPLE.mapsplice > mapsplice.log 2> mapsplice.log

	echo
	echo "3. Mapsplice"
	date '+%m/%d/%y %H:%M:%S'
	echo

	# not sure --fusion or --all-chromosomes-files is used with this version of Mapsplice
	# fusion is already run.
	if [ ! -e prep_2.fastq ]
	then
		python $MAPSPLICE_DIR/mapsplice_multi_thread.py \
			-1 prep_1.fastq \
			--chromosome-files-dir $REFERENCE_DIR/chromosomes \
			-o $SAMPLE.mapsplice \
			--Bowtieidx $REFERENCE_DIR/ebwt/humanchridx_M_rCRS \
			-Q fq \
			--threads $THREADS \
			--fusion --all-chromosomes-files $REFERENCE
	elif [ -e prep_2.fastq ]
	then
		python $MAPSPLICE_DIR/mapsplice_multi_thread.py \
			-1 prep_1.fastq -2 prep_2.fastq \
			--chromosome-files-dir $REFERENCE_DIR/$REFERENCE/chromosomes \
			-o $SAMPLE.mapsplice \
			--Bowtieidx $REFERENCE_DIR/$REFERENCE/bowtie_index/$REFERENCE\
			-Q fq \
			--threads $THREADS \
			--pairend \
			--fusion --all-chromosomes-files $REFERENCE_DIR/$REFERENCE/$REFERENCE.fa
	fi

	if [ $? -ne 0 ]
	then
		echo "ERROR: MapSplice failed"
		exit -1
	fi
	
	rm prep_1.fastq
	touch prep_1.fastq
	if [ -e prep_2.fastq ]; then
		rm prep_2.fastq
		touch prep_2.fastq
	fi
fi

##############################################################################
# 4. Add read groups
##############################################################################
if [ ! -e $SAMPLE.bam ]
then
	echo
	echo "4. Add read groups"
	date '+%m/%d/%y %H:%M:%S'
        echo

	## omitting these tags because i don't know them for my samples:
	java -Xmx4G -jar $PICARD_JAR AddOrReplaceReadGroups INPUT=$SAMPLE.mapsplice/alignments.bam OUTPUT=${SAMPLE}.bam RGSM=${SAMPLE} RGID=${SAMPLE} RGLB=TruSeq RGPL=illumina RGPU=barcode VALIDATION_STRINGENCY=SILENT  TMP_DIR=$TMP_DIR SORT_ORDER=coordinate CREATE_INDEX=true

	if [ $? -ne 0 ]
	then
		echo Error adding read group
		rm $SAMPLE.ba?
		exit -1
	fi
	
	rm $SAMPLE.mapsplice/alignments.bam
fi

##############################################################################
## omitting this step because it's stupid. and i aligned in phred33 to begin with
# 5. Convert back to phred33
##############################################################################
# java -Xmx512M -jar ubu.jar sam-convert --phred64to33 --in working/rg_alignments.bam .out working/phred33_alignments.bam > working/sam_convert.log 2> working/sam_convert.log

##############################################################################
# 6. Sort by coordinate (I think this step can be combined with adding read groups)
##############################################################################
# By commenting out Step 5, this is now done by Step 4.
#if [ ! -e $SAMPLE.genome.aln.bam ]
#then
#	ln -s $SAMPLE.rg.alignments.bam $SAMPLE.genome.aln.bam
#	ln -s $SAMPLE.rg.alignments.bai $SAMPLE.genome.aln.bai

#	echo "6. Sort by coordinate"
	## I'm gonna use picard instead of samtools
	# samtools sort ${SAMPLE}_rg_alignments.bam ${SAMPLE}.genome.aln
#	java -Xmx8g -jar -Djava.io.tmpdir=. $PICARD_JAR SortSam I=${SAMPLE}.rg.alignments.bam O=${SAMPLE}.genome.aln.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
#fi

##############################################################################
# 7. Flagstat
##############################################################################
if [ ! -e ${SAMPLE}.flagstat ]
then
	echo
	echo "7. Flagstat"
        date '+%m/%d/%y %H:%M:%S'
        echo

	$SAMTOOLS flagstat ${SAMPLE}.bam > ${SAMPLE}.flagstat
fi

##############################################################################
# 8. Index
##############################################################################
## don't need to do this because I've already indexed it with picard.
# samtools index ${SAMPLE}.genome.aln.bam

##############################################################################
# 9. Sort by chromosome, then read id
##############################################################################
if [ ! -e ${SAMPLE}.chromReadSorted.bam ]
then
	echo
	echo "Sort by chromosome, then read id"
        date '+%m/%d/%y %H:%M:%S'
        echo

	time perl $UBU_DIR/src/perl/sort_bam_by_reference_and_name.pl --input ${SAMPLE}.bam --output ${SAMPLE}.chromReadSorted.bam --temp-dir . --samtools $SAMTOOLS
fi

##############################################################################
# 10. Translate to transcriptome coords
##############################################################################
if [ ! -e ${SAMPLE}.transcriptome.bam ]
then
	echo
	echo "Translate to transcriptome coords"
        date '+%m/%d/%y %H:%M:%S'
        echo

	#time java -Xms8g -Xmx8g -jar $UBU_JAR sam-xlate --bed $REFERENCE_DIR/../unc_hg19.bed --in ${SAMPLE}.genome.chromReadSorted.bam --out ${SAMPLE}.transcriptome.bam --order $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse
	time java -Xms8g -Xmx8g -jar $UBU_JAR sam-xlate --bed $REFERENCE_DIR/$REFERENCE/annotation/tcga/unc_hg19.bed --in ${SAMPLE}.chromReadSorted.bam --out ${SAMPLE}.transcriptome.bam --order $REFERENCE_DIR/$REFERENCE/rsem_index/$REFERENCE.transcripts.fa --xgtags --reverse
fi

##############################################################################
# 11. Filter indels, large inserts, zero mapping quality from transcriptome bam
##############################################################################
if [ ! -e ${SAMPLE}.transcriptome.filtered.bam ]
then
	echo
	echo "Filter indels, large inserts, zero mapping quality from transcriptome bam"
        date '+%m/%d/%y %H:%M:%S'
        echo

	time java -Xmx1g -jar $UBU_JAR sam-filter --in ${SAMPLE}.transcriptome.bam --out ${SAMPLE}.transcriptome.filtered.bam --strip-indels --max-insert 10000 --mapq 1
fi

##############################################################################
# 11b. Strip read names of /1 and /2 
##############################################################################
if [ ! -e ${SAMPLE}.transcriptome.filtered.stripped.bam ]
then
	echo
	echo Stripping /1 and /2 from read names
        date '+%m/%d/%y %H:%M:%S'
        echo

	java -jar $RGTOOLS_DIR/StripPairInfoFromReadName.jar I=${SAMPLE}.transcriptome.filtered.bam O=${SAMPLE}.transcriptome.filtered.stripped.bam S="/1" S="/2"
fi

##############################################################################
# 12. RSEM
##############################################################################
if [ ! -e ${SAMPLE}.rsem.genes.results ]
then
	echo
	echo "RSEM"
        date '+%m/%d/%y %H:%M:%S'
        echo
	
	## i don't know why they have a --gcr-output-file flag in their command, but it is not a valid rsem parameter so I'm omitting it
	if [ ! -e prep_2.fastq ]
        then
		$RSEM_DIR/rsem-calculate-expression --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.filtered.stripped.bam $REFERENCE_DIR/rsem_index/$REFERENCE $SAMPLE.rsem
	else
		$RSEM_DIR/rsem-calculate-expression --paired-end --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.filtered.stripped.bam $REFERENCE_DIR/$REFERENCE/rsem_index/$REFERENCE $SAMPLE.rsem
	fi

	if [ $? -ne 0 ]
	then
		echo "Error running RSEM"
		exit -1
	fi

	echo Completed RSEM
	date '+%m/%d/%y %H:%M:%S'
        echo

	## example cmd
	# rsem-calculate-expression --gcr-output-file --paired-end --bam --estimate-rspd -p 8 working/transcriptome_alignments_filtered.bam /datastore/tier1data/nextgenseq/seqware-analysis/mapsplice_rsem/rsem_ref/hg19_M_rCRS_ref rsem > working/rsem.log 2> working/rsem.log

	# 13. Strip trailing tabs from rsem.isoforms.results
#	echo
#	echo "Strip trailing tabs from rsem.isoforms.results"
	#perl $UBU_DIR/src/perl/strip_trailing_tabs.pl --input $SAMPLE.rsem.isoforms.results --temp $SAMPLE.rsem.orig.isoforms.results 
#	mv $SAMPLE.rsem.isoforms.results $SAMPLE.rsem.orig.isoforms.results
#	sed 's/\\t\$//g' $SAMPLE.rsem.orig.isoforms.results > $SAMPLE.rsem.isoforms.results

	# 14. Prune isoforms from gene quant file
#	echo
#	echo "Prune isoforms from gene quant file"
#	mv $SAMPLE.rsem.genes.results $SAMPLE.rsem.orig.genes.results
#	sed /^uc0/d $SAMPLE.rsem.orig.genes.results > $SAMPLE.rsem.genes.results
fi

# 15. Normalize gene quant
#if [ ! -e $SAMPLE.rsem.genes.normalized_results ]
#then
#	echo
#	echo "Normalize gene quant"
#	perl $UBU_DIR/src/perl/quartile_norm.pl -c 2 -q 75 -t 1000 -o $SAMPLE.rsem.genes.normalized_results $SAMPLE.rsem.genes.results
#fi

# 16. Normalize isoform quant
#if [ ! -e $SAMPLE.rsem.isoforms.normalized_results ]
#then
#	echo
#	echo "Normalize isoform quant"
#	perl $UBU_DIR/src/perl/quartile_norm.pl -c 2 -q 75 -t 300 -o $SAMPLE.rsem.isoforms.normalized_results $SAMPLE.rsem.isoforms.results
#fi

##############################################################################
# 17. Junction counts
##############################################################################
if [ ! -e $SAMPLE.junction_quantification.txt ]
then
	echo
	echo Counting splice junctions
	java -Xmx512M -jar $UBU_JAR sam-junc --junctions $REFERENCE_DIR/$REFERENCE/annotation/tcga/splice_junctions.txt --in $SAMPLE.bam --out $SAMPLE.junction_quantification.txt 
	
	if [ $? -ne 0 ]
	then
		echo Error counting splice junctions
		rm SAMPLE.junction_quantification.txt
		exit -1
	fi
fi

##############################################################################
# 18. Exon counts
##############################################################################
if [ ! -e $SAMPLE.exon_quantification.txt ]
then
	echo
	echo Counting exons
	$BEDTOOLS/coverageBed -split -abam $SAMPLE.bam -b $REFERENCE_DIR/$REFERENCE/annotation/tcga/composite_exons.bed | normalizeBedToolsExonQuant.pl $SAMPLE.bam $REFERENCE_DIR/$REFERENCE/annotation/tcga/composite_exons.bed > $SAMPLE.exon_quantification.txt 

	if [ $? -ne 0 ]
	then
		echo Error counting exons
		rm $SAMPLE.exon_quantification.txt
		exit -1	
	fi
fi

##############################################################################
# 19. Cleanup large intermediate output
##############################################################################
if [ $AWS -eq 1 ]
then
	rm *.gz
fi

# Remove temp dir created
rm -rf *_`uname -n`_*
# Removing unused local files 
rm *.fastq *.chromReadSorted.bam *.rsem.transcript.bam *.rsem.transcript.sorted.bam* *.transcriptome.bam *.transcriptome.filtered.bam *.transcriptome.filtered.stripped.bam

##############################################################################
# 13.  Transfer data to final resting place
##############################################################################
cd ..
if [ $AWS == 0 ]
then
	echo "Moving ${SAMPLE_DIR} to ${OUTDIR}/${SAMPLE}"
	mv ${SAMPLE_DIR} ${OUTDIR}/${SAMPLE}
else
	aws s3 cp --recursive $SAMPLE_DIR ${BUCKET}/${SAMPLE}

	if [ $? -ne 0 ]
	then
		echo Error copying $SAMPLE_DIR to s3
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