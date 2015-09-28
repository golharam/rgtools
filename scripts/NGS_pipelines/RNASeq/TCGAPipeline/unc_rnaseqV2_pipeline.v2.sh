#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

# Pipeline is from: https://raw.githubusercontent.com/cc2qe/sandbox/master/unc_rnaseqV2_pipeline.sh
# RG Version 0.03
HELP=0
AWS=0
OUTDIR=`pwd`
RUN_FASTQC=0
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
                --fastqc)
                RUN_FASTQC=1
		;;

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

		-h|--help)
		HELP=1
		;;

                -o|--outdir)
                OUTDIR="$2"
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
	echo "Usage 1: $0 <options> [-s|--sample <sample name>] [-1|--fastq1 <path/to/fastq1>] [-2|--fastq2 <path/to/fastq2>] [-t|--threads <threads to use>]"
	echo "Usage 2: $0 <options> [-s|--sample <SRA ID>] [--sraftp <ftp location>] [-t|--threads <threads to use>]" 
	echo "Options:"
	echo "	-a|--aws"
	echo "  --delete-intermediate (delete intermediate files, not including source fq.gz) (not yet implemented)"
	echo "	-h|--help"
	echo "  -o|--outdir <output directory> [default=current working directory]"
	echo "  --fastqc (default: do not run fastqc)"
	exit 0
fi

if [ $AWS == 0 ]; then
	EXT_PKGS_DIR=/apps/sys/galaxy/external_packages
	REFERENCE_DIR=/ngs/ngs15/golharr/NGS/mRNAseq_TCGA/hg19_M_rCRS
else
	EXT_PKGS_DIR=/ngs/apps
	REFERENCE_DIR=/ngs/reference/mRNAseq_TCGA/hg19_M_rCRS
fi

if [ -z "$TMP_DIR" ]; then
        TMP_DIR=/scratch
fi

BEDTOOLS=$EXT_PKGS_DIR/bedtools-2.23.0/bin
FASTQC=$EXT_PKGS_DIR/FastQC-0.11.2/fastqc
MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.129/picard.jar
REFERENCE=$REFERENCE_DIR/hg19_M_rCRS.fa
RGTOOLS_DIR=$EXT_PKGS_DIR/rgtools
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
RSEM_DIR=$EXT_PKGS_DIR/RSEM-1.2.20
SAMTOOLS=$EXT_PKGS_DIR/samtools-0.1.19/samtools
SRA_DIR=$EXT_PKGS_DIR/sratoolkit-2.5.2/bin
UBU_DIR=$EXT_PKGS_DIR/ubu
UBU_JAR=$EXT_PKGS_DIR/ubu-1.2-jar-with-dependencies.jar

# Start Analysis
echo "Processing $SAMPLE"
date '+%m/%d/%y %H:%M:%S'
echo "Pipeline: $VERSION"
echo "Subsample: $SUBSAMPLE"
echo "TmpDir: $TMP_DIR"
echo "AWS: $AWS"
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

# If SRA, download from SRA
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

echo FASTQ1: $FASTQ1
if [ -n "$FASTQ2" ]
then
	echo FASTQ2: $FASTQ2
fi
date '+%m/%d/%y %H:%M:%S'
analysis_date_started=$(date +"%s")
echo

# 0a. Download from kraken
if [ $AWS == 1 ]
then
	echo "Running on AWS"
	BASEFQ1=`basename $FASTQ1`
	BASEFQ2=`basename $FASTQ2`
	if [ ! -e $BASEFQ1 ]
	then
		echo "Downloading $BASEFQ1"
		date '+%m/%d/%y %H:%M:%S'
		echo

		#S3LOC=`aws s3 ls --recursive s3://bmsrd-ngs/ngs/ngs18/ | grep $FQ1 | cut -d ' ' -f 4`
		#S3LOC="s3://bmsrd-ngs/$S3LOC"
		#echo aws s3 cp $S3LOC .
		#aws s3 cp $S3LOC .
		scp -q golharr@kraken.pri.bms.com:$FASTQ1 .

		if [ $? -ne 0 ]; then
			echo Error downloading $BASEFQ1
			exit -1
		fi
	fi

	if [ ! -e $BASEFQ2 ]
	then
		echo
		echo "Downloading $BASEFQ2"
		date '+%m/%d/%y %H:%M:%S'
		echo

		#S3LOC=`aws s3 ls --recursive s3://bmsrd-ngs/ngs/ngs18/ | grep $FQ2 | cut -d ' ' -f 4`
		#S3LOC="s3://bmsrd-ngs/$S3LOC"
		#echo aws s3 cp $S3LOC .
		#aws s3 cp $S3LOC .
		scp -q golharr@kraken.pri.bms.com:$FASTQ2 .

		if [ $? -ne 0 ]; then
			echo Error downloading $BASEFQ2
		fi
	fi
	FASTQ1=$BASEFQ1
	FASTQ2=$BASEFQ2
fi

# 0b. Unzip the fastqs
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

# Run FastQC
if [ "$RUN_FASTQC" -eq 1 ] 
then
	echo "Running FastQC"
	date '+%m/%d/%y %H:%M:%S'
	echo

	if [ -e ${SAMPLE}_2.fastq ]
	then
		$FASTQC -t 2 -o . ${SAMPLE}_1.fastq ${SAMPLE}_2.fastq
	else
		$FASTQC -o . ${SAMPLE}_1.fastq
	fi

	echo
fi

# 1. Format fastq 1 for Mapsplice
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

# 2. Format fastq 2 for Mapsplice
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

# 3. Mapsplice
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

	# Not sure the last line of parameters are used
	if [ ! -e prep_2.fastq ]
	then
		python $MAPSPLICE_DIR/mapsplice_multi_thread.py \
			-1 prep_1.fastq \
			--chromosome-files-dir $REFERENCE_DIR/chromosomes \
			-o $SAMPLE.mapsplice \
			--Bowtieidx $REFERENCE_DIR/ebwt/humanchridx_M_rCRS \
			-Q fq \
			-X $THREADS \
			--fusion --all-chromosomes-files $REFERENCE -p $THREADS
	elif [ -e prep_2.fastq ]
	then
		python $MAPSPLICE_DIR/mapsplice_multi_thread.py \
			-1 prep_1.fastq -2 prep_2.fastq \
			--chromosome-files-dir $REFERENCE_DIR/chromosomes \
			-o $SAMPLE.mapsplice \
			--Bowtieidx $REFERENCE_DIR/ebwt/humanchridx_M_rCRS \
			-Q fq \
			-X $THREADS \
			--pairend \
			--fusion --all-chromosomes-files $REFERENCE -p $THREADS
	fi

	if [ $? -ne 0 ]
	then
		echo "ERROR: MapSplice failed"
		exit -1
	fi
fi

# 4. Add read groups
if [ ! -e $SAMPLE.rg.alignments.bam ]
then
	echo
	echo "4. Add read groups"
	date '+%m/%d/%y %H:%M:%S'
        echo

	## omitting these tags because i don't know them for my samples:
	java -Xmx4G -jar $PICARD_JAR AddOrReplaceReadGroups INPUT=$SAMPLE.mapsplice/alignments.bam OUTPUT=${SAMPLE}.rg.alignments.bam RGSM=${SAMPLE} RGID=${SAMPLE} RGLB=TruSeq RGPL=illumina RGPU=barcode VALIDATION_STRINGENCY=SILENT  TMP_DIR=$TMP_DIR SORT_ORDER=coordinate CREATE_INDEX=true

	if [ $? -ne 0 ]
	then
		echo Error adding read group
		rm $SAMPLE.rg.alignments.bam
		exit -1
	fi
fi

## omitting this step because it's stupid. and i aligned in phred33 to begin with
# 5. Convert back to phred33
# java -Xmx512M -jar ubu.jar sam-convert --phred64to33 --in working/rg_alignments.bam .out working/phred33_alignments.bam > working/sam_convert.log 2> working/sam_convert.log

# 6. Sort by coordinate (I think this step can be combined with adding read groups)
if [ ! -e $SAMPLE.genome.aln.bam ]
then
	ln -s $SAMPLE.rg.alignments.bam $SAMPLE.genome.aln.bam
	ln -s $SAMPLE.rg.alignments.bai $SAMPLE.genome.aln.bai

#	echo "6. Sort by coordinate"
	## I'm gonna use picard instead of samtools
	# samtools sort ${SAMPLE}_rg_alignments.bam ${SAMPLE}.genome.aln
#	java -Xmx8g -jar -Djava.io.tmpdir=. $PICARD_JAR SortSam I=${SAMPLE}.rg.alignments.bam O=${SAMPLE}.genome.aln.bam SO=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
fi

# 7. Flagstat
if [ ! -e ${SAMPLE}.genome.aln.bam.flagstat ]
then
	echo
	echo "7. Flagstat"
        date '+%m/%d/%y %H:%M:%S'
        echo

	$SAMTOOLS flagstat ${SAMPLE}.genome.aln.bam > ${SAMPLE}.genome.aln.bam.flagstat
fi

# 8. Index
## don't need to do this because I've already indexed it with picard.
# samtools index ${SAMPLE}.genome.aln.bam

# 9. Sort by chromosome, then read id
if [ ! -e ${SAMPLE}.alignments.chromReadSorted.bam ]
then
	echo
	echo "Sort by chromosome, then read id"
        date '+%m/%d/%y %H:%M:%S'
        echo

	time perl $UBU_DIR/src/perl/sort_bam_by_reference_and_name.pl --input ${SAMPLE}.genome.aln.bam --output ${SAMPLE}.alignments.chromReadSorted.bam --temp-dir $TMP_DIR --samtools $SAMTOOLS
fi

# 10. Translate to transcriptome coords
if [ ! -e ${SAMPLE}.transcriptome.aln.bam ]
then
	echo
	echo "Translate to transcriptome coords"
        date '+%m/%d/%y %H:%M:%S'
        echo

	time java -Xms8g -Xmx8g -jar $UBU_JAR sam-xlate --bed $REFERENCE_DIR/../unc_hg19.bed --in ${SAMPLE}.alignments.chromReadSorted.bam --out ${SAMPLE}.transcriptome.aln.bam --order $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse
fi


# 11. Filter indels, large inserts, zero mapping quality from transcriptome bam
if [ ! -e ${SAMPLE}.transcriptome.aln.filtered.bam ]
then
	echo
	echo "Filter indels, large inserts, zero mapping quality from transcriptome bam"
        date '+%m/%d/%y %H:%M:%S'
        echo

	time java -Xmx1g -jar $UBU_JAR sam-filter --in ${SAMPLE}.transcriptome.aln.bam --out ${SAMPLE}.transcriptome.aln.filtered.bam --strip-indels --max-insert 10000 --mapq 1
fi

# Strip read names of /1 and /2 
if [ ! -e ${SAMPLE}.transcriptome.aln.filtered.stripped.bam ]
then
	echo
	echo Stripping /1 and /2 from read names
        date '+%m/%d/%y %H:%M:%S'
        echo

	java -jar $RGTOOLS_DIR/StripPairInfoFromReadName.jar I=${SAMPLE}.transcriptome.aln.filtered.bam O=${SAMPLE}.transcriptome.aln.filtered.stripped.bam S="/1" S="/2"
fi

# 12. RSEM
if [ ! -e ${SAMPLE}.rsem.genes.results ]
then
	echo
	echo "RSEM"
        date '+%m/%d/%y %H:%M:%S'
        echo
	
	## i don't know why they have a --gcr-output-file flag in their command, but it is not a valid rsem parameter so I'm omitting it
	if [ ! -e prep_2.fastq ]
        then
		$RSEM_DIR/rsem-calculate-expression --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.aln.filtered.stripped.bam $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref $SAMPLE.rsem
	else
		$RSEM_DIR/rsem-calculate-expression --paired-end --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.aln.filtered.stripped.bam $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref $SAMPLE.rsem
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

# 17. Junction counts
#if [ ! -e $SAMPLE.junction_quantification.txt ]
#then
#	echo
#	echo Counting splice junctions
#	java -Xmx512M -jar $UBU_JAR sam-junc --junctions $REFERENCE_DIR/../splice_junctions.txt --in $SAMPLE.genome.aln.bam --out $SAMPLE.junction_quantification.txt 
	
#	if [ $? -ne 0 ]
#	then
#		echo Error counting splice junctions
#		rm SAMPLE.junction_quantification.txt
#		exit -1
#	fi
#fi

# 18. Exon counts
#if [ ! -e $SAMPLE.bt.exon_quantification.txt ]
#then
#	echo
#	echo Counting exons
#	$BEDTOOLS/coverageBed -split -abam $SAMPLE.genome.aln.bam -b $REFERENCE_DIR/../composite_exons.bed | perl $REFERENCE_DIR/../normalizeBedToolsExonQuant.pl $SAMPLE.genome.aln.bam $REFERENCE_DIR/../composite_exons.bed > $SAMPLE.bt.exon_quantification.txt 

#	if [ $? -ne 0 ]
#	then
#		echo Error counting exons
#		rm $SAMPLE.bt.exon_quantification.txt
#		exit -1	
#	fi
#fi

# 19. Cleanup large intermediate output
if [ $AWS == 1 ]
then
	rm $FASTQ1 $FASTQ2
fi
rm -rf *.fastq 

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

