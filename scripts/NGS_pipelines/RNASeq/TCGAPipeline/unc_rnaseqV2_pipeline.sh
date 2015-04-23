#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

# From: https://raw.githubusercontent.com/cc2qe/sandbox/master/unc_rnaseqV2_pipeline.sh
if [ -z $SAMPLE ] || [ -z $FASTQ1 ] || [ -z $FASTQ2 ] || [ -z $THREADS ] || [ -z $AWS ]
then
	if [ $# -lt 5 ]
	then
		echo usage $0 [SAMPLE] [FASTQ1] [FASTQ2] [THREADS] [AWS]
		exit 1
	else
		SAMPLE=$1
		FASTQ1=$2
		FASTQ2=$3
		THREADS=$4
		AWS=$5
	fi
fi

if [ $AWS == 0 ]; then
	EXT_PKGS_DIR=/apps/sys/galaxy/external_packages
	REFERENCE_DIR=/ngs/ngs15/golharr/NGS/mRNAseq_TCGA/hg19_M_rCRS
	TMP_DIR=/dev/shm
	PROJECT_DIR=`pwd`
else
	EXT_PKGS_DIR=/ngs/apps
	REFERENCE_DIR=/ngs/reference/mRNAseq_TCGA/hg19_M_rCRS
	TMP_DIR=/scratch
	PROJECT_DIR=/ng20/golharr/M2GEN-TCGA
fi

MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multithreads_12_07/bin
#MAPSPLICE_DIR=$EXT_PKGS_DIR/MapSplice_multi_threads_2.0.1.9/bin
UBU_DIR=$EXT_PKGS_DIR/ubu
UBU_JAR=$EXT_PKGS_DIR/ubu-1.2-jar-with-dependencies.jar
PICARD_JAR=$EXT_PKGS_DIR/picard-tools-1.129/picard.jar
SAMTOOLS=$EXT_PKGS_DIR/samtools-0.1.19/samtools
REFERENCE=$REFERENCE_DIR/hg19_M_rCRS.fa
BEDTOOLS=$EXT_PKGS_DIR/bedtools-2.23.0/bin
#RSEM_DIR=$EXT_PKGS_DIR/rsem-1.1.13
RSEM_DIR=$EXT_PKGS_DIR/rsem-1.2.19

echo "Processing $SAMPLE from seq files:"
echo $FASTQ1
echo $FASTQ2
date '+%m/%d/%y %H:%M:%S'
analysis_date_started=$(date +"%s")
echo

cd $TMP_DIR
rm -rf $SAMPLE

if [ ! -d $SAMPLE ]
then
	mkdir $SAMPLE
fi
cd $SAMPLE
echo "Working in `uname -n`:`pwd`"
echo

# 0a. Download from S3
if [ $AWS == 1 ]
then
	echo "Running on AWS"
	FQ1=`basename $FASTQ1`
	FQ2=`basename $FASTQ2`
	if [ ! -e $FQ1 ]
	then
		echo "Downloading $FQ1"
		date '+%m/%d/%y %H:%M:%S'
		echo

		S3LOC=`aws s3 ls --recursive s3://bmsrd-ngs/ngs/ngs18/ | grep $FQ1 | cut -d ' ' -f 4`
		S3LOC="s3://bmsrd-ngs/$S3LOC"
		echo aws s3 cp $S3LOC .
		aws s3 cp $S3LOC .

		if [ $? -ne 0 ]; then
			echo Error downloading $S3LOC from S3
			exit -1
		fi
	fi

	if [ ! -e $FQ2 ]
	then
		echo
		echo "Downloading $FQ2"
		date '+%m/%d/%y %H:%M:%S'
		echo

		S3LOC=`aws s3 ls --recursive s3://bmsrd-ngs/ngs/ngs18/ | grep $FQ2 | cut -d ' ' -f 4`
		S3LOC="s3://bmsrd-ngs/$S3LOC"
		echo aws s3 cp $S3LOC .
		aws s3 cp $S3LOC .

                if [ $? -ne 0 ]; then
                        echo Error downloading $S3LOC from S3
                        exit -1
                fi
	fi
	FASTQ1=$FQ1
	FASTQ2=$FQ2
fi

# 0b. Unzip the fastqs
if [ ! -e ${SAMPLE}_1.fastq ]; then
	echo
	echo "Unzipping $FASTQ1 > ${SAMPLE}_1.fastq"
	date '+%m/%d/%y %H:%M:%S'
	echo

	gunzip -dc $FASTQ1 > ${SAMPLE}_1.fastq
                                
	if [ $? -ne 0 ]; then
		echo "Error Unzipping $FASTQ1 > ${SAMPLE}_1.fastq"
		rm ${SAMPLE}_1.fastq
		exit -1
	fi
fi

if [ ! -e ${SAMPLE}_2.fastq ]; then
        echo "Unzipping $FASTQ2 > ${SAMPLE}_2.fastq"
	date '+%m/%d/%y %H:%M:%S'
	echo

        gunzip -dc $FASTQ2 > ${SAMPLE}_2.fastq

        if [ $? -ne 0 ]; then
                echo "Error Unzipping $FASTQ2 > ${SAMPLE}_2.fastq"
                rm ${SAMPLE}_2.fastq
                exit -1
        fi
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
fi

# 2. Format fastq 2 for Mapsplice
if [ ! -e prep_2.fastq ] 
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

	python $MAPSPLICE_DIR/mapsplice_multi_thread.py --fusion --all-chromosomes-files $REFERENCE --pairend -p $THREADS -X 8 -Q fq --chromosome-files-dir $REFERENCE_DIR/chromosomes --Bowtieidx $REFERENCE_DIR/ebwt/humanchridx_M_rCRS -1 prep_1.fastq -2 prep_2.fastq -o $SAMPLE.mapsplice

	if [ $? -ne 0 ]
	then
		echo "ERROR: MapSplice failed"
		exit -1
	fi
fi

# 4. Add read groups
if [ ! -e $SAMPLE.rg.alignments.bam ]
then
	echo "4. Add read groups"
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
	echo "7. Flagstat"
	$SAMTOOLS flagstat ${SAMPLE}.genome.aln.bam > ${SAMPLE}.genome.aln.bam.flagstat
fi

# 8. Index
## don't need to do this because I've already indexed it with picard.
# samtools index ${SAMPLE}.genome.aln.bam

# 9. Sort by chromosome, then read id
if [ ! -e ${SAMPLE}.alignments.chromReadSorted.bam ]
then
	echo "Sort by chromosome, then read id"
	time perl $UBU_DIR/src/perl/sort_bam_by_reference_and_name.pl --input ${SAMPLE}.genome.aln.bam --output ${SAMPLE}.alignments.chromReadSorted.bam --temp-dir $TMP_DIR --samtools $SAMTOOLS
fi

# 10. Translate to transcriptome coords
if [ ! -e ${SAMPLE}.transcriptome.aln.bam ]
then
	echo "Translate to transcriptome coords"
	time java -Xms8g -Xmx8g -jar $UBU_JAR sam-xlate --bed $REFERENCE_DIR/../unc_hg19.bed --in ${SAMPLE}.alignments.chromReadSorted.bam --out ${SAMPLE}.transcriptome.aln.bam --order $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse
fi


# 11. Filter indels, large inserts, zero mapping quality from transcriptome bam
if [ ! -e ${SAMPLE}.transcriptome.aln.filtered.bam ]
then
	echo "Filter indels, large inserts, zero mapping quality from transcriptome bam"
	time java -Xmx1g -jar $UBU_JAR sam-filter --in ${SAMPLE}.transcriptome.aln.bam --out ${SAMPLE}.transcriptome.aln.filtered.bam --strip-indels --max-insert 10000 --mapq 1
fi

# 12. RSEM
if [ ! -e ${SAMPLE}.rsem.genes.results ]
then
	echo "RSEM"
	## i don't know why they have a --gcr-output-file flag in their command, but it is not a valid rsem parameter so I'm omitting it
	time $RSEM_DIR/rsem-calculate-expression --paired-end --bam --estimate-rspd -p $THREADS  ${SAMPLE}.transcriptome.aln.filtered.bam $REFERENCE_DIR/rsem_ref/hg19_M_rCRS_ref $SAMPLE.rsem
	## example cmd
	# rsem-calculate-expression --gcr-output-file --paired-end --bam --estimate-rspd -p 8 working/transcriptome_alignments_filtered.bam /datastore/tier1data/nextgenseq/seqware-analysis/mapsplice_rsem/rsem_ref/hg19_M_rCRS_ref rsem > working/rsem.log 2> working/rsem.log

	# 13. Strip trailing tabs from rsem.isoforms.results
	echo "Strip trailing tabs from rsem.isoforms.results"
	#perl $UBU_DIR/src/perl/strip_trailing_tabs.pl --input $SAMPLE.rsem.isoforms.results --temp $SAMPLE.rsem.orig.isoforms.results 
	mv $SAMPLE.rsem.isoforms.results $SAMPLE.rsem.orig.isoforms.results
	sed 's/\\t\$//g' $SAMPLE.rsem.orig.isoforms.results > $SAMPLE.rsem.isoforms.results

	# 14. Prune isoforms from gene quant file
	echo "Prune isoforms from gene quant file"
	mv $SAMPLE.rsem.genes.results $SAMPLE.rsem.orig.genes.results
	sed /^uc0/d $SAMPLE.rsem.orig.genes.results > $SAMPLE.rsem.genes.results
fi

# 15. Normalize gene quant
if [ ! -e $SAMPLE.rsem.genes.normalized_results ]
then
	echo "Normalize gene quant"
	perl $UBU_DIR/src/perl/quartile_norm.pl -c 2 -q 75 -t 1000 -o $SAMPLE.rsem.genes.normalized_results $SAMPLE.rsem.genes.results
fi

# 16. Normalize isoform quant
if [ ! -e $SAMPLE.rsem.isoforms.normalized_results ]
then
	echo "Normalize isoform quant"
	perl $UBU_DIR/src/perl/quartile_norm.pl -c 2 -q 75 -t 300 -o $SAMPLE.rsem.isoforms.normalized_results $SAMPLE.rsem.isoforms.results
fi

# 17. Junction counts
if [ ! -e $SAMPLE.junction_quantification.txt ]
then
	echo Counting splice junctions
	java -Xmx512M -jar $UBU_JAR sam-junc --junctions $REFERENCE_DIR/../splice_junctions.txt --in $SAMPLE.genome.aln.bam --out $SAMPLE.junction_quantification.txt 
	
	if [ $? -ne 0 ]
	then
		echo Error counting splice junctions
		rm SAMPLE.junction_quantification.txt
		exit -1
	fi
fi

# 18. Exon counts
if [ ! -e $SAMPLE.bt.exon_quantification.txt ]
then
	echo Counting exons
	$BEDTOOLS/coverageBed -split -abam $SAMPLE.genome.aln.bam -b $REFERENCE_DIR/../composite_exons.bed | perl $REFERENCE_DIR/../normalizeBedToolsExonQuant.pl $SAMPLE.genome.aln.bam $REFERENCE_DIR/../composite_exons.bed > $SAMPLE.bt.exon_quantification.txt 

	if [ $? -ne 0 ]
	then
		echo Error counting exons
		rm $SAMPLE.bt.exon_quantification.txt
		exit -1	
	fi
fi

# 19. Cleanup large intermediate output
if [ $AWS == 1 ]
then
	rm $FASTQ1 $FASTQ2
fi
#rm -rf *.fastq *.bam *.bai $SAMPLE.mapsplice/alignments.bam $SAMPLE.mapsplice/logs 
rm -rf *.fastq 

cd ..
if [ $AWS == 0 ]
then
	mv $SAMPLE ${PROJECT_DIR}/${SAMPLE}.local
else
	scp $SAMPLE golharr@kraken.pri.bms.com:$PROJECT_DIR
fi

echo
echo $SAMPLE Complete
date '+%m/%d/%y %H:%M:%S'
analysis_date_finished=$(date +"%s")
diff=$(($analysis_date_finished-$analysis_date_started))
echo "Total analysis took $(($diff / 60)) minutes and $(($diff % 60)) seconds."
