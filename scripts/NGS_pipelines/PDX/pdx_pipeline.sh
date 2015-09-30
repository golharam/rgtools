#!/bin/bash -e

#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -pe orte 8

VERSION=0.1

external_packages=/ngs/apps
reference_genomes=/ngs/reference

bwaApp="$external_packages/bwa-0.7.12/bwa"
bowtie2App="$external_packages/bowtie2-2.2.6/bowtie2"
javaApp="$external_packages/jdk1.6.0_45/bin/java"
picardJar="$external_packages/picard-tools-1.129/picard.jar"
samtoolsApp="$external_packages/samtools-0.1.19/samtools"
countReadPairsAppsJar="/home/golharr/workspace/rgtools/dist/CountMappedReadPairs.jar" 

graftGenome="hg19"
graftReference="$reference_genomes/$graftGenome/$graftGenome.fa"
graftbwaReference="$reference_genomes/$graftGenome/bwa_index/$graftGenome"
graftbowtie2Reference="$reference_genomes/$graftGenome/bowtie2_index/$graftGenome"

hostGenome="mm10"
hostReference="$reference_genomes/$hostGenome/$hostGenome.fa"
hostbwaReference="$reference_genomes/$hostGenome/bwa_index/$hostGenome"
hostbowtie2Reference="$reference_genomes/$hostGenome/bowtie2_index/$hostGenome"

graftSamFile="analysis/pdx/$sample.$graftGenome.sam"
graftMapped1="analysis/pdx/$sample.$graftGenome.aligned.fastq.1.gz"
graftMapped2="analysis/pdx/$sample.$graftGenome.aligned.fastq.2.gz"
graftUnmapped1="analysis/pdx/$sample.$graftGenome.unaligned.fastq.1.gz"
graftUnmapped2="analysis/pdx/$sample.$graftGenome.unaligned.fastq.2.gz"

graftAlignedHostSamFile="analysis/pdx/$sample.graftAligned.$hostGenome.sam"
graftUnalignedHostSamFile="analysis/pdx/$sample.graftUnaligned.$hostGenome.sam"

echo "PDX Pipeline $VERSION"
echo "Running $SAMPLE"
echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"

cd /scratch

aws s3 cp $FASTQ1 .
FASTQ1=`basename $FASTQ1`

aws s3 cp $FASTQ2 .
FASTQ2=`basename $FASTQ2`

# Map sample to hg19 (graft) to get graft mapped and unmapped
$bowtie2App --no-mixed --un-conc-gz analysis/pdx/$SAMPLE.$graftGenome.unaligned.fastq.gz --al-conc-gz analysis/pdx/$SAMPLE.$graftGenome.aligned.fastq.gz -p 8 -1 $FASTQ1 -2 $FASTQ2 \
		--no-unal --rg-id $sample --rg 'SM:$SAMPLE\tLB:$SAMPLE\tPL:illumina' -S $graftSamFile -x $graftbowtie2Reference
}

# Map graft mapped to mm10 to get graft.mapped.host.[mapped/unmapped]
$bowtie2App --no-mixed --un-conc-gz analysis/pdx/$sample.graftAlignedHostUnaligned.fastq.gz --al-conc-gz analysis/pdx/$sample.graftAlignedHostAligned.fastq.gz \
		-p 8 -1 $graftMapped_1 -2 $graftMapped_2 \
		--no-unal --rg-id $sample --rg 'SM:$sample\tLB:$sample\tPL:illumina' -S $graftAlignedHostSamFile -x $hostbowtie2Reference

# Map graft unmapped to mm10 to get graft.unmapped.host.[mapped/unmapped]
$bowtie2App --no-mixed --un-conc-gz analysis/pdx/$sample.graftUnalignedHostUnaligned.fastq.gz --al-conc-gz analysis/pdx/$sample.graftUnalignedHostAligned.fastq.gz \
		-p 8 -1 $graftUnmapped_1 -2 $graftUnmapped_2 \
		--no-unal --rg-id $sample --rg 'SM:$sample\tLB:$sample\tPL:illumina' -S $graftUnalignedHostSamFile -x $hostbowtie2Reference
