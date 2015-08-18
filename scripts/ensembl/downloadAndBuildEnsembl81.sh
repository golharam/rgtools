#!/bin/bash

ENSEMBL_RELEASE_NUM=81
ENSEMBL_RELEASE="release-$ENSEMBL_RELEASE_NUM"

SPECIES="Homo_sapiens"
ASSEMBLY="GRCh38"

#SPECIES="Mus_musculus"
#ASSEMBLY="GRCm38"

#SPECIES="Rattus_norvegicus"
#ASSEMBLY="Rnor_6.0"

LC_SPECIES=`echo $SPECIES | awk '{print tolower($0)}'`

if [ ! -e README ]; then
	echo "Downloading README"
	wget ftp://ftp.ensembl.org/pub/${ENSEMBL_RELEASE}/fasta/${LC_SPECIES}/dna/README
fi

if [ ! -e ${SPECIES}.${ASSEMBLY}.dna.primary_assembly.fa.gz ]; then
	echo "Downloading genome"
	wget ftp://ftp.ensembl.org/pub/${ENSEMBL_RELEASE}/fasta/${LC_SPECIES}/dna/${SPECIES}.${ASSEMBLY}.dna.primary_assembly.fa.gz
fi

if [ ! -e ${ASSEMBLY}.fa ]; then
	echo "Uncompressing genome"
	gunzip -dc ${SPECIES}.${ASSEMBLY}.dna.primary_assembly.fa.gz > ${ASSEMBLY}.fa
fi

if [ ! -e ${ASSEMBLY}.dict ]; then
	echo "Building genome dictionary"
	/apps/sys/galaxy/external_packages/bin/java \
		-jar /apps/sys/galaxy/external_packages/picard-tools-1.137/picard.jar \
		CreateSequenceDictionary \
		R=${ASSEMBLY}.fa O=${ASSEMBLY}.dict GENOME_ASSEMBLY=${ASSEMBLY}
fi

if [ ! -e ${ASSEMBLY}.fa.fai ]; then
	echo "Build genome index"
	/apps/sys/galaxy/external_packages/samtools-0.1.19/samtools faidx ${ASSEMBLY}.fa
fi

if [ ! -d bwa_index ]; then
	echo "Building BWA Index"
	mkdir bwa_index
	cd bwa_index
	echo "/apps/sys/galaxy/external_packages/bwa-0.7.12/bwa index -p ${ASSEMBLY} ../${ASSEMBLY}.fa"| qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bwa
	cd ..
fi

if [ ! -d bowtie_index ]; then
	echo "Building bowtie index"
	mkdir bowtie_index
	cd bowtie_index
	echo "/apps/sys/galaxy/external_packages/bowtie-1.1.1/bowtie-build ../${ASSEMBLY}.fa ./${ASSEMBLY}" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bowtie
	cd ..
fi

if [ ! -d bowtie2_index ]; then
        echo "Building bowtie2 index"
        mkdir bowtie2_index
        cd bowtie2_index
        echo "/apps/sys/galaxy/external_packages/bowtie2-2.2.6/bowtie2-build ../${ASSEMBLY}.fa ./${ASSEMBLY}" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bowtie2
        cd ..
fi

if [ ! -d annotation/ensembl-${ENSEMBL_RELEASE} ]; then
	echo "Downloading annotation"
	mkdir -p annotation/ensembl-${ENSEMBL_RELEASE}
	cd annotation/ensembl-${ENSEMBL_RELEASE}
	wget ftp://ftp.ensembl.org/pub/${ENSEMBL_RELEASE}/gtf/${LC_SPECIES}/${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gtf.gz
	wget ftp://ftp.ensembl.org/pub/${ENSEMBL_RELEASE}/gff3/${LC_SPECIES}/${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gff3.gz
	
	gunzip ${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gtf.gz
	gunzip ${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gff3.gz

	# This part is for Stefan's pipeline
	/home/golharr/workspace/rgtools/scripts/ensembl/buildGeneModelFromEnsembl${ENSEMBL_RELEASE_NUM}.pl ${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gtf
	##
	cd ../..
fi
