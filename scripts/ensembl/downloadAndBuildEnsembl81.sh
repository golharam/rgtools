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

if [ ! -e ${ASSEMBLY}.ercc.fa ]; then
	echo "Building genome with ERCC"
	cat ${ASSEMBLY}.fa ../ercc/ERCC92/ERCC92.fa > ${ASSEMBLY}.ercc.fa
fi

if [ ! -e ${ASSEMBLY}.dict ]; then
	echo "Building genome dictionary"
	/apps/sys/galaxy/external_packages/bin/java \
		-jar /apps/sys/galaxy/external_packages/picard-tools-1.137/picard.jar \
		CreateSequenceDictionary \
		R=${ASSEMBLY}.fa O=${ASSEMBLY}.dict GENOME_ASSEMBLY=${ASSEMBLY}
fi

if [ ! -e ${ASSEMBLY}.ercc.dict ]; then
        echo "Building genome dictionary with ERCC"
        /apps/sys/galaxy/external_packages/bin/java \
                -jar /apps/sys/galaxy/external_packages/picard-tools-1.137/picard.jar \
                CreateSequenceDictionary \
                R=${ASSEMBLY}.ercc.fa O=${ASSEMBLY}.ercc.dict GENOME_ASSEMBLY=${ASSEMBLY}
fi

if [ ! -e ${ASSEMBLY}.fa.fai ]; then
	echo "Build genome index"
	/apps/sys/galaxy/external_packages/samtools-0.1.19/samtools faidx ${ASSEMBLY}.fa
fi

if [ ! -e ${ASSEMBLY}.ercc.fa.fai ]; then
        echo "Build genome with ercc index"
        /apps/sys/galaxy/external_packages/samtools-0.1.19/samtools faidx ${ASSEMBLY}.ercc.fa
fi

if [ ! -d bwa_index ]; then
	echo "Building BWA Index"
	mkdir bwa_index
	cd bwa_index
	echo "/apps/sys/galaxy/external_packages/bwa-0.7.12/bwa index -p ${ASSEMBLY} ../${ASSEMBLY}.fa"| qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bwa
	cd ..
fi

if [ ! -e bwa_index/GRCh38.ercc.sa ]; then
	echo "Building BWA with ercc index"
	cd bwa_index
	echo "/apps/sys/galaxy/external_packages/bwa-0.7.12/bwa index -p ${ASSEMBLY}.ercc ../${ASSEMBLY}.ercc.fa"| qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-ercc-bwa
	cd ..
fi

if [ ! -d bowtie_index ]; then
	echo "Building bowtie index"
	mkdir bowtie_index
	cd bowtie_index
	echo "/apps/sys/galaxy/external_packages/bowtie-1.1.1/bowtie-build ../${ASSEMBLY}.fa ./${ASSEMBLY}" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bowtie
	cd ..
fi

if [ ! -e bowtie_index/GRCh38.ercc.1.ebwt ]; then
	echo "Building bowtie with ercc index"
	cd bowtie_index
	echo "/apps/sys/galaxy/external_packages/bowtie-1.1.1/bowtie-build ../${ASSEMBLY}.ercc.fa ./${ASSEMBLY}.ercc" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-ercc-bowtie
	cd ..
fi

if [ ! -d bowtie2_index ]; then
        echo "Building bowtie2 index"
        mkdir bowtie2_index
        cd bowtie2_index
        echo "/apps/sys/galaxy/external_packages/bowtie2-2.2.6/bowtie2-build ../${ASSEMBLY}.fa ./${ASSEMBLY}" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-bowtie2
        cd ..
fi

if [ ! -e bowtie2_index/GRCh38.ercc.1.bt2 ]; then
	echo "Building bowtie2 with ercc index"
	cd bowtie2_index
	echo "/apps/sys/galaxy/external_packages/bowtie2-2.2.6/bowtie2-build ../${ASSEMBLY}.ercc.fa ./${ASSEMBLY}.ercc" | qsub -cwd -j y -S /bin/bash -N ${ASSEMBLY}-ercc-bowtie2
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

if [ ! -d rsem_index.${ENSEMBL_RELEASE} ]; then
	echo "Building RSEM index for Ensembl ${ENSEMBL_RELEASE}"
	mkdir -p rsem_index.${ENSEMBL_RELEASE}
	cd rsem_index.${ENSEMBL_RELEASE}
	/apps/sys/galaxy/external_packages/rsem-1.2.22/rsem-prepare-reference \
		--gtf ../annotation/ensembl-${ENSEMBL_RELEASE}/${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gtf \
		--bowtie --bowtie-path /apps/sys/galaxy/external_packages/bowtie-1.1.1 \
		--bowtie2 --bowtie2-path /apps/sys/galaxy/external_packages/bowtie2-2.2.6 \
		--star --star-path /apps/sys/galaxy/external_packages/STAR_2.4.1d/bin/Linux_x86_64_static \
		-p 8 \
		../${ASSEMBLY}.fa ./${ASSEMBLY}
	cd ..
fi

if [ ! -d rsem_index.${ENSEMBL_RELEASE}.ercc ]; then
	echo "Building RSEM index with ERCC for Ensembl"
	mkdir -p rsem_index.${ENSEMBL_RELEASE}.ercc
	cd rsem_index.${ENSEMBL_RELEASE}.ercc
	cat ../annotation/ensembl-${ENSEMBL_RELEASE}/${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.gtf ../../ercc/ERCC92/ERCC92.gtf > ./${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.ercc.gtf
        /apps/sys/galaxy/external_packages/rsem-1.2.22/rsem-prepare-reference \
                --gtf ${SPECIES}.${ASSEMBLY}.${ENSEMBL_RELEASE_NUM}.ercc.gtf \
                --bowtie --bowtie-path /apps/sys/galaxy/external_packages/bowtie-1.1.1 \
                --bowtie2 --bowtie2-path /apps/sys/galaxy/external_packages/bowtie2-2.2.6 \
                --star --star-path /apps/sys/galaxy/external_packages/STAR_2.4.1d/bin/Linux_x86_64_static \
                -p 8 \
                ../${ASSEMBLY}.ercc.fa ./${ASSEMBLY}.ercc
        cd ..

fi
