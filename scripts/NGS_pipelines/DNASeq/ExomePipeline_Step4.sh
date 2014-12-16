#!/bin/bash

# SGE Grid Engine Parameters
#$ -pe orte 4
#$ -j y

# PBS Parameters
#PBS -l nodes=1:ppn=4
#PBS -j oe

# Example to run:
# FASTQ files are assumed to be in $BASEDIR as specified in exome_common_parameters.ini
# qsub -v PARAMETERS_FILE=`pwd`/exome_common_parameters.ini `pwd`/ExomePipeline_Step4.sh

VERSION="v0.1"

echo "Running Exome Analysis Pipeline ($VERSION) on `uname -n`"
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
##  Exome Analysis Pipeline 		                    ##
##############################################################

# We can define this as an environment variable.
if [ -z ${PARAMETERS_FILE} ]
then
#	Because each project may have different parameters, its better load 
#	a project specific file from the current directory, if its not
#       specified.
	PARAMETERS_FILE=./exome_common_parameters.ini
fi

if [ ! -e ${PARAMETERS_FILE} ] 
then
	echo "Unable to locate ${PARAMETERS_FILE}"
	exit -1
fi

echo "Using parameters file $PARAMETERS_FILE"
date '+%m/%d/%y %H:%M:%S' 
echo
. ${PARAMETERS_FILE}

##############################################################
##  Make analysis directory                                 ##
##############################################################

if [ -z "$BASE_DIR" ]
then
	echo "BASE_DIR is not set"
	exit 1
fi

if [ ! -d "${BASE_DIR}" ] 
then
	echo "$BASE_DIR does not exist"
	exit -1
fi

cd ${BASE_DIR}

echo Using Analysis Directory: `pwd` 
date '+%m/%d/%y %H:%M:%S'
echo

##############################################################
## Print out version of tools used in this analysis         ##
##############################################################

echo The following applications will be used:
echo fastqValidator: ${FASTQ_VALIDATOR}
echo
echo FastQC: ${FASTQC}
echo
echo BWA: ${BWA}
 ${BWA} 2>&1 | head -n 3 | tail -n 1
echo BWA Reference: ${BWA_REFERENCE}
echo
echo Java: ${JAVA}
 ${JAVA} -version 2>&1 | head -n 1
echo
echo PicardTools: ${PICARD_TOOLS}
echo
echo GATK: ${GATK_JAR}
echo GATK Version `${JAVA} -jar $GATK_JAR -version`
echo
echo SAMTools: ${SAM_TOOLS}
echo SAMTools `${SAM_TOOLS} 2>&1 | head -n 3 | tail -n 1`
echo

if [ ! -e "cohort.raw.vcf" ]
then
	echo Calling Variants using Haplotype Caller
	date '+%m/%d/%y %H:%M:%S'
	echo

	$JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T HaplotypeCaller  \
		-R $REFERENCE \
		-I bamfiles.recal.list \
		--dbsnp $GATK_KNOWN_SNPS \
		-o cohort.raw.vcf

	if [ $? -ne 0 ]
	then
		echo "Error running HaplotypeCaller"
		rm cohort.raw.vcf 2>/dev/null
		exit -1
	fi

        echo Finished running Haplotype Caller
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

# Need to apply VariantRecalibration here but without a good set of variants to calibrate to, this is impossible.

# Then we can evaluate the variant call set:
if [ ! -e "cohort.raw.eval.txt" ]
then
	echo Evalulating variant quality
	date '+%m/%d/%y %H:%M:%S'
        echo

	$JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
                -T VariantEval \
                --eval cohort.raw.vcf \
                -R $REFERENCE \
		--dbsnp $GATK_KNOWN_SNPS \
                -o cohort.raw.eval.txt

        if [ $? -ne 0 ]
        then
                echo "Error running VariantEval"
                rm cohort.raw.eval.txt 2>/dev/null
                exit -1
        fi
	
        echo Finished evalulating variant quality 
	date '+%m/%d/%y %H:%M:%S'
        echo
fi

# Until then, just filtered for variants with at least 10 reads:
if [ ! -e "cohort.raw.filtered.vcf" ]
then
	echo Filtering variant based on read depth (DP>=10)
        date '+%m/%d/%y %H:%M:%S'
        echo

	$JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
		-T SelectVariants \
		-R $REFERENCE \
		-V cohort.raw.vcf \
		-o cohort.filtered.vcf \
		-select "DP>=10" \
		-env \
		-ef
		
        if [ $? -ne 0 ]
        then
                echo "Error selecting variants"
                rm cohort.filtered.vcf 2>/dev/null
                exit -1
        fi

        echo Finished filtering variants
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

if [ ! -e "cohort.filtered.eval.txt" ]
then
        echo Evalulating variant quality
        date '+%m/%d/%y %H:%M:%S'
        echo

        $JAVA -jar $GATK_DIR/GenomeAnalysisTK.jar \
                -T VariantEval \
                --eval cohort.filtered.vcf \
                -R $REFERENCE \
                --dbsnp $GATK_KNOWN_SNPS \
                -o cohort.filtered.eval.txt

        if [ $? -ne 0 ]
        then
                echo "Error running VariantEval"
                rm cohort.filtered.eval.txt 2>/dev/null
                exit -1
        fi

        echo Finished evalulating variant quality 
        date '+%m/%d/%y %H:%M:%S'
        echo
fi

echo
echo Finished analyzing ${SAMPLE_NAME} 
date '+%m/%d/%y %H:%M:%S' 
echo
@
