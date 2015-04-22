#!/bin/bash

THREADS=8
IFS=$'\n'
#c=0
while read line
do
	SAMPLE=$(echo $line | cut -f 1)
	FQ1=$(echo $line | cut -f 2)
	FQ2=$(echo $line | cut -f 3)

	if [ ! -e $SAMPLE/$SAMPLE.bt.exon_quantification.txt ]; then
		#c=$(($c+1))
		#echo $c. $SAMPLE
		echo qsub -N j$SAMPLE -v SAMPLE=${SAMPLE},FASTQ1=${FQ1},FASTQ2=${FQ2},THREADS=$THREADS,AWS=1 ./unc_rnaseqV2_pipeline.aws.sh

		#if [ $c -ge 20 ]; then
		#	exit 0
		#fi
	fi
done < samples.txt

#./makeExporessionMatrix.sh
