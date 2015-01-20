#!/bin/sh

for f in `find /ngs/ngs18/m2gen-d? -name '*.clipped.fastq.gz'`
do 
  ./add_file_to_library.py $GALAXY_KEY $GALAXY_URL $f "M2GEN/Clipped_FastQs"
done

for f in `find /ngs/ngs18/m2gen-d? -name '*_1.fastq.gz'`
do
  ./add_file_to_library.py $GALAXY_KEY $GALAXY_URL $f "M2GEN/Raw_FastQs"
done

for f in `find /ngs/ngs18/m2gen-d? -name '*_2.fastq.gz'`
do
  ./add_file_to_library.py $GALAXY_KEY $GALAXY_URL $f "M2GEN/Raw_FastQs"
done

