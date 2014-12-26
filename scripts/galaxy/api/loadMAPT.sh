#!/bin/sh

for f in `find /home/golharr/ngs/work/BMS_MAPT/FASTQ/ -name '*.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f "MAPT"
done

