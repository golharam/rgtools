#!/bin/sh

for f in `find /ngs/ngs17/F14FTSUSAT0355-0383/Standard_Analysis -name '*.clean.fq.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f F14FTSUSAT0355-0383
done

