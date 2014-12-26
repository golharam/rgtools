#!/bin/sh

for f in `find /ngs/ngs07/benaroya2/Bio_Inf_Ext_Drv/compoundAnnoStudy -name '*.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f "Benaroya/Batch 2"
done

for f in `find /ngs/ngs07/benaroya2/Bio_Inf_Ext_Drv/LibsForFASTQtoBMS/ -name '*.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f "Benaroya/Batch 3"
done

