#!/bin/sh

for f in `find /home/golharr/ngs/work/PD1_FFPE_CancerPanel/source_data/BaseCalls -name '*.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f PD1_FFPE_CancerPanel/BaseCalls
done

for f in `find /home/golharr/ngs/work/PD1_FFPE_CancerPanel/source_data/OldCalls -name '*.gz'`
do 
  ./add_file_to_library.py $GALAXY_API_KEY $GALAXY_API_URL $f PD1_FFPE_CancerPanel/OldCalls
done

