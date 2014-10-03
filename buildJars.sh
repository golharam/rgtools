#!/bin/bash

for xml in `ls build_scripts/*.xml`
do
	ant -f $xml
done
