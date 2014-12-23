#!/bin/bash

# Specify the library id as the first parameter

curl http://kraken.pri.bms.com:18080/api/libraries/$1/contents?key=$GALAXY_API_KEY

