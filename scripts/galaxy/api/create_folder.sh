#!/bin/bash

folder=$1
parent_folder=$2
desc=$3

curl --data "encoded_parent_folder_id=$2&name=$1&param=$3" "http://kraken.pri.bms.com:18080/api/folders?key=2fd87eb290956bee983cc33e4c227c3e"
