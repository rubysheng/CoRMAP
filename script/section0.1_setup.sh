#!/bin/bash
#title          :section0.1_setup.sh
#description    :generate directories and subdirectories for each datasets.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section0.1_setup.sh list_of_projectnames
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

export ROOT=$(pwd)

# Generate $ROOT/data/ directory to hold raw data, intermediate data, and final reports for each dataset
while read line
do
  echo "Generating a folder in '"$ROOT"/data/' for the dataset - "${line}
  mkdir -p $ROOT/data/${line}
done < $1

ls -l ./data

echo
echo

# Generate $ROOT/download_input/ directory to hold download links of all datasets
echo "Generating a directory '"$ROOT"/download_input/' to hold input download links"
mkdir $ROOT/download_input
