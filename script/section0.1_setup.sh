#!/bin/bash
#title          :section0.1_setup.sh
#description    :generate directories and subdirectories for each datasets.
#author         :Ruby(Yiru) Sheng
#usage          :source ./script/section0.1_setup.sh list_of_projectnames
#bash_version   :4.4.19(1)-release
#============================================================================

# Generate ./data/ directory to hold data for each dataset
while read line
do
  echo "Generating a folder in './data/' for dataset - "${line}
  mkdir -p ./data/${line}
done < $1

ls -l ./data

echo
echo

# Generate ./download_input/ directory to hold download links of all datasets
echo "Generating a directory './download_input/' to hold input download links"
mkdir ./download_input
