#!/bin/bash
#title          :section3.1_rawdata_qualitycheck.sh
#description    :Quality check for raw data before quality control.
#author         :Ruby(Yiru) Sheng
#date           :20191018
#version        :1.2
#usage          :./section3.1_rawdata_qualitycheck.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================


############################
# raw data quality control #
############################
function raw_qc () {
  #mkdir ./raw_fastqc/
  # generate quality control with each untrimmed runs
  fastqc -o ./raw_fastqc/ -f fastq -t 5 --extract SRR*.fastq.gz
  # combine reports of all runs in a study to one quality control report
  conda activate multiqc  # activates environment
  multiqc ./raw_fastqc/ -o ./raw_fastqc/multiqc_output/
  conda deactivate         # deactivates

  touch ./raw_fastqc/multiqc_output/log_list.txt
  ls -l ./raw_fastqc/multiqc_output/multiqc_data/ > ./raw_fastqc/multiqc_output/log_list.txt

  # remove the separate fastqc files just keep the multiqc report
  #find ./raw_fastqc -type f -name "*_Replicate_*" -exec rm {} \;
  #find ./raw_fastqc -type d -name "*_Replicate_*" -exec rm -r {} \;
  find ./raw_fastqc -type f -name "SRR*" -exec rm {} \;
  find ./raw_fastqc -type d -name "SRR*" -exec rm -r {} \;
}


function check_rawqc () {
  for file in ./*
  do
    cd ${file}
    # run "raw data quality check" only in the unprocessed datasets
    if [ ! -d ./raw_qc/multiqc_output ]; then
      raw_qc
      echo "${file} is completed with raw data quality check"
    else
      echo "Raw data have been checked. Skip this step."
    fi
    cd ..
  done
  # `find . -type d -name "raw_qc" | awk -F'/' '{ print $2 }'`
}

check_rawqc
