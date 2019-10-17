#!/bin/bash
#title          :section3.1_rawdata_qualitycheck.sh
#description    :Quality check for raw data before quality control.
#author         :Ruby(Yiru) Sheng
#date           :20191016
#version        :1.1
#usage          :./section3.1_rawdata_qualitycheck.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

function check () {
  # find subdirectories that already have "raw_qc"
  find 
}
############################
# raw data quality control #
############################
function raw_qc () {
  # mkdir ./raw_fastqc/
  # generate quality control with each untrimmed runs
  fastqc -o ./raw_fastqc/ -f fastq -t 5 --extract *.fastq.gz #SRR*.fastq.gz
  # combine reports of all runs in a study to one quality control report
  conda activate multiqc  # activates environment
  multiqc ./raw_fastqc/ -o ./raw_fastqc/multiqc_output/
  conda deativate         # deactivates

  touch ./raw_fastqc/multiqc_output/log_list.txt
  ls -l ./raw_fastqc/multiqc_output/multiqc_data/ > ./raw_fastqc/multiqc_output/log_list.txt

  # remove the separate fastqc files just keep the multiqc report
  find ./raw_fastqc -type f -name "*_Replicate_*" -exec rm {} \;
  find ./raw_fastqc -type d -name "*_Replicate_*" -exec rm -r {} \;
  # find ./raw_fastqc -type f -name "SRR*" -exec rm {} \;
  # find ./raw_fastqc -type d -name "SRR*" -exec rm -r {} \;
}
raw_qc
