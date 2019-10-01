#!/bin/bash

function SR_ONLY () {
  #1. trim
  #trim_sr >& trim_report.log &
  #2. generate cleaned data fastqc report
  #multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ >& trimedqc_report.log &
  #3. prepare: rename files, and input the grouping table of samples
  #rename_sr >& rename_report.log &
  #generateruninf
  #4. de novo assembly
  #assembly_sr >& assembly_report.log &
  #5.alignment and transcript quantification
  cd ${PRJNA_PATH}
  mkdir ./transcripts_count/
  count_sr >& count_report.log &
  expressionmx >& expressionmx_report.log &
  #rerun_abundance_estimation
#  echo 1
}
function PE_ONLY () {
  #1. trim
  #trim_pe >& trim_report.log &
  #2. generate cleaned data fastqc report
  #multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ >& trimedqc_report.log &
  #3. prepare: rename files, and input the grouping table of samples
  #rename_pe >& rename_report.log &
  #generateruninf
  #4. de novo assembly
  #assembly_pe >& assembly_report.log &
  #5.alignment and transcript quantification
  cd ${PRJNA_PATH}
  mkdir ./transcripts_count/
  count_pe >& count_report.log &
  expressionmx >& expressionmx_report.log &
#  echo 2
}

function BOTH () {
  #1. trim
  #trim_pe >& trimpe_report.log &
  #trim_sr >& trimsr_report.log &
  #2. generate cleaned data fastqc report
  #multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ >& trimedqc_report.log &
  #3. prepare: rename files, combine two types of layout data into one file, and input the grouping table of samples
  #rename_pe >& renamepe_report.log &
  #rename_sr >& renamesr_report.log &
  #pretrinity_both >& pretrinity_report.log &
  #generateruninf
  #4. de novo assembly
  #assembly_bo >& assembly_report.log &
  #5.alignment and transcript quantification
  #cd ${PRJNA_PATH}
  #mkdir ./transcripts_count/
  #count_bo >& count_report.log &
  #expressionmx >& expressionmx_report.log &
#  echo 3
}

function define () {
  # test layout type:
  SE_N=$[`ls *.fastq.gz|grep -v "_"| wc -l|cut -f1 -d' '`]
  PE_L_N=$[`ls *.fastq.gz| grep -c "_1"|cut -f1 -d' '`]
  PE_R_N=$[`ls *.fastq.gz| grep -c "_2"|cut -f1 -d' '`]
  # only pe
  if [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -eq '0' ]; then
    echo "Only Paired-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    PE_ONLY

  # only sr
  elif [ ${PE_L_N} -eq '0' ] && [ ${SE_N} -gt '0' ]; then
    echo "Only Single-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    SR_ONLY

  # both
  elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -gt '0' ]; then
    echo "Both types of layout here"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    BOTH

  # error
  else
    echo "there is no data can be analyzed"
  fi
}

function define_finddir () {
  # find the directory under ./trim/
  DIR_SR="./trim/SR/"
  DIR_PE="./trim/PE/"
  if [ -d "$DIR_SR" ] && [ ! -d "$DIR_PE" ]; then
    echo "Only Single-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    SR_ONLY

  elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
    echo "Only Paired-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    PE_ONLY

  elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
    echo "Both types of layout here"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    BOTH
  # error
  else
    echo "there is no data can be analyzed"
  fi
}

#function mainflow () {
#  line=$1
#  PRJNA_PATH=$(pwd)
#  echo "$PRJNA_PATH"
#  # raw data quality check
#  #raw_qc
#  # define the type of data layout, and run various processes
#  #define
#  # generate the transcript count matrix
#  #expressionmx
#  #rerun_exp
#  #rerun_plotcount
#  # differential expression analysis
#  #difexpre
#  #rerun_dge
#  # assembly quanlity assessment
#  #TrinityStats.pl ./trinity_out_dir/Trinity.fasta > ./trinity_out_dir/N50_stats_output.txt
#  # annotate
#  #annotation
#  echo #############################################
#  echo "#the data of project ${line} is all finished#"
#  echo #############################################
#  cd ${ROOT_DIR}
#  echo "Back to ${ROOT_DIR}"
#  echo
#}
