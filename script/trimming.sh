#!/bin/bash
# this is only for trimming:

function SR_ONLY () {
  # #1. trim
  # trim_sr > trim_report.log
  # #2. generate cleaned data fastqc report
  # multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ > trimedqc_report.log
  # #3. prepare: rename files, and input the grouping table of samples
  # rename_sr > rename_report.log
  #generateruninf
  #4. de novo assembly
  assembly_sr > assembly_report.log
  #5.alignment and transcript quantification
  # cd ${PRJNA_PATH}
  # mkdir ./transcripts_count/
  # count_sr > count_report.log
  # expressionmx > expressionmx_report.log
  # rerun_abundance_estimation
#  echo 1
}
function PE_ONLY () {
  # #1. trim
  # trim_pe > trim_report.log
  # #2. generate cleaned data fastqc report
  # multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ > trimedqc_report.log
  # #3. prepare: rename files, and input the grouping table of samples
  # rename_pe > rename_report.log
  #generateruninf
  #4. de novo assembly
  assembly_pe > assembly_report.log
  #5.alignment and transcript quantification
  # cd ${PRJNA_PATH}
  # mkdir ./transcripts_count/
  # count_pe > count_report.log
  # expressionmx > expressionmx_report.log
#  echo 2
}

function BOTH () {
  #1. trim
  #trim_pe > trimpe_report.log
  #trim_sr > trimsr_report.log
  #2. generate cleaned data fastqc report
  #multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ > trimedqc_report.log
  #3. prepare: rename files, combine two types of layout data into one file, and input the grouping table of samples
  #rename_pe > renamepe_report.log
  #rename_sr > renamesr_report.log
  #pretrinity_both > pretrinity_report.log
  #generateruninf
  #4. de novo assembly
  assembly_bo > assembly_report.log
  #5.alignment and transcript quantification
  #cd ${PRJNA_PATH}
  #mkdir ./transcripts_count/
  #count_bo > count_report.log
  #expressionmx > expressionmx_report.log
#  echo 3
}
