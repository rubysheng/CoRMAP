#!/bin/bash
# this is for define the dataset layout types
# note: separate into 3 mainsteps: trimming, assembly, (quantify&annotations)

function def_trimming () {
  # test layout type:
  SE_N=$[`ls *.fastq.gz|grep -v "_"| wc -l|cut -f1 -d' '`]
  PE_L_N=$[`ls *.fastq.gz| grep -c "_1"|cut -f1 -d' '`]
  PE_R_N=$[`ls *.fastq.gz| grep -c "_2"|cut -f1 -d' '`]
#  PE_L_N=$[`ls *.fastq.gz| grep -c "_R1"|cut -f1 -d' '`]
#  PE_R_N=$[`ls *.fastq.gz| grep -c "_R2"|cut -f1 -d' '`]
  # only pe
  if [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -eq '0' ]; then
    echo "Only Paired-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    #1. trim
    trim_pe > trim_report.log
    #2. generate cleaned data fastqc report
    conda activate multiqc  # activates environment
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ #> trimedqc_report.log
    conda deactivate         # deactivates
    #3. prepare: rename files, and input the grouping table of samples
    rename_pe #> rename_report.log

  # only sr
  elif [ ${PE_L_N} -eq '0' ] && [ ${SE_N} -gt '0' ]; then
    echo "Only Single-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    #1. trim
    trim_sr #> trim_report.log
    #2. generate cleaned data fastqc report
    conda activate multiqc  # activates environment
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ #> trimedqc_report.log
    conda deactivate         # deactivates
    #3. prepare: rename files, and input the grouping table of samples
    rename_sr #> rename_report.log

  # both
  elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -gt '0' ]; then
    echo "Both types of layout here"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    #1. trim
    trim_pe #>> trimpe_report.log
    trim_sr #>> trimsr_report.log
    #2. generate cleaned data fastqc report
    conda activate multiqc  # activates environment
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ #>> trimedqc_report.log
    conda deactivate         # deactivates
    #3. prepare: rename files, combine two types of layout data into one file, and input the grouping table of samples
    pretrinity_both #>> pretrinity_report.log

  # error
  else
    echo "there is no data can be analyzed"
    exit 1
  fi
}

function def_finddir_assembly () {
  # find the directory under ./trim/
  DIR_SR="./trim/SR/"
  DIR_PE="./trim/PE/"
  echo $(pwd)
  # check if the group design table has been input
  if [ -e sample_file_?.txt ]; then
    # SR
    if [ -d "$DIR_SR" ] && [ ! -d "$DIR_PE" ]; then
      echo "Only Single-end reads"
      cd ./trim/SR/
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      assembly > assembly_report.log
      #assembly_sr > assembly_report.log
    # PE
    elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
      echo "Only Paired-end reads"
      cd ./trim/PE/
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      assembly > assembly_report.log
      #assembly_pe > assembly_report.log
    # BOTH
    elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
      echo "Both types of layout here"
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      #assembly > assembly_report.log
      assembly_bo > assembly_report.log
    # error
    else
      echo "there is no data can be analyzed"
      #exit 1
    fi
  else
    echo "no input file for the group design"
    #exit 1
  fi

}


function def_quantification () {
  echo ${PRJNA_PATH}

  # check if there is Trinity.fasta
  cd trinity_out_dir
  if [ -e Trinity.fasta ]; then
    cd ${PRJNA_PATH}
    count > quantify.log
  else
    echo "no Trinity.fasta as input file"
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
