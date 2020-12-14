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
    trim_pe 2>&1 | tee trim.log
    #2. generate cleaned data fastqc report
    conda activate multiqc  # activates environment
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ 2>&1 | tee trimedqc.log
    conda deactivate         # deactivates
    #3. prepare: rename files, and input the grouping table of samples
    rename_pe

  # only sr
  elif [ ${PE_L_N} -eq '0' ] && [ ${SE_N} -gt '0' ]; then
    echo "Only Single-end reads"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    #1. trim
    trim_sr 2>&1 | tee trim.log
    #2. generate cleaned data fastqc report
    conda activate multiqc
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ 2>&1 | tee trimedqc.log
    conda deactivate
    #3. prepare: rename files, and input the grouping table of samples
    rename_sr

  # both
  elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -gt '0' ]; then
    echo "Both types of layout here"
    #mkdir ${PRJNA_PATH}/trinity_out_dir/
    #1. trim
    trim_pe 2>&1 | tee trimpe.log
    trim_sr 2>&1 | tee trimsr.log
    #2. generate cleaned data fastqc report
    conda activate multiqc
    multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ 2>&1 | tee trimedqc.log
    conda deactivate
    #3. prepare: rename files, combine two types of layout data into one file, and input the grouping table of samples
    pretrinity_both 2>&1 | tee  pretrinityBOTH.log

  # error
  else
    echo "there is no data can be analyzed"
    exit 1
  fi
}



function def_Normalize () {
  find . -type f -name *_renamed.fq.gz -exec cp -v {} /home/lewis/Documents/Normalize/input_data/ \;
  cd /home/lewis/Documents/Normalize/input_data/
  # in ~/Documents/Normalize/input_data/
  # test layout type:
  SE_N=$[`ls -1 *_trimmed_renamed.fq.gz | wc -l | cut -f1 -d' '`]
  PE_L_N=$[`ls -1 *_1_val_1_renamed.fq.gz | wc -l | cut -f1 -d' '`]
  PE_R_N=$[`ls -1 *_2_val_2_renamed.fq.gz | wc -l | cut -f1 -d' '`]
  # SR
  if [ ${PE_L_N} -eq '0' ] && [ ${SE_N} -gt '0' ]; then
    echo "Only Single-end reads"
    normalize_sr 2>&1 | tee  ../normalize.log

  # PE
  elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -eq '0' ]; then
    echo "Only Paired-end reads"
    normalize_pe 2>&1 | tee ../normalize.log
  # BOTH
  elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -gt '0' ]; then
    echo "Both types of layout here"
    echo "Data has been normalized"
    # normalized
  # error
  else
    echo "there is no data can be analyzed"
  fi

}



function def_finddir_normalize () {
  # find the directory under ./trim/
  DIR_SR="./trim/SR/"
  DIR_PE="./trim/PE/"
  echo $(pwd)
  # check if the group design table has been input
  if [ -e sample_file_*.txt ]; then
    # SR
    if [ -d "$DIR_SR" ] && [ ! -d "$DIR_PE" ]; then
      echo "Only Single-end reads"
      cd ./trim/SR/
      # source section3.4_normalizeBOTH.sh
      # normalization
      normalize_sr 2>&1 | tee  ${PRJNA_PATH}/normalize.log
    # PE
    elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
      echo "Only Paired-end reads"
      cd ./trim/PE/
      # source section3.4_normalizeBOTH.sh
      # normalization
      normalize_pe 2>&1 | tee ${PRJNA_PATH}/normalize.log
    # BOTH
    elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
      echo "Both types of layout here"
      echo "Data has been normalized"
      # normalized
    # error
    else
      echo "there is no data can be analyzed"
    fi
  else
    echo "no input file for the group design"
  fi

}



function def_finddir_assembly () {
  # find the directory under ./trim/
  DIR_SR="./trim/SR/"
  DIR_PE="./trim/PE/"
  echo $(pwd)
  PRJNA_PATH=$(pwd)
  echo  "PRJNA_PATH="${PRJNA_PATH}
  # check if the group design table has been input
  if [ -e sample_file_*.txt ]; then
    # SR
    if [ -d "$DIR_SR" ] && [ ! -d "$DIR_PE" ]; then
      echo "Only Single-end reads"
      cd ./trim/SR/
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      assembly 2>&1 | tee ${PRJNA_PATH}/assembly_report.log
      #assembly_sr > assembly_report.log
    # PE
    elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
      echo "Only Paired-end reads"
      cd ./trim/PE/
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      assembly 2>&1 | tee ${PRJNA_PATH}/assembly_report.log
      #assembly_pe > assembly_report.log
    # BOTH
    elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
      echo "Both types of layout here"
      #mkdir ${PRJNA_PATH}/trinity_out_dir/
      #4. de novo assembly
      #assembly > assembly_report.log
      assembly_bo 2>&1 | tee ${PRJNA_PATH}/assembly_report.log
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


function def_quant_expmx () {
  echo ${PRJNA_PATH}
  # check if there is Trinity.fasta
  cd trinity_out_dir
  if [ -e Trinity.fasta ]; then
    cd ${PRJNA_PATH}
    count 2>&1 | tee ${PRJNA_PATH}/quantify.log

    # generate the expression matrix
    cd ${PRJNA_PATH}
    echo ====================================
    echo === convert to expression matrix ===
    echo ====================================

    expressionmx 2>&1 | tee ${PRJNA_PATH}/exprmx.log

    echo ===================================
    echo === check the expression matrix ===
    echo ===================================

  else
    echo "no Trinity.fasta as input file"
  fi
}


# function def_ortho () {
#   # under a Taxonomy_lst directory
#   # check if there is a text file named as "*.lst" providing the TAX_CODE and PATH_TO_DATADIR
#   if [ -e *.lst ]; then
#     # check if the input directory existed
#     input_dir="$(pwd)/input/"
#     if [ ! -e ${input_dir} ]; then
#       mkdir -v ${input_dir}
#     fi
#     echo ==========================
#     echo === ortholog searching ===
#     echo ==========================
#
#     run_ortho_pip 2>&1 | tee orthomcl_pip.log
#   fi
# }
# 
# ortho_sum 2>&1 | tee sum.log
