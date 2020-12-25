#!/bin/bash
# this is for define the dataset layout types
# note: separate into 3 mainsteps: trimming, assembly, (quantify&annotations)

# function def_trimming () {
#   source section3.2_trim.sh
#   # 1. trim
#     # only pe
#   if [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -eq '0' ]; then
#     echo "Only Paired-end reads"
#     trim_pe
#
#     # only sr
#   elif [ ${PE_L_N} -eq '0' ] && [ ${SE_N} -gt '0' ]; then
#     echo "Only Single-end reads"
#     trim_sr
#
#     # both
#   elif [ ${PE_L_N} -gt '0' ] && [ ${PE_L_N}=${PE_R_N} ] && [ ${SE_N} -gt '0' ]; then
#     echo "Both types of layout here"
#     trim_pe
#     trim_sr
#
#     # error
#   else
#     echo "there is no data can be analyzed"
#   fi
#
#   # 2. generate cleaned data fastqc report
#   conda activate multiqc  # activates environment
#   multiqc ./trim/ --outdir ./trim/trimmed_fastqc/ #> trimedqc_report.log
#   conda deactivate         # deactivates
#
# }

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


function def_finddir_normalize () {
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
      # source section3.4_normalizeBOTH.sh
      # normalization
      normalize_sr 2>&1 | tee normalize.log


    # PE
    elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
      echo "Only Paired-end reads"
      cd ./trim/PE/
      # source section3.4_normalizeBOTH.sh
      # normalization
      normalize_pe 2>&1 | tee normalize.log


    # BOTH
    elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
      echo "Both types of layout here"
      echo "Data has been normalized"
      # mkdir ${PRJNA_PATH}/normalization/
      # OUTDIR="${PRJNA_PATH}/normalization/"
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
  # source ${SCRIPT_LOC}/section4.1_assembly.sh
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
      assembly > assembly_report.log
    # PE
    elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
      echo "Only Paired-end reads"
      cd ./trim/PE/
      assembly > assembly_report.log
    # BOTH
    elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
      echo "Both types of layout here"
      assembly_bo > assembly_report.log
    # error
    else
      echo "there is no data can be analyzed"
    fi
  else
    echo "no input file for the group design"
  fi

}


function def_quantification () {
  source ${SCRIPT_LOC}/section4.2_quantify.sh

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

function def_expmx () {
  source ${SCRIPT_LOC}/section4.3_matrix.sh

  cd ${PRJNA_PATH}
  echo ===================================
  echo === change to expression matrix ===
  echo ===================================

  expressionmx > exprmx.log

  echo ===================================
  echo === check the expression matrix ===
  echo ===================================

}

function def_anno () {
  cd ${PRJNA_PATH}
  echo ========================
  echo === start annotation ===
  echo ========================

  # check if the local annotation database has been built
  DB="${MAIN}/annotation/Ruby_transpip.sqlite"
  if test -f "$DB"; then
    echo "$DB exist"
  else
    echo "start generate the local database first..."
    $TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Ruby_transpip
  fi

  # Select annotation type
  printf "%s\n%s\n%s\n\n" "What do you want after annotation?(input with number please)" \
         "1 for known sequence data" \
         "2 for known sequence data and protein prediction"
  read -r typechosen

  # Start to process data.
  case $typechosen in
    1) annotation &> anno.log
      ;;
    2) predict_pro &> anno.log
      ;;
    *) #error
      printf "%s\n%s\n\n" "I did not understand your selection." \
             "Quit."
      ;;
  esac

  echo ======================
  echo === end annotation ===
  echo ======================

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
