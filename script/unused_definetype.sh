#!/bin/bash

# ----- log : 20200214 ------------------------------------------ #
# the changed part in pip.sh (scripts below are the original version)
    # printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n"
           # "5 for make the local blast database." \
           # "6 for local blast." \
           # "7 for load to SQL database." \
  # 5) #make the local blast database
  # def_mkdb
  #   ;;
  # 6) #local blast
  # def_blast
  #   ;;
  # 7) #load to SQL database
  # def_loadsql
  #   ;;
# ----- log end : 20200214 -------------------------------------- #
function def_mkdb () {
  echo =====================================
  echo === make the local blast database ===
  echo =====================================
  mk_local_blastdb 2>&1 | tee ${PRJNA_PATH}/mkdb.log
}

function def_blast () {
  lcblast 2>&1 | tee ${PRJNA_PATH}/blast.log
}


function def_loadsql () {
  echo ${PRJNA_PATH}
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
  echo ">>> Loding to SQL database >>>"
  echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

  load_sql 2>&1 | tee ${PRJNA_PATH}/load_sql.log

  echo "============================="
  echo "=== Loded to SQL database ==="
  echo "============================="
}

###############################################
# function def_finddir_quantify () {
#   source ${scriptloc}/section4.2_quantify.sh
#
#   # check if there is Trinity.fasta
#   cd trinity_out_dir
#   if [ -e Trinity.fasta ]; then
#     cd ../
#
#     # find the directory under ./trim/
#     DIR_SR="./trim/SR/"
#     DIR_PE="./trim/PE/"
#     # SR
#     if [ -d "$DIR_SR" ] && [ ! -d "$DIR_PE" ]; then
#       echo "Only Single-end reads"
#       cd ./trim/SR/
#       count > quantify.log
#     # PE
#     elif [ -d "$DIR_PE" ] && [ ! -d "$DIR_SR" ]; then
#       echo "Only Paired-end reads"
#       cd ./trim/PE/
#       count > quantify.log
#     # BOTH
#     elif [ -d "$DIR_PE" ] && [ -d "$DIR_SR" ]; then
#       echo "Both types of layout here"
#       # cd
#       # ???
#       # count > quantify.log
#   else
#     echo "no Trinity.fasta as input file"
#   fi
# }
###############################################

function def_expmx () {
  cd ${PRJNA_PATH}
  echo ===================================
  echo === change to expression matrix ===
  echo ===================================

  expressionmx > exprmx.log

  echo ===================================
  echo === check the expression matrix ===
  echo ===================================

}

export MAIN=/media/lewis/Seagate_Backup_Plus_Drive/ruby

function def_anno () {
  cd ${PRJNA_PATH}
  echo ========================
  echo === start annotation ===
  echo ========================



  # check if the local annotation database has been built
  DB="${MAIN}/annotation/Ruby_transpip.sqlite"
  if [ -e $DB ] ; then
    echo "$DB exist"
  else
    echo "start generate the local database first..."
    cd ${MAIN}/annotation
    $TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Ruby_transpip
  fi

  # Select annotation type
  printf "%s\n%s\n%s\n\n" "What do you want after annotation?(input with number please)" \
         "1 for known sequence data" \
         "2 for known sequence data and protein prediction"
  read -r typechosen

  # Start to process data.
  case $typechosen in
    1) annotation > anno.log
      ;;
    2) predict_pro > anno.log
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
