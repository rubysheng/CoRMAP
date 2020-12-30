#!/bin/bash
#title          :section3.3_rename.sh
#description    :Change the formate of sequence headers for running trinity smoothly.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.2
#usage          :./section3.3_rename.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

##########
# rename #
##########

function rename_sr () {
  echo ==== START RECONSTRUCT ====
  cd ./trim/SR/
  if [[ ! -f  "*_renamed.fq.gz" ]]; then
    for var in `find . -type f -name "*_trimmed.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm *_trimmed.fq.gz
  fi
  echo 'All single-end sequencing data have renamed!'
  cd ../..
  echo 'Back to main directory'
  echo ==== END RECONSTRUCT ====
  echo
}

function rename_pe () {
  echo ==== START RECONSTRUCT ====
  cd ./trim/PE/
  if [[ ! -f  "*_renamed.fq.gz" ]]; then
    for var in `find . -type f -name "*_val_[12].fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm *_R[12]_val_[12].fq.gz
  fi
  echo 'All pair-end sequencing data have renamed!'
  cd ../..
  echo 'Back to main directory'
  echo ==== END RECONSTRUCT ====
  echo
}
