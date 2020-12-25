#!/bin/bash
#title          :movetrimreport.sh
#description    :move trimming report to a new directory
#author         :Ruby(Yiru) Sheng
#date           :20191030
#version        :2.0
#usage          :./movetrimreport.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

function mvtrimreport () {
  echo $1
  DIR=$1
  echo ${DIR}
  mkdir ${DIR}/trim_report_dir/
  OUTDIR="${DIR}/trim_report_dir/"
  cp -r ${DIR}/trim/trimmed_fastqc ${OUTDIR}
  rm -r ${DIR}/trim/trimmed_fastqc

  DIR_SR="${DIR}/trim/SR/"
  DIR_PE="${DIR}/trim/PE/"
  # SR
  if [ -d "$DIR_SR" ]; then
    cp ${DIR}/trim/SR/*_trimming_report.txt ${OUTDIR}
    rm ${DIR}/trim/SR/*_trimming_report.txt
    cp ${DIR}/trim/SR/*_fastqc.* ${OUTDIR}
    rm ${DIR}/trim/SR/*_fastqc.*
  # PE
  elif [ -d "$DIR_PE" ]; then
    cp ${DIR}/trim/PE/*_trimming_report.txt ${OUTDIR}
    rm ${DIR}/trim/PE/*_trimming_report.txt
    cp ${DIR}/trim/PE/*_fastqc.* ${OUTDIR}
    rm ${DIR}/trim/PE/*_fastqc.*
  fi

}

# find directories in "data", and add to a text file in "ruby"
ls -1  > ../path.txt

# go to each location directory by read lines of the "path.log" via a for loop
ROOT_DIR=$(pwd)
for line in `cat ../path.txt`; do
  echo "${line}"
  cd ${line}
  mvtrimreport ${line} #> test_log_file.txt || true
  # nohup mvtrimreport ${line} &
  cd ${ROOT_DIR}
done
