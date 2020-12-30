#!/bin/bash
#title          :mv_rawdataandreport.sh
#description    :move raw data and report to new backup directory
#author         :Ruby(Yiru) Sheng
#date           :20191106
#version        :1.0
#usage          :source $RUBY_SCRIPTS/mv_rawdataandreport.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

function mv_exclude () {
  echo $1
  OUTDIR="/media/lewis/New_Seagate_Drive_8TB/ruby/rawdata_reports/$1"
  mv !(trim|transcripts_count) $OUTDIR
}

function mv_include () {
  echo $1
  OUTDIR="/media/lewis/New_Seagate_Drive_8TB/ruby/rawdata_reports/$1"
  mv *fastq.gz $OUTDIR
}
# find directories in "data", and add to a text file in "ruby"
#ls -1  > ../path.txt

# go to each location directory by read lines of the "path.log" via a for loop
ROOT_DIR="/media/lewis/Seagate_Backup_Plus_Drive/ruby/data/"
for line in `cat /media/lewis/New_Seagate_Drive_8TB/ruby/path.txt`; do
  echo "${line}"
  cd ${line}
  mv_include ${line} #> test_log_file.txt || true
  cd ${ROOT_DIR}
done
