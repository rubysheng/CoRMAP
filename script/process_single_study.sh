#!/bin/bash
############################################################################
# This script is for process a study from the beginning (after downloaded) #
# created date: 2019-10-03 #
# author : Ruby Sheng #
# input files format: *.fastq.gz ... #
# output files: expression matrix, annotation table ... #
############################################################################

date
export RUBY_GITHUB=/media/lewis/Seagate_Backup_Plus_Drive/ruby/github/bombina/comparative-transcriptomic-analysis-pip

source $RUBY_GITHUB/script/process.sh # checked: works well
source $RUBY_GITHUB/script/define_type.sh

function mainflow() {
  #this function is to call specific processing step
  if [[ $1 == "trimming" ]]; then
    # check the quality of raw sequences
    #raw_qc > raw_qc.log
    # classify the layout types
    def_trimming
  elif [[ $1 == "assembly" ]]; then
    # classify the layout types
    def_finddir_assembly
  fi
}


read -p "Which step of processing will run?" ANSWER

mainflow $ANSWER
