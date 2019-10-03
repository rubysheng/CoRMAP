#!/bin/bash
############################################################################
# This script is for process a study from the beginning (after downloaded) #
# created date: 2019-10-03 #
# author : Ruby Sheng #
# input files format: *.fastq.gz ... #
# output files: expression matrix, annotation table ... #
############################################################################

date

source $RUBY_GITHUB/script/process.sh # checked: works well
source $RUBY_GITHUB/script/define_type.sh

function mainflow() {
  #this function is to call specific processing step
  if [[ '$step' == "trimming" ]]; then
    # check the quality of raw sequences
    raw_qc
    # classify the layout types
    def_trimming
  elif [[ '$step' == "assembly" ]]; then
    # classify the layout types
    def_finddir_assembly
  fi
}
