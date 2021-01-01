#!/bin/bash
############################################################################
# This script is for process one dataset from the beginning (after downloaded)
# created date: 2019-10-03
# author : Ruby Sheng
# input files format: *.fastq.gz ...
# output files: expression matrix, annotation table ...
############################################################################

date

export PRJNA_PATH=$(pwd)
source $SCRIPTS/process.sh # checked: works well
source $SCRIPTS/define_type.sh

_choosestep(){

  # Select steps
  printf "%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n" "Which step do you want to process? (input with number please)" \
         "1 for trim / quality control." \
         "2 for normalize." \
         "3 for assembly." \
         "4 for quantify & generate expression matrix." \
         "5 for annotate with predicted protein sequence & run orthomcl-pip." \
         "0 for all."
  read -r stepchosen

  # Start to process data.
  case $stepchosen in
    1) #trim
    # classify the layout types
    def_trimming
      ;;
    2) #normalize
    # need to find the "SR" or "PE" folder
    def_normalize
      ;;
    3) #assembly
    # need to find the "SR" or "PE" folder
    def_finddir_assembly
      ;;
    4) #quantify & generate expression matrix
    # need to find the "SR" or "PE" folder
    def_quant_expmx
      ;;
    5) #annotate with predicted protein sequence & run orthomcl-pip
    # need to prepare a list of TAX_CODEs and PATH_TO_DATADIRs
    def_ortho
      ;;
    0) #all
    def_trimming
    def_normalize
    def_finddir_assembly
    def_quant_expmx
    def_ortho
      ;;
    *) #error
      printf "%s\n%s\n\n" "I did not understand your selection." \
             "Press <Ctrl-c> to quit."
      _choosestep
      ;;
  esac
}

_choosestep
