#!/bin/bash
#title          :section1_pipeline_ui.sh
#description    :This script is for process one dataset from the beginning (after downloaded)
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section1_pipeline_ui.sh
#input files    :*.fastq.gz ...
#output files   :expression matrix, annotation table ...
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

date

export PRJNA_PATH=$(pwd)
source $CMRP_PATH/script/section1.1_process.sh # checked: works well
source $CMRP_PATH/script/section1.2_define_type.sh

_choosestep(){

  # Select steps
  printf "%s\n%s\n%s\n%s\n%s\n%s\n\n" "Which step do you want to process? (input with number please)" \
         "1 for trim / quality control." \
         "2 for normalize." \
         "3 for assembly." \
         "4 for quantify & generate expression matrix." \
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
    0) #all
    def_trimming
    def_normalize
    def_finddir_assembly
    def_quant_expmx
      ;;
    *) #error
      printf "%s\n%s\n\n" "I did not understand your selection." \
             "Press <Ctrl-c> to quit."
      _choosestep
      ;;
  esac
}

_choosestep
