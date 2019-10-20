#!/bin/bash
#title          :section3.4_normalizeBOTH.sh
#description    :Start Normalization Respectively for datasets containing both types pf layout.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.3
#usage          :./section3.4_normalizeBOTH.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

# SOURCE THE "PRJNA_PATH" !!!!!!!

##################
# normalization #
##################

#### usage guidence :   insilico_read_normalization.pl ####
#  --seqType <string>      :type of reads: ( 'fq' or 'fa')
#  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for
#                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
#  --max_cov <int>         :targeted maximum coverage for reads.
#
#  If paired reads:
#      --left  <string>    :left reads
#      --right <string>    :right reads
#
#  Or, if unpaired reads:
#      --single <string>   :single reads
#
#  Or, if you have read collections in different files you can use 'list' files, where each line in a list
#  file is the full path to an input file.  This saves you the time of combining them just so you can pass
#  a single file for each direction.
#      --left_list  <string> :left reads, one file path per line
#      --right_list <string> :right reads, one file path per line








################################
# both- prepare before trinity #
################################

function pretrinity_both () {
  echo ====Start Normalization Respectively====
  cd ./trim/SR/
  SRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize single-end sequencing data"
  insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
    --single ${SRFILES} --CPU 16
  echo "end of single-end sequencing data normalization"
  echo
  cd ../PE/
  LEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  RIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize paired-end sequencing data"
  insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
    --left ${LEFTFILES} --right ${RIGHTFILES} --CPU 16
  echo "end of single-end sequencing data normalization"
  cd ${PRJNA_PATH}
  echo
  echo ====End Normalization Respectively====
  echo
  echo "Start combining"
  zcat ${PRJNA_PATH}/trim/SR/SRR*.fq.gz ${PRJNA_PATH}/trim/PE/SRR*.fq.gz | gzip -c > ${PRJNA_PATH}/bigfile.fastq.gz
  echo "Combined fastq file created!"
  echo
}
