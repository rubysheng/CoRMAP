#!/bin/bash
#title          :section3.4_normalizeBOTH.sh
#description    :Normalize all reads through the dataset and generate normalized file(s)
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
###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( 'fq' or 'fa')
#  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for
#                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
#
#
#  --max_cov <int>         :targeted maximum coverage for reads.
#
#
#  If paired reads:
#      --left  <string>    :left reads   (if specifying multiple files, list them as comma-delimited. eg. leftA.fq,leftB.fq,...)
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
#
####################################
##  Misc:  #########################
#
#  --pairs_together                :process paired reads by averaging stats between pairs and retaining linking info.
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( "pwd" )
#
#  --CPU <int>                     :number of threads to use (default: = 2)
#  --PARALLEL_STATS                :generate read stats in parallel for paired reads
#
#  --KMER_SIZE <int>               :default 25
#
#  --max_CV <int>                   :maximum coeff of var (default: 10000)
#
#  --min_cov <int>                 :minimum kmer coverage for a read to be retained (default: 0)
#
#  --no_cleanup                    :leave intermediate files
#  --tmp_dir_name <string>         default("tmp_normalized_reads");
#
###############################################################################

##################
# normalization #
##################

function normalize_sr () {
  mkdir ${PRJNA_PATH}/normalization/
  OUTDIR="${PRJNA_PATH}/normalization/"
  echo ==== Trinity In silico Read Normalization START ====
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 50G --max_cov 50 \
                                 --single `ls -m *.fq* |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 16 --output ${OUTDIR} --PARALLEL_STATS
  echo ==== Trinity In silico Read Normalization END ====
}

function normalize_pe () {
  mkdir ${PRJNA_PATH}/normalization/
  OUTDIR="${PRJNA_PATH}/normalization/"
  echo ==== Trinity In silico Read Normalization START ====
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 50G --max_cov 50 \
                                 --left `ls -m *1_val_1_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --right `ls -m *2_val_2_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 16 --output ${OUTDIR} --PARALLEL_STATS
  echo ==== Trinity In silico Read Normalization END ====
}



################################
# both- prepare before trinity #
################################

function pretrinity_both () {
  echo ====Start Normalization Respectively====
  cd ./trim/SR/
  SRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize single-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
    --single ${SRFILES} --CPU 16
  echo "end of single-end sequencing data normalization"
  echo
  cd ../PE/
  LEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  RIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize paired-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
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
