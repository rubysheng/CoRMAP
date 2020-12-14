#!/bin/bash
#title          :section3.2_trim.sh
#description    :Quality control, and generate the combined quality check report.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.2
#usage          :./section3.2_trim.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

########
# trim #
########

# test layout type:
export SE_N=$[`ls *.fastq.gz|grep -v "_"| wc -l|cut -f1 -d' '`]
export PE_L_N=$[`ls *.fastq.gz| grep -c "_1"|cut -f1 -d' '`]
export PE_R_N=$[`ls *.fastq.gz| grep -c "_2"|cut -f1 -d' '`]



function trim_sr () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_fastqc/
  mkdir ./trim/SR/
  echo ===================================
  echo "there are ${SE_N} single-end reads"
  echo ===================================
  # trim
  $TRIMGALORE_HOME/trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/SR/ `ls *.fastq.gz | grep -v "_"` &> trim_sr.log
  echo " trim_galore ended as $(date)"
  echo ==== End Trimming ====
  echo
}


function trim_pe () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_fastqc/
  mkdir ./trim/PE/
  echo =====================================
  echo "there are ${PE_L_N} paired-end reads"
  echo =====================================
  # create the list of file names
  ls -1 *_1.fastq.gz | cat | sed 's/_1.fastq.gz//g' | cat > pairname.txt
  FILE=pairname.txt
  # trim each pair
  while IFS= read -r line; do
      Fq1n="${line}_1.fastq.gz"
      Fq2n="${line}_2.fastq.gz"
      echo ${Fq1n}
      echo ${Fq2n}
      $TRIMGALORE_HOME/trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/PE/ --paired ./${Fq1n} ./${Fq2n} &> trim_pe.log
      echo "${line} finished trimming"
  done < "$FILE"
  echo " trim_galore ended as $(date)"
  echo ==== End Trimming ====
  echo
}
