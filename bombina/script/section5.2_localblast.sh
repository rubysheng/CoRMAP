#!/bin/bash
#title          :section5.2_localblast.sh
#description    :Generate the local database for blast.
#author         :Ruby(Yiru) Sheng
#date           :20191118
#version        :1.1
#usage          :source $RUBY_SCRIPTS/section5.2_localblast.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

function local_blastdb() {
  PREFIX=`basename $(pwd)`
  mkdir blast_out_dir
  NAME="${PREFIX}_Trinity.fasta.RSEM.transcripts.fa"
  sed "s/TRINITY/${PREFIX}/g" ./transcripts_count/Trinity.fasta.RSEM.transcripts.fa > ./blast_out_dir/${NAME}
  
  # convert Fasta file to the format of a table
  # awk -f '$RUBY_SCRIPTS/FaToTb' *RSEM.transcripts.fa > transcripts_seq.tb
  # source a R script for visualization of transctipts counts
  # Rscript $RUBY_SCRIPTS/load_db.R
  # cd ..
  cd ./blast_out_dir 
  OUT_NAME="${PREFIX}_blast"
  makeblastdb -in ${NAME} -parse_seqids -hash_index -blastdb_version 5 -title "${PREFIX}" -dbtype nucl -out ${OUT_NAME}
  DBNAME="${OUT_NAME}"
  
  # if [ -e "/media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list" ]; then
  #   cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  # else
  #   find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -type f -name *_Trinity.fasta.RSEM.transcripts.fa > /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list
  #   cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  # fi
  
  # while read line
  # do
  #   echo "${line}"
  #   INPUT="${line}_Trinity.fasta.RSEM.transcripts.fa"
  
  #   OUTPUT="bln_q${line}_d${PREFIX}.outfmt6"
  #   time blastn -query ${INPUT} -db ${DBNAME} -evalue 1e-3 -outfmt 6 -num_threads 6 -out ${OUTPUT}
  #   wc -l ${OUTPUT}
  # 
  # done < query_file_loc.list


}

#local_blastdb
cd ..
PREFIX=`basename $(pwd)`
echo "${PREFIX}"
cd blast_out_dir
OUT_NAME="${PREFIX}_blast"
DBNAME="${OUT_NAME}"

if [ -e "../../query_file_loc.list" ]; then
    cp -v ../../query_file_loc.list .
else
    cd ../..
    find /media/lewis/New_Seagate_Drive_8TB/ruby/testset -type f -name *_Trinity.fasta.RSEM.transcripts.fa > query_file_loc.list
    echo "${PREFIX}"
    cd ${PREFIX}/blast_out_dir
    cp -v ../../query_file_loc.list .
    
fi

  # if [ -e "/media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list" ]; then
  #   cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  # else
  #   find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -type f -name *_Trinity.fasta.RSEM.transcripts.fa > /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list
  #   cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  # fi


while read line
do
  echo "${line}"
  INPUT="${line}" #path/to/*_Trinity.fasta.RSEM.transcripts.fa
  FILE=`basename ${line}`
  echo "${FILE}"
  
  OUTPUT="bln_q${FILE%_Trinity.fasta.RSEM.transcripts.fa}_d${PREFIX}.outfmt6"
  echo "${OUTPUT}"

  time blastn -query ${INPUT} -db ${DBNAME} -evalue 1e-3 -outfmt 6 -num_threads 6 -out ${OUTPUT}
  wc -l ${OUTPUT}

done < query_file_loc.list

