#!/bin/bash

#=======================================================================================================
##############
# annotation #
##############

function annotation () {
  cd ${PRJNA_PATH}
  PRJNA_NAME=`basename $(pwd)`
  SQL_NAME="${PRJNA_NAME}.sqlite"

  # create the outdir to hold results
  mkdir anno_result
  cd anno_result
  # copy the SQL file from basic annotation folder, and change name to be the current dataset Name
  cp ${MAIN}/annotation/Ruby_transpip.sqlite ./

  mv Ruby_transpip.sqlite ${SQL_NAME}

  TransDecoder.LongOrfs -t ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta

  # Search Trinity transcripts
  blastx -query ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta -db ${MAIN}/annotation/uniprot_sprot.pep -num_threads 32 \
    -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx.outfmt6

  # Combine to SQL
  $TRINOTATE_HOME/Trinotate ${SQL_NAME} init \
    --gene_trans_map ${PRJNA_PATH}/transcripts_count/Trinity.fasta.gene_trans_map \
    --transcript_fasta  ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta \
    --transdecoder_pep ./Trinity.fasta.transdecoder_dir/longest_orfs.pep

  # Loading BLAST homologies: Load transcript hits
  $TRINOTATE_HOME/Trinotate ${SQL_NAME} LOAD_swissprot_blastx blastx.outfmt6
  cd ${PRJNA_PATH}
}

function predict_pro () {
  annotation

  # Search Transdecoder-predicted proteins
  blastp -query ./Trinity.fasta.transdecoder_dir/longest_orfs.pep \
    -db ${MAIN}/annotation/uniprot_sprot.pep \
    -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
  # Run hmmscan
  hmmscan --cpu 32 --domtblout TrinotatePFAM.out \
    Pfam-A.hmm ./Trinity.fasta.transdecoder_dir/longest_orfs.pep  > pfam.log

  # Loading BLAST homologies: Load protein hits
  Trinotate ${SQL_NAME} LOAD_swissprot_blastp blastp.outfmt6
  # Load Pfam domain entries
  Trinotate ${SQL_NAME} LOAD_pfam TrinotatePFAM.out

}


function mk_local_blastdb() {
  # in the path/to/the/directory/of/PRJNA*(Bombina*)

  PREFIX=`basename $(pwd)`
  mkdir blast_out_dir
  NAME="${PREFIX}_Trinity.fasta.RSEM.transcripts.fa"
  sed "s/TRINITY/${PREFIX}/g" ./transcripts_count/Trinity.fasta.RSEM.transcripts.fa > ./blast_out_dir/${NAME}

  cd ./blast_out_dir
  OUT_NAME="${PREFIX}_blast"
  makeblastdb -in ${NAME} -parse_seqids -hash_index -blastdb_version 5 -title "${PREFIX}" -dbtype nucl -out ${OUT_NAME}
  DBNAME="${OUT_NAME}"
}

function lcblast() {
  # in the path/to/the/directory/of/PRJNA*(Bombina*)
  PREFIX=`basename $(pwd)`
  #echo "${PREFIX}"
  cd blast_out_dir
  DBNAME="${PREFIX}_blast"

  if [ -e "/media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list" ]; then
    cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  else
    cd ../..
    find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -type f -name *_Trinity.fasta.RSEM.transcripts.fa > /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list
    echo "${PREFIX}"
    cd ${PREFIX}/blast_out_dir
    cp -v /media/lewis/Seagate_Backup_Plus_Drive/ruby/annotation/query_file_loc.list .
  fi

  # set the parameter:
  printf "%s\n%s\n\n" "How is the e-value you want to set?" \
         "(input must be like 1e-3 or 1e-4)"
  read -r E_VALUE

  echo
  echo ==== start the local blast ====
  echo
  echo "The database for blast is : ${DBNAME}"

  while read line
  do
    echo
    echo "The Query file is : ${line}"
    echo "The database for blast is : ${DBNAME}"
    echo
    INPUT="${line}" #path/to/*_Trinity.fasta.RSEM.transcripts.fa
    FILE=`basename ${line}`
    OUTPUT="bln_q${FILE%_Trinity.fasta.RSEM.transcripts.fa}_d${PREFIX}.outfmt6"

    if [[ ${FILE%_Trinity.fasta.RSEM.transcripts.fa} != ${PREFIX} ]]; then
      time blastn -query ${INPUT} -db ${DBNAME} -evalue ${E_VALUE} -outfmt 6 -num_threads 6 -out ${OUTPUT}
      echo
      echo "The output of this local blasting is called ${OUTPUT}"
      echo "The counts of similar sequences:"
      wc -l ${OUTPUT}
    fi


  done < query_file_loc.list

  echo
  echo ==== end the local blast ====
  echo

  # move all blast results the query of which are using the current dataset to this directory
  cd ../..
  find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -type f -name "bln_q${PREFIX}_d*.outfmt6" -exec mv {} /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/${PREFIX}/blast_out_dir/ \;

  cd ${PREFIX}/blast_out_dir

}

# blast_outfmt6:
        # qseqid: query (e.g., gene) sequence id
        # sseqid: subject (e.g., reference genome) sequence id
        # pident: percentage of identical matches
        # length: alignment length
        # mismatch: number of mismatches
        # gapopen: number of gap openings
        # qstart: start of alignment in query
        # qend: end of alignment in query
        # sstart: start of alignment in subject
        # send: end of alignment in subject
        # evalue: expect value
        # bitscore: bit score



#=======================================================================================================
#####################
# make the database #
#####################

function load_sql() {

  cd transcripts_count

  # convert Fasta file to the format of a table
  awk -f '$RUBY_SCRIPTS/FaToTb' ../blast_out_dir/*RSEM.transcripts.fa > transcripts_seq.tb

  # load transcript infomation, expression matrix, and blast results to A SQL DATABASE
  Rscript $RUBY_SCRIPTS/load_db.R
  cd ..

}
