#!/bin/bash
#title          :section5.1_annotationdb.sh
#description    :Generate the local database for annotation.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.1
#usage          :./section5.1_annotationdb.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================


# SOURCE THE "PRJNA_PATH" !!!!!!!




##############
# annotation #
##############
  #Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
  #makeblastdb -in uniprot_sprot.pep -dbtype prot
  #gunzip Pfam-A.hmm.gz
  #hmmpress Pfam-A.hmm
  #/media/heyland-lab/"Seagate Backup Plus Drive"/ruby/annotation/

function annotation () {
  cd ${PRJNA_PATH}
  cp /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/annotation/* .
  TransDecoder.LongOrfs -t ./trinity_out_dir/Trinity.fasta
  #Search Trinity transcripts
  blastx -query ./trinity_out_dir/Trinity.fasta -db uniprot_sprot.pep -num_threads 8 \
    -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
  #Search Transdecoder-predicted proteins
  blastp -query ./Trinity.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.pep \
    -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
  #Run hmmscan
  hmmscan --cpu 5 --domtblout TrinotatePFAM.out \
    Pfam-A.hmm ./Trinity.fasta.transdecoder_dir/longest_orfs.pep  > pfam.log

}
