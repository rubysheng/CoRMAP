#!/bin/bash
#title          :section2.3_annotation_represent.sh
#description    :Generate the local database for annotation.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.3_annotation_represent.sh
##usage          : 1. source $CMRP_PATH/script/section2.3_annotation_represent.sh
#                  2. source /media/lewis/New_Seagate_Drive_8TB/ruby/CMRP/script/section2.3_annotation_represent.sh --source-only ; load_annodb 2>&1 | tee load_annodb.log
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

##############
# annotation #
##############

function load_annodb() {
  # in the rodent_lst/
  # go to annotation directory
  # in the rodent_lst/annotation
  if [ ! -d annotation/ ]; then
    mkdir annotation
  fi
  cd annotation

  # create the outdir to hold results
  if [ ! -d anno_output/ ]; then
    mkdir anno_output
  fi
  # in the rodent_lst/annotation/anno_output
  cd anno_output

  ## obtain annotation databases by Trinotate
  $TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  CMRP_transpip
  ## Output:
  ##     $CMRP_PATH/sample/orthologs/annotation/anno_output/CMRP_transpip.sqlite
  ##     $CMRP_PATH/sample/orthologs/annotation/anno_output/uniprot_sprot.pep
  ##     $CMRP_PATH/sample/orthologs/annotation/anno_output/Pfam-A.hmm.gz

  ## Prepare the protein database for blast searches
  #    (default) from Trinotate database => uniprot_sprot.pep
  makeblastdb -in uniprot_sprot.pep -dbtype prot -out swissprot_TRIN -parse_seqids -blastdb_version 5

  # other options
  # #    1. from UniProtKB => Reviewed (Swiss-Prot) database
  # #        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  # wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
  # gzip -d uniprot_sprot.fasta.gz
  # makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot_UPKB -parse_seqids -blastdb_version 5
  # #    2. from NCBI FTP server => BLAST-specific directory => Preformatted BLAST database files => swissprot :Protein sequences from the swiss-prot sequence database (last major update).
  # #        ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
  # # * "nr" is non-redundant protein, "nt" is non-redundant nucleotide
  # update_blastdb.pl --decompress swissprot
  # blastdbcmd -entry all -db swissprot -dbtype prot -out swissprot_NCBI.pep.fasta
  # makeblastdb -in swissprot_NCBI.pep.fasta -dbtype prot -out swissprot_NCBI -parse_seqids -blastdb_version 5

  #### Uncompress and prepare the Pfam database for use with 'hmmscan' like so:
  gunzip Pfam-A.hmm.gz
  hmmpress Pfam-A.hmm

  # back to rodent_lst/
  cd ../..
}



function annotation () {
  # in the rodent_lst/annotation/anno_output
  cd annotation/anno_output
  $TRINOTATE_HOME/Trinotate CMRP_transpip.sqlite init \
      --gene_trans_map ../RD_perg.gene_trans_map \
      --transcript_fasta  ../RD_uniq_dnaseq.fasta \
      --transdecoder_pep ../RD_uniq_mergepepseq.fasta
  blastx -query ../RD_uniq_dnaseq.fasta -db swissprot_TRIN -num_threads 14 \
      -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > allRD_blastx_TRIN.outfmt6
  $TRINOTATE_HOME/Trinotate CMRP_transpip.sqlite LOAD_swissprot_blastx allRD_blastx_TRIN.outfmt6
  blastp -query ../RD_uniq_mergepepseq.fasta -db swissprot_TRIN  -num_threads 16 \
      -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > allRD_blastp_TRIN.outfmt6
  hmmscan --cpu 16 --domtblout allRD_TrinotatePFAM.out \
    Pfam-A.hmm ../RD_uniq_mergepepseq.fasta > allRD_pfam.log
  $TRINOTATE_HOME/Trinotate CMRP_transpip.sqlite LOAD_swissprot_blastp allRD_blastp_TRIN.outfmt6
  $TRINOTATE_HOME/Trinotate CMRP_transpip.sqlite LOAD_pfam allRD_TrinotatePFAM.out
  $TRINOTATE_HOME/Trinotate CMRP_transpip.sqlite report  > allRD_trinotate_annotation_report.xls
}

function main() {
  load_annodb 2>&1 | tee load_annodb.log
  annotation 2>&1 | tee annotation.log
}

if [ "${1}" != "--source-only" ]; then
    main "${@}"
fi
