#!/bin/bash
#title          :section2.3_annotation_represent.sh
#description    :Generate the local database for annotation.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.3_annotation_represent.sh
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

##############
# annotation #
##############

## we are in ${MAIN}/annotation/
mkdir annotation/
## only gene annotation not include protein prediction

## obtain annotation databases by Trinotate
$TRINOTATE_HOME/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Ruby_transpip
## Output:
##     ${MAIN}/annotation/Ruby_transpip.sqlite
##     ${MAIN}/annotation/uniprot_sprot.pep
##     ${MAIN}/annotation/Pfam-A.hmm.gz

## Prepare the protein database for blast searches
#    1. from Trinotate database => uniprot_sprot.pep
makeblastdb -in uniprot_sprot.pep -dbtype prot -out swissprot_TRIN -parse_seqids -blastdb_version 5


#    2. from UniProtKB => Reviewed (Swiss-Prot) database
#        ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot_UPKB -parse_seqids -blastdb_version 5


#    3. from NCBI FTP server => BLAST-specific directory => Preformatted BLAST database files => swissprot :Protein sequences from the swiss-prot sequence database (last major update).
#        ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
# NOTE: "nr" is non-redundant protein, "nt" is non-redundant nucleotide
update_blastdb.pl --decompress swissprot
blastdbcmd -entry all -db swissprot -dbtype prot -out swissprot_NCBI.pep.fasta
makeblastdb -in swissprot_NCBI.pep.fasta -dbtype prot -out swissprot_NCBI -parse_seqids -blastdb_version 5


#### Uncompress and prepare the Pfam database for use with 'hmmscan' like so:
# gunzip ${MAIN}/annotation/Pfam-A.hmm.gz
# hmmpress ${MAIN}/annotation/Pfam-A.hmm


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
  blastx -query ../RDrat_perg.dnaseq.fasta -db swissprot_TRIN -num_threads 32 \
    -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_TRIN.outfmt6


  blastx -query ../RDrat_perg.dnaseq.fasta -db swissprot_UPKB -num_threads 32 \
    -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_UPKB.outfmt6


  blastx -query ../RDrat_perg.dnaseq.fasta -db swissprot_NCBI -num_threads 32 \
    -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_NCBI.outfmt6

  # blastx -query ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta -db swissprot_TRIN -num_threads 32 \
  #   -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_TRIN.outfmt6
  #
  # blastx -query ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta -db swissprot_UPKB -num_threads 32 \
  #   -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_UPKB.outfmt6
  #
  # blastx -query ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta -db swissprot_NCBI -num_threads 32 \
  #   -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > blastx_NCBI.outfmt6

  # generate gene_trans_map
  sed -ibak 's/[.]/_/g' ../RDrat_perg.dnaseq.fasta
  sed -ibak 's/[.]/_/g' ../RDrat_perg.pepseq.fasta
  $TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl ../RDrat_perg.dnaseq.fasta > ../RDrat_perg.gene_trans_map

  # Add the TransDecoder prediction details to the HEADER
  grep '>' ../RDrat_perg.pepseq.fasta | cut -d'_' -f3 | sort | uniq  # check which datasets are included
    # PRJNA390522
    # PRJNA451011
  grep 'PRJNA390522' ../RDrat_perg.pepseq.fasta > ../RD_PRJNA390522_pephead.lst
  grep 'PRJNA451011' ../RDrat_perg.pepseq.fasta > ../RD_PRJNA451011_pephead.lst
  # sed -ibak 's/^[^_]*_//' ../RD_PRJNA*.lst
  sed -ibak 's/_/\t/1; s/_/./g; s/[.]/_/1' ../RD_PRJNA*_pephead.lst

  # in each of the dataset
  # PE_NAME="${FA_NEWNAME}.transdecoder.pep"
  # PE_NEWNAME="${file}_pep.fasta"
  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
     /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA390522/PRJNA390522_RSEM.fasta.transdecoder.pep > ../tmp.fa
    # PRJNA390522_RSEM.fasta.transdecoder.pep > ./tmp.fa
  printf "\n" >> ../tmp.fa
  # sed -i 's/[.]/_/g' ../tmp.fa

  while IFS= read -r line; do
    group_num=`echo $line | awk '{print $1}'`
    pattern=`echo $line | awk '{print $2}'`
    echo ${group_num}
    echo ${pattern}
    awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line  # > ./RD_PRJNA390522.transdecoder.pep   #../tmp.fa >> ../RD_PRJNA390522.transdecoder.pep
    group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
    pattern=${pattern%.i*}
    pattern_new=`echo $pattern | sed 's/[.]/_/g'`
    c="${group_num}_${pattern_new}"
    sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
    cat ../tmp_line >> ../RD_PRJNA390522.transdecoder.pep
    rm ../tmp_line
  done < ../RD_PRJNA390522_pephead.lst ##RD_PRJNA390522_pephead.lst

  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
     /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA451011/PRJNA451011_RSEM.fasta.transdecoder.pep > ../tmp.fa
  printf "\n" >> ../tmp.fa

  while IFS= read -r line; do
    group_num=`echo $line | awk '{print $1}'`
    pattern=`echo $line | awk '{print $2}'`
    echo ${group_num}
    echo ${pattern}
    awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
    group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
    pattern=${pattern%.i*}
    pattern_new=`echo $pattern | sed 's/[.]/_/g'`
    c="${group_num}_${pattern_new}"
    sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
    cat ../tmp_line >> ../RD_PRJNA451011.transdecoder.pep
    rm ../tmp_line
  done < ../RD_PRJNA451011_pephead.lst

  cat ../RD_PRJNA*.transdecoder.pep > RDrat_perg.merge.transdecoder.pep
  # Combine to SQL
  $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite init \
    --gene_trans_map ../RDrat_perg.gene_trans_map \
    --transcript_fasta  ../RDrat_perg.dnaseq.fasta \
    --transdecoder_pep ../RDrat_perg.merge.transdecoder.pep


  # $TRINOTATE_HOME/Trinotate ${SQL_NAME} init \
  #   --gene_trans_map ${PRJNA_PATH}/transcripts_count/Trinity.fasta.gene_trans_map \
  #   --transcript_fasta  ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta \
  #   --transdecoder_pep ./Trinity.fasta.transdecoder_dir/longest_orfs.pep

  # Loading BLAST homologies: Load transcript hits
  # $TRINOTATE_HOME/Trinotate ${SQL_NAME} LOAD_swissprot_blastx blastx.outfmt6
  $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastx blastx_TRIN.outfmt6


  cd ${PRJNA_PATH}
}

function predict_pro () {
  annotation

  # Search Transdecoder-predicted proteins
  blastp -query ../RDrat_perg.pepseq.fasta \
    -db swissprot_TRIN \
    -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_TRIN.outfmt6

  blastp -query ../RDrat_perg.pepseq.fasta \
    -db swissprot_UPKB \
    -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_UPKB.outfmt6

  blastp -query ../RDrat_perg.pepseq.fasta \
    -db swissprot_NCBI \
    -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_NCBI.outfmt6


  # Run hmmscan
  hmmscan --cpu 16 --domtblout TrinotatePFAM.out \
    Pfam-A.hmm ./Trinity.fasta.transdecoder_dir/longest_orfs.pep  > pfam.log

  # Loading BLAST homologies: Load protein hits
  $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastp blastp_TRIN.outfmt6
  # Trinotate ${SQL_NAME} LOAD_swissprot_blastp blastp.outfmt6
  # Load Pfam domain entries
  $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_pfam TrinotatePFAM.out
  # Trinotate ${SQL_NAME} LOAD_pfam TrinotatePFAM.out

}



head -4 {TrinotatePFAM.out,blastp_TRIN.outfmt6,blastx_TRIN.outfmt6,../RDrat_perg.gene_trans_map,../RDrat_perg.dnaseq.fasta,../RDrat_perg.merge.transdecoder.pep}

sed -ibak 's/[.]D/_D/g; s/[.]c/_c/g; s/[.]g/_g/g; s/[.]i/_i/g; s/[.]p/_p/g' {TrinotatePFAM.out,blastp_TRIN.outfmt6,blastx_TRIN.outfmt6}

head -6 {TrinotatePFAM.out,blastp_TRIN.outfmt6,blastx_TRIN.outfmt6}

$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastp blastp_TRIN.outfmt6
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_pfam TrinotatePFAM.out
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastx blastx_TRIN.outfmt6
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite report  > trinotate_annotation_report.xls
