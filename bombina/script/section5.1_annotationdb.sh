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

########## makeblastdb #########################################################
################################################################################
# USAGE
#   makeblastdb [-h] [-help] [-in input_file] [-input_type type]
#     -dbtype molecule_type [-title database_title] [-parse_seqids]
#     [-hash_index] [-mask_data mask_data_files] [-mask_id mask_algo_ids]
#     [-mask_desc mask_algo_descriptions] [-gi_mask]
#     [-gi_mask_name gi_based_mask_names] [-out database_name]
#     [-blastdb_version version] [-max_file_sz number_of_bytes]
#     [-logfile File_Name] [-taxid TaxID] [-taxid_map TaxIDMapFile] [-version]
#
# DESCRIPTION
#    Application to create BLAST databases, version 2.9.0+
#
# REQUIRED ARGUMENTS
#  -dbtype <String, `nucl', `prot'>
#    Molecule type of target db
#
# OPTIONAL ARGUMENTS
#  -h
#    Print USAGE and DESCRIPTION;  ignore all other parameters
#  -help
#    Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters
#  -version
#    Print version number;  ignore other arguments
#
#  *** Input options
#  -in <File_In>
#    Input file/database name
#    Default = `-'
#  -input_type <String, `asn1_bin', `asn1_txt', `blastdb', `fasta'>
#    Type of the data specified in input_file
#    Default = `fasta'
#
#  *** Configuration options
#  -title <String>
#    Title for BLAST database
#    Default = input file name provided to -in argument
#  -parse_seqids
#    Option to parse seqid for FASTA input if set, for all other input types
#    seqids are parsed automatically
#  -hash_index
#    Create index of sequence hash values.
#
#  *** Sequence masking options
#  -mask_data <String>
#    Comma-separated list of input files containing masking data as produced by
#    NCBI masking applications (e.g. dustmasker, segmasker, windowmasker)
#  -mask_id <String>
#    Comma-separated list of strings to uniquely identify the masking algorithm
#     * Requires:  mask_data
#     * Incompatible with:  gi_mask
#  -mask_desc <String>
#    Comma-separated list of free form strings to describe the masking algorithm
#    details
#     * Requires:  mask_id
#  -gi_mask
#    Create GI indexed masking data.
#     * Requires:  parse_seqids
#     * Incompatible with:  mask_id
#  -gi_mask_name <String>
#    Comma-separated list of masking data output files.
#     * Requires:  mask_data, gi_mask
#
#  *** Output options
#  -out <String>
#    Name of BLAST database to be created
#    Default = input file name provided to -in argumentRequired if multiple
#    file(s)/database(s) are provided as input
#  -blastdb_version <Integer, 4..5>
#    Version of BLAST database to be created
#    Default = `4'
#  -max_file_sz <String>
#    Maximum file size for BLAST database files
#    Default = `1GB'
#  -logfile <File_Out>
#    File to which the program log should be redirected
#
#  *** Taxonomy options
#  -taxid <Integer, >=0>
#    Taxonomy ID to assign to all sequences
#     * Incompatible with:  taxid_map
#  -taxid_map <File_In>
#    Text file mapping sequence IDs to taxonomy IDs.
#    Format:<SequenceId> <TaxonomyId><newline>
#     * Requires:  parse_seqids
#     * Incompatible with:  taxid
## NOTE: which makeblastdb
#        /media/lewis/Workstation_Data2/Bombina_Analysis/NCBI_Database/ncbi-blast-2.9.0+/bin/makeblastdb



########## Trinotate #########################################################
################################################################################
# usage: ./Trinotate <sqlite.db> <command> <input> [...]
#
#      <commands>:
#
#          * Initial import of transcriptome and protein data:
#
#              init --gene_trans_map <file> --transcript_fasta <file> --transdecoder_pep <file>
#
#          * Transdecoder protein search results:
#
#              LOAD_swissprot_blastp <file>
#              LOAD_pfam <file>
#              LOAD_tmhmm <file>
#              LOAD_signalp <file>
#
#           * Trinity transcript search results:
#
#              LOAD_swissprot_blastx <file>
#              LOAD_rnammer <file>
#
#
#           * Load custom blast results using any searchable database
#
#
#              LOAD_custom_blast --outfmt6 <file> --prog <blastp|blastx> --dbtype <database_name>
#
#
#           * report generation:
#
#              report [ -E (default: 1e-5) ] [--pfam_cutoff DNC|DGC|DTC|SNC|SGC|STC (default: DNC=domain noise cutoff)]


##############
# annotation #
##############

## we are in ${MAIN}/annotation/

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

# http://131.104.148.35:3000/cgi-bin/index.cgi
$TRINOTATE_HOME/run_TrinotateWebserver.pl 3000
/usr/sbin/lighttpd -D -f /media/lewis/Seagate_Backup_Plus_Drive/ruby/applications/Trinotate-Trinotate-v3.2.1/TrinotateWeb.conf/lighttpd.conf.port3000
2020-04-14 21:03:11: (log.c.164) server started
Can't locate CGI.pm in @INC (you may need to install the CGI module) (@INC contains: /home/lewis/miniconda3/lib/site_perl/5.26.2/x86_64-linux-thread-multi /home/lewis/miniconda3/lib/site_perl/5.26.2 /home/lewis/miniconda3/lib/5.26.2/x86_64-linux-thread-multi /home/lewis/miniconda3/lib/5.26.2 .) at /media/lewis/Seagate_Backup_Plus_Drive/ruby/applications/Trinotate-Trinotate-v3.2.1/TrinotateWeb/cgi-bin/index.cgi line 7.
BEGIN failed--compilation aborted at /media/lewis/Seagate_Backup_Plus_Drive/ruby/applications/Trinotate-Trinotate-v3.2.1/TrinotateWeb/cgi-bin/index.cgi line 7.
