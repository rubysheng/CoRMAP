#!/bin/bash
#title          :section4.2_quantify.sh
#description    :Transcripts quantification.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.2
#usage          :./section4.2_quantify.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================


# SOURCE THE "PRJNA_PATH" !!!!!!!




##################
# quantification #
##################

#### usage guidence :   $TRINITY_HOME/util/align_and_estimate_abundance.pl ####
########################
#  Essential parameters:
########################
#
#  --transcripts <string>           transcript fasta file
#
#  --seqType <string>               fq|fa
#
#  If Paired-end:
#
#     --left <string>
#     --right <string>
#
#   or Single-end:
#
#      --single <string>
#
#   or (preferred):
#
#      --samples_file <string>    tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                      # if single-end instead of paired-end, then leave the 4th column above empty.
#
#
#
#  --est_method <string>           abundance estimation method.
#                                        alignment_based:  RSEM
#                                        alignment_free: kallisto|salmon
#
###################################
#  Potentially optional parameters:
###################################
#
# --output_dir <string>            write all files to output directory
#                                  (note, if using --samples_file, output_dir will be set automatically according to replicate name))
#
#
#  if alignment_based est_method:
#       --aln_method <string>            bowtie|bowtie2 alignment method.  (note: RSEM requires either bowtie or bowtie2)
#
###########
# Optional:
# #########
#
# --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
#
# --samples_idx <int>               restricte processing to sample entry (index starts at one)
#
#
# --thread_count                   number of threads to use (default = 4)
#
# --debug                          retain intermediate files
#
#  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
#     or
#  --trinity_mode                   Setting --trinity_mode will automatically generate the gene_trans_map and use it.
#
#
#  --prep_reference                 prep reference (builds target index)
#
#
########################################
#
#  Parameters for single-end reads:
#
#  --fragment_length <int>         specify RNA-Seq fragment length (default: 200)
#  --fragment_std <int>            fragment length standard deviation (defalt: 80)
#
########################################
#
#   bowtie-related parameters: (note, tool-specific settings are further below)
#
#  --max_ins_size <int>             maximum insert size (bowtie -X parameter, default: 800)
#  --coordsort_bam                  provide coord-sorted bam in addition to the default (unsorted) bam.
#
########################################
#  RSEM opts:
#
#  --bowtie_RSEM <string>          if using 'bowtie', default: "--all --best --strata -m 300 --chunkmbs 512"
#  --bowtie2_RSEM <string>         if using 'bowtie2', default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
#                                ** if you change the defaults, specify the full set of parameters to use! **
#
#  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
#  --rsem_add_opts <string>        additional parameters to pass on to rsem-calculate-expression
#
##########################################################################
#  kallisto opts:
#
#  --kallisto_add_opts <string>  default:
#
##########################################################################
#
#  salmon opts:
#
#  --salmon_idx_type <string>    quasi|fmd (defalt: quasi)
#  --salmon_add_opts <string>    default:
#
#
#  Example usage
#
#   ## Just prepare the reference for alignment and abundance estimation
#
#    /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
#
#   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
#
#    /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
#
##  ## prep the reference and run the alignment/estimation
#
#    /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir
#
#   ## Use a samples.txt file:
#
#    /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode --samples_file samples.txt --seqType fq
#
#########################################################################

## sample_file_*.txt :
##     MB-na�ve	MB-na�ve_rep1	./SRR7299147_trimmed_renamed.fq.gz
##   or :(bombina)
##     Control_min3to9H	Control_min3to9H_rep1	./SRR6320648_1_val_1_renamed.fq.gz	./SRR6320648_2_val_2_renamed.fq.gz

### now we are in ${PRJNA_PATH}

function count () {
  # create a transcripts quantification folder to hold the results
  mkdir ${PRJNA_PATH}/transcripts_count
  #cd ${PRJNA_PATH}/transcripts_count

  # move input files to this folder
  find . -type f -name *renamed.fq.gz -exec cp -v {} ./transcripts_count \;
  find . -name sample_file_*.txt -exec cp {} ./transcripts_count \;

  # go to ${PRJNA_PATH}/transcripts_count
  cd ${PRJNA_PATH}/transcripts_count

  $TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq \
      --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta \
      --est_method RSEM \
      --aln_method bowtie2 --trinity_mode --prep_reference \
      --samples_file sample_file_*.txt

  # create a folder to hold processed(trim and rename) data
  mkdir ${PRJNA_PATH}/processed_data
  # move renamed files to the folder
  mv *renamed.fq.gz ${PRJNA_PATH}/processed_data

  # move sample_file to the ${PRJNA_PATH}
  mv sample_file_*.txt ${PRJNA_PATH}

  cd ${PRJNA_PATH}
}
