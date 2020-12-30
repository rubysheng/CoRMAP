#!/bin/bash
#title          :section4.3_matrix.sh
#description    :Change the output format to be a readable matrix.
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.2
#usage          :./section4.3_matrix.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================


# SOURCE THE "PRJNA_PATH" !!!!!!!

####    ??????????????
# where is Trinity.fasta.gene_trans_map
#     now change the outdir for last step(quantification) to ${PRJNA_PATH}/transcripts_count/
#     so where should this map in?????
#     map not generate in the "trinity" step
# --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \?

####################################################################################
#
# Usage:  /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/abundance_estimates_to_matrix.pl --est_method <method>  sample1.results sample2.results ...
#
#      or  /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/abundance_estimates_to_matrix.pl --est_method <method> --quant_files file.listing_target_files.txt
#
#      Note, if only a single input file is given, it's expected to contain the paths to all the target abundance estimation files.
#
# Required:
#
#  --est_method <string>           RSEM|eXpress|kallisto|salmon  (needs to know what format to expect)
#
#  --gene_trans_map <string>           the gene-to-transcript mapping file. (if you don't want gene estimates, indicate 'none'.
#
#
# Options:
#
#  --cross_sample_norm <string>         TMM|UpperQuartile|none   (default: TMM)
#
#  --name_sample_by_basedir             name sample column by dirname instead of filename
#      --basedir_index <int>            default(-2)
#
#  --out_prefix <string>                default: value for --est_method
#
#  --quant_files <string>              file containing a list of all the target files.
#
######################################################################################


####################################################################################
#
# Usage: /home/lewis/Trinity/trinityrnaseq-Trinity-v2.8.5/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl tpm.matrix
#
######################################################################################



#####################
# expression matrix #
#####################

function expressionmx () {
  # create the expression matrix
  cd ${PRJNA_PATH}/transcripts_count/
  ls */RSEM.isoforms.results > isoforms.quant_files.txt
  $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method RSEM \
    --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \
    --cross_sample_norm TMM \
    --out_prefix rsem-gene \
    --name_sample_by_basedir \
    --quant_files isoforms.quant_files.txt
  ## output files (in ${PRJNA_PATH}/transcripts_count):
  #    rsem-gene.isoform.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
  #    rsem-gene.isoform.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
  #    rsem-gene.isoform.TMM.EXPR.matrix : a matrix of TMM-normalized expression values

  # count Numbers of Expressed Transcripts
  $TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
          rsem-gene.isoform.TPM.not_cross_norm | tee rsem-gene.isoform.TPM.not_cross_norm.counts_by_min_TPM

  # source a R script for visualization of transctipts counts
  #Rscript ../quantify.R
  cd ${PRJNA_PATH}

}
