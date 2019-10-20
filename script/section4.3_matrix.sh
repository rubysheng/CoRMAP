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




#####################
# expression matrix #
#####################

function expressionmx () {
  # create the expression matrix
  cd ${PRJNA_PATH}/transcripts_count/
  ls */RSEM.genes.results > genes.quant_files.txt
  abundance_estimates_to_matrix.pl --est_method RSEM \
    --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \
    --cross_sample_norm TMM \
    --out_prefix rsem-gene \
    --name_sample_by_basedir \
    --quant_files genes.quant_files.txt
  ## output files:
  #    rsem-gene.isoform.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
  #    rsem-gene.isoform.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
  #    rsem-gene.isoform.TMM.EXPR.matrix : a matrix of TMM-normalized expression values
  # count Numbers of Expressed Transcripts
  count_matrix_features_given_MIN_TPM_threshold.pl \
          rsem-gene.isoform.TPM.not_cross_norm | tee rsem-gene.isoform.TPM.not_cross_norm.counts_by_min_TPM
  # source a R script for visualization of transctipts counts
  #Rscript ../quantify.R
  cd ${PRJNA_PATH}
}
