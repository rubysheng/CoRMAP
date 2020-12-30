#!/bin/bash
###############################################################################
###############################################################################
# Generate basis statistics of the assembled transcriptome for all project
# after intra-taxa orthologs searching.
# Require using a Trinity utility script TrinityStats.pl.
###############################################################################
###############################################################################

while IFS= read -r project; do
  echo ${project}
  file_path="./"${project}"/trinity_out_dir/Trinity.fasta"
  $TRINITY_HOME/util/TrinityStats.pl  ${file_path}
done <  /home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/analyze/expre_value/pro_acc
