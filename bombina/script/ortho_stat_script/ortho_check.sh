#!/bin/bash
#=======================================================================================================
##########################
# orthomcl results check #
##########################

function ortho_chk_blast() {
  #statements
  # in the "*_lst/analyze/per_group" directory
  # build a local blast database per group
  for dnafasta in `ls -1 *.dnaseq.fasta`; do
    makeblastdb -in ${dnafasta} -dbtype nucl -hash_index -out dna_pergroup
    output_name=${dnafasta%dnaseq.fasta}"blastn.fmt7"
    blastn -query ${dnafasta} -db dna_pergroup -evalue 1e-5 -out  ${output_name} -outfmt "7" -num_alignments 15 -num_threads 15
  done

  for pepfasta in `ls -1 *.pepseq.fasta`; do
    makeblastdb -in ${pepfasta} -dbtype prot -hash_index -out pep_pergroup
    output_name=${pepfasta%pepseq.fasta}"blastp.fmt7"
    blastp -query ${pepfasta} -db pep_pergroup -evalue 1e-5 -out  ${output_name} -outfmt "7" -num_alignments 15 -num_threads 15
  done
  # -out output_file_name
  # -outfmt      
  #  0 = pairwise,
  #  1 = query-anchored showing identities,
  #  2 = query-anchored no identities,
  #  3 = flat query-anchored, show identities,
  #  4 = flat query-anchored, no identities,
  #  5 = XML Blast output,
  #  6 = tabular,
  #  7 = tabular with comment lines,
  #  8 = Text ASN.1,
  #  9 = Binary ASN.1
  # 10 = Comma-separated values
  rm dna_pergroup*
  rm pep_pergroup*
}

ortho_chk_blast 2>&1 | tee ortho_check.log
