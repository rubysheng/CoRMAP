#!/usr/bin/env bash
singularity exec -e Trinity.simg  Trinity \
          --seqType fq \
          --left `pwd`/reads.left.fq.gz  \
          --right `pwd`/reads.right.fq.gz \
          --max_memory 1G --CPU 4 \
          --output `pwd`/trinity_out_dir

singularity exec -e Trinity.simg /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq \
              --est_method RSEM \
              --aln_method bowtie --trinity_mode --prep_reference \
              --samples_file

singularity exec -e /data/Trinity.simg /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl \
  --transcripts `pwd`/../../test_DATA/Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode \
  --prep_reference --samples_file `pwd`/samples.txt --SS_lib_type R --seqType fq
touch test_RSEM



# ## In the rodent_lst/ directory
# study_lst=`ls -1 ./input/*_pep.fasta`
# for file in `ls -1 ./input/*_pep.fasta`; do
#   file=`basename $file`
#   sp_name=${file%.fasta}
#   echo "========================================================================="
#   echo ${sp_name%_pep}" is processing with extracting sequences from orthologs."
#   rm -v ./analyze/${sp_name%_pep}.dnaseq.fasta ./analyze/${sp_name%_pep}.pepseq.fasta

while IFS= read -r line; do
  group_num=`echo $line | cut -d' ' -f1`
  pepname=`echo $line | cut -d' ' -f3`
  awk -v gp="$group_num" -v pat="$pepname" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./PRJNA529794_pep.fasta >> ./PRJNA529794.pepseq.fasta
  isoform=${pepname%.p*}
  awk -v gp="$group_num" -v pat="$isoform" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./PRJNA529794_RSEM.fasta >> ./PRJNA529794.dnaseq.fasta
done < ./clustered_PRJNA529794.lst
#
#   echo ${sp_name%_pep}" is end with processing of extracting sequences from orthologs."
#   echo "========================================================================="
# done
