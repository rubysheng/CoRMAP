#!/bin/bash
#SBATCH --time=1-10:00:00
#SBATCH --account=def-heylanda
#SBATCH --job-name=extract_seq
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ysheng@uoguelph.ca
#SBATCH --ntasks=1
#SBATCH --output=extract_seq.log
date;hostname;pwd


while IFS= read -r line; do
   group_num=`echo $line | cut -d' ' -f1`
   pepname=`echo $line | cut -d' ' -f3`
   awk -v gp="$group_num" -v pat="$pepname" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./PRJNA529794_pep.fasta >> ./PRJNA529794.pepseq.fasta
   isoform=${pepname%.p*}
   awk -v gp="$group_num" -v pat="$isoform" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./PRJNA529794_RSEM.fasta >> ./PRJNA529794.dnaseq.fasta
done < ./clustered_PRJNA529794.lst

date
