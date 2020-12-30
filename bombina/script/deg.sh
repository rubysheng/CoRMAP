#!/usr/bin/env bash

scp -r ruby@graham.computecanada.ca:/home/ruby/projects/def-heylanda/ruby/all_taxa/expression/deg/ .

echo 'PRJNA252803
PRJNA287145
PRJNA287152
PRJNA302146
PRJNA316996
PRJNA381689
PRJNA390522
PRJNA419677
PRJNA422916
PRJNA450614
PRJNA451011
PRJNA475804
PRJNA529794
PRJNA541005' > lst


while IFS= read -r project; do
  input="./deg/count_matrix_"${project}".txt"
  outdir=${project}"_dir"
  sample="./deg/samples_described_"${project}".txt"

  $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
          --matrix ${input} \
          --method DESeq2 \
          --output ${outdir} \
          --samples_file ${sample}

  cd ${outdir}
  input="../deg/count_matrix_"${project}".txt"
  $TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ${input} -P 1e-1 -C 1
  cd ..
done <  lst

source ./dea.sh 2>&1 | tee DEA.log


$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
          --matrix ./deg/count_matrix_Bombina.txt \
          --method DESeq2 \
          --output Bombina_dir \
          --samples_file ./deg/samples_described_Bombina.txt

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl \
          --matrix ./deg/count_matrix_PRJNA541005.txt \
          --method DESeq2 \
          --output PRJNA541005_dir \
          --samples_file ./deg/samples_described_PRJNA541005.txt

          
while IFS= read -r line; do
  awk -v pat="${line}" '{ if ($0 ~ pat) {print} }' glmLRT_*.csv
done < uniq_dge_group.lst

while IFS= read -r line; do
  awk -v pat="${line}_" '{ if ($0 ~ pat) {print} }' glmLRT_*.csv > tmp
  BO_perg_ave=`grep "BO" tmp | cut -d"," -f 1,2 | sed "s/,/ /" | awk '{ total += $2; count++ } END { print total/count }'`
  FI_perg_ave=`grep "FI" tmp | cut -d"," -f 1,2 | sed "s/,/ /" | awk '{ total += $2; count++ } END { print total/count }'`
  IN_perg_ave=`grep "IN" tmp | cut -d"," -f 1,2 | sed "s/,/ /" | awk '{ total += $2; count++ } END { print total/count }'`
  RD_perg_ave=`grep "RD" tmp | cut -d"," -f 1,2 | sed "s/,/ /" | awk '{ total += $2; count++ } END { print total/count }'`
  RW_perg_ave=`grep "RW" tmp | cut -d"," -f 1,2 | sed "s/,/ /" | awk '{ total += $2; count++ } END { print total/count }'`
  echo "group:"${line}  "BO:"${BO_perg_ave}  "FI:"${FI_perg_ave}  "IN:"${IN_perg_ave}  "RD:"${RD_perg_ave}  "RW:"${RW_perg_ave} >> uniq_dge_group_ave_logfc.csv
done < head10.lst
# group_num=`echo $line | cut -d'_' -f 3 | cut -d'.' -f 1`
# input_file="allsp_group_"$group_num".pepseq.fasta"
# awk -v pat=">RD" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | sed "s/>RD/>g${group_num}_RD/g" >> ../RD_perg.pepseq.fasta
