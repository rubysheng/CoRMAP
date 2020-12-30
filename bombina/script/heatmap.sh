#!/usr/bin/env bash

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


sudo fallocate -l 1G /swapfile
