#!/bin/bash

#=======================================================================================================
############################
# raw data quality control #
############################
function raw_qc () {
  mkdir ./raw_qc/
  # generate quality control with each untrimmed runs
  fastqc -o ./raw_qc/ -f fastq -t 5 --extract *.fastq.gz #SRR*.fastq.gz
  # combine reports of all runs in a study to one quality control report
  multiqc ./raw_qc/ -o ./raw_qc/multiqc_output/

  touch ./raw_qc/multiqc_output/log_list.txt
  ls -l ./raw_qc/multiqc_output/multiqc_data/ > ./raw_qc/multiqc_output/log_list.txt

  # remove the separate fastqc files just keep the multiqc report
  find ./raw_qc -type f -name "SRR*" -exec rm {} \;
  find ./raw_qc -type d -name "SRR*" -exec rm -r {} \;
}
#=======================================================================================================
########
# trim #
########
function trim_sr () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_qc/
  mkdir ./trim/SR/
  echo ===================================
  echo "there are ${SE_N} single-end reads"
  echo ===================================
  # trim
  $TRIMGALORE_HOME/trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/SR/ `ls *.fastq.gz | grep -v "_"` > trim.log
  echo " trim_galore ended as $(date)"
  echo ==== End Trimming ====
  echo
}


function trim_pe () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_qc/
  mkdir ./trim/PE/
  echo =====================================
  echo "there are ${PE_L_N} paired-end reads"
  echo =====================================
  # trim
  ls *_1.fastq.gz | cat | sed 's/_1.fastq.gz//g' | cat > pairname.txt
  # trim each pair
  FILE=pairname.txt
  while IFS= read -r line; do
      Fq1n="${line}_1.fastq.gz"
      Fq2n="${line}_2.fastq.gz"
      echo ${Fq1n}
      echo ${Fq2n}
      $TRIMGALORE_HOME/trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/PE/ --paired ./${Fq1n} ./${Fq2n} &>> trim.log
      echo "${line} finished trimming"
  done < "$FILE"
  echo " trim_galore ended as $(date)"
  echo ==== End Trimming ====
  echo
}
#=======================================================================================================
##########
# rename #
##########

function rename_sr () {
  echo ==== START RECONSTRUCT ====
  cd ./trim/SR/
  if [[ ! -f  "*_renamed.fq.gz" ]]; then
    for var in `find . -type f -name "*_trimmed.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm *_trimmed.fq.gz
  fi
  echo 'All single-end sequencing data have renamed!'
  cd ../..
  echo 'Back to main directory'
  echo ==== END RECONSTRUCT ====
  echo
}

function rename_pe () {
  echo ==== START RECONSTRUCT ====
    for var in `find . -type f -name "*_1.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    for var in `find . -type f -name "*_2.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm *_val_?.fq.gz
  echo 'All pair-end sequencing data have renamed!'
  echo ==== END RECONSTRUCT ====
  echo
}
#=======================================================================================================
################################
# both- prepare before trinity #
################################

function pretrinity_both () {
  cd ${PRJNA_PATH}
  echo ====Start Normalization Respectively====
  mkdir ./normalize
  OUTDIR="../../normalize"

  cd ./trim/SR/
  echo "start to normalize single-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                                 --single `ls -m *.fq* |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 12 --output ${OUTDIR}
  echo "end of single-end sequencing data normalization"
  echo

  cd ../PE/
  echo "start to normalize paired-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                 --left `ls -m *1_val_1_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                 --right `ls -m *2_val_2_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                 --CPU 10 --output ${OUTDIR} --pairs_together --PARALLEL_STATS
  echo "end of single-end sequencing data normalization"

  cd ${PRJNA_PATH}
  echo
  echo ====End Normalization Respectively====
  echo

  echo "Start combining"
  cd ./normalize
  cat SRR*_ext_all_reads.normalized_K25_maxC50_minC0_pctSD10000.fq > bigfile.fq
  cd ${PRJNA_PATH}
  echo "Combined fastq file created!"
  echo
}

#=======================================================================================================
##################
# normalization #
##################


function normalize_sr () {
  OUTDIR="./normalize"
  echo ==== Trinity In silico Read Normalization START ====
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                                 --single `ls -m *.fq* |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 12 --output ${OUTDIR}
  echo ==== Trinity In silico Read Normalization END ====
}

function normalize_pe () {
  OUTDIR="./normalize"
  echo ==== Trinity In silico Read Normalization START ====
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                                 --left `ls -m *1_val_1_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --right `ls -m *2_val_2_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 10 --output ${OUTDIR} --pairs_together --PARALLEL_STATS
  echo ==== Trinity In silico Read Normalization END ====
}


#=======================================================================================================
####################
# de novo assembly #
####################


function assembly () {
  cd ${PRJNA_PATH}
  if [[ ! -e "runs_info.txt" ]]; then
    echo "There is no runs information, named 'runs_info.txt' for assembly under this directory. "
    return 1
  fi

  echo ==== De Novo Assembly START ====
  cd ./processed_data/
  mv -v ../runs_info.txt .
  $TRINITY_HOME/Trinity --seqType fq --max_memory 60G --CPU 10 \
    --samples_file runs_info.txt \
    --full_cleanup

  if [[ -e './trinity_out_dir/*Trinity.fasta' ]]; then
    mv -v ./trinity_out_dir/*Trinity.fasta ./trinity_out_dir/Trinity.fasta
    echo
    echo ==== De Novo Assembly END ====
    echo

    cd  ${PRJNA_PATH}
    cp -r -v ./processed_data/trinity_out_dir .
    rm -r -v ./processed_data/trinity_out_dir

  else
    echo === Assembly has some error! ===
    return 1
  fi
}


#=======================================================================================================
##################
# quantification #
##################


function count () {
  cd  ${PRJNA_PATH}

  cd processed_data

  # # create a transcripts quantification folder to hold the results
  # mkdir ./transcripts_count
  #
  # cd ${PRJNA_PATH}/transcripts_count

  $TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq \
      --transcripts ../trinity_out_dir/Trinity.fasta \
      --est_method RSEM --thread_count 16 \
      --aln_method bowtie2 --trinity_mode --prep_reference \
      --samples_file runs_info.txt

  cd  ${PRJNA_PATH}
  cp -r -v ./processed_data/trinity_out_dir .
  rm -r -v ./processed_data/trinity_out_dir


  # move sample_file to the ${PRJNA_PATH}
  mv sample_file_*.txt ${PRJNA_PATH}

  cd ${PRJNA_PATH}
}


function count_sr () {
  #cd ./trim/SR/
  cd ./trinity_out_dir

  # count transcripts of each sample sepreately
  #for TSRFILES in `find . -type f -name "SRR*_trimmed_renamed.fq.gz"`; do
  #  echo '====sample ${TSRFILES%_trimmed_renamed.fq.gz} starts to align===='
    #mkdir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}

  align_and_estimate_abundance.pl --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta --seqType fq \
      --est_method RSEM \
      --aln_method bowtie --trinity_mode --prep_reference \
      --samples_file sample_file_*.txt

      # --single ${TSRFILES}
      # --gene_trans_map ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta.gene_trans_map \
      # --output_dir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz} \


#    echo = preview result table =
#    head -5 ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}/RSEM.isoforms.results | column -t
#    echo
    cd ${PRJNA_PATH}
#  done
}

function rerun_abundance_estimation () {
  cd ./trim/SR/
  for TSRFILES in `find . -type f -name "SRR*_trimmed_renamed.fq.gz"`; do
    echo '====sample ${TSRFILES%_trimmed_renamed.fq.gz} starts to align===='
    mkdir -p ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}
    align_and_estimate_abundance.pl --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta --seqType fq \
      --est_method salmon --output_dir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz} \
      #--gene_trans_map ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta.gene_trans_map \
      --aln_method bowtie --trinity_mode --prep_reference \
      --single ${TSRFILES}
    echo = preview result table =
    head -5 ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}/RSEM.isoforms.results | column -t
    echo
    cd ${PRJNA_PATH}/trim/SR/
  done
}

function count_pe () {
  cd ./trim/PE/
  # count transcripts of each sample sepreately
  for TPEFILES in `find . -type f -name "SRR*_1_val_1_renamed.fq.gz"`; do
    echo '====sample ${TPEFILES%_1_val_1_renamed.fq.gz} starts to align===='
    mkdir ${PRJNA_PATH}/transcripts_count/${TPEFILES%_1_val_1_renamed.fq.gz}
    TLEFTFILES=TPEFILES
    TRIGHTFILES="${TPEFILES%_1_val_1_renamed.fq.gz}'_2_val_2_renamed.fq.gz'"
    align_and_estimate_abundance.pl --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta --seqType fq \
      --est_method RSEM --output_dir ${PRJNA_PATH}/transcripts_count/${TPEFILES%_1_val_1_renamed.fq.gz} \
      --gene_trans_map ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta.gene_trans_map \
      --aln_method bowtie --trinity_mode --prep_reference \
      --samples_file sample_file_?.txt
      #--left ${TLEFTFILES} --right ${TRIGHTFILES}
    echo = preview result table =
    head -5 ${PRJNA_PATH}/transcripts_count/${TPEFILES%_1_val_1_renamed.fq.gz}/RSEM.isoforms.results | column -t
    echo
    cd ${PRJNA_PATH}/trim/PE/
  done
}

function count_bo () {
  count_sr
  count_pe
}

#=======================================================================================================
#####################
# expression matrix #
#####################

function rerun_exp () {
  # create the expression matrix
  cd ${PRJNA_PATH}/transcripts_count/
  ls */quant.sf > genes.quant_files.txt
  abundance_estimates_to_matrix.pl --est_method salmon \
    --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \
    --cross_sample_norm TMM \
    --out_prefix salmon-gene \
    --name_sample_by_basedir \
    --quant_files genes.quant_files.txt
}

function rerun_plotcount () {
  cd ${PRJNA_PATH}/transcripts_count/
  # count Numbers of Expressed Transcripts
  count_matrix_features_given_MIN_TPM_threshold.pl \
          salmon-gene.isoform.TPM.not_cross_norm | tee salmon-gene.isoform.TPM.not_cross_norm.counts_by_min_TPM
  # source a R script for visualization of transctipts counts
  Rscript ../rerun_quantify.R
  cd ${PRJNA_PATH}
}

function expressionmx () {
  find . -name Trinity.fasta.gene_trans_map -exec mv {} ${PRJNA_PATH}/transcripts_count/ \;
  # create the expression matrix
  cd ${PRJNA_PATH}/transcripts_count/
  find . -name RSEM.isoforms.results > isoforms.quant_files.txt
  #find . -name RSEM.genes.results > genes.quant_files.txt
  $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method RSEM \
    --gene_trans_map Trinity.fasta.gene_trans_map \
    --cross_sample_norm TMM \
    --out_prefix rsem-gene \
    --name_sample_by_basedir \
    --quant_files isoforms.quant_files.txt
#    --quant_files genes.quant_files.txt
#    --gene_trans_map Trinity.fasta.gene_trans_map \
  ## output files:
  #    rsem-gene.isoform.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
  #    rsem-gene.isoform.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
  #    rsem-gene.isoform.TMM.EXPR.matrix : a matrix of TMM-normalized expression values

  # count Numbers of Expressed Transcripts
  $TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
          rsem-gene.gene.TPM.not_cross_norm | tee rsem-gene.gene.TPM.not_cross_norm.counts_by_min_TPM
#  $TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
#          rsem-gene.isoform.TPM.not_cross_norm | tee rsem-gene.isoform.TPM.not_cross_norm.counts_by_min_TPM

  # source a R script for visualization of transctipts counts
  #Rscript ../quantify.R
  cd ${PRJNA_PATH}
}


#=======================================================================================================
################
# DEG analysis #
################

function rerun_dge () {
  cd ${PRJNA_PATH}/transcripts_count/
  # differential expression analysis
  run_DE_analysis.pl \
    --matrix salmon-gene.isoform.counts.matrix \
    --method edgeR \
    --dispersion 0.1 \
    --output edgeR
  #>>
  echo `ls -l`
  echo '== go to edgeR =='
  echo `ls -l ./edgeR/`
  ls ./edgeR/ > edgeR.log
  head edgeR.log | column -t
  # select the transcripts the FDR<=0.05
  #for pair in `find ./edgeR/ -type f -name 'salmon-gene.isoform.counts.matrix.*.edgeR.DE_results'`; do
   # SIGRESULT="${pair}_FDR005"
   # touch ${SIGRESULT}
   # sed '1,1d' ${pair} | awk '{ if ($5 <= 0.05) print;}' > ${SIGRESULT}
   # echo `wc -l ${SIGRESULT}`
  #done
  # heatmap plot
  cd ${PRJNA_PATH}/transcripts_count/edgeR
  analyze_diff_expr.pl --matrix ../*.isoform.TMM.EXPR.matrix -P 1e-3 -C 2
  ##     -P: p-value cutoff; -C: log 2 fold change cutoff (-C 2 means 4 times)
  cd ${PRJNA_PATH}
}

function difexpre () {
  cd ${PRJNA_PATH}/transcripts_count/
  # differential expression analysis
  run_DE_analysis.pl \
    --matrix rsem-gene.isoform.counts.matrix \
    --method edgeR \
    --dispersion 0.1 \
    --output edgeR
  # generate interactive volcano and MA plots using Glimma
  #?????
  # heatmap plot
  cd ${PRJNA_PATH}/transcripts_count/edgeR
  analyze_diff_expr.pl --matrix ../*.isoform.TMM.EXPR.matrix -P 1e-3 -C 2
  ##     (can run GO enrichment analysis with GO annotations file)
  ##     -P: p-value cutoff; -C: log 2 fold change cutoff (-C 2 means 4 times)
  cd ${PRJNA_PATH}
}

#=======================================================================================================
###############################
# ortholog searching pipeline #
###############################


function run_ortho_pip() {
  source ${RUBY_SCRIPTS}/run_orthomcl.sh -l *.lst -o output
}
