#!/bin/bash

#=======================================================================================================
############################
# raw data quality control #
############################
function raw_qc () {
  mkdir ./raw_fastqc/
  # generate quality control with each untrimmed runs
  fastqc -o ./raw_fastqc/ -f fastq -t 5 --extract SRR*.fastq.gz
  # combine reports of all runs in a study to one quality control report
  multiqc ./raw_fastqc/ -o ./raw_fastqc/multiqc_output/

  touch ./raw_fastqc/multiqc_output/log_list.txt
  ls -l ./raw_fastqc/multiqc_output/multiqc_data/ > ./raw_fastqc/multiqc_output/log_list.txt

  # remove the separate fastqc files just keep the multiqc report
  find ./raw_fastqc -type f -name "SRR*" -exec rm {} \;
  find ./raw_fastqc -type d -name "SRR*" -exec rm -r {} \;
}
#=======================================================================================================
########
# trim #
########
function trim_sr () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_fastqc/
  mkdir ./trim/SR/
  echo ===================================
  echo "there are ${SE_N} single-end reads"
  echo ===================================
  # trim
  trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/SR/ `ls *.fastq.gz | grep -v "_"`
  echo " trim_galore ended as $(date)"
  echo ==== End Trimming ====
  echo
}


function trim_pe () {
  echo ==== Start Trimming ====
  echo "trim_galore cut adapters started at $(date)"
  mkdir ./trim
  mkdir ./trim/trimmed_fastqc/
  mkdir ./trim/PE/
  echo =====================================
  echo "there are ${PE_L_N} paired-end reads"
  echo =====================================
  # trim
  ls *_1.fastq.gz | cat | sed 's/_1.fastq.gz//g' | cat > pairname.txt
  # trim each pair
  FILE=pairname.txt
  while IFS= read -r line; do
      Fq1n="${line}'_1.fastq.gz'"
      Fq2n="${line}'_2.fastq.gz'"
      echo ${Fq1n}
      echo ${Fq2n}
      trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/PE/ --paired ./${Fq1n} ./${Fq2n}
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
    for var in `find . -type f -name "SRR*_trimmed.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm SRR*_trimmed.fq.gz
  fi
  echo 'All single-end sequencing data have renamed!'
  cd ../..
  echo 'Back to main directory'
  echo ==== END RECONSTRUCT ====
  echo
}

function rename_pe () {
  echo ==== START RECONSTRUCT ====
  cd ./trim/PE/
  if [[ ! -f  "*_renamed.fq.gz" ]]; then
    for var in `find . -type f -name "SRR*_val_[12].fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm SRR*_[12]_val_[12].fq.gz
  fi
  echo 'All pair-end sequencing data have renamed!'
  cd ../..
  echo 'Back to main directory'
  echo ==== END RECONSTRUCT ====
  echo
}
#=======================================================================================================
################################
# both- prepare before trinity #
################################

function pretrinity_both () {
  echo ====Start Normalization Respectively====
  cd ./trim/SR/
  SRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize single-end sequencing data"
  insilico_read_normalization.pl --seqType fq --JM 10G --max_cov 50 \
    --single ${SRFILES} --CPU 4
  echo "end of single-end sequencing data normalization"
  echo
  cd ../PE/
  LEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  RIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize paired-end sequencing data"
  insilico_read_normalization.pl --seqType fq --JM 10G --max_cov 50 \
    --left ${LEFTFILES} --right ${RIGHTFILES} --CPU 4
  echo "end of single-end sequencing data normalization"
  cd ${PRJNA_PATH}
  echo
  echo ====End Normalization Respectively====
  echo
  ## the output files include: ()
  #    ??????
  echo "Start combining"
  zcat ${PRJNA_PATH}/trim/SR/SRR*.fq.gz ${PRJNA_PATH}/trim/PE/SRR*.fq.gz | gzip -c > ${PRJNA_PATH}/bigfile.fastq.gz
  echo "Combined fastq file created!"
  echo
}

#=======================================================================================================
##################
# normalization #
##################
function normal_sr () {
  echo ==== Trinity In silico Read Normalization START ====
  insilico_read_normalization.pl --seqType fq --JM 50G --max_cov 50 \
                                 --single *.fq* --CPU 16
 #  --seqType <string>      :type of reads: ( 'fq' or 'fa')
 #  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for
 #                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
 #  --max_cov <int>         :targeted maximum coverage for reads.
 #
 #  If paired reads:
 #      --left  <string>    :left reads
 #      --right <string>    :right reads
 #
 #  Or, if unpaired reads:
 #      --single <string>   :single reads
 #
 #  Or, if you have read collections in different files you can use 'list' files, where each line in a list
 #  file is the full path to an input file.  This saves you the time of combining them just so you can pass
 #  a single file for each direction.
 #      --left_list  <string> :left reads, one file path per line
 #      --right_list <string> :right reads, one file path per line
  echo ==== Trinity In silico Read Normalization END ====
}

function normal_pe () {
  echo ==== Trinity In silico Read Normalization START ====

  insilico_read_normalization.pl --seqType fq --JM 50G --max_cov 50 \
                                 --left `ls -m *` --CPU 16
 #  --seqType <string>      :type of reads: ( 'fq' or 'fa')
 #  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for
 #                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
 #  --max_cov <int>         :targeted maximum coverage for reads.
 #
 #  If paired reads:
 #      --left  <string>    :left reads
 #      --right <string>    :right reads
 #
 #  Or, if unpaired reads:
 #      --single <string>   :single reads
 #
 #  Or, if you have read collections in different files you can use 'list' files, where each line in a list
 #  file is the full path to an input file.  This saves you the time of combining them just so you can pass
 #  a single file for each direction.
 #      --left_list  <string> :left reads, one file path per line
 #      --right_list <string> :right reads, one file path per line
  echo ==== Trinity In silico Read Normalization END ====
}

#=======================================================================================================
####################
# de novo assembly #
####################
#     ______  ____   ____  ____   ____  ______  __ __
#    |      ||    \ |    ||    \ |    ||      ||  |  |
#    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
#    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
#      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
#      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
#      |__|  |__|\_||____||__|__||____|  |__|  |____/
#
# Required:
#
#  --seqType <string>      :type of reads: ('fa' or 'fq')
#
#  --max_memory <string>   :suggested max memory to use by Trinity where limiting can be enabled. (jellyfish, sorting
#                               , etc)
#                          provided in Gb of RAM, ie.  '--max_memory 10G'
#
#  If paired reads:
#      --left  <string>    :left reads, one or more file names (separated by commas, no spaces)
#      --right <string>    :right reads, one or more file names (separated by commas, no spaces)
#
#  Or, if unpaired reads:
#      --single <string>   :single reads, one or more file names, comma-delimited (note, if single file contains pairs,
#      can use flag: --run_as_paired )
#
#  Or, --samples_file <string>      tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#
#                                   if single-end instead of paired-end, then leave the 4th column above empty.
#
##  Misc:
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#
#  --CPU <int>                     :number of CPUs to use, default: 2
#  --min_contig_length <int>       :minimum assembled contig length to report
#                                   (def=200)
#
#  --long_reads <string>           :fasta file containing error-corrected or circular consensus (CCS) pac bio reads
#                                   (** note: experimental parameter **, this functionality continues to be under
#                                   development)
#
#  --genome_guided_bam <string>    :genome guided mode, provide path to coordinate-sorted bam file.
#                                   (see genome-guided param section under --show_full_usage_info)
#
#  --jaccard_clip                  :option, set if you have paired reads and
#                                   you expect high gene density with UTR
#                                   overlap (use FASTQ input file format for reads).
#                                   (note: jaccard_clip is an expensive operation, so avoid using it unless
#                                   necessary due to finding excessive fusion transcripts w/o it.)
#
#  --trimmomatic                   :run Trimmomatic to quality trim reads
#                                        see '--quality_trimming_params' under full usage info for tailored settings.
#
#
#  --no_normalize_reads            :Do *not* run in silico normalization of reads. Defaults to max. read coverage of 50.
#                                       see '--normalize_max_read_cov' under full usage info for tailored settings.
#                                       (note, as of Sept 21, 2016, normalization is on by default)
#
#  --no_distributed_trinity_exec   :do not run Trinity phase 2 (assembly of partitioned reads), and stop after
#  generating command list.
#
#
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( your current working directory: "/media/heyland-lab/184EE9BC1B1C919A/
#                                   Ruby_extra/RawData/PRJNA475804/trim/SR/trinity_out_dir"
#                                    note: must include 'trinity' in the name as a safety precaution! )
#
#  --workdir <string>              :where Trinity phase-2 assembly computation takes place (defaults to
#                                   --output setting).
#                                  (can set this to a node-local drive or RAM disk)
#
#  --full_cleanup                  :only retain the Trinity fasta file, rename as ${output_dir}.Trinity.fasta
#
#  --cite                          :show the Trinity literature citation
#
#  --verbose                       :provide additional job status info during the run.
#
#  --version                       :reports Trinity version (Trinity-v2.6.6) and exits.
#
#  --show_full_usage_info          :show the many many more options available for running Trinity (expert usage).
#



function assembly_sr () {
  echo ==== De Novo Assembly START ====
  mv ./sample_file_?.txt ./trim/SR/
  cd ./trim/SR/
  #???????
  #TSRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  mkdir ./trinity_out_dir/
  Trinity --seqType fq --max_memory 100G --CPU 16 \
    --single single.norm.fq #--no_normalize_reads
    #--samples_file sample_file_?.txt   #--single ${TSRFILES}
  cd ../../
  mv ./trim/SR/trinity_out_dir ./
  #ls -l
  echo
  echo ==== De Novo Assembly END ====
  echo
}

function assembly_sr_nochwo () {
  echo ==== De Novo Assembly START ====
  mv ./sample_file_?.txt ./trim/SR/
  cd ./trim/SR/
  #???????
  #TSRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  mkdir ./trinity_out_dir/
  Trinity --seqType fq --max_memory 50G --CPU 16 \
    --samples_file sample_file_?.txt   #--single ${TSRFILES}
  cd ../../
  mv ./trim/SR/trinity_out_dir ./
  #ls -l
  echo
  echo ==== De Novo Assembly END ====
  echo

}


function assembly_pe () {
  echo ==== De Novo Assembly START ====
  mv ./sample_file_?.txt ./trim/PE/
  cd ./trim/PE/
  #????????
  #TLEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  #TRIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  mkdir ./trinity_out_dir/
  Trinity --seqType fq --max_memory 100G --CPU 16 \
     --samples_file sample_file_?.txt
    #--output "${PRJNA_PATH}/trinity_out_dir/" \
  #Trinity --seqType fq --max_memory 50G --CPU 6 \
  #  --output ${PRJNA_PATH}/trinity_out_dir/ \
  #  --left ${TLEFTFILES} --right ${TRIGHTFILES}
  cd ../../
  mv ./trim/PE/trinity_out_dir ./
  #ls -l
  echo
  echo ==== De Novo Assembly END ====
  echo
}

function assembly_bo () {
  echo ==== De Novo Assembly START ====
  Trinity --seqType fq --max_memory 100G --CPU 16 \
    --no_normalize_reads --run_as_paired \
    --single bigfile.fastq.gz
  echo
  echo ==== De Novo Assembly END ====
  echo
}
#=======================================================================================================
##################
# quantification #
##################

function count_sr () {
  cd ./trim/SR/
  # count transcripts of each sample sepreately
  for TSRFILES in `find . -type f -name "SRR*_trimmed_renamed.fq.gz"`; do
    echo '====sample ${TSRFILES%_trimmed_renamed.fq.gz} starts to align===='
    mkdir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}
    align_and_estimate_abundance.pl --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta --seqType fq \
      --est_method RSEM --output_dir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz} \
      --gene_trans_map ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta.gene_trans_map \
      --aln_method bowtie --trinity_mode --prep_reference \
      --samples_file sample_file_?.txt
      #--single ${TSRFILES}
    echo = preview result table =
    head -5 ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}/RSEM.isoforms.results | column -t
    echo
    cd ${PRJNA_PATH}/trim/SR/
  done
}

function rerun_abundance_estimation () {
  cd ./trim/SR/
  for TSRFILES in `find . -type f -name "SRR*_trimmed_renamed.fq.gz"`; do
    echo '====sample ${TSRFILES%_trimmed_renamed.fq.gz} starts to align===='
    mkdir -p ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz}
    align_and_estimate_abundance.pl --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta --seqType fq \
      --est_method salmon --output_dir ${PRJNA_PATH}/transcripts_count/${TSRFILES%_trimmed_renamed.fq.gz} \
      --gene_trans_map ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta.gene_trans_map \
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
##############
# annotation #
##############
  #Build_Trinotate_Boilerplate_SQLite_db.pl Trinotate
  #makeblastdb -in uniprot_sprot.pep -dbtype prot
  #gunzip Pfam-A.hmm.gz
  #hmmpress Pfam-A.hmm
  #/media/heyland-lab/"Seagate Backup Plus Drive"/ruby/annotation/

function annotation () {
  cd ${PRJNA_PATH}
  cp /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/annotation/* .
  TransDecoder.LongOrfs -t ./trinity_out_dir/Trinity.fasta
  #Search Trinity transcripts
  blastx -query ./trinity_out_dir/Trinity.fasta -db uniprot_sprot.pep -num_threads 8 \
    -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx.outfmt6
  #Search Transdecoder-predicted proteins
  blastp -query ./Trinity.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.pep \
    -num_threads 8 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp.outfmt6
  #Run hmmscan
  hmmscan --cpu 5 --domtblout TrinotatePFAM.out \
    Pfam-A.hmm ./Trinity.fasta.transdecoder_dir/longest_orfs.pep  > pfam.log

}
