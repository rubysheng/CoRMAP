#!/bin/bash

#=======================================================================================================
############################
# raw data quality control #
############################
function raw_qc () {
  mkdir ./raw_fastqc/
  # generate quality control with each untrimmed runs
  fastqc -o ./raw_fastqc/ -f fastq -t 5 --extract *.fastq.gz #SRR*.fastq.gz
  # combine reports of all runs in a study to one quality control report
  conda activate multiqc  # activates environment
  multiqc ./raw_fastqc/ -o ./raw_fastqc/multiqc_output/
  conda deativate         # deactivates

  touch ./raw_fastqc/multiqc_output/log_list.txt
  ls -l ./raw_fastqc/multiqc_output/multiqc_data/ > ./raw_fastqc/multiqc_output/log_list.txt

  # remove the separate fastqc files just keep the multiqc report
  #find ./raw_fastqc -type f -name "SRR*" -exec rm {} \;
  #find ./raw_fastqc -type d -name "SRR*" -exec rm -r {} \;
  find ./raw_fastqc -type f -name "*_Replicate_*" -exec rm {} \;
  find ./raw_fastqc -type d -name "*_Replicate_*" -exec rm -r {} \;
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
  $TRIMGALORE_HOME/trim_galore --phred33 --fastqc --gzip --trim-n --output_dir ./trim/SR/ `ls *.fastq.gz | grep -v "_"` > trim.log
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
  #ls *_R1.fastq.gz | cat | sed 's/_R1.fastq.gz//g' | cat > pairname.txt
  ls *_1.fastq.gz | cat | sed 's/_1.fastq.gz//g' | cat > pairname.txt
  # trim each pair
  FILE=pairname.txt
  while IFS= read -r line; do
      Fq1n="${line}_1.fastq.gz"
      Fq2n="${line}_2.fastq.gz"
      #Fq1n="${line}_R1.fastq.gz"
      #Fq2n="${line}_R2.fastq.gz"
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
      # zcat ${var} | awk '{ if (NR%4 == 1) { print $1"_"$2"/1" } else { print } }'  | gzip -c > ${NEWNAME}
    done
    for var in `find . -type f -name "*_2.fq.gz"`; do
      echo ${var}
      NEWNAME="${var%.fq.gz}_renamed.fq.gz"
      echo ${NEWNAME}
      zcat ${var} | awk '{{print (NR%4 == 1) ? substr($1,1,11) "_" ++i "/" substr($2,length($NF),1): $0}}'  | gzip -c > ${NEWNAME}
    done
    rm *_val_?.fq.gz
# mv ./*_renamed.fq.gz ./trim/PE/
#  for var in `find . -type f -name "2_Hour_Replicate_1_T2R1_E7_R1.fastq"`; do
#    echo ${var}
#    NEWNAME="${var%.fastq}_renamed.fq"
#    echo ${NEWNAME}
#    cat ${var} | awk '{ if (NR%4 == 1) { print $1"_"$2"/1" } else { print } }' > ${NEWNAME}
#  done
  echo 'All pair-end sequencing data have renamed!'
  echo ==== END RECONSTRUCT ====
  echo
}


#=======================================================================================================
################################
# both- prepare before trinity #
################################

function pretrinity_both () {
  echo ====Start Normalization Respectively====
  mkdir ./normalize_out_dir
  OUTDIR="../../normalize_out_dir"

  cd ./trim/SR/
  # SRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize single-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                                 --single `ls -m *.fq* |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 12 --output ${OUTDIR}
  # insilico_read_normalization.pl --seqType fq --JM 10G --max_cov 50 --single ${SRFILES} --CPU 4
  echo "end of single-end sequencing data normalization"
  echo

  cd ../PE/
  # LEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  # RIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  echo "start to normalize paired-end sequencing data"
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                 --left `ls -m *1_val_1_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                 --right `ls -m *2_val_2_renamed.fq.gz |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                 --CPU 10 --output ${OUTDIR} --pairs_together --PARALLEL_STATS
  # insilico_read_normalization.pl --seqType fq --JM 10G --max_cov 50 --left ${LEFTFILES} --right ${RIGHTFILES} --CPU 4
  echo "end of single-end sequencing data normalization"

  cd ${PRJNA_PATH}
  echo
  echo ====End Normalization Respectively====
  echo

  echo "Start combining"
  # zcat ${PRJNA_PATH}/trim/SR/SRR*.fq.gz ${PRJNA_PATH}/trim/PE/SRR*.fq.gz | gzip -c > ${PRJNA_PATH}/bigfile.fastq.gz
  cd ./normalize_out_dir
  cat SRR*_ext_all_reads.normalized_K25_maxC50_minC0_pctSD10000.fq > bigfile.fq
  echo "Combined fastq file created!"
  echo
}

#=======================================================================================================
##################
# normalization #
##################


function normalize_sr () {
  #mkdir ${PRJNA_PATH}/normalization/
  OUTDIR="../normalize_out_dir"
  #OUTDIR="../normalize_out_dir_1129"
  echo ==== Trinity In silico Read Normalization START ====
  $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq --JM 40G --max_cov 50 \
                                 --single `ls -m *.fq* |sed 's/ //g' | sed ":a;N;s/\n//g;ta"` \
                                 --CPU 12 --output ${OUTDIR}
  echo ==== Trinity In silico Read Normalization END ====
}

function normalize_pe () {
  #mkdir ${PRJNA_PATH}/normalization/
  OUTDIR="../normalize_out_dir"
  #OUTDIR="../normalize_out_dir_1129"
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


function assembly () {
  # now we are in ${PRJNA_PATH}/trim/PE( or SR)
  #PRJNA_PATH=$(pwd)/../..
  #echo  "PRJNA_PATH="${PRJNA_PATH}

  echo ==== De Novo Assembly START ====
  cp -v ${PRJNA_PATH}/sample_file_*.txt /home/lewis/Documents/DeNovoAssambly/

  cp -v *_renamed.fq.gz /home/lewis/Documents/DeNovoAssambly/


  #cd ./trim/SR/
  cd /home/lewis/Documents/DeNovoAssambly/

  #mkdir ./trinity_out_dir/

  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 60G --CPU 10 \
    --samples_file sample_file_*.txt \
    --full_cleanup
  conda deactivate
    #--single single.norm.fq --no_normalize_reads
    #--samples_file sample_file_?.txt   #--single ${TSRFILES}
    #--single *_ext_all_reads.normalized_K25_maxC50_minC0_pctSD10000.fq
    #--no_normalize_reads

  #ls -lh . # now in the directory : /home/lewis/Documents/DeNovoAssambly

  if [[ -e './*Trinity.fasta' ]]; then
    echo
    echo ==== De Novo Assembly END ====
    echo

    rm -v sample_file_*.txt
    rm -v SRR*

    cd  ${PRJNA_PATH}
    cd ../../
    mkdir trinity_out_dir
    mv -v ~/Documents/DeNovoAssambly/* ./trinity_out_dir/
    mv ./trinity_out_dir/trinity_out_dir.Trinity.fasta ./trinity_out_dir/Trinity.fasta
    mv ./trinity_out_dir/trinity_out_dir.Trinity.fasta.gene_trans_map ./trinity_out_dir/Trinity.fasta.gene_trans_map


  else
    echo === Assembly has some error! ===
  fi




}

function assembly_sr () {
  echo ==== De Novo Assembly START ====
  #mv ./sample_file_?.txt ./trim/SR/
  #cp ./trim/SR/*_ext_all_reads.normalized_K25_maxC50_minC0_pctSD10000.fq /home/lewis/Documents/DeNovoAssambly/
  #cp ./trim/SR/*
  #cp sample_file_?.txt /home/lewis/Documents/DeNovoAssambly/


  #cd ./trim/SR/
  cd /home/lewis/Documents/DeNovoAssambly/

  #TSRFILES=`ls -m *.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`

  #mkdir ./trinity_out_dir/

  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 100G --CPU 16 \
    --samples_file sample_file_?.txt --no_normalize_reads \
    --full_cleanup
  conda deactivate
    #--single single.norm.fq --no_normalize_reads
    #--samples_file sample_file_?.txt   #--single ${TSRFILES}
    #--single *_ext_all_reads.normalized_K25_maxC50_minC0_pctSD10000.fq

  cd  ${PRJNA_PATH}
  #cd ../../

  #mv ./trim/SR/trinity_out_dir ./

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
  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 50G --CPU 16 \
    --samples_file sample_file_?.txt \
    --full_cleanup
    #--single ${TSRFILES}
  conda deactivate
  cd ../../
  mv ./trim/SR/trinity_out_dir ./
  #ls -l
  echo
  echo ==== De Novo Assembly END ====
  echo

}


function assembly_pe () {
  echo ==== De Novo Assembly START ====
  mv ./sample_file_*.txt ./trim/PE/
  cd ./trim/PE/
  #mkdir ./trinity_out_dir/
  #TLEFTFILES=`ls -m SRR*_1_val_1_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  #TRIGHTFILES=`ls -m SRR*_2_val_2_renamed.fq.gz | sed 's/ //g' | sed ':t;N;s/\n//;b t'`
  #--output "${PRJNA_PATH}/trinity_out_dir/" \
  #Trinity --seqType fq --max_memory 50G --CPU 6 \
  #  --output ${PRJNA_PATH}/trinity_out_dir/ \
  #  --left ${TLEFTFILES} --right ${TRIGHTFILES}
  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 30G --CPU 16 \
    --samples_file sample_file_*.txt \
    --full_cleanup
  conda deactivate
  #cd ../../
  #mv ./trim/PE/trinity_out_dir ./
  if [[ -e './trinity_out_dir/Trinity.fasta' ]]; then
    echo
    echo ==== De Novo Assembly END ====
    echo
  else
    echo === Assembly has some error! ===
    #exit 1
  fi
}

function assembly_bo () {
  echo ==== De Novo Assembly START ====
  #mkdir ./trinity_out_dir/
  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 100G --CPU 16 \
    --no_normalize_reads --run_as_paired \
    --single bigfile.fastq.gz
    --full_cleanup
  conda deactivate
  echo
  echo ==== De Novo Assembly END ====
  echo
}

#=======================================================================================================
##################
# quantification #
##################


function count () {
  # create a transcripts quantification folder to hold the results
  mkdir ${PRJNA_PATH}/transcripts_count
  #cd ${PRJNA_PATH}/transcripts_count

  # move input files to this folder
  find . -type f -name *renamed.fq.gz -exec mv {} ${PRJNA_PATH}/transcripts_count \;
  find . -name sample_file_*.txt -exec mv {} ${PRJNA_PATH}/transcripts_count \;

  # go to ${PRJNA_PATH}/transcripts_count
  cd ${PRJNA_PATH}/transcripts_count

  $TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq \
      --transcripts ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta \
      --est_method RSEM --thread_count 16 \
      --aln_method bowtie2 --trinity_mode --prep_reference \
      --samples_file sample_file_*.txt

  # create a folder to hold processed(trim and rename) data
  mkdir ${PRJNA_PATH}/processed_data
  # move renamed files to the folder
  mv *renamed.fq.gz ${PRJNA_PATH}/processed_data
  #rm -r ${PRJNA_PATH}/trim/

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


#=======================================================================================================
############################
# orthomcl results summary #
############################

function ortho_pepsum() {
  #in the *_lst/ directory
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%.fasta}
    bioawk -c fastx '{ print "HEADER|"$name }' ./input/$sp_name.fasta > ./analyze/orgin_${sp_name%_pep}.lst
    sed -i "s/HEADER/${sp_name}/g" ./analyze/orgin_${sp_name%_pep}.lst
    cat ./output/groups/groups.txt | awk -v species="$sp_name" '{for(i=1; i<=NF; i++) if ($i ~ species) print $i}' > ./analyze/clustered_${sp_name%_pep}.lst
    cat ./analyze/orgin_${sp_name%_pep}.lst ./analyze/clustered_${sp_name%_pep}.lst | sort | uniq -u > ./analyze/unclustered_${sp_name%_pep}.lst
    #cat ./output/groups/unclustered_${sp_name%_pep}.lst
  done
  # | paste -s -d '\n' > ./analyze/unclustered.lst
}

function ortho_dnasum() {
  #in the *_lst/ directory
  cd ./analyze/
  for file in `ls -1 *.lst`; do
    echo $file
    dna_lst_name=${file%.lst}"_dnalst"
    echo $dna_lst_name
    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      echo ${trimmed%.p*}
      isoform=${trimmed%.p*}
      echo "$isoform" >> tmp_dnalst
    done < "$file"
    cat tmp_dnalst | sort | uniq > "$dna_lst_name"
    rm tmp_dnalst
  done
  cd ..
}

function ortho_dnastat() {
  #in the *_lst/ directory
  # generate a output file to hold the starts
  echo  "Species_Acc" "Count_Totaldna" "Count_Clustered" "Count_Unclustered" > ./analyze/cluster.dnastat  # HEADER
  i=1
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%_pep.fasta}
    echo "SPECIES "$i " Name : "${sp_name}
    echo
    # amount of all pep sequences (per species)
    SP_OR_CNT=`wc -l < ./analyze/orgin_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " total dna sequences : "${SP_OR_CNT}
    echo
    # amount of clustered pep sequences (per species)
    SP_CL_CNT=`wc -l < ./analyze/clustered_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " clustered dna sequences : "${SP_CL_CNT}
    echo
    # amount of unclustered pep sequences (per species)
    SP_UN_CNT=`wc -l < ./analyze/unclustered_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " unclustered dna sequences : "${SP_UN_CNT}
    echo
    # hold the stats to output file
    # awk -v name="${sp_name}" -v total="${SP_OR_CNT}" -v clustered="${SP_CL_CNT}" -v unclustered="${SP_UN_CNT}" '{ print name,total, clustered, unclustered }' >> cluster.stat
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >>./analyze/cluster.dnastat
    let i=i+1
  done
}

function ortho_pepstat() {
  # in the *_lst/ directory
  echo "Path : "$(pwd)
  file=`basename $(pwd)`
  echo "Taxonomy Group Name : "${file%_lst}
  echo

  ### step1 : count the total number of clustering groups in a Taxonomy group
  TOTAL_GP_COUNT=`grep -c "group" output/groups/groups.txt`
  echo "Total clustering groups : "${TOTAL_GP_COUNT}
  echo

  ### step2 : counts of proteins for each species in each cluster group
  ###########
  # cat output/groups/groups.txt | \
  # awk 'BEGIN {print "GROUP_NUM", "PRJNA?????1", "PRJNA?????2", "PRJNA?????3", "PRJNA?????4", "PRJNA?????5", "PRJNA?????6", "PRJNA?????7", "PRJNA?????8", "PRJNA?????9", "PRJNA????10"}
  # { sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
  # if ($i ~ "PRJNA?????1") sp1 += 1; \
  # else if ($i ~ "PRJNA?????2") sp2 += 1; \
  # else if ($i ~ "PRJNA?????3") sp3 += 1; \
  # else if ($i ~ "PRJNA?????4") sp4 += 1; \
  # else if ($i ~ "PRJNA?????5") sp5 += 1; \
  # else if ($i ~ "PRJNA?????6") sp6 += 1; \
  # else if ($i ~ "PRJNA?????7") sp7 += 1; \
  # else if ($i ~ "PRJNA?????8") sp8 += 1; \
  # else if ($i ~ "PRJNA?????9") sp9 += 1; \
  # else if ($i ~ "PRJNA????10") sp10 += 1; \
  # print $1, sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10 }' > ./analyze/groups.pepcounts.txt
  ###########

  ### step3 : count unclustered/clustered sequences per species
  study_all=`ls -1 input/*_pep.fasta | wc -l`
  echo "Amount Of Species : "${study_all}
  echo
  # generate a output file to hold the starts
  echo  "Species_Acc" "Count_Totalpep" "Count_Clustered" "Count_Unclustered" > ./analyze/cluster.pepstat  # HEADER
  i=1
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%_pep.fasta}
    echo "SPECIES "$i " Name : "${sp_name}
    echo
    # amount of all pep sequences (per species)
    SP_OR_CNT=`wc -l < ./analyze/orgin_${sp_name%_pep}.lst`
    echo "SPECIES "$i " total pep sequences : "${SP_OR_CNT}
    echo
    # amount of clustered pep sequences (per species)
    SP_CL_CNT=`wc -l < ./analyze/clustered_${sp_name%_pep}.lst`
    echo "SPECIES "$i " clustered pep sequences : "${SP_CL_CNT}
    echo
    # amount of unclustered pep sequences (per species)
    SP_UN_CNT=`wc -l < ./analyze/unclustered_${sp_name%_pep}.lst`
    echo "SPECIES "$i " unclustered pep sequences : "${SP_UN_CNT}
    echo
    # hold the stats to output file
    # awk -v name="${sp_name}" -v total="${SP_OR_CNT}" -v clustered="${SP_CL_CNT}" -v unclustered="${SP_UN_CNT}" '{ print name,total, clustered, unclustered }' >> cluster.stat
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >>./analyze/cluster.pepstat
    let i=i+1
  done
}

function ortho_allsp_g() {
  #in the *_lst/ directory
  cd analyze/stat/
  echo `head -1 groups.counts.stat` SUM > allsp_g.counts.stat

  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print $0,sum}' >> allsp_g.counts.stat

  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{ print $1 }' > tmp.allsp_g.lst

  while IFS= read -r line; do
    awk -v pat="$line" '{ if ($0 ~ pat) {print} }' ../../output/groups/groups.txt >> ../allsp_groups.txt
  done < tmp.allsp_g.lst
  rm tmp.allsp_g.lst
  cd ../../
}


function ortho_perg_seq() {
  #in the *_lst/ directory
  # check if the fasta files prepared
  if [ ! -d fasta_dir/ ]; then
    mkdir fasta_dir
    for file in `ls -1 ./input/*_pep.fasta`; do
      file=`basename $file`
      file=${file%_pep.fasta}"_RSEM.fasta"
      echo ${file}
      find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -name ${file} -exec cp -v {} fasta_dir/ \;
    done
  fi

  # randomly select 10 (test with 2) groups to test the similarity of sequences in groups
  cd analyze/
  mkdir per_group
  cd per_group
  shuf -n 2 ../../output/groups/groups.txt > 10g_groups.txt  # test all groups
  while IFS= read -r line; do
    echo ${line} > 10g_pergroup.txt
    group_num=`awk '{print $1}' 10g_pergroup.txt`
    group_num=${group_num%:}
    echo ${group_num}
    sed "s/ P/\nP/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname

    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      echo ${trimmed%.p*}
      isoform=${trimmed%.p*}
      awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 10g_pergroup.danseq.fasta
    done < 10g_pergroup.pepname

    output=$group_num".danseq.fasta"
    mv 10g_pergroup.danseq.fasta ${output}
    rm 10g_pergroup*
  done < 10g_groups.txt


  shuf -n 2 ../allsp_groups.txt > allsp_10g_groups.txt  # only test those all-species-included groups
  while IFS= read -r line; do
    echo ${line} > 10g_pergroup.txt
    group_num=`awk '{print $1}' 10g_pergroup.txt`
    group_num=${group_num%:}
    echo ${group_num}
    sed "s/ P/\nP/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname

    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      echo ${trimmed%.p*}
      isoform=${trimmed%.p*}
      awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 10g_pergroup.danseq.fasta
    done < 10g_pergroup.pepname

    output="allsp_"$group_num".danseq.fasta"
    mv 10g_pergroup.danseq.fasta ${output}
    rm 10g_pergroup*
  done < allsp_10g_groups.txt
  cd ../../

}

#=======================================================================================================
##################################
# orthomcl results check - blast #
##################################
function ortho_chk_blast() {
  #statements
  # in the "*_lst/analyze/per_group" directory
  # build a local blast database per group
  for dnafasta in `ls -1 *.danseq.fasta`; do
    makeblastdb -in ${dnafasta} -dbtype nucl -hash_index -out dna_pergroup
    output_name=${dnafasta%danseq.fasta}"blast.fm6"
    blastn -query ${dnafasta} -db dna_pergroup -evalue 1e-5 -out  ${output_name} -outfmt "6" -num_alignments 10 -num_threads 10
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
}
