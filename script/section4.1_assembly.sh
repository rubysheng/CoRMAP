#!/bin/bash
#title          :section4.1_assembly.sh
#description    :de novo assembly
#author         :Ruby(Yiru) Sheng
#date           :20191020
#version        :1.2
#usage          :./section4.1_assembly.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

# SOURCE THE "PRJNA_PATH" !!!!!!!




####################
# de novo assembly #
####################

#### usage guidence :   $TRINITY_HOME/Trinity ####
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


function assembly () {
  # now we are in ${PRJNA_PATH}/trim/PE( or SR)
  # go to this location in def_finddir_assembly!!!!

  echo ==== De Novo Assembly START ====

  # move renamed files and group inoformation to internal Drive
  cp *_renamed.fq.gz /home/lewis/Documents/DeNovoAssambly/
  cp ${PRJNA_PATH}/sample_file_?.txt /home/lewis/Documents/DeNovoAssambly/

  cd /home/lewis/Documents/DeNovoAssambly/
  # create the output directory to hold results
  mkdir ./trinity_out_dir/

  conda activate salmon
  $TRINITY_HOME/Trinity --seqType fq --max_memory 100G --CPU 16 \
    --samples_file sample_file_?.txt --no_normalize_reads \
    --full_cleanup
  conda deactivate

  # move output back to the dataset directory and free the space
  mv trinity_out_dir ${PRJNA_PATH}/
  rm *

  cd  ${PRJNA_PATH}
  echo
  echo ==== De Novo Assembly END ====
  echo
}
