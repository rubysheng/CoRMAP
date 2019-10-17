#!/bin/bash
#SBATCH --time=1-00:00:00                   # Time limit hrs:min:sec
#SBATCH --account=def-lukens
#SBATCH --job-name=bombina
#SBATCH --mem=70G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=bombina.normalization.log     # Standard output and error log
date;hostname;pwd

module load trinity
module load samtools
module load jellyfish
module load gcc/5.4.0
module load trinity
module load salmon/0.14.0
module load jellyfish/2.2.6


/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc5.4/trinity/2.8.5/trinityrnaseq-Trinity-v2.8
.5/util/insilico_read_normalization.pl \
  --seqType fq --JM 50G  --max_cov 50 --min_cov 1 --CPU 16 \
  --output /project/6001295/ruby/transcriptional_profiles/trinity_out_dir/insilico_read_normalization \
  --left /project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_1_T24R1_E4_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_2_T24R2_E5_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_3_T24R3_E6_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_1_T2R1_E7_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_2_T2R2_E9_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_3_T2R3_E10_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_1_T4R1_E1_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_2_T4R2_E2_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_3_T4R3_E3_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_1_C1_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_2_C2_R1_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_3_C3_R1_renamed.fq.gz \
  --right /project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_1_T24R1_E4_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_2_T24R2_E5_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/24_Hour_Replicate_3_T24R3_E6_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_1_T2R1_E7_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_2_T2R2_E9_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/2_Hour_Replicate_3_T2R3_E10_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_1_T4R1_E1_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_2_T4R2_E2_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/4_Hour_Replicate_3_T4R3_E3_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_1_C1_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_2_C2_R2_renamed.fq.gz,/project/6001295/ruby/transcriptional_profiles/Control_Replicate_3_C3_R2_renamed.fq.gz \
   --pairs_together --PARALLEL_STATS

date


Trinity --seqType fq --max_memory 60G --samples_file sample_file_bombina.txt --SS_lib_type FR --trimmomatic --CPU 14

Trinity --seqType fq --max_memory 60G --left left.norm.fq --right right.norm.fq --SS_lib_type FR --CPU 14\
  --no_normalize_reads