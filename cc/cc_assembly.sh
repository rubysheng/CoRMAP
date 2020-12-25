#!/bin/bash
#SBATCH --time=1-10:00:00                   # Time limit hrs:min:sec
#SBATCH --account=def-heylanda
#SBATCH --job-name=assembly_directly
#SBATCH --mem=150G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=assembly.log             # Standard output and error log
date;hostname;pwd

#module load trinity
module load samtools
module load jellyfish
module load gcc/5.4.0
#module load trinity
module load salmon/0.14.0
module load jellyfish/2.2.6
module load bowtie2



function assembly () {
  mkdir ./trinity_out_dir/
  echo ==== Trinity Assembly START ====
  OUTDIR="./trinity_out_dir/"

  /home/ruby/trinityrnaseq-2.8.5/Trinity --seqType fq --max_memory 120G --CPU 30 \
    --samples_file sample_file_*.txt  \
    --full_cleanup  --output ${OUTDIR} --bflyCalculateCPU
  echo ==== De Novo Assembly END ====
}

assembly


date
