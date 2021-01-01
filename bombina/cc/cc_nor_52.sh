#!/bin/bash
#SBATCH --time=2-20:00:00                   # Time limit hrs:min:sec
#SBATCH --account=def-heylanda
#SBATCH --job-name=normalize
#SBATCH --mem=150G
#SBATCH --cpus-per-task=18
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=normalize.log             # Standard output and error log
date;hostname;pwd

module load trinity
module load samtools
module load jellyfish
module load gcc/5.4.0
module load trinity
module load salmon/0.14.0
module load jellyfish/2.2.6


# mkdir normalization_dir
OUTDIR="./normalization_dir/"
echo ==== Trinity In silico Read Normalization START ====
$EBROOTTRINITY/trinityrnaseq-Trinity-v2.8.5/util/insilico_read_normalization.pl --seqType fq --JM 150G --max_cov 50 \
                               --single SRR8809606_trimmed_renamed.fq.gz,SRR8809607_trimmed_renamed.fq.gz,SRR8809608_trimmed_renamed.fq.gz,SRR8809609_trimmed_renamed.fq.gz,SRR8809610_trimmed_renamed.fq.gz,SRR8809611_trimmed_renamed.fq.gz,SRR8809612_trimmed_renamed.fq.gz,SRR8809613_trimmed_renamed.fq.gz,SRR8809614_trimmed_renamed.fq.gz,SRR8809615_trimmed_renamed.fq.gz,SRR8809616_trimmed_renamed.fq.gz,SRR8809617_trimmed_renamed.fq.gz,SRR8809618_trimmed_renamed.fq.gz,SRR8809619_trimmed_renamed.fq.gz,SRR8809620_trimmed_renamed.fq.gz,SRR8809621_trimmed_renamed.fq.gz,SRR8809622_trimmed_renamed.fq.gz,SRR8809623_trimmed_renamed.fq.gz,SRR8809624_trimmed_renamed.fq.gz,SRR8809625_trimmed_renamed.fq.gz,SRR8809626_trimmed_renamed.fq.gz,SRR8809627_trimmed_renamed.fq.gz,SRR8809628_trimmed_renamed.fq.gz,SRR8809629_trimmed_renamed.fq.gz,SRR8809630_trimmed_renamed.fq.gz,SRR8809631_trimmed_renamed.fq.gz,SRR8809632_trimmed_renamed.fq.gz,SRR8809633_trimmed_renamed.fq.gz,SRR8809634_trimmed_renamed.fq.gz,SRR8809635_trimmed_renamed.fq.gz,SRR8809636_trimmed_renamed.fq.gz,SRR8809637_trimmed_renamed.fq.gz,SRR8809638_trimmed_renamed.fq.gz,SRR8809639_trimmed_renamed.fq.gz,SRR8809640_trimmed_renamed.fq.gz,SRR8809641_trimmed_renamed.fq.gz,SRR8809642_trimmed_renamed.fq.gz,SRR8809643_trimmed_renamed.fq.gz \
                               --CPU 16 --output ${OUTDIR}
echo ==== Trinity In silico Read Normalization END ====

date
