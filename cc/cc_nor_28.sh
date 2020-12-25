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
$EBROOTTRINITY/trinityrnaseq-Trinity-v2.8.5/util/insilico_read_normalization.pl --seqType fq --JM 100G --max_cov 50 \
                                 --left SRR2105773_1_val_1_renamed.fq.gz,SRR2105774_1_val_1_renamed.fq.gz,SRR2105775_1_val_1_renamed.fq.gz,SRR2105776_1_val_1_renamed.fq.gz,SRR2105777_1_val_1_renamed.fq.gz,SRR2105778_1_val_1_renamed.fq.gz,SRR2105779_1_val_1_renamed.fq.gz,SRR2105780_1_val_1_renamed.fq.gz,SRR2105781_1_val_1_renamed.fq.gz,SRR2105782_1_val_1_renamed.fq.gz,SRR2105783_1_val_1_renamed.fq.gz,SRR2105784_1_val_1_renamed.fq.gz,SRR2105785_1_val_1_renamed.fq.gz,SRR2105786_1_val_1_renamed.fq.gz,SRR2105787_1_val_1_renamed.fq.gz,SRR2105788_1_val_1_renamed.fq.gz,SRR2105789_1_val_1_renamed.fq.gz,SRR2105790_1_val_1_renamed.fq.gz,SRR2105791_1_val_1_renamed.fq.gz,SRR2105792_1_val_1_renamed.fq.gz,SRR2105793_1_val_1_renamed.fq.gz,SRR2105794_1_val_1_renamed.fq.gz,SRR2105795_1_val_1_renamed.fq.gz,SRR2105796_1_val_1_renamed.fq.gz,SRR2105797_1_val_1_renamed.fq.gz,SRR2105798_1_val_1_renamed.fq.gz,SRR2105799_1_val_1_renamed.fq.gz,SRR2105800_1_val_1_renamed.fq.gz,SRR2105801_1_val_1_renamed.fq.gz,SRR2105802_1_val_1_renamed.fq.gz \
                                 --right SRR2105773_2_val_2_renamed.fq.gz,SRR2105774_2_val_2_renamed.fq.gz,SRR2105775_2_val_2_renamed.fq.gz,SRR2105776_2_val_2_renamed.fq.gz,SRR2105777_2_val_2_renamed.fq.gz,SRR2105778_2_val_2_renamed.fq.gz,SRR2105779_2_val_2_renamed.fq.gz,SRR2105780_2_val_2_renamed.fq.gz,SRR2105781_2_val_2_renamed.fq.gz,SRR2105782_2_val_2_renamed.fq.gz,SRR2105783_2_val_2_renamed.fq.gz,SRR2105784_2_val_2_renamed.fq.gz,SRR2105785_2_val_2_renamed.fq.gz,SRR2105786_2_val_2_renamed.fq.gz,SRR2105787_2_val_2_renamed.fq.gz,SRR2105788_2_val_2_renamed.fq.gz,SRR2105789_2_val_2_renamed.fq.gz,SRR2105790_2_val_2_renamed.fq.gz,SRR2105791_2_val_2_renamed.fq.gz,SRR2105792_2_val_2_renamed.fq.gz,SRR2105793_2_val_2_renamed.fq.gz,SRR2105794_2_val_2_renamed.fq.gz,SRR2105795_2_val_2_renamed.fq.gz,SRR2105796_2_val_2_renamed.fq.gz,SRR2105797_2_val_2_renamed.fq.gz,SRR2105798_2_val_2_renamed.fq.gz,SRR2105799_2_val_2_renamed.fq.gz,SRR2105800_2_val_2_renamed.fq.gz,SRR2105801_2_val_2_renamed.fq.gz,SRR2105802_2_val_2_renamed.fq.gz \
                                 --CPU 16 --output ${OUTDIR} --pairs_together
echo ==== Trinity In silico Read Normalization END ====
date
