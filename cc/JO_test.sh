#!/bin/bash
#SBATCH --time=20:00:00                   # Time limit hrs:min:sec
#SBATCH --account=def-heylanda
#SBATCH --job-name=ortho
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=ortho_JO.log             # Standard output and error log
date;hostname;pwd

python2.7 ~/Documents/test_OSA/JustOrthologs/wrapper.py -g1 ./PRJNA240970/orthosearch_out_dir/Trinity.fasta.transdecoder_dir/cds.gff3 -g2 ./PRJNA241010/orthosearch_out_dir/Trinity.fasta.transdecoder_dir/cds.gff3 -r1 ./PRJNA240970/PRJNA240970.Trinity.fasta -r2 ./PRJNA241010/PRJNA241010.Trinity.fasta -all -t 16 -o JO_240VS241_out

date
