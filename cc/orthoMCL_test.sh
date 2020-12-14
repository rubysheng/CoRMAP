#!/bin/sh
#SBATCH --time=0-20:00:00                   # Time limit day-hrs:min:sec
#SBATCH --account=def-heylanda
#SBATCH --job-name=test_orthomcl
#SBATCH --mem=100G
#SBATCH --cpus-per-task=500
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=normalize.log             # Standard output and error log
date;hostname;pwd


export VERSION=1.13.4 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

# https://dl.google.com/go/go1.13.4.linux-amd64.tar.gz
# /usr/local/go

echo 'export GOPATH=/usr/local/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

module load trinity
module load samtools
module load jellyfish
module load gcc/5.4.0
module load trinity
module load salmon/0.14.0
module load jellyfish/2.2.6

module load singularity/2.5
