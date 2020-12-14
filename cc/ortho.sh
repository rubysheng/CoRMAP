#!/bin/bash
#SBATCH --account=def-heylanda
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ysheng@uoguelph.ca
#SBATCH --job-name=orthomcl
#SBATCH --output=orthomcl.log

date;hostname;pwd

# Load modules:
module --force purge
module load nixpkgs/16.09 intel/2017.5
module load icc/.2016.4.258 impi/2017.4.239
# module load OrthoMCL/2.0.9-custom-Perl-5.24.0
module load orthomcl/2.0.9
module load gcc/7.3.0 blast+/2.4.0
module load gcccore/.5.4.0 libxml2/2.9.4

echo "Working directory: `pwd`"
echo "Hostname: `hostname`"
echo "Starting run at: `date`"

# Download source and configure
cd ~
git clone https://github.com/apetkau/orthomcl-pipeline.git
cd orthomcl-pipeline

perl scripts/orthomcl-pipeline-setup.pl # set paths to dependencies
cat etc/orthomcl-pipeline.conf # parameters in this file can be adjusted; consult the instruction linked above

# Testing
export PATH=~/orthomcl-pipeline/bin:~/orthomcl-pipeline/scripts:$PATH
perl t/test_pipeline.pl -m /home/ruby/scratch/orthomcl.conf \
    -s fork -t ~/tmp # replace the path to orthomcl.config with your own




echo "Program finished with exit code $? at: `date`"

exit



    # ---
    # blast:
    #   F: 'm S'
    #   b: '100000'
    #   e: '1e-5'
    #   v: '100000'
    # filter:
    #   max_percent_stop: '20'
    #   min_length: '10'
    # mcl:
    #   inflation: '1.5'
    # path:
    #   blastall: /opt/software/BLAST/2.2.26-Linux_x86_64/bin/blastall
    #   formatdb: /opt/software/BLAST/2.2.26-Linux_x86_64/bin/formatdb
    #   mcl: /opt/software/MCL/14.137-intel-2016b/bin/mcl
    #   orthomcl: /opt/software/OrthoMCL/orthomclsoftware-custom/bin
    # scheduler: fork
    # split: '4'


# ---------------------------------------------------------------------------
# singularity pull --name orthomcl.simg shub://ISU-HPC/orthomcl
#
# mkdir -p original complaintFasta
# mv input/* original/
# singularity shell orthomcl.simg
# cd complaintFasta
# for fasta in ../original/*.pep; do
# orthomclAdjustFasta $(basename ${fasta%.*}) ${fasta} 1
# done
#
# cd ..
# orthomclFilterFasta complaintFasta 10 20
