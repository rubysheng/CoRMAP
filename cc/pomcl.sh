#!/bin/bash
#SBATCH --account=def-heylanda
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --time=15:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ysheng@uoguelph.ca
#SBATCH --job-name=portho
#SBATCH --output=porthomcl-0211.log

date;hostname;pwd

# Load modules:

module load nixpkgs/16.09
module load gcc/7.3.0
module load blast+/2.10.0
module load mcl

# Run scripts

source /home/ruby/projects/def-heylanda/ruby/porthomcl_test/PorthoMCL/porthomcl.sh /home/ruby/projects/def-heylanda/ruby/porthomcl_test/container/

echo "Program finished with exit code $? at: `date`"

exit
