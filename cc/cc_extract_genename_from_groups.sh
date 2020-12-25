#!/bin/bash
#SBATCH --time=02:00:00                   # Time limit hrs:min:sec
#SBATCH --account=def-heylanda
#SBATCH --job-name=extract_genename
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=ALL                   # Mail events (NONE, BEGIN, END, FAIL..)
#SBATCH --mail-user=ysheng@uoguelph.ca    # Where to send mail
#SBATCH --ntasks=1                        # Run a single task
#SBATCH --output=extract_genename.log             # Standard output and error log
date;hostname;pwd


echo ========= Parameter 1 =============
echo $1
echo ========= Parameter 2 =============
echo $2
echo ========= Parameter 3 =============
echo $3
echo ========= Parameter 4 =============
echo $4
echo ========= Start =============

for line in `cat $1`; do
    awk -v line="$line" '$0 ~ line {print}' $2 | tr " " "\n" > tmp
    awk -v pattern="$3" '$0 ~ pattern {print}' tmp >> $4
done

rm tmp


# sbatch --output=PRJNA316996.exGN.out \
#     cc_extract_genename_from_groups.sh \
#     allsp_g.mouse_23.groups.counts.stat \
#     raw_g.rodent.groups.txt \
#     PRJNA316996 \
#     allsp_g.PRJNA316996.genenames
#
# fly
# insect
  # PRJNA302146 -
  # PRJNA419677 -
  # PRJNA475804 -
#
# mouse
# rodent
  # PRJNA529794 -
  # PRJNA316996
