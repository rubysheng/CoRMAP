#!/usr/bin/env bash

while IFS= read -r project; do
  file_path="./"${project}"/trinity_out_dir/Trinity.fasta"
  $TRINITY_HOME/util/TrinityStats.pl  ${file_path}
done <  /home/lewis/Documents/test_OSA/orthomcl_out_dir/5_taxa/analyze/expre_value/pro_acc


# source $RUBY_SCRIPTS/stats.sh 2>&1 | tee /media/lewis/New_Seagate_Drive_8TB/ruby/results/assembly/EN50.stats.log

# project="Bombina"; echo `grep ${project} orgin_bom.lst | wc -l` ` grep  ${project}  clustered_bom.lst | wc -l` ` grep  ${project}  unclustered_bom.lst | wc -l`
#
# echo `head -1 ${file} | awk 'BEGIN {IFS=" "} END {print NF}'`

while IFS= read -r file; do
  echo ${file} `head -1 ${file} | awk 'BEGIN {IFS=" "} END {print NF}'`
done < ../lst

# lst
# ls -lh ./*/processed_data/*
# mkdir /media/lewis/New_Seagate_Drive_8TB/ruby/data/
# while IFS= read -r project; do
#   folder_name="/media/lewis/New_Seagate_Drive_8TB/ruby/data/"${project}"/"
#   mkdir ${folder_name}
#   cd ${project}
#   filename=${project}"_*"
#   cp -r -v ${filename} ${folder_name}
#   cp -v transcripts_count/* ${folder_name}
#   cp -v trinity_out_dir/* ${folder_name}
#   cd ..
# done < lst

# rsync -av --progress ./ ${folder_name} --exclude processed_data
#
# project="PRJNA252803"; folder_name="/media/lewis/New_Seagate_Drive_8TB/ruby/data/"${project}"/"; mkdir ${folder_name}; cd ${project}; filename=${project}"_*" ;cp -r -v ${filename} ${folder_name}; cp -v transcripts_count/* ${folder_name}; cd ..
