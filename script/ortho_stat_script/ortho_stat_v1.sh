#!/bin/bash

function ortho_sum() {
  #in the *_lst/ directory
  study_lst=`ls -1 ./input/`
  for file in `ls -1 ./input`; do
    sp_name=${file%.fasta}
    bioawk -c fastx '{ print "HEADER|"$name }' ./input/$sp_name.fasta > ./output/groups/orgin_${sp_name%_pep}.lst
    sed -i "s/HEADER/${sp_name}/g" ./output/groups/orgin_${sp_name%_pep}.lst
    cat ./output/groups/groups.txt | awk -v species="$sp_name" '{for(i=1; i<=NF; i++) if ($i ~ species) print $i}' > ./output/groups/clustered_${sp_name%_pep}.lst
    cat ./output/groups/orgin_${sp_name%_pep}.lst ./output/groups/clustered_${sp_name%_pep}.lst | sort | uniq -u 
  done | 
  paste -s -d '\n' > ./output/groups/unclustered.lst

  #cat <(cat ./output/groups/groups.txt | awk -v species="$sp_name" '{ for(i=1; i<=NF; i++) if ($i ~ species) print $i }') \
  #  <(bioawk -c fastx '{ print $sp_name$name }' ./input/$sp_name.fasta) | sort | uniq -u;
  #done | \
  #paste -s -d ' ' | \
  #paste -d ' ' <(echo "unclustered:") - | \
  #cat ./output/groups/groups.txt - > ./output/groups/groups.all.txt
}

ortho_sum 2>&1 | tee sum.log

