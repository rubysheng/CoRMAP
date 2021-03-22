#!/usr/bin/env bash

while IFS= read -r line; do
  echo ${line} > 550g_pergroup_14.txt
  group_num=`awk '{print $1}' 550g_pergroup_14.txt`
  group_num=${group_num%:}
  echo ${group_num}
  sed -e "s/ fis/\nfis/g" -e "s/ ins/\nins/g" -e "s/ bom/\nbom/g" -e "s/ mam/\nmam/g" -e "s/ worm/\nworm/g" 550g_pergroup_14.txt | tail -n +2 > 550g_pergroup_14.pepname

  while IFS= read -r line; do
    trimmed=`echo $line | cut -d'|' -f 2`
    awk -v pat="$trimmed" '{ if ($0 ~ pat) {print; getline; print;} }' ../../input/* >> 550g_pergroup_14.pepseq.fasta
    echo ${trimmed%.p*}
    isoform=${trimmed%.p*}
    awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 550g_pergroup_14.dnaseq.fasta
  done < 550g_pergroup_14.pepname

  output="allsp_"$group_num".dnaseq.fasta"
  mv 550g_pergroup_14.dnaseq.fasta ${output}
  output="allsp_"$group_num".pepseq.fasta"
  mv 550g_pergroup_14.pepseq.fasta ${output}
  rm 550g_pergroup_14*
done < allsp_203g_14






# scp ./{allsp_group_6899.pepseq.fasta,file2,file3,file4} /home/usr/destination/
# awk -v pat="mam|RDrat" '{ if ($0 ~ pat) {print; getline; print;} }' allsp_group_6899.pepseq.fasta | head -2

while IFS= read -r line; do
  group_num=`echo $line | cut -d'_' -f 3 | cut -d'.' -f 1`
  input_file="allsp_group_"$group_num".pepseq.fasta"
  awk -v pat="mam|RDrat" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | head -2 | sed "s/>RD/>g${group_num}_RD/g" >> RDrat_perg.pepseq.fasta
done < lst

cd ..

# mv ../RD_perg.pepseq.fasta ../RDmus_missedg.pepseq.fasta
# rm ../RDmus_perg.pepseq.fasta
while IFS= read -r line; do
  group_num=`echo $line | cut -d'_' -f 3 | cut -d'.' -f 1`
  input_file="allsp_group_"$group_num".pepseq.fasta"
  awk -v pat=">RD" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | sed "s/>RD/>g${group_num}_RD/g" >> ../RD_perg.pepseq.fasta
done < lst
while IFS= read -r line; do
  group_num=`echo $line | cut -d'_' -f 3 | cut -d'.' -f 1`
  input_file="allsp_group_"$group_num".dnaseq.fasta"
  awk -v pat=">RD" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | sed "s/>RD/>g${group_num}_RD/g" >> ../RD_perg.dnaseq.fasta
done < lst





while IFS= read -r line; do
  group_num=`echo $line`
  input_file="allsp_group_"$group_num".dnaseq.fasta"
  awk -v pat="RDmus" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | head -2 | sed "s/>RD/>g${group_num}_RD/g" >> ../RDmus_perg.dnaseq.fasta
done < missed_gp



ls -1 allsp_group_*.dnaseq.fasta > lst_dna
# ls -1 | head -n-1 > lst
while IFS= read -r line; do
  group_num=`echo $line | cut -d'_' -f 3 | cut -d'.' -f 1`
  input_file="allsp_group_"$group_num".dnaseq.fasta"
  awk -v pat="mam|RDrat" '{ if ($0 ~ pat) {print; getline; print;} }' $input_file | head -2 | sed "s/>RD/>g${group_num}_RD/g" >> RDrat_perg.dnaseq.fasta
done < lst_dna

# awk -v pat="mam|RDrat" '{ if ($0 ~ pat) {print; getline; print;} }' allsp_group_6899.pepseq.fasta | head -2 | sed 's/>RD/>g6899_RD/g' | sed -e "s/[.]/_/g"
# >g6899_RDrat_PRJNA390522_DN24978_c0_g1_i1_p1
# MASSASVRGLGLRVLACSPELPCAWRALHTSAVCAKNRAARVRVAKGDKPVSYEEAHAPHYIAHRKGWLSQHTGNLDGEDHAAERTLEDVFLRKFMLGTFPGCLADQIVLKRRANQVDICALVLRQLPAHKFYFLVGYSETLLSHFYKCPVRLHLQTVPSKVVYKYI*



# cp ../allsp_groups.spname.txt allsp_10g_groups.txt  # extract all-species-included groups
tail -3045 allsp_10g_groups.txt > allsp_3045g_groups.txt
split -d -l 203 allsp_10g_groups.txt allsp_203g_
# shuf -n 2 ../allsp_groups.txt > allsp_10g_groups.txt  # only test those all-species-included groups

while IFS= read -r line; do
  echo ${line} > 10g_pergroup.txt
  group_num=`awk '{print $1}' 10g_pergroup.txt`
  group_num=${group_num%:}
  echo ${group_num}
  sed -e "s/ fis/\nfis/g" -e "s/ ins/\nins/g" -e "s/ bom/\nbom/g" -e "s/ mam/\nmam/g" -e "s/ worm/\nworm/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname

  while IFS= read -r line; do
    trimmed=`echo $line | cut -d'|' -f 2`
    awk -v pat="$trimmed" '{ if ($0 ~ pat) {print; getline; print;} }' ../../input/* >> 10g_pergroup.pepseq.fasta
    echo ${trimmed%.p*}
    isoform=${trimmed%.p*}
    awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 10g_pergroup.dnaseq.fasta
  done < 10g_pergroup.pepname

  output="allsp_"$group_num".dnaseq.fasta"
  mv 10g_pergroup.dnaseq.fasta ${output}
  output="allsp_"$group_num".pepseq.fasta"
  mv 10g_pergroup.pepseq.fasta ${output}
  rm 10g_pergroup*
done < allsp_10g_groups.txt
cd ../../
