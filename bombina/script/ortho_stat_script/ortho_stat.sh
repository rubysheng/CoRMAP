#!/bin/bash
#=======================================================================================================
############################
# orthomcl results summary #
############################

mkdir ./analyze/
function ortho_sum() {
  #in the *_lst/ directory
  study_lst=`ls -1 ./input/*_pep.fasta`
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%.fasta}
    bioawk -c fastx '{ print "HEADER|"$name }' ./input/$sp_name.fasta > ./analyze/orgin_${sp_name%_pep}.lst
    sed -i "s/HEADER/${sp_name}/g" ./analyze/orgin_${sp_name%_pep}.lst
    cat ./output/groups/groups.txt | awk -v species="$sp_name" '{for(i=1; i<=NF; i++) if ($i ~ species) print $i}' > ./analyze/clustered_${sp_name%_pep}.lst
    cat ./analyze/orgin_${sp_name%_pep}.lst ./analyze/clustered_${sp_name%_pep}.lst | sort | uniq -u > ./analyze/unclustered_${sp_name%_pep}.lst
    #cat ./analyze/unclustered_${sp_name%_pep}.lst
  done
  # | paste -s -d '\n' > ./analyze/unclustered.lst
}

function ortho_stat() {
  # in the *_lst/ directory
  echo "Path : "$(pwd)
  file=`basename $(pwd)`
  echo "Taxonomy Group Name : "${file%_lst}
  echo

  ### step1 : count the total number of clustering groups in a Taxonomy group
  TOTAL_GP_COUNT=`grep -c "group" output/groups/groups.txt`
  echo "Total clustering groups : "${TOTAL_GP_COUNT}
  echo

  ### step2 : counts of proteins for each species in each cluster group
###########
# cat output/groups/groups.txt | \
# awk 'BEGIN {print "GROUP_NUM", "PRJNA?????1", "PRJNA?????2", "PRJNA?????3", "PRJNA?????4", "PRJNA?????5", "PRJNA?????6", "PRJNA?????7", "PRJNA?????8", "PRJNA?????9", "PRJNA????10"}
# { sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
# if ($i ~ "PRJNA?????1") sp1 += 1; \
# else if ($i ~ "PRJNA?????2") sp2 += 1; \
# else if ($i ~ "PRJNA?????3") sp3 += 1; \
# else if ($i ~ "PRJNA?????4") sp4 += 1; \
# else if ($i ~ "PRJNA?????5") sp5 += 1; \
# else if ($i ~ "PRJNA?????6") sp6 += 1; \
# else if ($i ~ "PRJNA?????7") sp7 += 1; \
# else if ($i ~ "PRJNA?????8") sp8 += 1; \
# else if ($i ~ "PRJNA?????9") sp9 += 1; \
# else if ($i ~ "PRJNA????10") sp10 += 1; \
# print $1, sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10 }' > output/groups/groups.counts.txt
###########

  ### step3 : count unclustered/clustered sequences per species
  study_all=`ls -1 input/*_pep.fasta | wc -l`
  echo "Amount Of Species : "${study_all}
  echo
  # generate a output file to hold the starts
  echo  "Species_Acc" "Count_Totalpep" "Count_Clustered" "Count_Unclustered" > output/groups/cluster.stat  # HEADER
  i=1
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%_pep.fasta}
    echo "SPECIES "$i " Name : "${sp_name}
    echo
    # amount of all pep sequences (per species)
    SP_OR_CNT=`wc -l < ./analyze/orgin_${sp_name%_pep}.lst`
    echo "SPECIES "$i " total pep sequences : "${SP_OR_CNT}
    echo
    # amount of clustered pep sequences (per species)
    SP_CL_CNT=`wc -l < ./analyze/clustered_${sp_name%_pep}.lst`
    echo "SPECIES "$i " clustered pep sequences : "${SP_CL_CNT}
    echo
    # amount of unclustered pep sequences (per species)
    SP_UN_CNT=`wc -l < ./analyze/unclustered_${sp_name%_pep}.lst`
    echo "SPECIES "$i " unclustered pep sequences : "${SP_UN_CNT}
    echo
    # hold the stats to output file
    # awk -v name="${sp_name}" -v total="${SP_OR_CNT}" -v clustered="${SP_CL_CNT}" -v unclustered="${SP_UN_CNT}" '{ print name,total, clustered, unclustered }' >> cluster.stat
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >> ./analyze/cluster.stat
    let i=i+1
  done
}


function ortho_dnasum() {
  #in the *_lst/ directory
  cd ./analyze/
  for file in `ls -1 *.lst`; do
    echo $file
    dna_lst_name=${file%.lst}"_dnalst"
    echo $dna_lst_name
    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      echo ${trimmed%.p*}
      isoform=${trimmed%.p*}
      echo "$isoform" >> tmp_dnalst
    done < "$file"
    cat tmp_dnalst | sort | uniq > "$dna_lst_name"
    rm tmp_dnalst
  done
  cd ..
}

function ortho_dnastat() {
  #in the *_lst/ directory
  # generate a output file to hold the starts
  echo  "Species_Acc" "Count_Totaldna" "Count_Clustered" "Count_Unclustered" > ./analyze/cluster.dnastat  # HEADER
  i=1
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%_pep.fasta}
    echo "SPECIES "$i " Name : "${sp_name}
    echo
    # amount of all dna sequences (per species)
    SP_OR_CNT=`wc -l < ./analyze/orgin_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " total dna sequences : "${SP_OR_CNT}
    echo
    # amount of clustered dna sequences (per species)
    SP_CL_CNT=`wc -l < ./analyze/clustered_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " clustered dna sequences : "${SP_CL_CNT}
    echo
    # amount of unclustered dna sequences (per species)
    SP_UN_CNT=`wc -l < ./analyze/unclustered_${sp_name%_pep}_dnalst`
    echo "SPECIES "$i " unclustered dna sequences : "${SP_UN_CNT}
    echo
    # hold the stats to output file
    # awk -v name="${sp_name}" -v total="${SP_OR_CNT}" -v clustered="${SP_CL_CNT}" -v unclustered="${SP_UN_CNT}" '{ print name,total, clustered, unclustered }' >> cluster.stat
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >>./analyze/cluster.dnastat
    let i=i+1
  done
}

function ortho_allsp_g() {
  #in the *_lst/ directory
  cd ./analyze/stat/
  
  echo `head -1 groups.spname.counts.stat` SUM > allsp_g.spname.counts.stat

  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.spname.counts.stat | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print $0,sum}' >> allsp_g.spname.counts.stat

  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.spname.counts.stat | awk '{ print $1 }' > tmp.allsp_g.spname.lst

  while IFS= read -r line; do
    awk -v pat="$line" '{ if ($0 ~ pat) {print} }' ../../output/groups/groups.txt >> ../allsp_groups.spname.txt
  done < tmp.allsp_g.spname.lst
  rm tmp.allsp_g.spname.lst
  cd ../../
  # echo `head -1 groups.counts.stat` SUM > allsp_g.counts.stat
  # 
  # awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print $0,sum}' >> allsp_g.counts.stat
  # 
  # awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{ print $1 }' > tmp.allsp_g.lst
  # 
  # while IFS= read -r line; do
  #   awk -v pat="$line" '{ if ($0 ~ pat) {print} }' ../../output/groups/groups.txt >> ../allsp_groups.txt
  # done < tmp.allsp_g.lst
  # rm tmp.allsp_g.lst
  # cd ../../
}


function ortho_perg_seq() {
  #in the *_lst/ directory
  # check if the fasta files prepared
  if [ ! -d fasta_dir/ ]; then
    mkdir fasta_dir
    for file in `ls -1 ./input/*_pep.fasta`; do
      file=`basename $file`
      file=${file%_pep.fasta}"_RSEM.fasta"
      echo ${file}
      find /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ -name ${file} -exec cp -v {} fasta_dir/ \;
    done
  fi

  # randomly select 10 (test with 2) groups to test the similarity of sequences in groups
  cd analyze/
  if [ ! -d per_group/ ]; then
    mkdir per_group
  fi

  cd per_group
  # shuf -n 10 ../../output/groups/groups.txt > 10g_groups.txt  # test all groups
  # while IFS= read -r line; do
  #   echo ${line} > 10g_pergroup.txt
  #   group_num=`awk '{print $1}' 10g_pergroup.txt`
  #   group_num=${group_num%:}
  #   echo ${group_num}
  #   sed "s/ P/\nP/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname
  # 
  #   while IFS= read -r line; do
  #     trimmed=`echo $line | cut -d'|' -f 2`
  #     awk -v pat="$trimmed" '{ if ($0 ~ pat) {print; getline; print;} }' ../../input/* >> 10g_pergroup.pepseq.fasta
  #     echo ${trimmed%.p*}
  #     isoform=${trimmed%.p*}
  #     awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 10g_pergroup.dnaseq.fasta
  #   done < 10g_pergroup.pepname
  # 
  #   output=$group_num".dnaseq.fasta"
  #   mv 10g_pergroup.dnaseq.fasta ${output}
  #   output=$group_num".pepseq.fasta"
  #   mv 10g_pergroup.pepseq.fasta ${output}
  #   rm 10g_pergroup*
  # done < 10g_groups.txt

  cp ../allsp_groups.spname.txt allsp_10g_groups.txt  # extract all-species-included groups
  # shuf -n 10 ../allsp_groups.txt > allsp_10g_groups.txt  # only test those all-species-included groups
  while IFS= read -r line; do
    echo ${line} > 10g_pergroup.txt
    group_num=`awk '{print $1}' 10g_pergroup.txt`
    group_num=${group_num%:}
    echo ${group_num}
    sed "s/ P/\nP/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname

    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      awk -v pat="$trimmed" '{ if ($0 ~ pat) {print; getline; print;} }' ../../input/* >> 10g_pergroup.pepseq.fasta
      # echo ${trimmed%.p*}
      # isoform=${trimmed%.p*}
      # awk -v pat="$isoform" '{ if ($0 ~ pat) {print; getline; print;} }' ../../fasta_dir/* >> 10g_pergroup.dnaseq.fasta
    done < 10g_pergroup.pepname

    # output="allsp_"$group_num".dnaseq.fasta"
    # mv 10g_pergroup.dnaseq.fasta ${output}
    output="allsp_"$group_num".pepseq.fasta"
    mv 10g_pergroup.pepseq.fasta ${output}
    rm 10g_pergroup*
  done < allsp_10g_groups.txt
  cd ../../

}

# step1
#sudo chmod 777 output/groups/*
#ortho_sum 2>&1 | tee ortho_sum.log
#ortho_stat 2>&1 | tee ortho_stat.log
#mkdir ./analyze/stat/
## mv ./analyze/groups.counts.txt ./analyze/stat/groups.counts.stat
## mv ./analyze/cluster.stat ./analyze/stat/cluster.stat
# ortho_dnasum
# ortho_dnastat

# step2
ortho_allsp_g

# step3
#ortho_perg_seq

# step4
# cd ./analyze/per_group
# source /media/lewis/New_Seagate_Drive_8TB/ruby/github/bombina/comparative-transcriptomic-analysis-pip/script/ortho_check.sh
# cd ../..

