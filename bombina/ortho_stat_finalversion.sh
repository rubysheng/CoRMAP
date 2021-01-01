#!/bin/bash
####################################################################################################################
####################################################################################################################
##### orthomcl results summary ######
####################################################################################################################
####################################################################################################################

function ortho_sum() {
  echo
  echo "================================ Start : ortho_sum() ================================"
  echo "\
  ####################################################################################################################
  # ortho_sum()
  #   Generate (3) lists of total, all clustered, and never clustered peptides' header for each project used in this orthomcl searching.
  #   Input: './output/groups/groups.txt'; './input/*_pep.fasta'
  #   Require: have had a folder under the current directory- './analyze'
  #   Output: 'orgin_\${sp_name%_pep}.lst'
  ####################################################################################################################"

  # in the *_lst/ directory

  study_lst=`ls -1 ./input/*_pep.fasta`
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%.fasta}
    bioawk -c fastx '{ print "HEADER|"$name }' ./input/$sp_name.fasta > ./analyze/orgin_${sp_name%_pep}.lst
    sed -i "s/HEADER/${sp_name}/g" ./analyze/orgin_${sp_name%_pep}.lst
    cat ./output/groups/groups.txt | awk -v species="$sp_name" '{for(i=1; i<=NF; i++) if ($i ~ species) print $i}' > ./analyze/clustered_${sp_name%_pep}.lst
    cat ./analyze/orgin_${sp_name%_pep}.lst ./analyze/clustered_${sp_name%_pep}.lst | sort | uniq -u > ./analyze/unclustered_${sp_name%_pep}.lst
  done

  echo
  echo "================================ End : ortho_sum() ================================"
  echo

}


function ortho_stat() {
  echo
  echo "================================ Start : ortho_stat() ================================"
  echo "\
  ####################################################################################################################
  # ortho_stat()
  #   Returns some stats : 1. amount of all groups;
  #                        2. (run mannually) a table counting peptides from each project in each group
  #                        3. a table of unclustered/clustered peptides sequences for each project
  #   Input: './output/groups/groups.txt'; './input/*_pep.fasta'; './analyze/orgin_\${sp_name%_pep}.lst'; './analyze/clustered_\${sp_name%_pep}.lst'; './analyze/unclustered_\${sp_name%_pep}.lst'
  #   Require: mannually run with 'output_dir_groups_counts-generator.sh'
  #   Output: './output/groups/groups.counts.txt', './analyze/cluster.stat'
  ####################################################################################################################"

  # in the *_lst/ directory

  echo "Path : "$(pwd)
  file=`basename $(pwd)`
  echo "Taxonomy Group Name : "${file%_lst}
  echo

  ####################################################################################################################
  ### step1 : count the total number of clustering groups in a Taxonomy group
  ####################################################################################################################
  TOTAL_GP_COUNT=`grep -c "group" ./output/groups/groups.txt`
  echo "Total clustering groups : "${TOTAL_GP_COUNT}
  echo

  ####################################################################################################################
  ### step2 : counts of proteins for each species in each cluster group
  ## Use the script named 'output_dir_groups_counts-generator.sh'
  ####################################################################################################################
  ## cat ./output/groups/groups.txt | \
    ##   awk 'BEGIN {print "GROUP_NUM", "PRJNA?????1", "PRJNA?????2", "PRJNA?????3", "PRJNA?????4", "PRJNA?????5", "PRJNA?????6", "PRJNA?????7", "PRJNA?????8", "PRJNA?????9", "PRJNA????10"}
  ## { sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
    ##   if ($i ~ "PRJNA?????1") sp1 += 1; \
    ##   else if ($i ~ "PRJNA?????2") sp2 += 1; \
    ##   else if ($i ~ "PRJNA?????3") sp3 += 1; \
    ##   else if ($i ~ "PRJNA?????4") sp4 += 1; \
    ##   else if ($i ~ "PRJNA?????5") sp5 += 1; \
    ##   else if ($i ~ "PRJNA?????6") sp6 += 1; \
    ##   else if ($i ~ "PRJNA?????7") sp7 += 1; \
    ##   else if ($i ~ "PRJNA?????8") sp8 += 1; \
    ##   else if ($i ~ "PRJNA?????9") sp9 += 1; \
    ##   else if ($i ~ "PRJNA????10") sp10 += 1; \
    ##   print $1, sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8, sp9, sp10 }' > output/groups/groups.counts.txt


  ####################################################################################################################
  ### step3 : count unclustered/clustered sequences per project
  ####################################################################################################################

  study_all=`ls -1 input/*_pep.fasta | wc -l`
  echo "Amount Of Species : "${study_all}
  echo
  # generate a output file to hold the starts
  echo  "Species_Acc" "Count_Totalpep" "Count_Clustered" "Count_Unclustered" > ./analyze/cluster.stat  # HEADER
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
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >> ./analyze/cluster.stat
    let i=i+1
  done

  echo
  echo "================================ End : ortho_stat() ================================"
  echo

}


function ortho_dnasum() {
  echo
  echo "================================ Start : ortho_dnasum() ================================"
  echo "\
  ####################################################################################################################
  # ortho_dnasum()
  #   Generate (3) lists of total, all clustered, and never clustered dna sequences' header for each project used in this orthomcl searching.
  #   Input: './analyze/*.lst'
  #   Require: complete previous steps
  #   Output: './analyze/\${file%.lst}\"_dnalst\"'
  ####################################################################################################################"

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

  echo
  echo "================================ End : ortho_dnasum() ================================"
  echo

}


function ortho_dnastat() {
  echo
  echo "================================ Start : ortho_dnastat() ================================"
  echo "\
  ####################################################################################################################
  # ortho_dnastat()
  #   \Returns some stats : a table of unclustered/clustered DNA sequences for each project
  #   \Input: './input/*_pep.fasta'; './analyze/orgin_\${sp_name%_pep}_dnalst'; './analyze/clustered_\${sp_name%_pep}_dnalst'; './analyze/unclustered_\${sp_name%_pep}_dnalst'
  #   \Require: complete previous steps
  #   \Output: './analyze/cluster.dnastat', './analyze/cluster.stat'
  ####################################################################################################################"

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
    echo ${sp_name} ${SP_OR_CNT} ${SP_CL_CNT} ${SP_UN_CNT} >>./analyze/cluster.dnastat
    let i=i+1
  done

  echo
  echo "================================ End : ortho_dnastat() ================================"
  echo

}


function ortho_allsp_g() {
  echo
  echo "================================ Start : ortho_allsp_g() ================================"
  echo "\
  ####################################################################################################################
  # ortho_allsp_g()
  #   Extract the groups that contain peptides from all projects from the peptides counting table, then use these group numbers to find the lines containing peptides headers in these groups from orginal groups output files.
  #   Input: './analyze/stat/groups.counts.stat', './output/groups/groups.txt'
  #   Require: have had a folder under the current directory- './analyze/stat/'; move 'groups.counts.stat' to './analyze/stat/'
  #   Output: './analyze/stat/allsp_g.counts.stat', './analyze/allsp_groups.txt'
  ####################################################################################################################"

  #in the *_lst/ directory
  cd ./analyze/stat/

  # Subset lines of all-project-included groups from the table of all groups peptides counts per project (from 'ortho_stat' function)
  echo `head -1 groups.counts.stat` SUM > allsp_g.counts.stat   # HEADER
  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print $0,sum}' >> allsp_g.counts.stat

  # use the group numbers to find the lines containing peptides headers in 'groups.txt'
  awk '{ for(i=2;i<=NF;i++) { if($i+0 < 1) next }} 1' groups.counts.stat | awk '{ print $1 }' > tmp.allsp_g.lst
  while IFS= read -r line; do
    awk -v pat="$line" '{ if ($0 ~ pat) {print} }' ../../output/groups/groups.txt >> ../allsp_groups.txt
  done < tmp.allsp_g.lst
  rm tmp.allsp_g.lst

  # return to *_lst/
  cd ../../

  echo
  echo "================================ End : ortho_allsp_g() ================================"
  echo

}


function ortho_perg_seq() {
  echo
  echo "================================ Start : ortho_perg_seq() ================================"
  echo "\
  ####################################################################################################################
  # ortho_perg_seq()
  #
  #   Input: ''
  #   Require: have had a folder under the current directory- './analyze/stat/'
  #   Output: './analyze/per_group/'
  ####################################################################################################################"

  #in the *_lst/ directory

  cd analyze/
  if [ ! -d per_group/ ]; then
    mkdir per_group
  fi

  cd per_group
  cp ../allsp_groups.txt all_per_groups.txt  # extract all-species-included groups
  while IFS= read -r line; do
    echo ${line} > 10g_pergroup.txt
    group_num=`awk '{print $1}' 10g_pergroup.txt`
    group_num=${group_num%:}
    echo ${group_num}
    sed "s/ P/\nP/g" 10g_pergroup.txt | tail -n +2 > 10g_pergroup.pepname

    while IFS= read -r line; do
      trimmed=`echo $line | cut -d'|' -f 2`
      awk -v pat="$trimmed" '{ if ($0 ~ pat) {print; getline; print;} }' ../../input/* >> 10g_pergroup.pepseq.fasta
    done < 10g_pergroup.pepname

    output="allsp_"$group_num".pepseq.fasta"
    mv 10g_pergroup.pepseq.fasta ${output}
    rm 10g_pergroup*
  done < all_per_groups.txt
  cd ../../

  echo
  echo "================================ End : ortho_perg_seq() ================================"
  echo

}


######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################

function ortho_selected_stat() {
  echo
  echo "================================ Start : ortho_selected_stat() ================================"
  echo "\
  ####################################################################################################################
  # ortho_selected_stat()
  #   Returns some stats : (run mannually) a table counting peptides from each project in each group, then pick up specific projects by columns from the table, finally filter the groups containing at least one seq from each selected group
  #   Input: './output/groups/groups.txt'
  #   Require: mannually run with 'output_dir_groups_counts-generator.sh'
  #   Output: './output/groups/groups.counts.txt'
  ####################################################################################################################"

  ## In the rodent_lst/ directory

  echo "Path : "$(pwd)
  file=`basename $(pwd)`
  echo "Taxonomy Group Name : "${file%_lst}
  echo

  # source output_dir_groups_counts-generator.sh

  cd ./analyze/stat/
  ## Run mannually with the following code block by change the numbers of columns
  # cat groups.counts.stat | cut -d' ' -f 1,4,5,6 > mus_13.groups.counts.stat
  mv mus_13.groups.counts.stat ../../../mus_1_3/analyze/stat/groups.counts.stat

  ## Go to the new subset folder
  cd ../../../mus_1_3/analyze/stat/
  echo `head -1 groups.counts.stat` SUM > allsp_g.counts.stat
  awk '{for(i=2;i<=NF;i++){if($i+0 < 1) next}} 1' groups.counts.stat | awk '{sum=0; for (i=2; i<=NF; i++) { sum+= $i } print $0,sum}' >> allsp_g.counts.stat

  ## Use the group numbers to find the lines containing peptides headers in 'groups.txt'
  awk '{ for(i=2;i<=NF;i++) { if($i+0 < 1) next }} 1' groups.counts.stat | awk '{ print $1 }' > tmp.allsp_g.lst
  while IFS= read -r line; do
    awk -v pat="$line" '{ if ($0 ~ pat) {print} }' ../../output/groups/groups.txt >> ../allsp_groups.txt
  done < tmp.allsp_g.lst
  rm tmp.allsp_g.lst

  ## Return to mus_1_3/
  cd ../../

  echo
  echo "================================ End : ortho_selected_stat() ================================"
  echo

}

function ortho_selected_sum() {
  echo
  echo "================================ Start : ortho_selected_sum() ================================"
  echo "\
  ####################################################################################################################
  # ortho_selected_sum()
  #   For those selected projects, after obtain the list of all clustered peptides' group number, header in 'allsp_groups.txt', and pepname for each selected project used in this orthomcl searching.
  #   Input: './analyze/allsp_groups.txt'; './input/*_pep.fasta'
  #   Require: have had a folder under the current directory- './analyze'
  #   Output: './analyze/clustered_\${sp_name%_pep}.lst'
  ####################################################################################################################

  # [Modify] in the mus_1_3/ directory
  # mkdir input
  # mv -v ../rodent_lst/input/PRJNA252803_pep.fasta ./input
  # mv -v ../rodent_lst/input/PRJNA529794_pep.fasta ./input"


  # study_lst=`ls -1 ./input/*_pep.fasta`
  for line in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $line`
    sp_name=${file%.fasta}
    cat ./analyze/allsp_groups.txt | awk -v species="$sp_name" '{for(i=1; i<=NF; i++) if ($i ~ species) print $1, $i}'   > tmp.lst
    cut -d'|' -f2 tmp.lst > tmp.pepname
    paste -d " " tmp.lst tmp.pepname > ./analyze/clustered_${sp_name%_pep}.lst
    rm tmp.lst tmp.pepname
  done

  echo
  echo "================================ End : ortho_selected_sum() ================================"
  echo


}

# ortho_selected_sum

function ortho_extract_seq() {
  echo
  echo "================================ Start : ortho_extract_seq() ================================"
  echo "\
  ####################################################################################################################
  # ortho_extract_seq()
  #   Returns two fasta files of peptides sequences and dna sequences
  #   Input: './analyze/clustered_\${sp_name%_pep}.lst'
  #   Require:
  #   Output: './output/groups/groups.counts.txt'
  ####################################################################################################################"

  ## In the rodent_lst/ directory
  study_lst=`ls -1 ./input/*_pep.fasta`
  for file in `ls -1 ./input/*_pep.fasta`; do
    file=`basename $file`
    sp_name=${file%.fasta}
    echo "========================================================================="
    echo ${sp_name%_pep}" is processing with extracting sequences from orthologs."
    rm -v ./analyze/${sp_name%_pep}.dnaseq.fasta ./analyze/${sp_name%_pep}.pepseq.fasta

    while IFS= read -r line; do
      group_num=`echo $line | cut -d' ' -f1`
      pepname=`echo $line | cut -d' ' -f3`
      awk -v gp="$group_num" -v pat="$pepname" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./input/${file} >> ./analyze/${sp_name%_pep}.pepseq.fasta
      isoform=${pepname%.p*}
      awk -v gp="$group_num" -v pat="$isoform" '{ if ($0 ~ pat) {print gp,$0 ; getline; print;} }' ./fasta_dir/${file%_pep.fasta}_RSEM.fasta >> ./analyze/${sp_name%_pep}.dnaseq.fasta
    done < ./analyze/clustered_${sp_name%_pep}.lst

    echo ${sp_name%_pep}" is end with processing of extracting sequences from orthologs."
    echo "========================================================================="
  done

  echo
  echo "================================ End : ortho_extract_seq() ================================"
  echo

}

# ortho_extract_seq























function main() {
  mkdir ./analyze/

  # check if the fasta files prepared
  if [ ! -d fasta_dir/ ]; then
    mkdir fasta_dir
    for file in `ls -1 ./input/*_pep.fasta`; do
      file=`basename $file`
      file=${file%_pep.fasta}"_RSEM.fasta"
      echo ${file}
      find $path_to_external_drive/Seagate_Backup_Plus_Drive/ruby/data/ -name ${file} -exec cp -v {} fasta_dir/ \;
    done
  fi

  # step1
  #sudo chmod 777 output/groups/*
  #ortho_sum 2>&1 | tee ortho_sum.log
  #ortho_stat 2>&1 | tee ortho_stat.log

  #mkdir ./analyze/stat/
  # mv ./analyze/groups.counts.txt ./analyze/stat/groups.counts.stat
  # mv ./analyze/cluster.stat ./analyze/stat/cluster.stat

  # ortho_dnasum
  # ortho_dnastat

  # step2
  # ortho_allsp_g

  # step3
  #ortho_perg_seq

  # step4
  # cd ./analyze/per_group
  # source /media/lewis/New_Seagate_Drive_8TB/ruby/github/bombina/comparative-transcriptomic-analysis-pip/script/ortho_check.sh
  # cd ../..
}

# if [ "${1}" != "--source-only" ]; then
#     main "${@}"
# fi
