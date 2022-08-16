#!/bin/bash
#title          :section0.3_download.sh
#description    :download raw RNA-Seq data from SRA by links and store by dataset names
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section0.2_download.sh list_of_projectnames
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

# check if exist folder ./download_input/
if [[ ! -d "./download_input" ]]; then
  source $SCRIPT/section0.1_setup.sh $1
fi

while read line
do
  echo "${line}"
  file="./download_input/"${line}"_download_links"
  echo ${file}
  if [[ -f "$FILE" ]]; then
    ascp -v -QT -l 400m -P 33001 \
      -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
      --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp \
      --file-list ${file} ./data/${line}
  else
    echo "There is no link file for download named as "${file}
  fi
done < $1
