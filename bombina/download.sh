#!/bin/bash

while read line
do
  echo "${line}"
  file="${line}.txt"
  echo ${file}
  ascp -v -QT -l 400m -P 33001 \
    -k 1 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
    --mode recv --host fasp.sra.ebi.ac.uk --user era-fasp \
    --file-list ${file} ./${line}
  ls -l ./${line} >> $2
done < $1
