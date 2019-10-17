#!/bin/bash

while read line
do
  echo "${line}"
  file="${line%.gz?}.gz"
  echo ${file}
  mv ${line} ${file}
done < $1