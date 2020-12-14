#!/bin/sh

LIST="$1"
INPUT="$2"

while IFS= read -r line
do
  #echo "line="$line
  grep $line ${INPUT}
  #awk '{ if ($1 == $line) print $0;}' ${INPUT}
  #awk '/"$line"/ {print $0;}' ${INPUT}
  #awk '{ if ($2 == $line) print $0;}' ${INPUT}
done < "${LIST}"
