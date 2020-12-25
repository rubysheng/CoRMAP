#!/bin/bash
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
