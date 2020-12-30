#!/bin/bash

while read line
do
  echo "${line}"
  mkdir ./${line}
done < $1
ls -l
