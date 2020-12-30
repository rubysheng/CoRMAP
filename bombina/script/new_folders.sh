#!/bin/bash

while read line
do
  echo "${line}"
  mkdir ./${line}
done < names.txt
ls -l

#cat /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.log | while read line; do cd $line; done

