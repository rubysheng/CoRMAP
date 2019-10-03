#!/bin/bash

# load needed files
source ./process.sh # checked: works well
source ./define_type.sh

# find directories path of all projects, and add to a text file
#find /media/heyland-lab/ -type d -name "PRJNA?*" > /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.txt
sed '/Trash/d' /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.txt > /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.log
Seagatelocation="Seagate\ Backup\ Plus\ Drive"
sed -i 's/Seagate Backup Plus Drive/$Seagatelocation/' /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.log


# go to each location directory by read lines of the "path.log" via a for loop
# cd /media/heyland-lab
# ROOT_DIR=$(pwd)
# for line in `cat /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.log`; do
#   echo "${line}"
#   cd ${line}
#   #mainflow ${line} > test_log_file.txt || true
#   # nohup mainflow ${line} &
# done
