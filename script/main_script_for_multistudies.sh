#!/bin/bash

# load needed files
export RUBY_GITHUB=/media/lewis/Seagate_Backup_Plus_Drive/ruby/github/comparative-transcriptomic-analysis-pip
source $RUBY_GITHUB/script/process.sh
source $RUBY_GITHUB/script/define_type.sh


# find directories path of all projects, and add to a text file
find /media/lewis/Seagate_Backup_Plus_Drive/ruby -type d -name "PRJNA?*" > /media/lewis/Seagate_Backup_Plus_Drive/ruby/path.txt
sed '/Trash/d' /media/lewis/Seagate_Backup_Plus_Drive/ruby/path.txt > /media/lewis/Seagate_Backup_Plus_Drive/ruby/path.log

# go to each location directory by read lines of the "path.log" via a for loop
# cd /media/lewis/Seagate_Backup_Plus_Drive/ruby
# ROOT_DIR=$(pwd)
# for line in `cat /media/lewis/Seagate_Backup_Plus_Drive/ruby/path.log`; do
#   echo "${line}"
#   cd ${line}
#   #mainflow ${line} > test_log_file.txt || true
#   # nohup mainflow ${line} &
# done
