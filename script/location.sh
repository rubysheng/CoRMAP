#!/bin/bash

sed '/Trash/d' /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.txt > /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.log
sed -i 's/\s/"\ "/g' /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.log
#sed -i 's/^/"/g;s/$/"/g' /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.log
#sed -i 's/\\s/\s/g' /media/heyland-lab/Seagate\ Backup\ Plus\ Drive/ruby/path.log

while read line
do
  echo "${line}"
  #cd ${line}
done < /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.log

cat /media/heyland-lab/Seagate_Backup_Plus_Drive/ruby/path.log | while read line; do cd $line; done
