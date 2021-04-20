#!/bin/bash
#title          :section2.2.1_counts-generator.sh
#description    :Generate every project's orthologous pep counts for each group.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.2.1_counts-generator.sh
#bash_version   :4.4.19(1)-release
#=======================================================================================================================


cat output/groups/groups.txt | \
awk 'BEGIN {print "GROUP_NUM", "PRJNA111111", "PRJNA222222", "PRJNA333333", "PRJNA444444", "PRJNA555555", "PRJNA666666"}
{ sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
if ($i ~ "PRJNA111111") sp1 += 1; \
else if ($i ~ "PRJNA222222") sp2 += 1; \
else if ($i ~ "PRJNA333333") sp3 += 1; \
else if ($i ~ "PRJNA444444") sp4 += 1; \
else if ($i ~ "PRJNA555555") sp5 += 1; \
else if ($i ~ "PRJNA666666") sp6 += 1; \
print $1, sp1, sp2, sp3, sp4, sp5, sp6}' > ./analyze/groups.counts.txt
