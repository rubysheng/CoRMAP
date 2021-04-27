#!/bin/bash
#title          :section2.2.1_counts-generator.sh
#description    :Generate every project's orthologous pep counts for each group.
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.2.1_counts-generator.sh
#bash_version   :4.4.19(1)-release
#=======================================================================================================================

mkdir -p -v ./analyze/
touch ./analyze/groups.counts.txt

cat output/groups/groups.txt | \
awk 'BEGIN {print "GROUP_NUM", "PRJNA252803", "PRJNA529794"}
{ sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
if ($i ~ "PRJNA252803") sp1 += 1; \
else if ($i ~ "PRJNA529794") sp2 += 1; \
print $1, sp1, sp2}' > ./analyze/groups.counts.txt
