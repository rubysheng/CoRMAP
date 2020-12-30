cat output/groups/groups.txt | \
awk 'BEGIN {print "GROUP_NUM", "PRJNA287145", "PRJNA287152", "PRJNA302146", "PRJNA419677", "PRJNA475804", "PRJNA541005"}
{ sp1=sp2=sp3=sp4=sp5=sp6=sp7=sp8=sp9=sp10=0; for(i=2; i<=NF; i++) \
if ($i ~ "PRJNA287145") sp1 += 1; \
else if ($i ~ "PRJNA287152") sp2 += 1; \
else if ($i ~ "PRJNA302146") sp3 += 1; \
else if ($i ~ "PRJNA419677") sp4 += 1; \
else if ($i ~ "PRJNA475804") sp5 += 1; \
else if ($i ~ "PRJNA541005") sp6 += 1; \
print $1, sp1, sp2, sp3, sp4, sp5, sp6}' > ./analyze/groups.counts.txt
