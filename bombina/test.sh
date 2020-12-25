#!/bin/bash
# How many ortholog clusters have a single gene in each species (full 1:1 orthologs)?
cat groups.counts.stat | \
awk '{ if ($2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 && $6 == 1  && $7 == 1 && $8 == 1 && $9 == 1 && $10 == 1 && $11 == 1) print $0 }' | \
wc -l
# How many have a single gene in at least two of the species (1:1 orthologs)?
cat groups.counts.stat | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0)  && ($5 == 1 || $5 == 0)  && ($6 == 1 || $6 == 0) && ($7 == 1 || $7 == 0) && ($8 == 1 || $8 == 0) && ($9 == 1 || $9 == 0) && ($10 == 1 || $10 == 0) && ($11 == 1 || $11 == 0) ) print $0 }' | \
wc -l
# How about how many of the above 1:1 orthologs are found in each genome? (The 1st column is the column number and the 2nd is the count.)
cat groups.counts.stat | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0)  && ($5 == 1 || $5 == 0)  && ($6 == 1 || $6 == 0) && ($7 == 1 || $7 == 0) && ($8 == 1 || $8 == 0) && ($9 == 1 || $9 == 0) && ($10 == 1 || $10 == 0) && ($11 == 1 || $11 == 0) ) print $0 }' | \
awk '{ for (i=2; i<=NF; i++) sum[i] += $i } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1
# Now let’s figure out how many groups contain at least two paralogs that are found in only 1 of the species. (The 1st column is the column number and the 2nd is the count.)
cat groups.counts.stat | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count == NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) sum[i] += $i } END {for (i in sum) print i, sum[i] }' | \
sort -k1,1
# How many groups are expanded (size > 1) in each species? (Note we have to factor out the last row of unclustered genes. The 1st column is the column number and the 2nd is the count.)
cat groups.counts.stat | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count < NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) if ($i > 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] - 1 }' | \
sort -k1,1
# How many groups are expanded (size = 1) in each species? This encompasses 1:1 orthologs from above. (The 1st column is the column number and the 2nd is the count.)
cat groups.counts.stat | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count < NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) if ($i == 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1



# How many ortholog clusters have a single gene in each species (full 1:1 orthologs)?
cat allsp_g.counts.stat | \
awk '{ if ($2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 && $6 == 1  && $7 == 1 && $8 == 1 && $9 == 1 && $10 == 1 && $11 == 1) print $0 }' | \
wc -l
# How many have a single gene in at least two of the species (1:1 orthologs)?
cat allsp_g.counts.stat | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0)  && ($5 == 1 || $5 == 0)  && ($6 == 1 || $6 == 0) && ($7 == 1 || $7 == 0) && ($8 == 1 || $8 == 0) && ($9 == 1 || $9 == 0) && ($10 == 1 || $10 == 0) && ($11 == 1 || $11 == 0) ) print $0 }' | \
wc -l
# How about how many of the above 1:1 orthologs are found in each genome? (The 1st column is the column number and the 2nd is the count.)
cat allsp_g.counts.stat | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0)  && ($5 == 1 || $5 == 0)  && ($6 == 1 || $6 == 0) && ($7 == 1 || $7 == 0) && ($8 == 1 || $8 == 0) && ($9 == 1 || $9 == 0) && ($10 == 1 || $10 == 0) && ($11 == 1 || $11 == 0) ) print $0 }' | \
awk '{ for (i=2; i<=(NF-1); i++) sum[i] += $i } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1
# Now let’s figure out how many groups contain at least two paralogs that are found in only 1 of the species. (The 1st column is the column number and the 2nd is the count.)
cat allsp_g.counts.stat | \
awk '{ count=0; for (i=2; i<=(NF-1); i++) if ($i == 0) count += 1; if (count == NF - 3) print $0 }' |  \
awk '{ for (i=2; i<=(NF-1); i++) sum[i] += $i } END {for (i in sum) print i, sum[i] }' | \
sort -k1,1
# How many groups are expanded (size > 1) in each species? (Note we have to factor out the last row of unclustered genes. The 1st column is the column number and the 2nd is the count.)
cat allsp_g.counts.stat | \
awk '{ count=0; for (i=2; i<=(NF-1); i++) if ($i == 0) count += 1; if (count < NF - 3) print $0 }' |  \
awk '{ for (i=2; i<=(NF-1); i++) if ($i > 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] - 1 }' | \
sort -k1,1
# How many groups are expanded (size = 1) in each species? This encompasses 1:1 orthologs from above. (The 1st column is the column number and the 2nd is the count.)
cat allsp_g.counts.stat | \
awk '{ count=0; for (i=2; i<=(NF-1); i++) if ($i == 0) count += 1; if (count < NF - 3) print $0 }' |  \
awk '{ for (i=2; i<=(NF-1); i++) if ($i == 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1



# count not all species inclued groups
cat groups.counts.stat | awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count == 1) print $0 }' | wc -l # if (count == 1 to 9) print $0
# count a species appearing among groups
