#!/bin/bash


sed -i 's/[.]/_/g' ../RD_perg.dnaseq.fasta
sed -i 's/[.]/_/g' ../RD_perg.pepseq.fasta

# generate gene_trans_map
$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl ../RD_perg.dnaseq.fasta > ../RD_perg.gene_trans_map

head ../RD_perg.gene_trans_map
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i1
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i2
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i3
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i5
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i6
# g1001_RDmus_PRJNA316996_DN121767_c0_g2	g1001_RDmus_PRJNA316996_DN121767_c0_g2_i8
# g1001_RDmus_PRJNA529794_DN100192_c1_g1	g1001_RDmus_PRJNA529794_DN100192_c1_g1_i3
# g1001_RDrat_PRJNA390522_DN25_c0_g1	g1001_RDrat_PRJNA390522_DN25_c0_g1_i1
# g1001_RDrat_PRJNA451011_DN59327_c1_g2	g1001_RDrat_PRJNA451011_DN59327_c1_g2_i1
# g1001_RDrat_PRJNA451011_DN59327_c1_g2	g1001_RDrat_PRJNA451011_DN59327_c1_g2_i1


# Add the TransDecoder prediction details to the HEADER
grep '>' ../RD_perg.pepseq.fasta | cut -d'_' -f3 | sort | uniq  # check which datasets are included


# PRJNA252803
# PRJNA316996
# PRJNA529794
# PRJNA390522
# PRJNA451011



grep 'PRJNA252803' ../RD_perg.pepseq.fasta > ../RD_PRJNA252803_pephead.lst
grep 'PRJNA316996' ../RD_perg.pepseq.fasta > ../RD_PRJNA316996_pephead.lst
grep 'PRJNA529794' ../RD_perg.pepseq.fasta > ../RD_PRJNA529794_pephead.lst
grep 'PRJNA390522' ../RD_perg.pepseq.fasta > ../RD_PRJNA390522_pephead.lst
grep 'PRJNA451011' ../RD_perg.pepseq.fasta > ../RD_PRJNA451011_pephead.lst


# sed -ibak 's/^[^_]*_//' ../RD_PRJNA*.lst
sed -ibak 's/_/\t/1; s/_/./g; s/[.]/_/1' ../RD_PRJNA*_pephead.lst

head ../RD_PRJNA*_pephead.lst
# ==> ../RD_PRJNA252803_pephead.lst <==
# >g1003	RDmus_PRJNA252803.DN11082.c0.g1.i1.p1
# >g10042	RDmus_PRJNA252803.DN19434.c0.g1.i1.p1
# >g10055	RDmus_PRJNA252803.DN10596.c0.g1.i1.p1
# >g10087	RDmus_PRJNA252803.DN18522.c0.g1.i1.p1
# >g1008	RDmus_PRJNA252803.DN29214.c0.g1.i1.p1
# >g10111	RDmus_PRJNA252803.DN2617.c0.g1.i1.p1
# >g1011	RDmus_PRJNA252803.DN27710.c0.g1.i1.p1
# >g1013	RDmus_PRJNA252803.DN11486.c0.g1.i1.p1
# >g1018	RDmus_PRJNA252803.DN1865.c0.g1.i1.p1
# >g10227	RDmus_PRJNA252803.DN17072.c0.g1.i1.p1
#
# ==> ../RD_PRJNA252803_pephead.lstbak <==
# >g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1_p1
# >g10042_RDmus_PRJNA252803_DN19434_c0_g1_i1_p1
# >g10055_RDmus_PRJNA252803_DN10596_c0_g1_i1_p1
# >g10087_RDmus_PRJNA252803_DN18522_c0_g1_i1_p1
# >g1008_RDmus_PRJNA252803_DN29214_c0_g1_i1_p1
# >g10111_RDmus_PRJNA252803_DN2617_c0_g1_i1_p1
# >g1011_RDmus_PRJNA252803_DN27710_c0_g1_i1_p1
# >g1013_RDmus_PRJNA252803_DN11486_c0_g1_i1_p1
# >g1018_RDmus_PRJNA252803_DN1865_c0_g1_i1_p1
# >g10227_RDmus_PRJNA252803_DN17072_c0_g1_i1_p1

head -2  /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA252803/PRJNA252803_RSEM.fasta.transdecoder.pep
# >RDmus_PRJNA252803.DN0.c0.g1.i1.p1 RDmus_PRJNA252803.DN0.c0.g1~~RDmus_PRJNA252803.DN0.c0.g1.i1.p1  ORF type:complete len:505 (+),score=151.88 RDmus_PRJNA252803.DN0.c0.g1.i1:77-1591(+)
# MSLICSISNEVPEHPCVSPVSNHVYERRLIEKYIAENGTDPINNQPLSEEQLIDIKVAHP

  # in each of the dataset
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
     /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA252803/PRJNA252803_RSEM.fasta.transdecoder.pep > ../tmp.fa
printf "\n" >> ../tmp.fa

while IFS= read -r line; do
  group_num=`echo $line | awk '{print $1}'`
  pattern=`echo $line | awk '{print $2}'`
  echo ${group_num}
  echo ${pattern}
  awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
  pattern=${pattern%.i*}
  pattern_new=`echo $pattern | sed 's/[.]/_/g'`
  c="${group_num}_${pattern_new}"
  sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
  cat ../tmp_line >> ../RD_PRJNA252803.transdecoder.pep
  rm ../tmp_line
done < ../RD_PRJNA252803_pephead.lst

head ../RD_PRJNA252803.transdecoder.pep
# >>g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1_p1 >g1003_RDmus_PRJNA252803_DN11082_c0_g1~~>g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1_p1  ORF type:complete len:973 (-),score=333.13 >g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1:924-3641(-)
# MTSATSPIILKWDPKSLEIRTLTVERLLEPLVTQVTTLVNTSNKGPSGKKKGRSKKAHVLAASVEQATQNFLEKGEQIAKESQDLKEELVAAVEDVRKQGETMRIASSEFADDPCSSVKRGTMVRAARALLSAVTRLLILADMADVMRLLSHLKIVEEALEAVKNATNEQDLANRFKEFGKEMVKLNYVAARRQQELKDPHCRDEMAAARGALKKNATMLYTASQAFLRHPDVAATRANRDYVFKQVQEAIAGISSAAQATSPTDEAKGHTGIGELAAALNEFDNKIILDPMTFSEARFRPSLEERLESIISGAALMADSSCTRDDRRERIVAECNAVRQALQDLLSEYMNNTGRKEKGDPLNIAIDKMTKKTRDLRRQLRKAVMDHISDSFLETNVPLLVLIEAAKSGNEKEVKEYAQVFREHANKLVEVANLACSISNNEEGVKLVRMAATQIDSLCPQVINAALTLAARPQSKVAQDNMDVFKDQWEKQVRVLTEAVDDITSVDDFLSVSENHILEDVNKCVIALQEGDVDTLDRTAGAIRGRAARVIHIINAEMENYEAGVYTEKVLEATKLLSETVMPRFAEQVEVAIEALSANVPQPFEENEFIDASRLVYDGVRDIRKAVLMIRTPEELEDDSDFEQEDYDVRSRTSVQTEDDQLIAGQSARAIMAQLPQEEKAKIAEQVEIFHQEKSKLDAEVAKWDDSGNDIIVLAKQMCMIMMEMTDFTRGKGPLKNTSDVINAAKKIAEAGSRMDKLARAVADQCPDSACKQDLLAYLQRIALYCHQLNICSKVKAEVQNLGGELIVSGLDSATSLIQAAKNLMNAVVLTVKASYVASTKYQKVYGTAAVNSPVVSWKMKAPEKKPLVKREKPEEFQTRVRRGSQKKHISPVQALSEFKAMDSF*

sed -i 's/>g/g/g' ../RD_PRJNA252803.transdecoder.pep

head ../RD_PRJNA252803.transdecoder.pep
# >g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1_p1 g1003_RDmus_PRJNA252803_DN11082_c0_g1~~g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1_p1  ORF type:complete len:973 (-),score=333.13 g1003_RDmus_PRJNA252803_DN11082_c0_g1_i1:924-3641(-)
#MTSATSPIILKWDPKSLEIRTLTVERLLEPLVTQVTTLVNTSNKGPSGKKKGRSKKAHVLAASVEQATQNFLEKGEQIAKESQDLKEELVAAVEDVRKQGETMRIASSEFADDPCSSVKRGTMVRAARALLSAVTRLLILADMADVMRLLSHLKIVEEALEAVKNATNEQDLANRFKEFGKEMVKLNYVAARRQQELKDPHCRDEMAAARGALKKNATMLYTASQAFLRHPDVAATRANRDYVFKQVQEAIAGISSAAQATSPTDEAKGHTGIGELAAALNEFDNKIILDPMTFSEARFRPSLEERLESIISGAALMADSSCTRDDRRERIVAECNAVRQALQDLLSEYMNNTGRKEKGDPLNIAIDKMTKKTRDLRRQLRKAVMDHISDSFLETNVPLLVLIEAAKSGNEKEVKEYAQVFREHANKLVEVANLACSISNNEEGVKLVRMAATQIDSLCPQVINAALTLAARPQSKVAQDNMDVFKDQWEKQVRVLTEAVDDITSVDDFLSVSENHILEDVNKCVIALQEGDVDTLDRTAGAIRGRAARVIHIINAEMENYEAGVYTEKVLEATKLLSETVMPRFAEQVEVAIEALSANVPQPFEENEFIDASRLVYDGVRDIRKAVLMIRTPEELEDDSDFEQEDYDVRSRTSVQTEDDQLIAGQSARAIMAQLPQEEKAKIAEQVEIFHQEKSKLDAEVAKWDDSGNDIIVLAKQMCMIMMEMTDFTRGKGPLKNTSDVINAAKKIAEAGSRMDKLARAVADQCPDSACKQDLLAYLQRIALYCHQLNICSKVKAEVQNLGGELIVSGLDSATSLIQAAKNLMNAVVLTVKASYVASTKYQKVYGTAAVNSPVVSWKMKAPEKKPLVKREKPEEFQTRVRRGSQKKHISPVQALSEFKAMDSF*

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
   /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA316996/PRJNA316996_RSEM.fasta.transdecoder.pep > ../tmp.fa
printf "\n" >> ../tmp.fa

while IFS= read -r line; do
  group_num=`echo $line | awk '{print $1}'`
  pattern=`echo $line | awk '{print $2}'`
  echo ${group_num}
  echo ${pattern}
  awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
  group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
  pattern=${pattern%.i*}
  pattern_new=`echo $pattern | sed 's/[.]/_/g'`
  c="${group_num}_${pattern_new}"
  sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
  cat ../tmp_line >> ../RD_PRJNA316996.transdecoder.pep
  rm ../tmp_line
done < ../RD_PRJNA316996_pephead.lst

sed -i 's/>g/g/g' ../RD_PRJNA316996.transdecoder.pep



awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
   /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA529794/PRJNA529794_RSEM.fasta.transdecoder.pep > ../tmp.fa
printf "\n" >> ../tmp.fa

while IFS= read -r line; do
  group_num=`echo $line | awk '{print $1}'`
  pattern=`echo $line | awk '{print $2}'`
  echo ${group_num}
  echo ${pattern}
  awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
  group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
  pattern=${pattern%.i*}
  pattern_new=`echo $pattern | sed 's/[.]/_/g'`
  c="${group_num}_${pattern_new}"
  sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
  cat ../tmp_line >> ../RD_PRJNA529794.transdecoder.pep
  rm ../tmp_line
done < ../RD_PRJNA529794_pephead.lst

sed -i 's/>g/g/g' ../RD_PRJNA529794.transdecoder.pep

# grep 'PRJNA390522' ../RD_perg.pepseq.fasta > ../RD_PRJNA390522_pephead.lst
# grep 'PRJNA451011' ../RD_perg.pepseq.fasta > ../RD_PRJNA451011_pephead.lst



awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
   /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA390522/PRJNA390522_RSEM.fasta.transdecoder.pep > ./tmp.fa
printf "\n" >> ./tmp.fa

while IFS= read -r line; do
  group_num=`echo $line | awk '{print $1}'`
  pattern=`echo $line | awk '{print $2}'`
  echo ${group_num}
  echo ${pattern}
  awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
  group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
  pattern=${pattern%.i*}
  pattern_new=`echo $pattern | sed 's/[.]/_/g'`
  c="${group_num}_${pattern_new}"
  sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
  cat ../tmp_line >> ../RD_PRJNA390522.transdecoder.pep
  rm ../tmp_line
done < ../RD_PRJNA390522_pephead.lst

sed -i 's/>g/g/g' ../RD_PRJNA390522.transdecoder.pep



awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' \
   /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA451011/PRJNA451011_RSEM.fasta.transdecoder.pep > ../tmp.fa
printf "\n" >> ../tmp.fa

while IFS= read -r line; do
  group_num=`echo $line | awk '{print $1}'`
  pattern=`echo $line | awk '{print $2}'`
  echo ${group_num}
  echo ${pattern}
  awk -v pat="$pattern" '{ if ($0 ~ pat) { print; getline; print;} }' ../tmp.fa > ../tmp_line
  group_num=`echo $line | awk '{print $1}' | sed 's/>//'`
  pattern=${pattern%.i*}
  pattern_new=`echo $pattern | sed 's/[.]/_/g'`
  c="${group_num}_${pattern_new}"
  sed -i -e "s/${pattern}/${c}/g; s/[.]i/_i/g; s/[.]p/_p/g" ../tmp_line
  cat ../tmp_line >> ../RD_PRJNA451011.transdecoder.pep
  rm ../tmp_line
done < ../RD_PRJNA451011_pephead.lst

sed -i 's/>g/g/g' ../RD_PRJNA451011.transdecoder.pep ????????????????????????????
# sed -i 's/^g/>&/' ../RD_PRJNA451011.transdecoder.pep



cat ../RD_PRJNA*.transdecoder.pep > ../allRD_perg.merge.transdecoder.pep

grep '>' ../RD_perg.dnaseq.fasta | sort | uniq > ../RD_uniq_dna.lst
grep '>' ../RD_perg.pepseq.fasta | sort | uniq > ../RD_uniq_pep.lst
grep '>' ../allRD_perg.merge.transdecoder.pep | sort | uniq > ../RD_uniq_merge.lst

while IFS= read -r line; do
  awk -v pat="$line" '{ if ($0 ~ pat) {print; getline; print;} }' ../RD_perg.dnaseq.fasta | head -2 >> ../RD_uniq_dnaseq.fasta
done < ../RD_uniq_dna.lst
while IFS= read -r line; do
  awk -v pat="$line" '{ if ($0 ~ pat) {print; getline; print;} }' ../RD_perg.pepseq.fasta | head -2 >> ../RD_uniq_pepseq.fasta
  awk -v pat="$line" '{ if ($0 ~ pat) {print; getline; print;} }' ../allRD_perg.merge.transdecoder.pep | head -2 >> ../RD_uniq_mergepepseq.fasta
done < ../RD_uniq_pep.lst
# while IFS= read -r line; do
#   awk -v pat="$line" '{ if ($0 ~ pat) {print; getline; print;} }' ../allRD_perg.merge.transdecoder.pep  >> ../RD_uniq_mergepepseq.fasta
# done < ../RD_uniq_merge.lst


# awk -v pat=">g963_RDrat_PRJNA451011_DN61761_c0_g1_i5_p11" '{ if ($0 ~ pat) {print; getline; print;} }' ../allRD_perg.merge.transdecoder.pep



$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite init --gene_trans_map ../RD_perg.gene_trans_map --transcript_fasta  ../RD_uniq_dnaseq.fasta --transdecoder_pep ../RD_uniq_mergepepseq.fasta   #../RD_perg.dnaseq.fasta \#../allRD_perg.merge.transdecoder.pep
blastx -query ../RD_uniq_dnaseq.fasta -db swissprot_TRIN -num_threads 14 \
    -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > allRD_blastx_TRIN.outfmt6
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastx allRD_blastx_TRIN.outfmt6
blastp -query ../RD_uniq_mergepepseq.fasta -db swissprot_TRIN  -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > allRD_blastp_TRIN.outfmt6     #../RD_perg.pepseq.fasta \
hmmscan --cpu 16 --domtblout allRD_TrinotatePFAM.out \
  Pfam-A.hmm ../RD_uniq_mergepepseq.fasta > allRD_pfam.log  #../RD_perg.pepseq.fasta
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastp allRD_blastp_TRIN.outfmt6
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_pfam allRD_TrinotatePFAM.out
# head -4 {allRD_TrinotatePFAM.out,allRD_blastp_TRIN.outfmt6,allRD_blastx_TRIN.outfmt6,../RD_perg.gene_trans_map,../RD_perg.dnaseq.fasta,../allRD_perg.merge.transdecoder.pep}
$TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite report  > allRD_trinotate_annotation_report.xls

####### Follow steps ##########

    cat ../{RD_PRJNA252803.transdecoder.pep,RD_PRJNA316996.transdecoder.pep,RD_PRJNA529794.transdecoder.pep} > ../RD_perg.merge.transdecoder.pep

    Combine to SQL
    $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite init \
      --gene_trans_map ../RD_perg.gene_trans_map \
      --transcript_fasta  ../RD_perg.dnaseq.fasta \
      --transdecoder_pep ../RD_perg.merge.transdecoder.pep

    Loading BLAST homologies: Load transcript hits
    blastx -query ../RD_perg.dnaseq.fasta -db swissprot_TRIN -num_threads 16 \
        -max_target_seqs 5 -outfmt 6 -evalue 1e-3 > RD_blastx_TRIN.outfmt6
    $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastx RD_blastx_TRIN.outfmt6

    blastp -query ../RD_perg.pepseq.fasta \
      -db swissprot_TRIN \
      -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_TRIN.outfmt6

    # Search Transdecoder-predicted proteins
    blastp -query ../RD_perg.merge.transdecoder.pep \
      -db swissprot_TRIN \
      -num_threads 16 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp_TRIN.outfmt6

    Run hmmscan
    hmmscan --cpu 16 --domtblout TrinotatePFAM.out \
      Pfam-A.hmm ../RD_perg.pepseq.fasta  > RD_pfam.log

    Loading BLAST homologies: Load protein hits
    $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_swissprot_blastp blastp_TRIN.outfmt6
    Load Pfam domain entries
    $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite LOAD_pfam TrinotatePFAM.out

    # Run hmmscan
    hmmscan --cpu 16 --domtblout TrinotatePFAM.out \
      Pfam-A.hmm ../RD_perg.merge.transdecoder.pep  > RD_pfam.log



  head -4 {TrinotatePFAM.out,blastp_TRIN.outfmt6,blastx_TRIN.outfmt6,../RD_perg.gene_trans_map,../RD_perg.dnaseq.fasta,../RD_perg.merge.transdecoder.pep}

  $TRINOTATE_HOME/Trinotate Ruby_transpip.sqlite report  > RD_trinotate_annotation_report.xls
