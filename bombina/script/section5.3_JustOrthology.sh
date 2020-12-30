#!/bin/bash
#title          :section5.3_JustOrthology.sh
#description    :Generate the local database for blast.
#author         :Ruby(Yiru) Sheng
#date           :20191118
#version        :1.1
#usage          :source $RUBY_SCRIPTS/section5.3_JustOrthology.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================

PRJNA_PATH=$(pwd)
mkdir orthosearch_out_dir
cd orthosearch_out_dir
TransDecoder.LongOrfs -t ${PRJNA_PATH}/trinity_out_dir/Trinity.fasta
# TransDecoder.LongOrfs -t ./Trinity.fasta

# ##########################
# # JustOrthologs  Wrapper #
# ##########################
# This wrapper is provided as a courtesy, and may or may not be the fastest or most memory efficient way to extract
# all CDS regions from 2 gff3 and 2 fasta files. Multithreading is only used for justOrthologs,
# and large fasta files are read into memory for the extracting phase. Typical runtime to extract CDS regions from a
# complete genome is under an hour, but it may run longer depending on the size of the genome/gff3 file and the
# processing capabilities of the machine you are running.
#
# ##########################
#
# ARGUMENT OPTIONS:
#
# The wrapper has the following options:
#  -h, --help    show this help message and exit
#  -g1 GFF3_ONE  1st GFF3 (gzip allowed with .gz)
#  -g2 GFF3_TWO  2nd GFF3 fasta file (gzip allowed with .gz)
#  -r1 REF_ONE   1st Reference Genome (gzip allowed with .gz)
#  -r2 REF_TWO   2nd Reference Genome (gzip allowed with .gz)
#  -fa1 FASTA1   1st Fasta file (only used without --e)
#  -fa2 FASTA2   2nd Fasta file (only used without --e)
#  -e            Extract CDS regions from genomes
#  -r            Run JustOrthologs
#  -f            Filters genes based on annotations
#  -s            Sort FASTA file for running JustOrthologs
#  -k            Keep All Temporary Files
#  -d            For Distantly Related Species (only with --r)
#  -c            Combine Both Algorithms In JustOrthologs For Best Accuracy
#  -o OUTPUT     Output File for --r
#  -t THREADS    Number of Cores (only affects -r option)
#  -all          Run --e, --f, --s, and --r
# ##########################
# REQUIREMENTS:
# wrapper.py uses Python version 2.7
#
# The wrapper needs to be used in the same directory as
# getNoException.py, justOrthologs.py, revised_gff3_parser.py, and sortFastaBySeqLen.sh
#
# ##########################
# Input Files:
#
# There are many ways to run this program, and the wrapper will tell you where your output files are created.
# Descriptions of what the algorithms can do are in the argument options.
#
# Typical usage requires 2 fasta files that may or may not be gzipped, and 2 gff3 files that may or may not be gzipped.
#
# ##########################
# USAGE:
#
# The default number of threads is 16. If you want to change that, use the -t option.
# However, the -t option will only affect justOrthologs.py
#
# Example of typical usage:
#
# python2.7 wrapper.py -g1 smallTest/wrapperTest/small_pan.gff3 -g2 smallTest/wrapperTest/small_human.gff3 -r1 smallTest/wrapperTest/small_pan.fasta.gz -r2 smallTest/wrapperTest/small_human.fasta.gz -all -o output
#
# Running the above command will produce a single output file called output in the current directory.
# Temporary files are created and removed, so don't worry if you see other files before the job is complete.
# Running the above command takes about 13 seconds of real time and less than 1 minute of user time on our hardware.
#
# ##########################

file=`basename $(pwd)`
sed -i 's+TRINITY+${file}+g' trinity_out_dir/Trinity.fasta.RSEM.transcripts.fa


TransDecoder.LongOrfs -t PRJNA240970_Trinity.fasta.RSEM.transcripts.fa
TransDecoder.LongOrfs -t PRJNA241010_Trinity.fasta.RSEM.transcripts.fa
mv PRJNA240970_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir PRJNA240970.rsem.transdecoder
mv PRJNA240970_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir.__checkpoints_longorfs PRJNA240970.rsem.transdecoder/transdecoder_dir.__checkpoints_longorfs
mv PRJNA241010_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir PRJNA241010.rsem.transdecoder
mv PRJNA241010_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir.__checkpoints_longorfs PRJNA241010.rsem.transdecoder/transdecoder_dir.__checkpoints_longorfs

cd PRJNA240970
sed -i '/^$/d' PRJNA24*.rsem.transdecoder/longest_orfs.gff3
python2.7 ../JustOrthologs/gff3_parser.py -g PRJNA24*.rsem.transdecoder/longest_orfs.gff3 \
  -f PRJNA24*_Trinity.fasta.RSEM.transcripts.fa -o extract
python2.7 ../JustOrthologs/getNoException.py -i extract -o filter
source ../JustOrthologs/sortFastaBySeqLen.sh filter sort


cd PRJNA241010
sed -i '/^$/d' PRJNA24*.rsem.transdecoder/longest_orfs.gff3
python2.7 ../JustOrthologs/gff3_parser.py -g PRJNA24*.rsem.transdecoder/longest_orfs.gff3 \
  -f PRJNA24*_Trinity.fasta.RSEM.transcripts.fa -o extract
python2.7 ../JustOrthologs/getNoException.py -i extract -o filter
source ../JustOrthologs/sortFastaBySeqLen.sh filter sort

python2.7 JustOrthologs/justOrthologs.py -q PRJNA240970/sort -s PRJNA241010/sort \
  -o PRJNA240970_241010 -c -t 16



# export JO_PATH=~/Documents/test_OSA/JustOrthologs/
# python2.7 ~/Documents/test_OSA/JustOrthologs/wrapper.py \
#   -g1 ./PRJNA240970/orthosearch_out_dir/Trinity.fasta.transdecoder_dir/cds.gff3 \
#   -g2 ./PRJNA241010/orthosearch_out_dir/Trinity.fasta.transdecoder_dir/cds.gff3 \
#   -r1 ./PRJNA240970/PRJNA240970.Trinity.fasta \
#   -r2 ./PRJNA241010/PRJNA241010.Trinity.fasta \
#   -all -t 4 -o JO_240VS241_out
#
# python2.7 ../application/JustOrthologs/wrapper.py \
#   -g1 ./PRJNA240970/cds.gff3 \
#   -g2 ./PRJNA241010/cds.gff3 \
#   -r1 ./PRJNA240970/PRJNA240970.Trinity.fasta \
#   -r2 ./PRJNA241010/PRJNA241010.Trinity.fasta \
#   -f -s -r -t 4 -o JO_240VS241_out
# python2.7 wrapper.py -g1 smallTest/wrapperTest/small_pan.gff3 -g2 smallTest/wrapperTest/small_human.gff3 \
#   -r1 smallTest/wrapperTest/small_pan.fasta.gz -r2 smallTest/wrapperTest/small_human.fasta.gz -all -o output

# ##########################
# #     JustOrthologs      #
# ##########################
# ARGUMENT OPTIONS:
#
#  -h, --help    							show this help message and exit
#  -q	Query Fasta File					A fasta file that has been divided into CDS regions and sorted based on the number of CDS regions.
#  -s	Subject Fasta File					A fasta file that has been divided into CDS regions and sorted based on the number of CDS regions.
#  -o	Output File							A path to an output file that will be created from this program.
#  -t	Number of Cores						The number of Cores available
#  -d	For More Distantly Related Species	A flag that implements the -d algorithm, which gives better accuracy for more distantly related species.
#  -c	Combine Both Algorithms				Combines both the normal and -d algorithms for the best overall accuracy.
#  -r Correlation value					An optional value for changing the required correlation value between orthologs.
#
# ##########################
# REQUIREMENTS:
#
# JustOrthologs uses Python version 2.7
#
# ##########################
# Input Files:
#
# This algorithm requires two fasta files that have header/sequence alternating lines with an asterisk (*)
# between every CDS region and at the end of each sequence. To create this file, look at the example files in
# the directory smallTest, or use the included wrapper to extract all CDS regions from a reference genome and
# a gff3 file.
#
# ##########################
# USAGE:
#
# Typically, the -c option should be used to get the best overall precision and recall.
#
# The default number of threads is the number of threads available. If you want to change that, use the -t option.
#
# The query, subject, and output files  must always be supplied
#
# Example usage:
#
# python justOrthologs.py -q smallTest/orthologTest/bonobo.fa -s smallTest/orthologTest/human.fa -o output -c -t 16
#
# Running the above command will produce a single output file called output in the current directory. For us, this test took
# approximately 4 seconds of real time and 41 seconds of user time.
#
# ##########################


# awk 'NR%2==1{print}NR%2==0{print $0"*"}' test.fa
#
#
# python2.7 justOrthologs.py -q ../PRJNA240970_JO.input.fa -s ../PRJNA241010_JO.input.fa -o output -c -t 6
#
# blastn -query PRJNA240970/blast_out_dir/PRJNA240970*RSEM*.fa -db PRJNA241010_blast -evalue 1e-3 -outfmt 6 -num_threads 6 -out output_blastn -max_target_seqs 1
