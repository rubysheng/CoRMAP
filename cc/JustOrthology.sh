#!/bin/bash
#title          :JustOrthology.sh
#description    :Group orthologs.
#author         :Ruby(Yiru) Sheng
#date           :20191202
#version        :1.2
#usage          :source $RUBY_SCRIPTS/JustOrthology.sh
#notes          :
#bash_version   :4.4.19(1)-release
#============================================================================
OPTIND=1         # Reset in case getopts has been used previously in the shell.
input_dir=""
output_dir=""
end_output=""

while getopts "h?:l:e:o:" opt; do
    case "$opt" in
    h)
        echo "-l is a required list with input directories "
        echo "-o is a required output directory for pairwise JustOrtholog comparisons"
        echo "-e is a required end output file for the combined analysis"
        echo "-h shows this help message"
        echo "press Ctrl+C to exit"
          exit 0
        ;;
    l) list=$OPTARG
        ;;
    e) end_output=$OPTARG
        ;;
    o)  output_dir=$OPTARG
        if [[ "${output_dir: -1}" != '/' ]]; then
            output_dir=${output_dir}"/"
        fi
        mkdir -p ${output_dir}
    esac
done
shift $(( OPTIND-1 ))



function JO_sort () {
  file=`basename $(pwd)`
  FA_NEWNAME="${file}_RSEM.fasta"
  sed "s+TRINITY+${file}+g" trinity_out_dir/Trinity.fasta.RSEM.transcripts.fa > ${FA_NEWNAME}
  TransDecoder.LongOrfs -t ${FA_NEWNAME} -O rsem.transdecoder
  GF_NEWNAME="${file}_RSEM..gff3"
  sed '/^$/d' rsem.transdecoder/longest_orfs.gff3 > ${GF_NEWNAME}
  echo
  echo ==Extract==
  python2.7 $JO_PATH/gff3_parser.py -g ${GF_NEWNAME} \
    -f ${FA_NEWNAME} -o extract
  echo
  echo ==Filter==
  python2.7 $JO_PATH/getNoException.py -i extract -o filter
  echo
  echo ==Sort==
  source $JO_PATH/sortFastaBySeqLen.sh filter sort
  echo
  echo == Output check list ==
  find . -name extract -exec ls -lh {} \;
  find . -name filter -exec ls -lh {} \;
  find . -name sort -exec ls -lh {} \;
}

set -- `cat ${list} | uniq`
for a; do
    cd ${a}
    JO_sort
    file_a="${a}/sort"
    cd ..
    shift
    for b; do
        cd ${b}
        JO_sort
        file_b="${b}/sort"
        cd ..
        echo
        echo "JO seaching input:"
        echo "    "${file_a}
        echo "    "${file_b}
        echo
        echo ==JO searching==
        mkdir -v ${output_dir}/${a}_${b}
        time python2.7 $JO_PATH/justOrthologs.py -q ${file_a} -s ${file_b} \
          -o ${output_dir}/${a}_${b} -c -t 16
        echo
        echo "Check outputs in "${output_dir}/${a}_${b}
    done
done

echo
echo ==Combine Ortho Groups==
python2.7 ${JO_PATH}/combineOrthoGroups.py -id ${output_dir} -o ${end_output}
echo
echo "Go to "${end_output}" and check the Combined Orthologs Groups"
cd ${end_output}
# TransDecoder.LongOrfs -t PRJNA240970_Trinity.fasta.RSEM.transcripts.fa
# TransDecoder.LongOrfs -t PRJNA241010_Trinity.fasta.RSEM.transcripts.fa
# mv PRJNA240970_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir PRJNA240970.rsem.transdecoder
# mv PRJNA240970_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir.__checkpoints_longorfs PRJNA240970.rsem.transdecoder/transdecoder_dir.__checkpoints_longorfs
# mv PRJNA241010_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir PRJNA241010.rsem.transdecoder
# mv PRJNA241010_Trinity.fasta.RSEM.transcripts.fa.transdecoder_dir.__checkpoints_longorfs PRJNA241010.rsem.transdecoder/transdecoder_dir.__checkpoints_longorfs
#
# cd PRJNA240970
# sed -i '/^$/d' PRJNA24*.rsem.transdecoder/longest_orfs.gff3
# python2.7 ../JustOrthologs/gff3_parser.py -g PRJNA24*.rsem.transdecoder/longest_orfs.gff3 \
#   -f PRJNA24*_Trinity.fasta.RSEM.transcripts.fa -o extract
# python2.7 ../JustOrthologs/getNoException.py -i extract -o filter
# source ../JustOrthologs/sortFastaBySeqLen.sh filter sort
#
#
# cd PRJNA241010
# sed -i '/^$/d' PRJNA24*.rsem.transdecoder/longest_orfs.gff3
# python2.7 ../JustOrthologs/gff3_parser.py -g PRJNA24*.rsem.transdecoder/longest_orfs.gff3 \
#   -f PRJNA24*_Trinity.fasta.RSEM.transcripts.fa -o extract
# python2.7 ../JustOrthologs/getNoException.py -i extract -o filter
# source ../JustOrthologs/sortFastaBySeqLen.sh filter sort

# python2.7 JustOrthologs/justOrthologs.py -q PRJNA240970/sort -s PRJNA241010/sort \
#   -o PRJNA240970_241010 -c -t 16



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
