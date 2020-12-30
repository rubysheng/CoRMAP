#!/bin/bash
#title          :JustOrthology.sh
#description    :Group orthologs.
#author         :Ruby(Yiru) Sheng
#date           :20191202
#version        :1.1
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

#echo
#echo ==Combine Ortho Groups==
#python2.7 ${JO_PATH}/combineOrthoGroups.py -id ${output_dir} -o ${end_output}
#echo
#echo "Go to "${end_output}" and check the Combined Orthologs Groups"
#cd ${end_output}
