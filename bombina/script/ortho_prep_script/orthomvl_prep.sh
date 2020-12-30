#!/usr/bin/bash

# input:
#   a list of target directories with absolute paths which can find *fasta.RSEM.transcripts.fa
#   example:
#       INdme /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA475804
#       INcgl /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA289731
#       INngi /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA241010
#       INnvi /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA240970

OPTIND=1         # Reset in case getopts has been used previously in the shell.
input_dir=""
output_dir=""
list=""

while getopts "h?:l:o:" opt; do
    case "$opt" in
    h)
        echo "-l is a required list with input directories (to find the raw FASTA files) "
        echo "-o is a required output directory"
        echo "-h shows this help message"
        echo "press Ctrl+C to exit"
          # exit 0
        ;;
    l) list=$OPTARG
        ;;
    # e) end_output=$OPTARG
    #     ;;
    o) output_dir=$OPTARG
        if [[ "${output_dir: -1}" != '/' ]]; then
            output_dir=${output_dir}"/"
        fi
        #mkdir -p ${output_dir}
    esac
done
shift $(( OPTIND-1 ))



function preprocessFA() {
  TAX_CODE="$1"
  OUT_DIR="$2"
  file=`basename $(pwd)`
  FA_NEWNAME="${file}_RSEM.fasta"
  if [ ! -e ${FA_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Renaming Trinity.fasta.RSEM.transcripts.fa to "${FA_NEWNAME}
    echo
    sed "s+TRINITY+${file}+g" trinity_out_dir/Trinity.fasta.RSEM.transcripts.fa > ${FA_NEWNAME}
    sed -i "s/_/./g" ${FA_NEWNAME}
  fi
  #echo "Adding the taxon code"
  #sed -i "s/>/>${TAX_CODE}_/g" ${FA_NEWNAME}
  #echo "Translate to protein sequence"
  #echo
  #TransDecoder.LongOrfs -t ${FA_NEWNAME}
  #TransDecoder.Predict -t ${FA_NEWNAME}
  echo "Moving to the input directory"
  PE_NAME=${FA_NEWNAME}".transdecoder.pep"
  awk '/>/ { print $1; getline; print $0 }' ${PE_NAME} > tmp.pep
  mv tmp.pep ${PE_NAME}
  cp -v ${PE_NAME} ${OUT_DIR}
}

if [ ! "${list}" = "" ]; then # have a input list
  input_dir="$(pwd)/input/"
  base_dir="$(pwd)"
  checkarray=()

  if [ ! -e ${input_dir} ]; then
    echo "Creatng the input directory"
    mkdir -v ${input_dir}
  fi

  while IFS=" " read -r f1 f2
  do
    printf 'Taxonomy_code: %s, Dir: %s\n' "$f1" "$f2"
    echo
    SP_CODE="$f1"
    DIR="$f2"
    cd ${DIR}
    preprocessFA ${SP_CODE} ${input_dir}
    echo
    echo "Check if the file is well - prepared in "${input_dir}
    cd ${base_dir}
    mvd_pep="input/"`basename ${DIR}`"_RSEM.fasta.transdecoder.pep"
    if [ ! -e ${mvd_pep} ]; then  # test if the pep file is moved to the destination
      echo "No file named "${mvd_pep}" moved to input directory"
    else
      find ${input_dir} -name ${mvd_pep} -exec ls -l {} \;
      acc=`basename ${DIR}`
      checkarray+=("${acc}")
    fi
  done < "${list}"
  count_l=`wc -l < ${list}`
fi

#if [ ! "${output_dir}" = "" ] && [ "${#checkarray[@]}" = "${count_l}" ]; then
  # if [ -d "${output_dir}" ]； then
  #   sudo rm -r ${output_dir}
  # fi
 # printf "pluteus123\n" | sudo -S perl5.22.1 $ORTHMCL_PIP/scripts/orthomcl-pipeline.pl \
  #  -i input/ -o ${output_dir} -m $ORTHMCL_PIP/orthomcl.conf --nocompliant
#fi

#sudo -S perl5.22.1 $ORTHMCL_PIP/scripts/orthomcl-pipeline.pl -i input/ -o output/ -m $ORTHMCL_PIP/orthomcl.conf --nocompliant

#head -n 1 *.fasta |grep "^>"

#sudo -S perl5.22.1 /home/lewis/Documents/test_OSA/orthomcl-pipeline/scripts/orthomcl-pipeline_nomysql.pl -i input/ -o output/ -m orthomcl.conf -c orthomcl-pip.conf --nocompliant

