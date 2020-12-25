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
  echo "Adding the taxon code"
  sed -i "s/>/>${TAX_CODE}_/g" ${FA_NEWNAME}
  echo "Translate to protein sequence"
  echo
  TransDecoder.LongOrfs -t ${FA_NEWNAME}
  TransDecoder.Predict -t ${FA_NEWNAME}
  echo "Moving to the input directory"
  PE_NAME=${FA_NEWNAME}".transdecoder.pep"
  awk '/>/ { print $1; getline; print $0 }' ${PE_NAME} > tmp.pep
  mv tmp.pep ${PE_NAME}
  FA_NEWNAME="${file}_RSEM_pep.fasta"
  awk '!/^>/ { printf "%s", $0; n = "\n" }
        /^>/ { print n $0; n = "" }
       END { printf "%s", n }
      '  ${PE_NAME} > ${FA_NEWNAME}
  mv -v ${FA_NEWNAME} ${OUT_DIR}
}

if [ ! "${list}" = "" ]; then # have a input list
  input_dir="$(pwd)/input"
  base_dir="$(pwd)"
  checkarray=()
  while IFS=" " read -r f1 f2
  do
    printf 'Taxonomy_code: %s, Dir: %s\n' "$f1" "$f2"
    echo
    SP_CODE="$f1"
    DIR="$f2"
    cd ${DIR}
    #preprocessFA ${SP_CODE} ${input_dir}
    echo
    echo "Check if the file is well - prepared in "${input_dir}
    cd ${base_dir}
    mvd_pep="input/"`basename ${DIR}`"_RSEM_pep.fasta"
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

if [ ! "${output_dir}" = "" ] && [ "${#checkarray[@]}" = "${count_l}" ]; then
  # if [ -d "${output_dir}" ]； then
  #   sudo rm -r ${output_dir}
  # fi
  # orthomcl-pipeline -i input/ -o output/ -m orthomcl.config -c orthomcl-pipeline.conf
	# 	Runs orthomcl using the given input/output directories.  Overrides parameters (blast, etc)
	# 	from file orthomcl-pipeline.conf.
  printf "pluteus123\n" | sudo -S perl5.22.1 $ORTHMCL_PIP/scripts/orthomcl-pipeline.pl \
    -i input/ -o ${output_dir} -m $ORTHMCL_PIP/orthomcl.conf --nocompliant
fi

#sudo -S perl5.22.1 $ORTHMCL_PIP/scripts/orthomcl-pipeline.pl -i input/ -o output/ -m $ORTHMCL_PIP/orthomcl.conf --nocompliant

#head -n 1 *.fasta |grep "^>"
# awk 'FNR==NR{a[$1];next}($1 in a){print}' file2 file1

#-----------------------------------------------------------------------------------------------------------------------------------------
# cat orthomcl_run.sh
# # output:
# # orthomcl-pipeline --no-cleanup --compliant -s 10 -m orthomcl.conf -i processed -o orthomcl_out
#
# mkdir orthomcl_out
# # Now we can construct our command and put it into a shell script. We’ll use settings that tell orthomcl-pipeline that our data is compliant with the required formatting (see above), that keeps all intermediate data (in case we want it), that splits our data into 10 chunks for running BLAST (I have 12 cores available, so I set this to 10), and we’ll pass in the orthomcl.conf file we created earlier with the names of the input and output directories.
#
# cat orthomcl_run.sh
# # # output:
# # # orthomcl-pipeline --no-cleanup --compliant -s 10 -m orthomcl.conf -i processed -o orthomcl_out
# # We can now run the pipeline, using screen to manage the terminal and writing a log file.
#
# screen -S orthomcl
# bash ./orthomcl_run.sh 2>&1 | tee orthomcl_run.log
# # That’s it! With these data and my system, this only took about a day to run. Adding more data (i.e., species) will increase the run time.
