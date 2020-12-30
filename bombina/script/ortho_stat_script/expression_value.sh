#!/usr/bin/bash

# input:
#   a list of target directories with absolute paths which can find *fasta.RSEM.transcripts.fa
#   example:
#       RWcel /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/PRJNA451011
#         *delimited by a space

# direction:
#   In a Taxonomy_directory with a list table as shown above,
#   create a output folder to hold the pep fasta quant_files
#       mkdir output
#   then run the command:
#       source /media/lewis/New_Seagate_Drive_8TB/ruby/github/bombina/comparative-transcriptomic-analysis-pip/script/expression_value.sh -l lst -o output 2>&1 | tee changename.log


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
    o) output_dir=$OPTARG
        if [[ "${output_dir: -1}" != '/' ]]; then
            output_dir=${output_dir}"/"
        fi
        #mkdir -p ${output_dir}
    esac
done
shift $(( OPTIND-1 ))

function preprocessTPM() {
  TAX_CODE="$1"
  OUT_DIR="$2"
  file=`basename $(pwd)`   # PRJNAXXXXXX
  TPM_NEWNAME=${file}".gene.TPM"
  Complete_name=${TAX_CODE}"_"${file}

  # change fasta file and gene_trans_map file format
  if [ ! -e ${TPM_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Reformating rsem-gene.gene.TPM.not_cross_norm to "${TPM_NEWNAME}
    echo
    sed "s+TRINITY+${Complete_name}+g" transcripts_count/rsem-gene.gene.TPM.not_cross_norm > ${TPM_NEWNAME}
  fi

    echo "Moving to the input directory"
    cp -v ${TPM_NEWNAME} ${OUT_DIR}
}

if [ ! "${list}" = "" ]; then # have a input list
  input_dir="$(pwd)/TPM_input/"
  base_dir="$(pwd)"
  checkarray=()

  if [ ! -e ${input_dir} ]; then
    echo "Creating the input directory to hold all protein fasta files"
    mkdir -v ${input_dir}
  fi

  while IFS=" " read -r f1 f2
  do
    printf 'Taxonomy_code: %s, Dir: %s\n' "$f1" "$f2"
    echo
    SP_CODE="$f1"
    DIR="$f2"
    cd ${DIR}
    preprocessTPM ${SP_CODE} ${input_dir}
    echo
    echo "Check if the file is well - prepared in "${input_dir}
    cd ${input_dir}
    mvd_tpm=`basename ${DIR}`".gene.TPM"
    if [ ! -e ${mvd_tpm} ]; then  # test if the tpm file is moved to the destination
      echo "No file named "${mvd_tpm}" moved to input directory"
    else
      cd ..
      find ${input_dir} -name ${mvd_tpm} -exec ls -l {} \;
      acc=`basename ${DIR}`
      checkarray+=("${acc}")
    fi
  done < "${list}"
  count_l=`wc -l < ${list}`
fi


# ################################################################################
# ### backup to cluster computer
# ls -1 /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ > pro_acc
# mkdir upload
# while IFS=" " read -r line
# do
#  File="/media/lewis/Seagate_Backup_Plus_Drive/ruby/data/"${line}"/transcripts_count/rsem-gene.gene.TPM.not_cross_norm"
#  NewName=${line}"_rsem-gene.gene.TPM.not_cross_norm"
#  cp -v ${File} upload/${NewName}  
# done < pro_acc
# cp -v /media/lewis/New_Seagate_Drive_8TB/Bombina/transcripts_count/rsem-gene.gene.TPM.not_cross_norm upload/Bombina_rsem-gene.gene.TPM.not_cross_norm
# scp upload/*  ruby@graham.computecanada.ca:/home/ruby/projects/def-heylanda/ruby/all_taxa/expression/
# rm -r upload/
# rm pro_acc
# ################################################################################





