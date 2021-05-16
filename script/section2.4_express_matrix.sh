#!/bin/bash
#title          :section2.4_express_matrix.sh
#description    : .
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.4_express_matrix.sh -l input_dir_list -o output_dir 2>&1 | tee express_matrix.log
#input files    :a list of target directories with absolute paths which can find *.gene.TPM.not_cross_norm
#output files   :renamed expression matrix of clustered orthologous genes
#bash_version   :4.4.19(1)-release
#=======================================================================================================================



OPTIND=1         # Reset in case getopts has been used previously in the shell.
input_dir=""
output_dir=""
list=""

while getopts "h?:l:o:" opt; do
    case "$opt" in
    h)
        echo "-l is a required list with input directories (to find the raw TPM expression matrix) "
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

  # change "gene TPM expression matrix (rsem-gene.gene.TPM.not_cross_norm)" format
  if [ ! -e ${TPM_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Reformating rsem-gene.gene.TPM.not_cross_norm to "${TPM_NEWNAME}
    echo
    find . -name "rsem-gene.gene.TPM.not_cross_norm" -exec cp -v {} ./tmp.gene_TPM.matrix \;
    sed "s+TRINITY+${Complete_name}+g" tmp.gene_TPM.matrix > ${TPM_NEWNAME}
    rm -v tmp.gene_TPM.matrix
  fi

    echo "Moving to the input directory"
    cp -v ${TPM_NEWNAME} ${OUT_DIR}
}

if [ ! "${list}" = "" ]; then # have a input list
  input_dir="$(pwd)/TPM_input/"
  base_dir="$(pwd)"
  checkarray=()

  if [ ! -e ${input_dir} ]; then
    echo "Creating the input directory to hold all gene expression matrix files"
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
