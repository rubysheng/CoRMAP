#!/bin/bash
#title          :section2.4_express_matrix.sh
#description    : .
#author         :Ruby(Yiru) Sheng
#usage          :source $CMRP_PATH/script/section2.4_express_matrix.sh -l input_dir_list 2>&1 | tee express_matrix.log
#input files    :a list of target directories with absolute paths which can find 'rsem-gene.gene.TPM.not_cross_norm'
#output files   :renamed expression matrix of all transcript isoforms
#bash_version   :4.4.19(1)-release
#=======================================================================================================================



OPTIND=1         # Reset in case getopts has been used previously in the shell.
input_dir=""
list=""

while getopts "h?:l:" opt; do
    case "$opt" in
    h)
        echo "-l is a required list with input directories (to find the raw TPM expression matrix) "
        echo "-h shows this help message"
        echo "press Ctrl+C to exit"
          # exit 0
        ;;
    l) list=$OPTARG
        ;;
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

function preprocessTMM() {
  TAX_CODE="$1"
  OUT_DIR="$2"
  file=`basename $(pwd)`   # PRJNAXXXXXX
  TMM_NEWNAME=${file}".gene.TMM"
  Complete_name=${TAX_CODE}"_"${file}

  # change "gene TMM expression matrix (rsem-gene.gene.TMM.EXPR.matrix)" format
  if [ ! -e ${TMM_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Reformating rsem-gene.gene.TMM.EXPR.matrix to "${TMM_NEWNAME}
    echo
    find . -name "rsem-gene.gene.TMM.EXPR.matrix" -exec cp -v {} ./tmp.gene_TMM.matrix \;
    sed "s+TRINITY+${Complete_name}+g" tmp.gene_TMM.matrix > ${TMM_NEWNAME}
    rm -v tmp.gene_TMM.matrix
  fi

    echo "Moving to the input directory"
    cp -v ${TMM_NEWNAME} ${OUT_DIR}
}

function preprocessRC() {
  TAX_CODE="$1"
  OUT_DIR="$2"
  file=`basename $(pwd)`   # PRJNAXXXXXX
  RC_NEWNAME=${file}".gene.RC"
  Complete_name=${TAX_CODE}"_"${file}

  # change "gene Raw Counts expression matrix (rsem-gene.gene.counts.matrix)" format
  if [ ! -e ${RC_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Reformating rsem-gene.gene.counts.matrix to "${RC_NEWNAME}
    echo
    find . -name "rsem-gene.gene.counts.matrix" -exec cp -v {} ./tmp.gene_RC.matrix \;
    sed "s+TRINITY+${Complete_name}+g" tmp.gene_RC.matrix > ${RC_NEWNAME}
    rm -v tmp.gene_RC.matrix
  fi

    echo "Moving to the input directory"
    cp -v ${RC_NEWNAME} ${OUT_DIR}
}


if [ ! "${list}" = "" ]; then # have a input list
  input_tpm_dir="$(pwd)/TPM_input/"
  input_tmm_dir="$(pwd)/TMM_input/"
  input_rc_dir="$(pwd)/RC_input/"
  base_dir="$(pwd)"
  checkarray=()

  if [ ! -e ${input_tpm_dir} ]; then
    echo "Creating the input directory to hold all gene expression matrix files (TPM)"
    mkdir -v ${input_tpm_dir}
  fi
  if [ ! -e ${input_tmm_dir} ]; then
    echo "Creating the input directory to hold all gene expression matrix files (TMM)"
    mkdir -v ${input_tmm_dir}
  fi
  if [ ! -e ${input_rc_dir} ]; then
    echo "Creating the input directory to hold all gene expression matrix files (Raw counts)"
    mkdir -v ${input_rc_dir}
  fi


  while IFS=" " read -r f1 f2
  do
    printf 'Taxonomy_code: %s, Dir: %s\n' "$f1" "$f2"
    echo
    SP_CODE="$f1"
    DIR="$f2"
    cd ${DIR}
    preprocessTPM ${SP_CODE} ${input_tpm_dir}
    preprocessTMM ${SP_CODE} ${input_tmm_dir}
    preprocessRC ${SP_CODE} ${input_rc_dir}
    echo
    echo "Only check if the TPM file is well - prepared in "${input_tpm_dir}
    cd ${input_tpm_dir}
    mvd_tpm=`basename ${DIR}`".gene.TPM"
    if [ ! -e ${mvd_tpm} ]; then  # test if the tpm file is moved to the destination
      echo "No file named "${mvd_tpm}" moved to input directory"
    else
      cd ..
      find ${input_tpm_dir} -name ${mvd_tpm} -exec ls -l {} \;
      acc=`basename ${DIR}`
      checkarray+=("${acc}")
    fi
  done < "${list}"


  count_l=`wc -l < ${list}`
fi
