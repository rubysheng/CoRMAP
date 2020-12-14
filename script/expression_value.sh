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
        mkdir -p ${output_dir}
    esac
done
shift $(( OPTIND-1 ))

function preprocessTPM() {
  TAX_CODE="$1"
  OUT_DIR="$2"
  file=`basename $(pwd)`   # PRJNAXXXXXX
  TPM_NEWNAME=${file}".gene.TPM"


  # change fasta file and gene_trans_map file format
  if [ ! -e ${TPM_NEWNAME} ]; then
    echo "No renamed file"
    echo "    Reformating rsem-gene.gene.TPM.not_cross_norm to "${TPM_NEWNAME}
    echo
    sed "s+TRINITY+${file}+g" transcripts_count/rsem-gene.gene.TPM.not_cross_norm > ${TPM_NEWNAME}
    sed -i "s/_/./g" ${TPM_NEWNAME}
    echo "Adding the taxon code"
    sed -i "s/^/${TAX_CODE}_/" ${TPM_NEWNAME}
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
    fi
  done < "${list}"
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

sed 's/worm|//; s/[.]/_/g' ../clustered_worm.lst > worm_clu.lst
sed 's/mam|//; s/[.]/_/g' ../clustered_mam.lst > mam_clu.lst
sed 's/bom|//; s/[.]/_/g' ../clustered_bom.lst > bom_clu.lst
sed 's/ins|//; s/[.]/_/g' ../clustered_ins.lst > ins_clu.lst
sed 's/fis|//; s/[.]/_/g' ../clustered_fis.lst > fis_clu.lst
# extract the Clustered genes' TPM
ls -1 /media/lewis/Seagate_Backup_Plus_Drive/ruby/data/ > pro_acc
# echo 'Bombina' >> pro_acc
mkdir TPM_output
# echo 'PRJNA422916' > re_acc
# echo 'PRJNA381689' >> re_acc
# echo 'PRJNA316996' >> re_acc
while IFS= read -r project; do
  genelist=${project}"_clugene"
  grep "$project" *.lst --no-filename | sed 's/_i.*//' | uniq > ${genelist}
  TPM_File="./TPM_input/"${project}".gene.TPM"
  Filtered_File="./TPM_output/"${project}".clugene.TPM"
  head -1 ${TPM_File} > ${Filtered_File}
  while IFS= read -r gene; do
    pro_acc=${project}
    echo ${pro_acc}
    echo ${gene}
    awk -v pat="$gene" '{ if ($0 ~ pat) {print} }' ${TPM_File} >> ${Filtered_File}
  done < ${genelist}
done <  re_acc # bom_acc ##pro_acc

# add group number to each gene
while IFS= read -r project; do
  grep "$project" ../per_group/allsp_group_*.dnaseq.fasta | sed 's/[.]i.*//; s/[.][.][/]per_group[/]allsp_group_//; s/[.]dnaseq[.]fasta:>/ /; s/[.]/_/g' | sort | uniq > gp_gene_map
  input="./TPM_output/"${project}".clugene.TPM"
  cp ${input} tmp
  while IFS=" " read -r f1 f2
  do
    group_num="$f1"
    gene_name="$f2"
    c="g"${group_num}"_"${gene_name}
    sed -i "s/${gene_name}/${c}/" tmp
  done < gp_gene_map
  rm gp_gene_map
  head -1 tmp > finalfile
  grep '^g' tmp >> finalfile
  rm tmp
  Filtered_File="./TPM_output/"${project}".allsp_clugene.TPM"
  mv finalfile ${Filtered_File}
done <  re_acc ##  pro_acc ##acc



# remove some duplicates entries and combine some multiple-group entries' group code
ls -1 *.allsp_clugene.TPMbak > download_file_lst
while IFS= read -r File; do
  output_file="../"${File%bak}
  echo ${output_file}
  sort ${File} | uniq > ${output_file}
  # tail -n+2 ${File} | awk '{print $1}' | cut -d'_' -f 2 | sort | uniq
  # echo
  # tail -n+2 ${File} | awk '{print $1}' | cut -d'_' -f 3 | sort | uniq -d
done < download_file_lst


head -1 PRJNA287145.allsp_clugene.TPM  > tmp
tail -n+2 PRJNA287145.allsp_clugene.TPM | sort | uniq >> tmp
mv tmp PRJNA287145.allsp_clugene.TPM
# echo PRJNA287145.allsp_clugene.TPM
# tail -n+2 PRJNA287145.allsp_clugene.TPM | awk '{print $1}' | cut -d'_' -f 2 | sort | uniq
grep "^g.*_g.*_PRJNA......_DN.*_c.*_g.*" PRJNA287145.allsp_clugene.TPM
awk -F="_PRJNA" '{gsub(/_g/,"/",$1)}1' OFS== PRJNA287145.allsp_clugene.TPM
sed -e ':b; s/^\([^_PRJNA]*\)*\_g/\1[/]/; tb;' PRJNA287145.allsp_clugene.TPM







tail -n+2 PRJNA252803.allsp_clugene.TPM | awk '{print $1}' | cut -d'_' -f 2 | sort | uniq
# RDmus
# g59
# g711
grep -E "_g59_|_g711_" PRJNA252803.allsp_clugene.TPM
sed -ibak 's+g2180_g59_+g2180/59_+; s+g251_g711_+g251/711_+' PRJNA252803.allsp_clugene.TPM

# after the R script (convert_to_longtb.R)
sed -ibak "s/ /\t/" long_format.txt
