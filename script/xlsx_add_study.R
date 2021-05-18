if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager","argparse","readxl","writexl")
BiocManager::install()

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()
parser$add_argument("--name_list", help="a file containing a list of dataset names", required=TRUE)
parser$add_argument("--input", help="the folder holding input long formated xlsx files (long_format_*_mulgp.xlsx)", required=TRUE)
parser$add_argument("--output", help="path to save the modified output xlsx files (*_mulgp_clugene_sty.xlsx)", required=TRUE)
parser$add_argument("--studytype", help="a xlsx file showing the type of learning category", required=TRUE)
args = parser$parse_args()

name_list_file = args$name_list
input_path = args$input
output_path = args$output
studytype_file = args$studytype

library(readxl)
library(writexl)


files <- scan(name_list_file, character(), quote = "")
studytype <- read_xlsx(studytype_file)

for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0(input_path,"long_format_",project,"_mulgp.xlsx")
  long_df <- read_xlsx(file_path)
  pro_sty <- studytype[which(studytype$`BioProject accession`==project),]
  for ( col in 1:ncol(pro_sty)) {
    # print(col)
    # print(as.character(pro_sty[1,col]))
    # print(sapply(pro_sty,class)[col])
    if (sapply(pro_sty,class)[col]=="numeric") {
      pro_sty[2:nrow(long_df),col] <- rep(as.numeric(pro_sty[1,col]), nrow(long_df)-1)
    } else {
      pro_sty[2:nrow(long_df),col] <- rep(as.character(pro_sty[1,col]), nrow(long_df)-1)
    }
  }

  output_tb <- cbind(long_df,pro_sty)
  output_tb$Taxa_group <- substr(output_tb$Species,1,2)
  # View(output_tb)
  # write.table(output_tb, "~/Documents/Documents/MSc/MSc_thesis/RNA-seq/expression/PRJNA252803_clugene_sty.txt", append = FALSE, sep = "\t", dec = ".", row.names = TRUE, col.names = TRUE, quote = FALSE)

  output_file_path <- paste0(output_path,project,"_mulgp_clugene_sty.xlsx")
  write_xlsx(output_tb,output_file_path)
}
