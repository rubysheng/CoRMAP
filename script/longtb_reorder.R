###########################################################################################
## Step 2:                                                                               ##
## 1) Convert a one column (TPM_Value) text file to a long format excel table            ##
## 2) Remove the repelicates of group number in gene name                                ##
## 3) Add more columns to split details for further analysis                             ##
###########################################################################################


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager","argparse","readxl","writexl")
BiocManager::install()

suppressPackageStartupMessages(library("argparse"))
parser = ArgumentParser()
parser$add_argument("--name_list", help="a file containing a list of dataset names", required=TRUE)
parser$add_argument("--input", help="the folder holding input long formated tables (long_format_*_mulgp.txt)", required=TRUE)
parser$add_argument("--output", help="path to save the output xlsx files", required=TRUE)
args = parser$parse_args()

name_list_file = args$name_list
input_path = args$input
output_path = args$output


library(readxl)
library(writexl)


files <- scan(name_list_file, character(), quote = "")


for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0(input_path,"long_format_",project,"_mulgp.txt")
  long_DF <- read.delim2(file_path, row.names=NULL)
  colnames(long_DF)[1] <- c("x")
  sp1_df <- data.frame(do.call('rbind', strsplit(as.character(long_DF$x),' ',fixed=TRUE)))
  sp2_df <- data.frame(do.call('rbind', strsplit(as.character(sp1_df$X2),'_',fixed=TRUE)))
  sp2_df[,ncol(sp2_df)+1] <- gsub('g','',sp2_df[,1])

  sp3_df <- data.frame(do.call('rbind', strsplit(as.character(sp2_df$V7),'/',fixed=TRUE)))
  sp2_df$V7 <- apply(sp3_df, 1,function(x) paste(unique(as.character.numeric_version(x)), collapse = "/"))



  # uniq_gn <- function(sp3_df,row1) {
  #   paste(unique(as.character.numeric_version(sp3_df)))
  #   # r_uniq_gn <- unique(as.character.numeric_version(sp3_df[row1]))
  #   # paste(r_uniq_gn, collapse = "/")
  # }
  # test <- sp3_df[c("1","2","526","1345","1346","1347","1353"),]
  # paste(unique(as.character.numeric_version(test["1345",])), collapse = "/")
  # test$uniq_gp <- apply(test, 1,function(x) paste(unique(as.character.numeric_version(x)), collapse = "/"))
  #
  # test$uniq_gp
  # r_uniq_gn <- unique(as.character.numeric_version(sp3_df["1345",]))
  # paste(r_uniq_gn, collapse = "/")
  # r_uniq_gn <- unique(as.character.numeric_version(sp3_df["1",]))
  # paste(r_uniq_gn, collapse = "/")
  #
  #
  # head(unique(sp3_df))


  sp2_df[,1] <- sp2_df[,ncol(sp2_df)]

  sp2_df$gene <- apply(sp2_df,1, function(x) paste0("g",x[1],"_",x[2],"_",x[3],"_",x[4],"_",x[5],"_",x[6]))

  colnames(sp1_df) <- c("treatment","gene")
  sp1_df <- sp1_df[ , c("gene", "treatment")]
  colnames(sp2_df)[1:3] <- c("Group_Num","Species","Project")
  out_long <- cbind(sp2_df[,c(1,2,3,8)],sp1_df[,2],long_DF[,2])
  colnames(out_long) <- c("Group_Num","Species","Project","gene","treatment","TPM_Value")
  output_file_path <- paste0(output_path,"long_format_",project,"_mulgp.xlsx")
  write_xlsx(out_long,output_file_path)
}
