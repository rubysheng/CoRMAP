# install.packages("readxl")
# install.packages("writexl")
library(readxl)
library(writexl)
library(stringr)
setwd("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/")

files<-c("Bombina","PRJNA252803","PRJNA287145","PRJNA287152","PRJNA302146","PRJNA316996","PRJNA381689","PRJNA390522","PRJNA419677","PRJNA422916","PRJNA450614","PRJNA451011","PRJNA475804","PRJNA529794","PRJNA541005")


duplicate_rows <- function(Group_Num,Species,Project,gene,treatment,TPM_Value,multi) {
  expanded_samples <- paste0("rep", 1:multi)
  Group_Nums <- unlist(str_split(Group_Num,"/"))
  repeated_rows <- data.frame("Group_Num"=Group_Nums,"Species"=Species,"Project"=Project,"gene"=gene,"treatment"=treatment,"TPM_Value"=TPM_Value,"multi"=expanded_samples)
  repeated_rows
}

rename_gene <- function(new_df2,col1,col2) {
  gene <- as.character(new_df2[col2])
  loc <- str_locate(gene,"_")[1,1]
  paste0("g",new_df2[col1],str_sub(gene,loc,nchar(gene)))
}

count_pt <- function(mul,col) { str_count(mul[col],"/")+1 }

for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/Newfolder/long_format_",project,"_mulgp.xlsx")
  long_df <- read_xlsx(file_path)
  slash_loc <- which(grepl(pattern="/", x=long_df$Group_Num))
  single <- data.frame(long_df[-slash_loc,])
  rownames(single) <- 1:nrow(single)
  mul <- data.frame(long_df[grepl(pattern="/", x=long_df$Group_Num),])
  rownames(mul) <- 1:nrow(mul)
  mul[,"multi"] <- apply(mul,1,count_pt, col="Group_Num")
  # mul
  expanded_rows <- Map(f = duplicate_rows,mul$Group_Num,mul$Species,mul$Project,mul$gene,mul$treatment,mul$TPM_Value,mul$multi)
  new_mul <- do.call(rbind, expanded_rows)
  new_mul$gene <- apply(new_mul, 1, rename_gene, col1="Group_Num",col2="gene")
  rownames(new_mul) <- 1:nrow(new_mul)
  mul <- new_mul[1:6]
  out_long <- rbind(mul,single)
  output_file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/Newfolder/long_format_",project,"_mulgp.xlsx")
  write_xlsx(out_long,output_file_path)
}
