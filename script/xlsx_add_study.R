# install.packages("readxl")
# install.packages("writexl")
library(readxl)
library(writexl)
setwd("~/Documents/Documents/MSc/MSc_thesis/RNA-seq/expression/")


files<-c("Bombina","PRJNA252803","PRJNA287145","PRJNA287152","PRJNA302146","PRJNA316996","PRJNA381689","PRJNA390522","PRJNA419677","PRJNA422916","PRJNA450614","PRJNA451011","PRJNA475804","PRJNA529794","PRJNA541005")

studytype <- read_xlsx("~/Documents/Documents/MSc/MSc_thesis/RNA-seq/expression/17_studytype.xlsx")

for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0("~/Documents/Documents/MSc/MSc_thesis/RNA-seq/expression/output/Newfolder/long_format_",project,"_mulgp.xlsx")
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

  output_file_path <- paste0("~/Documents/Documents/MSc/MSc_thesis/RNA-seq/expression/output/",project,"_mulgp_clugene_sty.xlsx")
  write_xlsx(output_tb,output_file_path)
}

