###########################################################################################
## Step 1:                                                                               ##
## Convert the TPM matrix to a long format text file that only includes one column       ##
## (TPM_Value).                                                                          ## 
## The rowname contains with condition and the gene name (the group number still have    ##
## some replicates), seperated by a space.                                               ##
###########################################################################################

files <- c("Bombina","PRJNA252803","PRJNA287145","PRJNA287152","PRJNA302146","PRJNA316996","PRJNA381689","PRJNA390522","PRJNA419677","PRJNA422916","PRJNA450614","PRJNA451011","PRJNA475804","PRJNA529794","PRJNA541005")

i=1

for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/input/",project,".allsp_clugene.TPM")
  # TPM_DF <- read.delim2("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/PRJNA252803.allsp_clugene.TPM",row.names = 1)
  TPM_DF <- read.delim2(file_path, row.names = 1)
  TPM_matrix <- as.matrix(TPM_DF)
  
  long_m <- matrix(TPM_matrix, dimnames=list(t(outer(colnames(TPM_matrix), rownames(TPM_matrix), FUN=paste)), NULL))
  
  out_long <- data.frame(long_m)
  colnames(out_long) <- "TPM_Value"
  output_file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/Newfolder/long_format_",project,"_mulgp.txt")
  write.table(out_long, output_file_path , append = FALSE, sep = "\t", dec = ",", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

