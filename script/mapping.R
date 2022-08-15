# getwd(): "/Users/ruby/Documents/Master/GitHub/CoRMAP"
setwd("/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/mapping/output/")
# in the "/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/mapping/output/" 


name_list_file = "/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/output/R_analysis/input_dataset_name.txt"
output_path = "./"


# create blank summary table for all runs from all datasets --------------------
mp_summary_tables <- data.frame(row.names =c("Dataset","Status","Assigned","Unassigned_Unmapped","Unassigned_Read_Type","Unassigned_Singleton","Unassigned_MappingQuality","Unassigned_Chimera","Unassigned_FragmentLength","Unassigned_Duplicate","Unassigned_MultiMapping","Unassigned_Secondary","Unassigned_NonSplit","Unassigned_NoFeatures","Unassigned_Overlapping_Length","Unassigned_Ambiguity") ) 

# summarize all counts.summary reports from the selected datasets --------------
files <- scan(name_list_file, character(), quote = "")

for (i in 1:length(files)){
  project <- files[i]
  temp <- list.files(path = project,pattern="*.summary$",full.names =TRUE)
  summary_files <- do.call(cbind,lapply(temp, FUN = function(x) read.table(file = x,header = FALSE,row.names = 1)))
  summary_files <- rbind(project,summary_files)
  colnames(summary_files) <- sub("_refseq.sam","",summary_files[2,])
  rownames(summary_files)[1] <- "Dataset"
  mp_summary_tables <- cbind(mp_summary_tables,summary_files)
}
# save files
save(mp_summary_tables, file = paste0(output_path,"counts_summaries.RData" ))
write.csv(mp_summary_tables,paste0(output_path,"counts_summaries.csv"), row.names = TRUE)

# generate the count matrix for each selected dataset --------------------------

for (i in 1:length(files)){
  project <- files[i]
  temp <- list.files(path = project,pattern="*.txt$",full.names =TRUE)
  files_list <- lapply(temp, FUN = function(x) read.table(file = x,skip = 1,header = 2)[,c(1,7)])
  
  count_matrix <- Reduce(function(x,y) merge(x,y,all=TRUE),files_list)
  colnames(count_matrix)[-1] <- sub("_refseq.sam","",colnames(count_matrix)[-1])
  rownames(count_matrix) <- count_matrix[,1]
  count_matrix <- count_matrix[,-1]
  
  save(count_matrix, file = paste0(output_path,project,"/count_matrix_all.RData"))
  write.csv(count_matrix,paste0(output_path,project,"/count_matrix_all.csv"), row.names = TRUE)
  
  
  count_matrix_cleaned <- count_matrix[rowSums(count_matrix[])>0,]
  save(count_matrix_cleaned, file = paste0(output_path,project,"/count_matrix_cleaned.RData"))
  write.csv(count_matrix_cleaned,paste0(output_path,project,"/count_matrix_cleaned.csv"), row.names = TRUE)
  
}

# load meta table for design groups --------------------------------------------


