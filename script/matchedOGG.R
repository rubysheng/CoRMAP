name_list_file="input_dataset_name.txt"
type="TMM"
output_path="./"
# annotation -------------------------------------------------------------------
OGG2GeneName <- read.delim("./Unique_OGG_GENENAME.txt", header = T, sep = "\t")
OGG2GeneName <- unique(OGG2GeneName)
OGG2GeneName.lst <- OGG2GeneName %>%
  split(f = as.factor(.$Group)) %>%
  lapply(., "[", , c("Gene_name"))
all_genenames.tb <- data.frame(sapply(OGG2GeneName.lst, FUN = function(values) paste(values,collapse = "|")))
colnames(all_genenames.tb) <- "gene_names"

# load -------------------------------------------------------------------------
files <- scan(name_list_file, character(), quote = "")
for (i in 1:length(files)){
  project <- files[i]
  load(paste0(project,'_avgTMM_DEOGG.RData'))
  assign(project,sigDownReg)
  assign(paste0(project,"_group"),group)
  rm(sigDownReg);rm(group);rm(project)
}

matched.df <- merge(PRJNA252803[-6],PRJNA529794[-6], by = "Group")
rownames(matched.df) <- matched.df$Group
matched.df <- merge(all_genenames.tb,matched.df, by = "row.names")
matched.df <- matched.df[-1]

output_file_path <- paste0(output_path,"matchedOGG_",type,"_",files[1],"_vs_",files[2],".csv")
write.csv(matched.df,row.names = FALSE, output_file_path, quote = FALSE)
