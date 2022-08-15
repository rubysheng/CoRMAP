PRJNA252803_deogg <- read.csv("PRJNA252803_mapping_TMM_DEOGG.csv",row.names = 1)
PRJNA529794_deogg <- read.csv("PRJNA529794_mapping_TMM_DEOGG.csv",row.names = 1)



# PRJNA252803_deogg <- read.csv("PRJNA252803_mapping_DESeq2_DEOGG.csv",row.names = 1)
# PRJNA529794_deogg <- read.csv("PRJNA529794_mapping_DESeq2_DEOGG.csv",row.names = 1)

colnames(PRJNA252803_deogg) <- paste("PRJNA252803",colnames(PRJNA252803_deogg),sep = "_")
colnames(PRJNA529794_deogg) <- paste("PRJNA529794",colnames(PRJNA529794_deogg),sep = "_")

common_deogg <- merge(PRJNA252803_deogg,PRJNA529794_deogg,by=0,all=FALSE)
dim(common_deogg)

# common_deogg <- common_deogg[,c("Row.names","PRJNA252803_log2FoldChange","PRJNA252803_pvalue", "PRJNA252803_padj","PRJNA529794_log2FoldChange","PRJNA529794_pvalue", "PRJNA529794_padj")]
common_deogg <- common_deogg[,c("Row.names","PRJNA252803_logFC","PRJNA252803_PValue", "PRJNA529794_logFC","PRJNA529794_PValue")]
common_deogg$`Gene Name` <- toupper(common_deogg$`Row.names`)

# common_deogg
output_path= "./"
output_file_path <- paste0(output_path,"common_mapping_TMM_DEOGG.csv")
write.csv(common_deogg,row.names = FALSE, output_file_path, quote = FALSE)

cormap_common_deogg <- read.csv("cormap_T3.csv",row.names = 1)

rownames(common_deogg) <- common_deogg$`Gene Name`
merge(common_deogg[,-1],cormap_common_deogg,by=0,all=FALSE)
