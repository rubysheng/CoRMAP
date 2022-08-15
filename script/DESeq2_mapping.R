setwd("/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/mapping/output/")
# in the "/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/mapping/output/" 


name_list_file = "/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/output/R_analysis/input_dataset_name.txt"
output_path = "./"
# project = "PRJNA529794"
project = "PRJNA252803"

pvalue = 0.05
adjmethod = "none"

# DEG
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ condition)
# dds <- DESeqDataSetFromMatrix(countData = mm252803_cts,
#                               colData = mm252803_meta_tb,
#                               design = ~ treatment)
mm529794_meta_tb$treatment_cleaned <- mm529794_meta_tb$treatment
mm529794_meta_tb$treatment_cleaned <- gsub("HomeCage","control",mm529794_meta_tb$treatment_cleaned)
mm529794_meta_tb$treatment_cleaned <- gsub("Non_Shock","control",mm529794_meta_tb$treatment_cleaned)
mm529794_meta_tb$treatment_cleaned <- factor(mm529794_meta_tb$treatment_cleaned)
dds <- DESeqDataSetFromMatrix(countData = mm529794_cts,
                              colData = mm529794_meta_tb,
                              design = ~ treatment_cleaned)
dds <- DESeqDataSetFromMatrix(countData = mm252803_cts,
                              colData = mm252803_meta_tb,
                              design = ~ treatment)




# dds

# not use LRT test (for time course experiments) ; but default Wald test
dds <- DESeq(dds)

# FDR(alpha=0.5), log2 fold change =1 ,"greaterAbs"
res50 <- results(dds,alpha=0.5,pAdjustMethod = "fdr",lfcThreshold=1,altHypothesis="greaterAbs")
res50Ordered <- res50[order(res50$pvalue),]
summary(res50)
sum(res50$padj < 0.5, na.rm=TRUE)

# default FDR(alpha=0.1), log2 fold change =1 ,"greaterAbs"
res10 <- results(dds,alpha=0.1,pAdjustMethod = "fdr",lfcThreshold=1,altHypothesis="greaterAbs")
res10Ordered <- res10[order(res10$pvalue),]
summary(res10)
sum(res10$padj < 0.1, na.rm=TRUE)

# FDR(alpha=0.05), log2 fold change =1 ,"greaterAbs"
res5 <- results(dds,alpha=0.05,pAdjustMethod = "fdr",lfcThreshold=1,altHypothesis="greaterAbs")
res5Ordered <- res5[order(res5$pvalue),]
summary(res5)
sum(res5$padj < 0.05, na.rm=TRUE)

# FDR(alpha=0.01), log2 fold change =1 ,"greaterAbs"
res1 <- results(dds,alpha=0.01,pAdjustMethod = "fdr",lfcThreshold=1,altHypothesis="greaterAbs")
res1Ordered <- res1[order(res1$pvalue),]
summary(res1)
sum(res1$padj < 0.01, na.rm=TRUE)

# FDR(alpha=0.001), log2 fold change =1 ,"greaterAbs"
res.1 <- results(dds,alpha=0.001,pAdjustMethod = "fdr",lfcThreshold=1,altHypothesis="greaterAbs")
res.1Ordered <- res.1[order(res.1$pvalue),]
summary(res.1)
sum(res.1$padj < 0.001, na.rm=TRUE)


# "none", pvalue>0.01
res1_none <- results(dds,pAdjustMethod = "none")
res1_noneOrdered <- res1_none[order(res1_none$pvalue),]
summary(res1_none)
sum(res1_none$pvalue < 0.01, na.rm=TRUE)

# "none", pvalue>0.05
res5_none <- results(dds,pAdjustMethod = "none")
res5_noneOrdered <- res5_none[order(res5_none$pvalue),]
summary(res5_none)
sum(res5_none$pvalue < 0.05, na.rm=TRUE)

# "none", pvalue>0.1
res10_none <- results(dds,pAdjustMethod = "none")
res10_noneOrdered <- res10_none[order(res10_none$pvalue),]
summary(res10_none)
sum(res10_none$pvalue < 0.1, na.rm=TRUE)


# ---------
# deg.1 <- res.1Ordered[1:sum(res.1$padj < 0.001, na.rm=TRUE),]
# deg1 <- res1Ordered[1:sum(res1$padj < 0.001, na.rm=TRUE),]
# deg5 <- res5Ordered[1:sum(res5$padj < 0.001, na.rm=TRUE),]
# deg10 <- res10Ordered[1:sum(res10$padj < 0.001, na.rm=TRUE),]
# deg50 <- res50Ordered[1:sum(res50$padj < 0.001, na.rm=TRUE),]
# 
# deg.1 <- res.1Ordered[1:sum(res.1$padj < 0.001, na.rm=TRUE),]
# deg1 <- res1Ordered[1:sum(res1$padj < 0.001, na.rm=TRUE),]
# deg5 <- res5Ordered[1:sum(res5$padj < 0.001, na.rm=TRUE),]
# deg10 <- res10Ordered[1:sum(res10$padj < 0.001, na.rm=TRUE),]
# deg50 <- res50Ordered[1:sum(res50$padj < 0.001, na.rm=TRUE),]

deg5_none <- res5_noneOrdered[1:sum(res5_none$pvalue < 0.05, na.rm=TRUE),]


# save files 
output_file_path <- paste0(output_path,project,"_mapping_DESeq2_DEOGG.RData")
save(deg5_none,file = output_file_path)
output_file_path <- paste0(output_path,project,"_mapping_DESeq2_DEOGG.csv")
write.csv(deg5_none,row.names = TRUE, output_file_path, quote = FALSE)




# DEG analysis edgeR -----------------------------------------------------------
group_rdata_path <- "/Users/ruby/Documents/Master/GitHub/CoRMAP/sample/output/R_analysis/"
# project <- "PRJNA252803"
project <- "PRJNA529794"
type <- "TMM"
adj.method <- "none"
pv <- 0.05

# mm252803_cts
# mm252803_meta_tb
# mm529794_cts
# mm529794_meta_tb

# countMatrix <- mm252803_cts
countMatrix <- mm529794_cts
# group <- factor(mm252803_meta_tb$treatment)
group <- factor(mm529794_meta_tb$treatment_cleaned)
# subGroup <- factor(rep(1:4,2))
subGroup <- factor(c(1:sum(mm529794_meta_tb$treatment_cleaned=="Fear_Conditioned"),1:sum(mm529794_meta_tb$treatment_cleaned=="control")))


y <- DGEList(counts = countMatrix, group = group)

keep <- rowSums(cpm(y)>1) >= 2
y_filter0 <- y[keep, , keep.lib.sizes=FALSE]

y_normed <- calcNormFactors(y_filter0 ,method = 'TMM')

# subGroup <- factor(substring(colnames(countMatrix), regexpr("_rep", colnames(countMatrix))+4))
design <- model.matrix(~ subGroup+group)
rownames(design) <- colnames(y_normed)

dge <- estimateDisp(y_normed, design, robust=TRUE)

#negative binomial generalized log-linear model 
fit <- glmFit(dge, design, robust = TRUE)     
lrt <- glmLRT(fit)   

resNoFilt <- topTags(lrt,n=nrow(lrt$table),adjust.method = adj.method,p.value = pv)
sigDownReg <- resNoFilt$table[resNoFilt$table$PValue<pv,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]

# save files 
output_file_path <- paste0(output_path,project,"_mapping_TMM_DEOGG.RData")
save(sigDownReg,file = output_file_path)
output_file_path <- paste0(output_path,project,"_mapping_TMM_DEOGG.csv")
write.csv(sigDownReg,row.names = TRUE, output_file_path, quote = FALSE)




# load(paste0(group_rdata_path,project,"_avg",type,"_DEOGG.RData"))
# load(paste0(group_rdata_path,project,"_avg",type,"_matrix_bygp.RData"))
# rm(sigDownReg)
# countMatrix
# group
# str(countMatrix)
# str(group)
# str(mm252803_cts)
# str(mm252803_meta_tb)
