setwd("/Users/ruby/Documents/GitHub/comparative-transcriptomic-analysis-pip/sample/output/R_analysis/")
cv_dataset <- "PRJNA529794"
output_path <- "./"
# load("./PRJNA252803_countMX.RData")
load("./PRJNA529794_countMX.RData")
adj.method <- 'fdr'
# library packages --------------------------------------------------------


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("ensembldb")
# library(DESeq2)
library(stringr)
library(tidyverse)
library(biomaRt)




# deg analysis ------------------------------------------------------------

library(edgeR)

y <- DGEList(counts = countMatrix, group = group)

keep <- rowSums(cpm(y)>1) >= 2
y_filter0 <- y[keep, , keep.lib.sizes=FALSE]

y_normed <- calcNormFactors(y_filter0 ,method = 'TMM')

subGroup <- factor(substring(colnames(countMatrix), regexpr("_rep", colnames(countMatrix))+4))
design <- model.matrix(~ subGroup+group)
rownames(design) <- colnames(y_normed)

dge <- estimateDisp(y_normed, design, robust=TRUE)

#negative binomial generalized log-linear model 
fit <- glmFit(dge, design, robust = TRUE)     
lrt <- glmLRT(fit)   


resNoFilt <- topTags(lrt,n=nrow(lrt$table))
sum(resNoFilt$table$FDR < 0.05)
sigDownReg <- resNoFilt$table[resNoFilt$table$FDR<0.05,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
# head(sigDownReg)

dge_de <- decideTestsDGE(lrt, adjust.method = adj.method, p.value = 0.05)  
isDE <- as.logical(dge_de)
DEnames <- rownames(dge)[isDE]



# biomartR ----------------------------------------------------------------

mart104.mm <- useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl")
# mart75.mm <- useMart("ENSEMBL_MART_ENSEMBL", "mmusculus_gene_ensembl", host = "feb2014.archive.ensembl.org")

# attributes <- listAttributes(mart75.mm)
# View(attributes)
# 
# filters <- listFilters(mart75.mm)
# View(filters)

attributes <- listAttributes(mart104.mm)
# View(attributes)

filters <- listFilters(mart104.mm)
# View(filters)


# load "OGG-GeneName" list  -----------------------------------------------

OGG2GeneName <- read.delim("./Unique_OGG_GENENAME.txt", header = T, sep = "\t")
OGG2GeneName <- unique(OGG2GeneName)
OGG2GeneName.lst <- OGG2GeneName %>%
  split(f = as.factor(.$Group)) %>%
  lapply(., "[", , c("Gene_name"))


# load the DE gene list ----------------------------------------------------------

# annotation_desp <- read.delim("~/Desktop/publication/Full_anno_by_group.txt", header = F, sep = "\t")
# annotation_desp$V1 <- str_replace(annotation_desp$V1," ","")
# annotation_desp$V2 <- sub("\\_.*", "", annotation_desp$V2)
# annotation_desp <- unique(annotation_desp)
# 
# annotation_desp_ls <- annotation_desp %>%
#   split(f = as.factor(.$V1))
# all_genenames <- lapply(annotation_desp_ls, "[", , c("V2"))

DEgp_names <- paste0("group_",DEnames)
DEgp_genenames <- OGG2GeneName.lst[which(names(OGG2GeneName.lst)%in%DEgp_names)]

# DEgp_genenames <- all_genenames[which(names(all_genenames)%in%DEgp_names)]


# degs.mm104 <- getBM(attributes = c("uniprot_gn_symbol", "ensembl_gene_id"), filters = "uniprot_gn_symbol", values = unique(unlist(DEgp_genenames)), bmHeader = T, mart = mart104.mm)
# colnames(degs.mm104) <- c("uniprot_gn_symbol","ensembl_gene_id")
# degs.mm104$ensembl_gene_id <- as.character(degs.mm104$ensembl_gene_id)

# deg <- scan("sample/output/R_analysis/ruby_deg.txt", character(), quote = "")
# deg <- deg[which(deg!="#N/A")]



# load the original gene list ---------------------------------------------

original_deg_path <- paste0("../../data/",cv_dataset,"/original_deg.txt")
ori <- scan(original_deg_path, character(), quote = "")

# ori.mm75 <- getBM(attributes = c("uniprot_genename", "ensembl_gene_id"), filters = "uniprot_genename", values = ori, bmHeader = T, mart = mart75.mm)
# colnames(ori.mm75) <- c("uniprot_genename","ensembl_gene_id")
# ori.mm75$ensembl_gene_id <- as.character(ori.mm75$ensembl_gene_id)

# ori.mm104 <- getBM(attributes = c("uniprot_gn_symbol", "ensembl_gene_id"), filters = "uniprot_gn_symbol", values = ori, bmHeader = T, mart = mart104.mm)
# colnames(ori.mm104) <- c("uniprot_gn_symbol","ensembl_gene_id")
# ori.mm104$ensembl_gene_id <- as.character(ori.mm104$ensembl_gene_id)

ori.mm104 <- getBM(attributes = c("uniprot_gn_symbol", "external_synonym"), filters = "uniprot_gn_symbol", values = ori, bmHeader = T, mart = mart104.mm)
colnames(ori.mm104) <- c("uniprot_gn_symbol","external_synonym")

oris.mm104 <- c(ori,ori.mm104$uniprot_gn_symbol,ori.mm104$external_synonym)
oris.mm104 <- unique(oris.mm104)

# length(unique(ori.mm75$uniprot_genename))
# length(unique(ori.mm75$uniprot_gn_symbol))
# # ori.mm75
# deg.mm104
# length(deg)
# unique(ori.mm75$uniprot_gn_symbol)


# compare two lists -------------------------------------------------------

sort(toupper(ori))


toupper(degs.mm104$uniprot_gn_symbol)

degs.mm104[which(degs.mm104$ensembl_gene_id %in% ori.mm104$ensembl_gene_id),]

dim(degs.mm104);dim(ori.mm104)

getBM(attributes = c("uniprot_gn_symbol", "external_synonym"), filters = "uniprot_gn_symbol", values = "POP4", bmHeader = T, mart = mart104.mm)
getBM(attributes = c("uniprot_gn_symbol", "uniprotswissprot","external_synonym"), filters = "uniprot_gn_symbol", values = "RPP29", bmHeader = T, mart = mart104.mm)
# entrezgene_id


unique(unlist(DEgp_genenames))[which(unique(unlist(DEgp_genenames))%in%toupper(oris.mm104))]
length(unique(unlist(DEgp_genenames))[which(unique(unlist(DEgp_genenames))%in%toupper(oris.mm104))]
)


# search in all DE gene list ----------------------------------------------

library(rlist)
same_genes <- unique(unlist(OGG2GeneName.lst))[which(unique(unlist(OGG2GeneName.lst))%in%toupper(oris.mm104))]
# same_genes <- unique(unlist(all_genenames))[which(unique(unlist(all_genenames))%in%toupper(oris.mm104))]

length(same_genes)

matched.df <- data.frame(matrix(nrow = 0,ncol = 2))
for (i in 1:length(same_genes)) {
  mat <- list.search(OGG2GeneName.lst, any(. == same_genes[i]))
  # mat <- list.search(all_genenames, any(. == same_genes[i]))
  if (length(names(mat))!=0) {
    # cat(sprintf("\"%s\" \"%s\" \"%s\"\n", i, names(mat),same_genes[i]))
    matched <- data.frame(sapply(mat,FUN = function(values) paste(values,collapse = "|")))
    matched$X1 <- rownames(matched)
    colnames(matched) <- c("X2","X1")
    matched <- matched[,c("X1","X2")]
    matched.df <- rbind(matched.df,matched)
  }
}

colnames(matched.df) <- c("matched_gp_nbr","matched_ori_gene")

matched.df


all_genenames.tb <- data.frame(sapply(OGG2GeneName.lst, FUN = function(values) paste(values,collapse = "|")))
# all_genenames.tb <- data.frame(sapply(all_genenames, FUN = function(values) paste(values,collapse = "|")))
colnames(all_genenames.tb) <- "gene_names"
# head(all_genenames.tb)



# find ori matched gp expression ------------------------------------------

countMatrix[which(rownames(countMatrix)%in%str_remove(matched.df$matched_gp_nbr,"group_")),]

matched_ori_exp <- resNoFilt$table[which(rownames(resNoFilt$table)%in%str_remove(matched.df$matched_gp_nbr,"group_")),]
de_exp <- resNoFilt$table[which(rownames(resNoFilt$table)%in%DEnames),]

# matched_ori_exp <- lrt$table[which(rownames(lrt$table)%in%str_remove(matched.df$matched_gp_nbr,"group_")),]
# de_exp <- lrt$table[which(rownames(lrt$table)%in%DEnames),]

rownames(matched_ori_exp) <- paste0("group_",rownames(matched_ori_exp) )
rownames(de_exp) <- paste0("group_",rownames(de_exp) )

matched_ori_exp <- merge(matched_ori_exp,all_genenames.tb, by = "row.names")
de_exp <- merge(de_exp,all_genenames.tb, by = "row.names")


output_file_path <- paste0(output_path,cv_dataset,"_matched_ori_exp.csv")
write.csv(matched_ori_exp,row.names = FALSE, output_file_path, quote = FALSE)
output_file_path <- paste0(output_path,cv_dataset,"_de_exp.csv")
write.csv(de_exp,row.names = FALSE, output_file_path, quote = FALSE)


exp_tb$Group_Num
raw_matched_ori_TPM <- exp_tb[which(exp_tb$Group_Num%in%str_remove(matched.df$matched_gp_nbr,"group_")),]

