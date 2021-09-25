#!/usr/bin/env Rscript

###########################################################################################
## Step _:                                                                               ##
## cross validation of sampled subgroups in specific dataset                             ##
## (TPM_Value).                                                                          ##
##     ##
##                                               ##
###########################################################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager","argparse")
# BiocManager::install()
# 
# suppressPackageStartupMessages(library(argparse))
# 
# parser = ArgumentParser()
# parser$add_argument("--name_list", help="a file containing a list of dataset names", required=TRUE)
# parser$add_argument("--input", help="the folder holding input orthologs expression matrices (end with *.allsp_clugene.TPM)", required=TRUE)
# parser$add_argument("--output", help="path to save the output tables", required=TRUE)
# args = parser$parse_args()
# 
# name_list_file = args$name_list
# mx_path = args$input
# output_path = args$output
# 
# files <- scan(name_list_file, character(), quote = "")



# name_list_file = "./sample/output/R_analysis/input_dataset_name.txt"
# input_path = "./sample/output/R_analysis/"
output_path = "./"
# treatcode_path = "treat_code.txt"
# projectname = "Project_name_by_species.txt"
cv_dataset = "PRJNA529794"
input_Expr_rdata = "./Expr_output.RData"


# load data ---------------------------------------------------------------

# files <- scan(name_list_file, character(), quote = "")
load(input_Expr_rdata)
# load(paste0(input_path,"fc_table.RData"))
exp_tb <- Expr_output

# functions ---------------------------------------------------------------


# subset all data of the selected dataset ---------------------------------
exp_tb <- exp_tb[which(exp_tb$Project==cv_dataset),]
exp_tb$Sample <- paste(exp_tb$Project,exp_tb$Sample,sep = "-")


# average TPM of genes in the same group for each dataset -----------------
avgtpm_gp <- data.frame(Group_Num=character(),
                        TPM_Value=numeric(),
                        Sample=character(),
                        Treatment=character(),
                        Project=character())

for (smp in unique(exp_tb$Sample)) {
  sg_smp <- subset(exp_tb,Sample==smp)
  avgTPM_by_gpnb <- cbind(aggregate(sg_smp[,"TPM_Value"], list(sg_smp$Group_Num), mean),
                          smp,unique(sg_smp$Treatment),unique(sg_smp$Project))
  colnames(avgTPM_by_gpnb) <- colnames(sg_smp)[c(1,9,6,7,4)]
  avgtpm_gp <- rbind(avgtpm_gp,avgTPM_by_gpnb)
}
rm(sg_smp,smp,avgTPM_by_gpnb)


if (cv_dataset == "PRJNA529794") {
  # generating the sample/treatment description tb --------------------------
  avgtpm_tb <- cbind(avgtpm_gp, Sample_spli<-data.frame(do.call('rbind', strsplit(as.character(avgtpm_gp$Sample), '_', fixed=TRUE))))
  # head(avgtpm_tb)
  for (cond in unique(avgtpm_tb$X1)) {
    max_pos_rep <- max(as.numeric(substring(unique(avgtpm_tb[which(avgtpm_tb$X1==cond & avgtpm_tb$X3=="pos"),"X4"]), 4)))
    avgtpm_tb[which(avgtpm_tb$X1==cond & avgtpm_tb$X3=="neg"),"X4"] <- paste0("rep",as.numeric(substring(avgtpm_tb[which(avgtpm_tb$X1==cond & avgtpm_tb$X3=="neg"),"X4"], 4))+max_pos_rep)
  }
  rm(max_pos_rep,cond)
  
  avgtpm_tb$X2 <- paste(avgtpm_tb$X1,avgtpm_tb$X4,sep="_")
  sample_cat <- unique(avgtpm_tb[,c("X1","X2")])
  colnames(sample_cat) <- c("Treatment","Sample")
  
  # sample subgroups for "Samples" ------------------------------------------
  # set.seed(1)
  smp_colnames <- c()
  for (tret in unique(sample_cat$Treatment)) {
    # if (tret %in% c("PRJNA529794-HomeCage","PRJNA529794-NonShock")) {
    #   smp_cols <- sample(sample_cat$Sample[which(sample_cat$Treatment==tret)],2)
    # } else {
    #   smp_cols <- sample(sample_cat$Sample[which(sample_cat$Treatment==tret)],4)
    # }
    smp_cols <- sample_cat$Sample[which(sample_cat$Treatment==tret)]
    smp_colnames <- c(smp_colnames,smp_cols)
  }
  rm(tret,smp_cols)
  # smp_colnames
  exp_colnames <- smp_colnames
  
  max_hc_rep <- max(as.numeric(substring(exp_colnames[str_detect(exp_colnames,"HomeCage")], regexpr("_rep", exp_colnames[str_detect(exp_colnames,"HomeCage")])+4)))
  exp_colnames[str_detect(exp_colnames,"NonShock")] <- paste0(substring(exp_colnames[str_detect(exp_colnames,"NonShock")], 0,regexpr("_rep", exp_colnames[str_detect(exp_colnames,"NonShock")])+3), as.numeric(substring(exp_colnames[str_detect(exp_colnames,"NonShock")], regexpr("_rep", exp_colnames[str_detect(exp_colnames,"NonShock")])+4))+max_hc_rep)
  
  exp_colnames <- str_replace_all(exp_colnames,pattern = c("HomeCage|NonShock"), replacement = "Control")
  
  # pats <- unique(substring(exp_colnames, regexpr("-", exp_colnames)+1,regexpr("_rep", exp_colnames)-1))
  # for (pat in pats) {
  #   exp_colnames[str_detect(exp_colnames,pattern = pat)] <- paste0(substring(exp_colnames[str_detect(exp_colnames,pattern = pat)], 0, regexpr("_rep", exp_colnames[str_detect(exp_colnames,pattern = pat)]) - 1),"-sub",1:4)
  # }
  # subgroup_info <- data.frame(Sample_name=smp_colnames,Subgroup_name=exp_colnames)
  
  # generate subgroups' expression matrix -----------------------------------
  # long to wide 
  library(tidyr)
  
  # head(avgtpm_tb)
  # str(avgtpm_tb)
  avgtpm_tb$Sample_name <- avgtpm_tb$X2
  exp_mx <- spread(avgtpm_tb[,c("Group_Num","TPM_Value","Sample_name")], Sample_name, TPM_Value)
  rownames(exp_mx) <- exp_mx$Group_Num
  countMatrix <- exp_mx[,smp_colnames]
  colnames(countMatrix) <- exp_colnames
  # head(countMatrix)
  group <- factor(substring(colnames(countMatrix), 0, regexpr("_rep", colnames(countMatrix)) - 1))
  # group <- factor(substring(colnames(countMatrix), 0, regexpr("-sub", colnames(countMatrix)) - 1))
} else {
  exp_mx <- spread(avgtpm_gp[,c("Group_Num","TPM_Value","Sample")], Sample, TPM_Value)
  rownames(exp_mx) <- exp_mx$Group_Num
  countMatrix <- exp_mx[,unique(avgtpm_gp$Sample)]
  # head(countMatrix)
  group <- factor(substring(colnames(countMatrix), 0, regexpr("_rep", colnames(countMatrix)) - 1))
}
output_file_path <- paste0(output_path,cv_dataset,"_countMX.RData")
save(countMatrix,group,file = output_file_path)



# DEG analysis ------------------------------------------------------------
library(edgeR)

y <- DGEList(counts = countMatrix, group = group)
# y

keep <- rowSums(cpm(y)>1) >= 2
y_filter0 <- y[keep, , keep.lib.sizes=FALSE]
# y_filter0 

y_normed <- calcNormFactors(y_filter0 ,method = 'TMM')
# y_normed

plotMDS(y_normed, col = rep(c('red', 'blue'), each = 5), dim = c(1, 2))

subGroup <- factor(substring(colnames(countMatrix), regexpr("_rep", colnames(countMatrix))+4))
# subGroup <- factor(substring(colnames(countMatrix), regexpr("-sub", colnames(countMatrix))+4))
design <- model.matrix(~ subGroup+group)
rownames(design) <- colnames(y_normed)
# design

dge <- estimateDisp(y_normed, design, robust=TRUE)
dge$common.dispersion


# plot
plotBCV(dge)

#negative binomial generalized log-linear model 
fit <- glmFit(dge, design, robust = TRUE)     
lrt <- glmLRT(fit)   

topTags(lrt)
# output_file_path <- paste0(output_path,cv_dataset,"_lrt.csv")
# write.csv(topTags(lrt, n = nrow(dge$counts)), output_file_path, quote = FALSE)
# output_file_path <- paste0(output_path,cv_dataset,"_lrt.RData")
# save(lrt,file = output_file_path)

# plot --------------------------------------------------------------------
dge_de <- decideTestsDGE(lrt, adjust.method = 'none', p.value = 0.05)  
summary(dge_de)

plotMD(lrt, status = dge_de, values = c(1, -1), col = c('blue', 'red'))     
abline(h = c(-1, 1), col = 'gray', lty = 2)


isDE <- as.logical(dge_de)
DEnames <- rownames(dge)[isDE]
head(DEnames)

plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

# MA plot -----------------------------------------------------------------

library(geneplotter)

plotMA(res, main="DESeq2", ylim=c(-2,2))


# HEATMAP diagram ---------------------------------------------------------

sum(res$padj < 0.1, na.rm=TRUE)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:1000]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("name","condition")])
# pdf('heatmap1000.pdf',width = 6, height = 7)
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
# dev.off()

