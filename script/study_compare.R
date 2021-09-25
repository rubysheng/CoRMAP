setwd("/Users/ruby/Documents/GitHub/comparative-transcriptomic-analysis-pip/sample/output/R_analysis/")
output_path <- "./"
adj.method <- 'none'
pv <- 0.05
cv_datasets <- c("PRJNA252803","PRJNA529794")

library(edgeR)

# annotation --------------------------------------------------------------
OGG2GeneName <- read.delim("./Unique_OGG_GENENAME.txt", header = T, sep = "\t")
OGG2GeneName <- unique(OGG2GeneName)
OGG2GeneName.lst <- OGG2GeneName %>%
  split(f = as.factor(.$Group)) %>%
  lapply(., "[", , c("Gene_name"))
all_genenames.tb <- data.frame(sapply(OGG2GeneName.lst, FUN = function(values) paste(values,collapse = "|")))
colnames(all_genenames.tb) <- "gene_names"


# deg analysis - STUDY 1 (mouse_1) ------------------------------------------------------------


cv_dataset <- cv_datasets[1]
load(paste0("./",cv_dataset,"_countMX.RData"))

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


resNoFilt <- topTags(lrt,n=nrow(lrt$table),adjust.method = adj.method,p.value = pv)
# sum(resNoFilt$table$FDR < 0.05)
sigDownReg <- resNoFilt$table[resNoFilt$table$PValue<pv,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]

rownames(sigDownReg) <- paste0("group_",rownames(sigDownReg) )
sigDownReg <- merge(sigDownReg,all_genenames.tb, by = "row.names")

assign(cv_dataset,sigDownReg)



# deg analysis - STUDY 2 (mouse_3) ------------------------------------------------------------


cv_dataset <- cv_datasets[2]
load(paste0("./",cv_dataset,"_countMX.RData"))

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


resNoFilt <- topTags(lrt,n=nrow(lrt$table),adjust.method = adj.method,p.value = pv)
# sum(resNoFilt$table$FDR < 0.05)
sigDownReg <- resNoFilt$table[resNoFilt$table$PValue<pv,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]

rownames(sigDownReg) <- paste0("group_",rownames(sigDownReg) )
sigDownReg <- merge(sigDownReg,all_genenames.tb, by = "row.names")

assign(cv_dataset,sigDownReg)


# match between two studies -------------------------------------------------------------------

cv_dataset <- cv_datasets[1]
colnames(PRJNA252803)[-1] <- paste0(cv_dataset,"_",colnames(PRJNA252803)[-1])
cv_dataset <- cv_datasets[2]
colnames(PRJNA529794)[-1] <- paste0(cv_dataset,"_",colnames(PRJNA529794)[-1])

colnames(PRJNA252803)[1] <- "Group"
colnames(PRJNA529794)[1] <- "Group" 
# head(PRJNA529794)
# dim(PRJNA529794)
matched.df <- merge(PRJNA252803[-6],PRJNA529794[-6], by = "Group")
rownames(matched.df) <- matched.df$Group
matched.df <- merge(all_genenames.tb,matched.df, by = "row.names")
matched.df <- matched.df[-1]

output_file_path <- paste0(output_path,"matchedOGG_",cv_datasets[1],"exp_vs_",cv_datasets[2],"exp.csv")
write.csv(matched.df,row.names = FALSE, output_file_path, quote = FALSE)

dim(matched.df)

# plot correlation --------------------------------------------------------

library(ggplot2)
library(plotly)

# Basic scatter plot

p <- matched.df %>%
  ggplot(aes(x=PRJNA252803_logFC, y=PRJNA529794_logFC,label=gene_names)) + 
  geom_point() +
  ggtitle("Correlation of expression values for matched OGGs") +
  xlab("Mouse_1 Log2(Fold-change)") + ylab("Mouse_3 Log2(Fold-change)") +
  geom_abline(intercept = 0, slope = 1, color="blue", linetype="dashed", size=.5) +
  geom_abline(intercept = 0, slope = -1, color="blue", linetype="dashed", size=.5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlim(-4, 4) + ylim(-4,4)


ggplotly(p)
