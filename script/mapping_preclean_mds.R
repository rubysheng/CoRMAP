library(Glimma)
library(edgeR)
library(DESeq2)

# ----------------
# library("pasilla")
# pasCts <- system.file("extdata",
#                       "pasilla_gene_counts.tsv",
#                       package="pasilla", mustWork=TRUE)
# pasAnno <- system.file("extdata",
#                        "pasilla_sample_annotation.csv",
#                        package="pasilla", mustWork=TRUE)
# cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
# coldata <- read.csv(pasAnno, row.names=1)
# coldata <- coldata[,c("condition","type")]
# coldata$condition <- factor(coldata$condition)
# coldata$type <- factor(coldata$type)



mm252803_cts <- as.matrix(read.csv("PRJNA252803/count_matrix_cleaned.csv",row.names = 1))
mm529794_cts <- as.matrix(read.csv("PRJNA529794/count_matrix_cleaned.csv",row.names = 1))

mm252803_meta_tb <- read.csv("PRJNA252803/meta_tb.csv",row.names = 1)
mm252803_meta_tb$treatment <- gsub(" ","_",mm252803_meta_tb$treatment)
mm252803_meta_tb$treatment <- factor(mm252803_meta_tb$treatment)
mm529794_meta_tb <- read.csv("PRJNA529794/meta_tb.csv",row.names = 1)
mm529794_meta_tb <- mm529794_meta_tb[order(rownames(mm529794_meta_tb)),]
mm529794_meta_tb$Cell_type <- factor(mm529794_meta_tb$Cell_type)
mm529794_meta_tb$Mouse_ID <- gsub(" ","_",mm529794_meta_tb$Mouse_ID)
mm529794_meta_tb$treatment <- gsub(" ","_",mm529794_meta_tb$treatment)
mm529794_meta_tb$Mouse_ID <- factor(mm529794_meta_tb$Mouse_ID)
mm529794_meta_tb$treatment <- factor(mm529794_meta_tb$treatment)


# ------------------
# dds <- DESeqDataSetFromMatrix(countData = cts,
#                               colData = coldata,
#                               design = ~ condition)
# dds

mm252803_dds <- DESeqDataSetFromMatrix(countData = mm252803_cts,
                                       colData = mm252803_meta_tb,
                                       design = ~ treatment)

mm529794_dds <- DESeqDataSetFromMatrix(countData = mm529794_cts,
                                       colData = mm529794_meta_tb,
                                       design = ~ treatment)

# ------------------
# data transform using vst (for sample size more than 20)
# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)

mm252803_vsd <- vst(mm252803_dds, blind=FALSE)
mm252803_rld <- rlog(mm252803_dds, blind=FALSE)
mm529794_vsd <- vst(mm529794_dds, blind=FALSE)


# ------------------
# # Data quality assessment by sample clustering and visualization
# sampleDists <- dist(t(assay(vsd)))
# 
library("RColorBrewer")
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
# colnames(sampleDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDists,
#          clustering_distance_cols=sampleDists,
#          col=colors)




mm252803_sampleDists <- dist(t(assay(mm252803_vsd)))
# mm252803_sampleDists <- dist(t(assay(mm252803_rld)))
mm252803_sampleDistMatrix <- as.matrix(mm252803_sampleDists)
rownames(mm252803_sampleDistMatrix) <- paste(rownames(mm252803_sampleDistMatrix),mm252803_vsd$treatment,sep = "-")
# rownames(mm252803_sampleDistMatrix) <- paste(rownames(mm252803_sampleDistMatrix),mm252803_rld$treatment,sep = "-")
colnames(mm252803_sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(mm252803_sampleDistMatrix,
         clustering_distance_rows=mm252803_sampleDists,
         clustering_distance_cols=mm252803_sampleDists,
         col=colors)


mm529794_sampleDists <- dist(t(assay(mm529794_vsd)))
# mm529794_sampleDists <- dist(t(assay(mm529794_rld)))
mm529794_sampleDistMatrix <- as.matrix(mm529794_sampleDists)
rownames(mm529794_sampleDistMatrix) <- paste(mm529794_vsd$treatment,
                                             mm529794_vsd$Cell_type,
                                             mm529794_vsd$Mouse_ID,sep="-")
# rownames(mm529794_sampleDistMatrix) <- paste(mm529794_rld$treatment)
colnames(mm529794_sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(mm529794_sampleDistMatrix,
         clustering_distance_rows=mm529794_sampleDists,
         clustering_distance_cols=mm529794_sampleDists,
         col=colors)

# ------------------
# # MDS:1
# glMDSPlot(dds)
# glMDSPlot(dds, 
#           labels = rownames(coldata), 
#           groups = coldata[,c("condition", "type")], 
#           folder = "mds")

glMDSPlot(mm252803_dds)
glMDSPlot(mm529794_dds)

# # MDS:2
# library(ggplot2)
# sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(rownames(sampleDistMatrix),vsd$condition, vsd$type, sep="-")
# colnames(sampleDistMatrix) <- NULL
# mds <- data.frame(cmdscale(sampleDistMatrix))
# mds <- cbind(mds, colData(vsd))
# df <- data.frame(Sample=rownames(mds), mds)
# qplot(X1,X2,color=condition,shape=type,data=as.data.frame(mds))
# ggplot(df, aes(x=X1, y=X2, label=Sample)) +
#   geom_mark_ellipse(aes(color = condition,
#                         label=condition),
#                     expand = unit(0.5,"mm"),
#                     label.buffer = unit(-5, 'mm'))+
#   geom_point(aes(color=condition))+
#   theme(legend.position = "none")


library(ggplot2)
mm529794_sampleDists <- dist(t(assay(mm529794_vsd)))
mm529794_sampleDistMatrix <- as.matrix(mm529794_sampleDists)
rownames(mm529794_sampleDistMatrix) <- paste(rownames(mm529794_sampleDistMatrix),
                                             mm529794_vsd$treatment,
                                             mm529794_vsd$Cell_type,
                                             mm529794_vsd$Mouse_ID,sep="-")
colnames(mm529794_sampleDistMatrix) <- NULL
mm529794_mds <- data.frame(cmdscale(mm529794_sampleDistMatrix))
mm529794_mds <- cbind(mm529794_mds, colData(mm529794_vsd))
mm529794_df <- data.frame(Sample=rownames(mm529794_mds), mm529794_mds)
qplot(X1,X2,color=treatment,shape=Cell_type,data=as.data.frame(mm529794_mds))
ggplot(mm529794_df, aes(x=X1, y=X2, label=Sample)) +
  geom_mark_ellipse(aes(color = treatment,
                        label=treatment),
                    expand = unit(0.5,"mm"),
                    label.buffer = unit(-5, 'mm'))+
  geom_point(aes(color=treatment))+
  theme(legend.position = "none")

# ------------------
# select samples for DEG Analysis
# PRJNA252803
# 48,49 & 52,53
mm252803_cts <- as.matrix(read.csv("PRJNA252803/count_matrix_cleaned.csv",row.names = 1))
cleaned_cts <- mm252803_cts[,which(colnames(mm252803_cts)%in%c("SRR1408848","SRR1408849","SRR1408852","SRR1408853"))]
mm252803_meta_tb <- as.matrix(read.csv("PRJNA252803/meta_tb.csv",row.names = 1))
cleaned_meta_tb <-  data.frame(treatment=mm252803_meta_tb[which(rownames(mm252803_meta_tb)%in%c("SRR1408848","SRR1408849","SRR1408852","SRR1408853")),])
cleaned_meta_tb$treatment <- gsub(" ","_",cleaned_meta_tb$treatment)
cleaned_meta_tb$treatment <- factor(cleaned_meta_tb$treatment)

cleaned_dds <- DESeqDataSetFromMatrix(countData = cleaned_cts,
                                       colData = cleaned_meta_tb,
                                       design = ~ treatment)

cleaned_vsd <- vst(cleaned_dds, blind=FALSE)

cleaned_sampleDists <- dist(t(assay(cleaned_vsd)))
cleaned_sampleDistMatrix <- as.matrix(cleaned_sampleDists)
rownames(cleaned_sampleDistMatrix) <- paste(rownames(cleaned_sampleDistMatrix),cleaned_vsd$treatment,sep = "-")
colnames(cleaned_sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(cleaned_sampleDistMatrix,
         clustering_distance_rows=cleaned_sampleDists,
         clustering_distance_cols=cleaned_sampleDists,
         col=colors)

glMDSPlot(cleaned_dds)

# select samples for DEG Analysis
# PRJNA529794
# combine homecage and non_shock together as "control"

# -----------
# DEG
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
# not use LRT test (for time course experiments) ; but default Wald test
# FDR 0.01
dds <- DESeq(dds)
res <- results(dds,
               alpha=0.01)
# res <- results(dds,
#                alpha=0.05,
#                contrast=c("condition","treated","untreated"),
#                name="condition_treated_vs_untreated")  #alpha -> FDR
# filter by log2 fold change 
res <- res[!is.na(res$log2FoldChange),]
res <- res[res$log2FoldChange >= 1, ]
res

mm252803_dds <- DESeqDataSetFromMatrix(countData = mm252803_cts,
                                       colData = mm252803_meta_tb,
                                       design = ~ treatment)

mm529794_meta_tb$treatment_cleaned <- mm529794_meta_tb$treatment
mm529794_meta_tb$treatment_cleaned <- gsub("HomeCage","control",mm529794_meta_tb$treatment_cleaned)
mm529794_meta_tb$treatment_cleaned <- gsub("Non_Shock","control",mm529794_meta_tb$treatment_cleaned)
mm529794_meta_tb$treatment_cleaned <- factor(mm529794_meta_tb$treatment_cleaned)
mm529794_dds <- DESeqDataSetFromMatrix(countData = mm529794_cts,
                                       colData = mm529794_meta_tb,
                                       design = ~ treatment_cleaned)
# not use LRT test (for time course experiments) ; but default Wald test
# FDR 0.01
mm252803_dds <- DESeq(mm252803_dds)
mm252803_res <- results(mm252803_dds,
               alpha=0.01)
mm252803_res <- mm252803_res[!is.na(mm252803_res$log2FoldChange),]
mm252803_res <- mm252803_res[mm252803_res$log2FoldChange >= 1, ]

mm529794_dds <- DESeq(mm529794_dds)
mm529794_res <- results(mm529794_dds,
                        alpha=0.01)
mm529794_res <- mm529794_res[!is.na(mm529794_res$log2FoldChange),]
mm529794_res <- mm529794_res[mm529794_res$log2FoldChange >= 1, ]

# ---------------------
mm252803_res
mm529794_res



