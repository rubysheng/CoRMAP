setwd("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/")
# install.packages("ggpubr")
# install.packages("PerformanceAnalytics")
library(ggpubr)
library(ggplot2)
library(tidyr)
library(plyr)
library(PerformanceAnalytics)


files <- c("Bombina","PRJNA252803","PRJNA287145","PRJNA287152","PRJNA302146","PRJNA316996","PRJNA381689","PRJNA390522","PRJNA419677","PRJNA422916","PRJNA450614","PRJNA451011","PRJNA475804","PRJNA529794","PRJNA541005")

# To load the data again
load("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/analysis/aveTPM_output.RData")

i=1
for (i in 1:length(files)){
  project <- files[i]
  sub_df <- aveTPM_output[which(aveTPM_output$project==project),]
  long_sub <- sub_df[,c(2,1,5,6,3,4)]
  long_sub$gene <- factor(long_sub$gene)
  
  wide_sub <- spread(long_sub, treatment,TPM)
  output_file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/wide_format_",project,"_mulgp.csv")
  write.csv(wide_sub,output_file_path, row.names = FALSE)
  # head(wide_sub)
  
  # mu <- ddply(long_sub, "treatment", summarise, grp.mean=mean(log(TPM)), grp.iqr=IQR(log(TPM)))
  # lim <- (max(long_sub$TPM)- max(mu$grp.mean))*0.0001
  # lim_p <- ggplot(long_sub, aes(x=TPM, color=treatment)) +
  #           geom_histogram(binwidth=0.001,fill="white") +
  #           xlim(0,lim) +
  #           geom_vline(data=mu, aes(xintercept=grp.iqr, color=treatment),
  #              linetype="dotdash") +
  #           theme(legend.position="top")
  # lim_p
  
  # p <- ggplot(long_sub, aes(x=TPM, color=treatment)) +
  #           geom_histogram(binwidth=1,fill="white") +
  #           geom_vline(data=mu, aes(xintercept=grp.mean, color=treatment),
  #            linetype="dashed") + 
  #           theme(legend.position="top")
  # p
  
  # ggarrange(p, lim_p, labels = c("A", "B"), ncol = 2, nrow = 1)
  
  wide_log_sub <- wide_sub[,5:ncol(wide_sub)]
  wide_log_sub$gene <- wide_sub$gene
  output_file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/wide_format_",project,"_mulgp.txt")
  write.table(wide_log_sub, file = output_file_path, sep = "\t", row.names = FALSE)
  # # cols <- colnames(wide_log_sub)
  # # wide_log_sub[cols] <- log(wide_log_sub[cols])
  # # distribution : density and correlation plot
  # # Export the image file
  # file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/analysis/logTPM_trt_",project,".bmp")
  # bmp(file_path, width=1000, height=1000, bg="white", res=300)
  # chart.Correlation(wide_log_sub, histogram=TRUE, pch=19)
  # # ggplot(long_sub, aes(x=log(TPM))) + geom_density(fill="white") + 
  # #   facet_grid(treatment ~ .)
  # dev.off()
  
  # dp <- ggdotplot(long_sub, x = "treatment", y = "TPM", 
  #                 color = "treatment", fill = "treatment",
  #                 dotsize = 0.05, binwidth = 1, width = 6)
  # dp
  
  # # Export the image file
  # file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/analysis/logTPM_gn_",project,".bmp")
  # bmp(file_path, width=2500, height=2500, bg="white", res=300)
  # ggplot(long_sub, aes(x=group_num, y=log(TPM))) + geom_point() + geom_rug()
  # dev.off()
  
}

# GN <- as.vector(unique(aveTPM_output$group_num))
# i=1
# for (i in 1:length(GN)){
#   gp_num <- GN[i]
#   sub_df <- aveTPM_output[which(aveTPM_output$group_num==gp_num),]
#   long_sub <- sub_df[,c(2,1,5,6,3,4)]
#   long_sub$gene <- factor(long_sub$gene)
#   
#   # wide_sub <- spread(long_sub, species, TPM)
#   # wide_sub
#   
#   # mu <- ddply(long_sub, "treatment", summarise, grp.mean=mean(log(TPM)), grp.iqr=IQR(log(TPM)))
#   # lim <- (max(long_sub$TPM)- max(mu$grp.mean))*0.0001
#   # lim_p <- ggplot(long_sub, aes(x=TPM, color=treatment)) +
#   #           geom_histogram(binwidth=0.001,fill="white") +
#   #           xlim(0,lim) +
#   #           geom_vline(data=mu, aes(xintercept=grp.iqr, color=treatment),
#   #              linetype="dotdash") +
#   #           theme(legend.position="top")
#   # lim_p
#   
#   # p <- ggplot(long_sub, aes(x=TPM, color=treatment)) +
#   #           geom_histogram(binwidth=1,fill="white") +
#   #           geom_vline(data=mu, aes(xintercept=grp.mean, color=treatment),
#   #            linetype="dashed") + 
#   #           theme(legend.position="top")
#   # p
#   
#   # ggarrange(p, lim_p, labels = c("A", "B"), ncol = 2, nrow = 1)
#   
#   # wide_log_sub <- wide_sub[,5:ncol(wide_sub)]
#   # cols <- colnames(wide_log_sub)
#   # wide_log_sub[cols] <- log(wide_log_sub[cols])
#   # # distribution : density and correlation plot
#   # # Export the image file
#   # file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/analysis/logTPM_trt_",project,".bmp")
#   # bmp(file_path, width=1000, height=1000, bg="white", res=300)
#   # chart.Correlation(wide_log_sub, histogram=TRUE, pch=19)
#   # # ggplot(long_sub, aes(x=log(TPM))) + geom_density(fill="white") + 
#   # #   facet_grid(treatment ~ .)
#   # dev.off()
#   
#   # dp <- ggdotplot(long_sub, x = "species", y = "TPM",
#   #                 color = "species", fill = "species",
#   #                 dotsize = 0.05, binwidth = 1, width = 6)
#   # dp
#   
#   # Export the image file
#   file_path <- paste0("C:/Users/tts_p/Documents/MSc/MSc_thesis/RNA-seq/expression/output/analysis/logTPM_gn_g",gp_num,".bmp")
#   bmp(file_path, width=2500, height=2500, bg="white", res=300)
#   ggplot(long_sub, aes(x=species, y=log(TPM))) + geom_point() + geom_rug()
#   dev.off()
#   
# }

