# install.packages("ggpubr")
# install.packages("PerformanceAnalytics")
library(ggpubr)
library(ggplot2)
library(tidyr)
library(plyr)
library(PerformanceAnalytics)


# To load the data again
load("aveTPM_output.RData")

TPM_output <- data.frame("project"=NULL,"Min."=NULL,"1st Qu."=NULL,"Median"=NULL,"Mean"=NULL,"3rd Qu."=NULL,"Max."=NULL,"IQR"=NULL,"stdev"=NULL)

GN <- as.vector(unique(aveTPM_output$group_num))
 
for (i in 1:length(GN)){
  gp_num <- GN[i]
  sub_df <- aveTPM_output[which(aveTPM_output$group_num==gp_num),]
  long_sub <- sub_df[,c(2,1,5,6,3,4)]
  long_sub$gene <- factor(long_sub$gene)
  
  # general TPM distribution stats
  stats_info <- summary(long_sub$TPM)
  stats_df <- data.frame("GROUP_NUM"=gp_num,t(matrix(stats_info)),IQR(long_sub$TPM),sd(long_sub$TPM))
  colnames(stats_df)[-1]<- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","IQR","stdev")
  TPM_output <- rbind(TPM_output,stats_df)
  
  # Export the image file
  file_path <- paste0("logTPM_gn_g",gp_num,".png")
  p <- ggplot(long_sub, aes(x=species, y=log(TPM))) + geom_point() + geom_rug() + theme(axis.text.x = element_text(face="bold", color="black", size=9, angle=45))
  ggsave(p,filename = file_path, width=3, height=3)

  
}

# save files
write.csv(TPM_output,"TPM_output_gp.csv", row.names = FALSE)


