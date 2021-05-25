library(readxl)
library(dplyr)
library(reshape2)
library(sva)
library(pamr)
library(limma)
library(rlist)
library(tidyverse)


name_list_file = "input_dataset_name.txt"
input_path = "./"
output_path = "./"
treatcode_path = "treat_code.txt"
projectname = "Project_name_by_species.txt"





# functions ---------------------------------------------------------------

# generate the expression matrix
tag <- function(input_df, col_in1, col_in2, col_out, pattern) {
  input_df[col_out] <- paste(input_df[col_in1],input_df[col_in2],sep = pattern)
}
SVA_Edata_generater <- function(data = Expr_output, group1 = "Group_Num", group2, group3){
  group_var <- c(group1, group2, group3)
  caled_df <- data %>%
    dplyr::group_by(across({{ group_var }})) %>%
    dplyr::summarize(TPM_Value = mean(as.numeric(TPM_Value)))
  caled_df$col_name <- apply(caled_df,1,tag, col_in1 = group2, col_in2 = group3, col_out = "col_name", pattern = "__")
  output_table <-  data.frame(caled_df)
  
  return(output_table)
}

div <- function(line){
  line[3] <- line[2]/ line[1]
  line[4] <- line[1]/ line[2]
  if (line[3]>=line[4]){ line[3] }else{ -line[4] }
}



# generate the a master table to save all TPM values ------------------------------------------------------------
files <- scan(name_list_file, character(), quote = "")

Expr_output <- data.frame("Group_Num"=NULL,"Taxa_Group"=NULL,"Species"=NULL,
                          "Project"=NULL,"gene"=NULL,"Sample"=NULL,"Treatment"=NULL,"TPM_value"=NULL)

for (i in 1:length(files)){
  project <- files[i]
  file_path <- paste0(input_path,project,"_mulgp_clugene_sty.xlsx")
  input_df <- read_xlsx(file_path)
  colnames(input_df)[c(5,10)] <- c("treatment","study_design")
  input_df <- as.data.frame(input_df)
  input_df$TPM_Value <- as.numeric(input_df$TPM_Value)
  treatmnt <- data.frame(do.call('rbind', strsplit(as.character(input_df$treatment),'_rep',fixed=TRUE)))
  input_df$treatment2 <- treatmnt$X1
  rm(treatmnt)
  
  expr_df <- input_df[,c(1,18,2,3,4,5,19,6)]
  colnames(expr_df)[6:7] <- c("Sample","Treatment")
  Expr_output <- rbind(Expr_output,expr_df)
  
}


treatcod <- read.delim(treatcode_path)
a <- unique(Expr_output$Treatment)
b <- treatcod$treatment
if (length(a[which(!a %in% b)])!=0) {
  add <- data.frame(treatment=a[which(!a %in% b)], TREAT_Code=c(1,2,1,2,2,2,1,1,1,1))
  rm(a);rm(b)
  treatcod <- rbind(treatcod,add)
  rm(add)
}
colnames(treatcod)[1] <- "Treatment"


Expr_output <- inner_join(Expr_output,treatcod)
Expr_output <- Expr_output[,c(1:7,9,8)]
save(Expr_output, file = paste0(output_path,"Expr_output.RData" ))
write.csv(Expr_output,paste0(output_path,"Expr_output.csv"), row.names = FALSE)


# average expression values by project and treat_code ------------------------
pj_treat <- SVA_Edata_generater(Expr_output, group2 = "Project", group3 = "TREAT_Code")

# generate the Expression data matrix
edata <- dcast(pj_treat[,c(1,5,4)], Group_Num ~ col_name, value.var = "TPM_Value")

df <- edata
rownames(df) <- df[,1]
df <- df[,-1]
df[is.na(df)] <- 0
test_df <- df
test_df[] <- lapply(test_df, function(x) as.numeric(as.character(x)))
# test_df <- as.matrix(test_df)


df_fc <- test_df
for (i in 1:length(files)){
  project <- files[i]
  project_cols <- test_df[,which(str_detect(colnames(test_df), project))]
  
  df_fc <- cbind(df_fc, round(apply(project_cols, 1, FUN=function(x2) div(x2)),2))
  colnames(df_fc)[length(colnames(df_fc))] <- paste0(project,"_fc")
}

fc_table <- df_fc[,c((ncol(test_df)+1):ncol(df_fc))]

colnames(fc_table) <- str_replace(colnames(fc_table),"_fc","")
fc_table$group_num <- paste0("group_",rownames(fc_table))

# interactive scatter plot ------------------------------------------------

library(ggplot2)
library(plotly)

cate_class <- function(line){
  x <- as.numeric(line[1])
  y <- as.numeric(line[2])
  if      (y>x & x>0)         { return("PRJNA529794>PRJNA252803 & both up")}
  else if (y>x & y>0 & 0>x )   { return("PRJNA529794>PRJNA252803 & PRJNA529794-up PRJNA252803-down")}
  else if (0>y & y>x  )   { return("PRJNA529794>PRJNA252803 & both down")}
  else if (y<x & x<0  )   { return("PRJNA529794<PRJNA252803 & both down")}
  else if (y<0 & 0<x & y<x)   { return("PRJNA529794<PRJNA252803 & PRJNA529794-down PRJNA252803-up")}
  else if (0<y & y<x)   { return("PRJNA529794<PRJNA252803 & both up")}
}

cate_class <- function(line){
  x <- as.numeric(line[1])
  y <- as.numeric(line[2])
  if      (y>x & x>0)         { return(1)}
  else if (y>x & y>0 & 0>x )   { return(2)}
  else if (0>y & y>x  )   { return(3)}
  else if (y<x & x<0  )   { return(4)}
  else if (y<0 & 0<x & y<x)   { return(5)}
  else if (0<y & y<x)   { return(6)}
}

fc_table$cate <- as.factor(apply(fc_table, 1, FUN=function(x2) cate_class(x2)))


p <- fc_table %>%
  ggplot( aes(PRJNA252803, PRJNA529794, text=group_num, color=cate)) +
  geom_point() +
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)

ggplotly(p)



write.csv(df_fc,paste0(output_path,"Expr_fc_output.csv"), row.names = FALSE)




