library(VennDiagram)
library(ggplot2)
library(plotly)

name_list_file = "input_dataset_name.txt"
input_path = "./"
output_path = "./"
treatcode_path = "treat_code.txt"
projectname = "Project_name_by_species.txt"


# load data ---------------------------------------------------------------

files <- scan(name_list_file, character(), quote = "")
load(paste0(input_path,"fc_table.RData"))

# functions ---------------------------------------------------------------
cate_class <- function(line){
  Categories <- c("PRJNA529794>PRJNA252803 & both up",
                  "PRJNA529794>PRJNA252803 & PRJNA529794-up PRJNA252803-down",
                  "PRJNA529794>PRJNA252803 & both down",
                  "PRJNA529794<PRJNA252803 & both down",
                  "PRJNA529794<PRJNA252803 & PRJNA529794-down PRJNA252803-up",
                  "PRJNA529794<PRJNA252803 & both up",
                  "same")
  x <- as.numeric(line[1])
  y <- as.numeric(line[2])
  if      (y>x & x>0)         {Categories[1]}
  else if (y>x & y>0 & 0>x )  {Categories[2]}
  else if (0>y & y>x  )       {Categories[3]}
  else if (y<x & x<0  )       {Categories[4]}
  else if (y<0 & 0<x & y<x)   {Categories[5]}
  else if (0<y & y<x)         {Categories[6]}
  else                        {Categories[7]}
}


# interactive scatter plot ------------------------------------------------
fc_table$cate <- unlist(apply(fc_table, 1, FUN=function(x2) cate_class(x2)))
p <- fc_table %>%
  ggplot( aes(PRJNA252803_log2fc, PRJNA529794_log2fc, text=group_num, color=cate)) +
  geom_point() +
  theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)

ggplotly(p)



# venn diagram ------------------------------------------------------------
# Generate 3 sets of 200 words
set1_up <- fc_table$group_num[which(fc_table[,1]>0)]
set1_down <- fc_table$group_num[which(fc_table[,1]<0)]
set2_up <- fc_table$group_num[which(fc_table[,2]>0)]
set2_down <- fc_table$group_num[which(fc_table[,2]<0)]

# |log2FC| > 1
set1_deg <- fc_table$group_num[which(abs(fc_table[,4])>1)]
set2_deg <- fc_table$group_num[which(abs(fc_table[,5])>1)]


# Chart
venn.diagram(
  x = list(set1_up,set2_up),
  category.names = paste0(files,rep(c("_upregulate"),2)),
  filename = paste0(output_path,'venn_up.png'),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500 , 
  width = 1500 , 
  resolution = 300,
  compression = "lzw",
  margin = .1,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("slateblue2","plum2"),
  
  # Numbers
  cex = 1,
  fontface = "bold",
  label.col = "black",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

venn.diagram(
  x = list(set1_down,set2_down),
  category.names = paste0(files,rep(c("_downregulate"),2)),
  filename = paste0(output_path,'venn_down.png'),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500 , 
  width = 1500 , 
  resolution = 300,
  compression = "lzw",
  margin = .2,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("deepskyblue","dodgerblue3"),
  
  # Numbers
  cex = 1,
  fontface = "bold",
  label.col = "black",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

venn.diagram(
  x = list(set1_deg, set2_deg),
  category.names = paste0(files,rep(c("_deg"),2)),
  filename = paste0(output_path,'venn_deg.png'),
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 1500 , 
  width = 1500 , 
  resolution = 300,
  compression = "lzw",
  margin = .1,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("tomato4","orange"),
  
  # Numbers
  cex = 1,
  fontface = "bold",
  label.col = "black",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer"
)


