library(ggalluvial)
library(ggplot2)
library(dplyr)
library(scPDtools)
# rt1 <- read.table("data/dgidb_export_2023-09-01.tsv",sep = "\t",header = T)
# rt1 <- rt1[,c(5,4)]
# rt2=read.table("data/pathway.csv",sep = ",",header = T)
# save(rt1,rt2,file="data/sankey.rda")
data("sankey")
# 合并两个数据框
merged_data <- merge(rt1, rt2, all = TRUE)
merged_data
newData=merged_data[,c(2,1,3)]
corLodes=to_lodes_form(newData, axes = 1:3, id = "Cohort")
pdf(file="Figure7/ggalluvial_pcos.pdf",width=12,height=6)


# define function
sankeyplot(corLodes)
dev.off()
