data("immune")
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)
library("IOBR")
rt = cbind(ID = rownames(rt),rt)
res_cibersort <- cell_bar_plot(rt)
## There are seven categories you can choose: box, continue2, continue, random, heatmap, heatmap3, tidyheatmap
## >>>>=== Palette option for random: 1: palette1; 2: palette2; 3: palette3;  4: palette4
cibersort_long <- res_cibersort$data
macaron = c("#A0C2E7",
            "#6894B9",
            "#8798A6",
            "#E0D9E0",
            "#EDBAA7",
            "#FADB7F",
            "#F3B646",
            "#EF9749",
            "#B27466",
            "#646F3F",
            "#899678",
            "#C2BC9A",
            "#868A63",
            "#C4C3BE",
            "#DFA0A6",
            "#98B3D9",
            "#E4BE92",
            "#CB6B7A",
            "#D5CBDA",
            "#f1707d", # 马卡龙草莓奶霜
            "#f15536",
            "#ef5767",
            "#ae716e",
            "#cb8e85",
            "#cf8878",
            "#c86f67",
            "#f1ccb8",
            "#f2debd",
            "#b8d38f",
            "#ddff95",
            "#ff9b6a",
            "#f1b8f1",
            "#d9b8f1",
            "#f1ccb8",
            "#f1f1b8",
            "#b8f1ed",
            "#e7dbca",
            "#e26538",
            "#f3d751",
            "#fd803a",
            "#fe997b",
            "#c490a0"
)
p1 <- cibersort_long %>%
  ggplot(aes(ID,fraction))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = macaron,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
p1

# 有顺序的箱线图
library(forcats)

p2 <- ggplot(cibersort_long,aes(fct_reorder(cell_type, fraction),fraction,fill = cell_type)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2,aes(color=cell_type))+
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = palette4)
p2


data("immune")
# remove the colomn without expression
# violin
rt <- rt[,c(1:6,11:14,16:22)]

## heatmap
pdf(file="Figure5/heatmap.pdf", width=7, height=4)
heatmap_immune(rt=rt,
              group_list = group_list,
              Control = "CT",
              Treat = "PC")
dev.off()

library(bulkPDtools)
## violin plot
library(vioplot)
library(stringr)
colnames(rt) = str_replace(colnames(rt),pattern = "_CIBERSORT",replacement = "")
pdf(file="./Figure5/vioplot.pdf", width=13, height=8)
violin_immune(rt=rt,
               group_list = group_list,
               Control = "CT",
               Treat = "PC")
dev.off()

## corplot
library(corrplot)
pdf("test/Figure5/corHeatmap.pdf", width=12, height=12)
corplot_immune(rt=rt,
              group_list = group_list,
              Control = "CT",
              Treat = "PC")
dev.off()
## lolliplot
# inputFile="data/LPIN1cor.immune.result.txt"
# data = read.table(inputFile, header=T, sep="\t", check.names=F)
# save(data,file = "data/LPIN1.rda")
data("LPIN1")
data$Cell = str_replace(data$Cell,pattern = "_CIBERSORT",replacement = "")
pdf(file="./Figure5/LPIN1_Lollipop.pdf", width=9.5, height=7)
lolliplot_immune(data = data)
dev.off()

# inputFile="data/NR4A1cor.immune.result.txt"
# data = read.table(inputFile, header=T, sep="\t", check.names=F)
# save(data,file = "data/NR4A1.rda")
data("NR4A1")
data$Cell = str_replace(data$Cell,pattern = "_CIBERSORT",replacement = "")
pdf(file="./Figure5/NR4A1_Lollipop.pdf", width=9.5, height=7)
lolliplot_immune(data = data)
dev.off()


data("gene_cor_heatmap")

library(psych)
library(ggplot2)
library(dplyr)
colnames(immune) = str_replace(colnames(immune),pattern = "_CIBERSORT",replacement = "")
pdf("Figure5/hubgeneheatmap.pdf",width = 7,height = 7)
gene_cor_plot(immune = immune,rt = rt)
dev.off()
