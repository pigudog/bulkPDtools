library(RColorBrewer)
# heatmap
color = c("#B3242A","#4048A4","#1B64A1","#A61C5D","#9496C4","#929FD6","#484EAA","#7668BA","#4569BB","#EE3432","#F47368")
color = c("#2874C5","white","#FC4E07")
color = c("#00AFBB","white","#FC4E07")

red = c("#F4A582","#FC8D62","#F47368","#FC4E07","#EE3432","#B3242A","#b03d26")
blue = c("#b7e1e9","#9ccfe6","#8DA0CB","#00AFBB","#2874C5","#1B64A1","#4068B2","#4048A4","#223271")


color = c("#00AFBB","#FC4E07",'#99c355',"#e6b707","#2874C5","#f87669","#868686","#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#999999")
# UMAP AND TSNE
color = c("#51C4C2","#0D8A8C","#4583B3","#C63596","#BE86BA","#8B66B8","#F78E26","#F172AD","#F7AFB9","#4068B2","#512A93","#223271")
color = c("#f7ae55","#c4dfa2","#70c17f","#f9e9ab","#b7e1e9","#7ca9cc")
color = c("#b03d26","#005f81","#9ccfe6","#e0897e","#a5a7ab")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
mycolors <- getPalette(25)
mycolors = c("#A6CEE3", "#68A6CD", "#2A7FB7", "#569EA4",
             "#99CD91", "#8CCC6E", "#52AF43" ,"#5C9E42",
             "#B89B74", "#F88A89","#ED4F50", "#E4201F",
             "#F06C45", "#FBB86B" ,"#FDA440", "#FE870D", "#ED8F47",
             "#D5A7A9", "#B294C7", "#865FAB","#825D99", "#C7B699", "#F8F18F" ,"#D4A55B" ,"#B15928")

# macaron
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
           "#f1707d", 
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

###################################################################
# draw_pca()
library(bulkPDtools)
####################################################################
# Shows the effect of removing batch effects on data structures
data("DESeq2_raw")
data("combined_raw")
# dds is the result of DESeq2 on the 3 data combined with count matrix
library(DESeq2)
## varianceStabilizingTransformation (VST)
vsd <- vst(dds, blind=FALSE)
dat <- assay(vsd)
## define group
Group <- factor(coldata$type,levels=c('GSE138518','GSE155489','GSE193123'))
# the dat has been log2 with DESeq2
mycolors = c( "#2A7FB7", "#B3242A",'#99c355')

draw_pca(exp=dat,
         group_list = Group,
         addEllipses = TRUE,
         style = "default", # ggplot2 default 3D
         color.label = "Group",
         title = "")+
  scale_color_manual(values = mycolors[1:3])
pdf("Figure1/scatterplot.pdf",width = 6,height = 5)
pca.plot = draw_pca(dat,Group)+
  scale_color_manual(values = mycolors[1:3])
pca.plot
dev.off()

# remove batch effect
library(sva)
assay(vsd)=ComBat(assay(vsd), batchType, par.prior=TRUE)
dat <- assay(vsd)

pdf("Figure1/ellipse.pdf",width = 6,height = 5)
pca.plot = draw_pca(dat,Group)+
  scale_color_manual(values = mycolors[1:3])
pca.plot
dev.off()

######################################################################
# we use the ordered matrix for volcano
data("key_train_exprSet")
data("DESeq2_raw")
DEG$g=ifelse(DEG$pvalue>0.05,'stable',
             ifelse( DEG$log2FoldChange >1,'up',
                     ifelse( DEG$log2FoldChange < -1,'down','stable') )
)

DEG <- cbind(symbol = rownames(DEG), DEG)
# DEG <- rename(DEG, regulate = change) #修改列名 将"change"列名修改为"regulate"
DEG <- add_regulate(DEG)
color <- c("#00AFBB","#999999","#fc4e07")
color <- c("#2fa1dd", "#999999", "#f87669")
color <- c("#2A7FB7", "#999999", "#B3242A")
library(ggrepel)
p <- ggvolcano(data = DEG,x="log2FoldChange",y="pvalue",output = F,label = "symbol",
               fills = color,
               colors = color,
               x_lab = "log2FC",
               y_lab = "-log10P.Value",
               legend_position = "UR",
               log2FC_cut = 1, FDR_cut = 0.05)
p
pdf(file = "Figure1/volcano.pdf",width = 6,height = 5)
p
dev.off()

########################################################################
# we use the DEG to plot heatmap
library(pheatmap)
dat <- rt
### extract the name of deg
DEG_genes <- DEG[DEG$pvalue<0.05&abs(DEG$log2FoldChange)>1,]
DEG_sorted <- DEG_genes[order(DEG_genes$log2FoldChange,decreasing = T), ]
tem1 <- head(rownames(DEG_sorted),50)
tem2 <- tail(rownames(DEG_sorted),50)
# we need to order our genes
DEG_gene_expr <- dat[c(tem1,tem2),]
DEG_gene_expr <- t(scale(t(DEG_gene_expr)))
DEG_gene_expr <- na.omit(DEG_gene_expr)
ac=data.frame(g=anno$group)
rownames(ac)=colnames(DEG_gene_expr)
pheatmap(DEG_gene_expr,
         color=colorRampPalette(c("#2A7FB7","white","#E4201F"))(100),
         scale="row",
         border_color=NA,
         fontsize=10,
         show_rownames=F,
         show_colnames = T,
         annotation_col=ac,
         cluster_cols = F,
         filename = './Figure1/heatmap.pdf')

