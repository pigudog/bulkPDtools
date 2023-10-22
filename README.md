# bulkPDtools
 a package for visualization of bulk RNA transcriptome analysis
 - test code in `./test`

# download
```r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("pigudog/bulkPDtools")
```
# 1. check data
## draw pca plot
```r
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
pca.plot = draw_pca(dat,Group)+
  scale_color_manual(values = mycolors[1:3])
pca.plot
```
![](README/Pasted%20image%2020231022200857.png)

```r
# remove batch effect
library(sva)
assay(vsd)=ComBat(assay(vsd), batchType, par.prior=TRUE)
dat <- assay(vsd)

pca.plot = draw_pca(dat,Group)+
  scale_color_manual(values = mycolors[1:3])
pca.plot
```
![](README/Pasted%20image%2020231022200913.png)


## volcano plot
```r
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

```
![](README/Pasted%20image%2020231022200928.png)
## heatmap
```r
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
```
![](README/Pasted%20image%2020231022200941.png)


# 2. GO,KEGG AND GSEA
![](README/GSEA_0030595.pdf)
# 3. ML
![](README/Pasted%20image%2020231022200952.png)
# 4. immune
![](README/Pasted%20image%2020231022201003.png)


# 5. sankey
![](README/Pasted%20image%2020231022201013.png)
