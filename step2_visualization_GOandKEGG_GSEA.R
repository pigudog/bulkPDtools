library(bulkPDtools)
library(clusterProfiler)
library(ggplot2)
source("./DotPlot.R")
# set the parameters
pvalueFilter=0.05
qvalueFilter=0.05
# GO UP
data("GO_up")
library(ggplot2)
# bubble plot

pdf(file="Figure2/bubble_go_up.pdf", width=9, height=7)
GO_plot(kk)
dev.off()

# GO DOWN
data("GO_down")
GO_plot(kk)
dev.off()


# KEGG UP
data("KEGG_up")
pdf(file="Figure2/bubble_kegg_up.pdf", width=9, height=7)
KEGG_plot(kk,showNum = 15)
dev.off()

# KEGG down
data("KEGG_down")
pdf(file="Figure2/bubble_kegg_down.pdf", width=9, height=7)
KEGG_plot(kk,showNum = 15)
dev.off()

###############################################
# GSEA
data("GSEA_result")
GO_result <- data.frame(GO)
KEGG_result <- data.frame(KEGG)
head(GO_result[,1:3],n = 20)
head(KEGG_result[,1:3],n = 20)

dotplotGsea(data = GO,topn = 10,
            order.by = 'NES',
            add.seg = T,
            scales = 'free')

p = gseaNb(object = KEGG,
       geneSetID = "hsa04060",
       subPlot = 2,
       addGene = T,
       kegg = T,
       markTopgene = T)
p

ID <- c("GO:0030595","GO:0097530","GO:0060326","GO:0097529","GO:0050727","GO:0070665")
for (i in ID){

  p = gseaNb(object = GO,
         geneSetID = i,
         subPlot = 2,
         addGene = T,
         kegg = T,
         addPval = T,
         markTopgene = T)
  ggsave(p,filename = paste0("Figure2/GSEA_",str_split(i,pattern = ":")[[1]][2],".pdf"),width = 5,height = 4)
}

ID <- c("hsa04060","hsa04062","hsa04668","hsa04657","hsa04064","hsa05417")
for (i in ID){

  p = gseaNb(object = KEGG,
             geneSetID = i,
             subPlot = 2,
             addGene = T,
             kegg = T,
             addPval = T,
             markTopgene = T)
  ggsave(p,filename = paste0("Figure2/GSEA_",i,".pdf"),width = 5,height = 4)
}
