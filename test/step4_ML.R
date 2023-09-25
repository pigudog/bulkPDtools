data("ML")
library(bulkPDtools)
library(pROC)
for (i in intersect_gene) {
  pdf(paste0("Figure4/train",i,'.pdf'),width = 4,height = 4)
  ROC_plot(X=as.numeric(X_train[,i]),
           Y=Y_train)
  dev.off()
}

for (i in intersect_gene) {
  pdf(paste0("Figure4/test",i,'.pdf'),width = 4,height = 4)
  ROC_plot(X=as.numeric(X_test[,i]),
           Y=Y_test)
  dev.off()
}


library(ggplot2)
library(ggpubr)
gene_test = cbind(data.frame(X_test),group=group_list_test)
box_plot(gene_test,X = "group",Y = 'LPIN1')
box_plot(gene_test,X = "group",Y = 'ACSS2')

# t.test is permitted in lognormalized data
gene_test = cbind(data.frame(X_test),group=group_list_test)
for (i in intersect_gene) {
  pdf(paste0("./test/Figure4/test_expr_",i,'.pdf'),width = 5,height = 5)
  p = box_plot(gene_test,X = "group",Y = i)
  print(p)
  dev.off()
}


gene_train = cbind(data.frame(X_train),group=group_list_train)
for (i in intersect_gene) {
  pdf(paste0("test/Figure4/train_expr_",i,'.pdf'),width = 5,height = 5)
  p <-box_plot(gene_train,X = "group",Y = i)
  print(p)
  dev.off()
}

# combined with xgboost
data("xgboost")
pROC::plot.roc(as.factor(Y_train),,print.auc=T)
pROC::plot.roc(as.factor(Y_test),XGB_test_Predictions,print.auc=T)
pdf(paste0("test/Figure4/train_combined.pdf"),width = 4,height = 4)
ROC_plot(X = XGB_train_Predictions,Y = as.factor(Y_train))
dev.off()
pdf(paste0("test/Figure4/test_combined.pdf"),width = 4,height = 4)
ROC_plot(X = XGB_test_Predictions,Y = as.factor(Y_test))
dev.off()
