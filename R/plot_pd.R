#' Title
#'
#' @param immune
#' @param rt
#'
#' @return
#' @export
#'
#' @examples
gene_cor_plot <- function(immune=immune,rt=rt) {
  pp <- corr.test(rt,immune,method="spearman",adjust = "fdr")
  cor <- pp$r
  pvalue <- pp$p
  heatmap <- melt(cor)
  colnames(heatmap)=c('sample','gene','cor')
  heatmap=mutate(heatmap,pvalue=melt(pvalue)[,3]) %>%
    mutate(signif = sapply(pvalue, function(x) myfun(x)))
  p = ggplot(heatmap,aes(sample,gene,col=cor))+
    geom_tile(color="grey70",fill="white",size=1)+
    geom_point(aes(size = abs(cor)),shape=15) +
    geom_text(aes(label=signif),size=6,color="black",
              hjust=0.5,vjust=0.7)+
    labs(x = NULL,y = NULL,color=NULL) +
    # scale_color_viridis_c(option = "inferno")+
    # scale_color_tableau()+
    # scale_color_continuous(type = "gradient")+
    # scale_color_gradient2()+
    scale_color_gradient2(low = "#2874C5",mid = "white",high = "#EE3432")+
    # scale_color_carto_c(palette = "BurgYl")+
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    theme(text=element_text(),
          axis.ticks.x = element_blank(),axis.text.x = element_text( angle =45),
          axis.ticks.y=element_blank(),
          panel.border = element_rect(fill=NA,color="grey70",
                                      size=2, linetype="solid")) +
    scale_size(range=c(1,10),guide=NULL)+
    guides(color = guide_colorbar(direction = "vertical",
                                  reverse = F,barwidth = unit(.5, "cm"),
                                  barheight = unit(15, "cm")))
  return(p)
}


#' Title
#'
#' @param pval
#'
#' @return
#' @export
#'
#' @examples
myfun <- function(pval) {
  stars = ""
  if(pval <= 0.001)
    stars = "***"
  if(pval > 0.001 & pval <= 0.01)
    stars = "**"
  if(pval > 0.01 & pval <= 0.05)
    stars = "*"
  if(pval > 0.05 & pval <= 0.1)
    stars = ""
  stars
}


#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
lolliplot_immune <- function(data){
  p.col = c('gold','pink','orange','LimeGreen','darkgreen')
  p.col = c('#00AFBB','#A6D854','gold','orange','#FC4E07')
  p.col = c("#F4A582","#FC8D62","#FC4E07","#EE3432","#B3242A")
  fcolor = function(x,p.col){
    color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],
                                                               ifelse(x>0.2,p.col[4], p.col[5])
    )))
    return(color)
  }

  p.cex = seq(2.5, 5.5, length=5)
  fcex = function(x){
    x=abs(x)
    cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],
                                                             ifelse(x<0.4,p.cex[4],p.cex[5]))))
    return(cex)
  }
  points.color = fcolor(x=data$pvalue,p.col=p.col)

  data$points.color = points.color
  points.cex = fcex(x=data$cor)
  data$points.cex = points.cex
  data=data[order(data$cor),]

  xlim = ceiling(max(abs(data$cor))*10)/10
  layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
  par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
  plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Correlation Coefficient",ylab="",yaxt="n",yaxs="i",axes=F)
  rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
  grid(ny=nrow(data),col="white",lty=1,lwd=2)
  #????ͼ?ε??߶?
  segments(x0=data$cor,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
  #????ͼ?ε?ԲȦ
  points(x=data$cor,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
  #??ͼ????չʾ????ϸ????????
  text(par('usr')[1],1:nrow(data),data$Cell,adj=1,xpd=T,cex=1.5)
  #չʾ?????Լ?????pvalue
  pvalue.text=ifelse(data$pvalue<0.001,'<0.001',sprintf("%.03f",data$pvalue))
  redcutoff_cor=0
  redcutoff_pvalue=0.05
  text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$cor)>redcutoff_cor & data$pvalue<redcutoff_pvalue,"red","black"),cex=1.5)
  axis(1,tick=F)

  #????ԲȦ??С??ͼ??
  par(mar=c(0,4,3,4))
  plot(1,type="n",axes=F,xlab="",ylab="")
  legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(cor)")

  #????ԲȦ??ɫ??ͼ??
  par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
  barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
  axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)

}


#' Title
#'
#' @param rt
#' @param group_list
#' @param Control
#' @param Treat
#' @param f
#'
#' @return
#' @export
#'
#' @examples
corplot_immune <- function(rt,group_list,
                          Control='CT',
                          Treat='PC'){
  con <- group_list==Control
  treat <- group_list==Treat
  conData=rt[con,]
  treatData=rt[treat,]
  conNum=nrow(conData)
  treatNum=nrow(treatData)
  data=t(rbind(conData,treatData))
  data=data[apply(data,1,sd)>0,]
  Type=c(rep("Control",conNum), rep("Treat",treatNum))
  names(Type)=colnames(data)
  Type=as.data.frame(Type)
  treatData=treatData[,apply(treatData,2,sd)>0]
  p = corrplot(corr=cor(treatData),
           method = "color",          #ͼ?ε?չʾ??ʽ
           order = "hclust",          #????ϸ??????????ʽ
           tl.col="black",            #??????ɫ
           number.cex = 0.8,          #????ϵ????????С
           addCoef.col = "black",     #????ϵ????????ɫ
           col=colorRampPalette(c("#2874C5","white","#B3242A"))(100),     #ͼ?ε???ɫ
  )
  return(p)
  }


#' Title
#'
#' @param rt
#' @param group_list
#' @param Control
#' @param Treat
#' @param f
#'
#' @return
#' @export
#'
#' @examples
violin_immune <- function(rt,group_list,
                           Control='CT',
                           Treat='PC',
                           f = "./Figure5/vioplot.pdf" ){
  con <- group_list==Control
  treat <- group_list==Treat
  conData=rt[con,]
  treatData=rt[treat,]
  conNum=nrow(conData)
  treatNum=nrow(treatData)
  rt=rbind(conData,treatData)
  outTab=data.frame()
  # pdf(file=f, width=13, height=8)
  par(las=1,mar=c(10,6,3,3))
  x=c(1:ncol(rt))
  y=c(1:ncol(rt))
  plot(x, y,
       xlim=c(0,(3*ncol(rt)-3)), ylim=c(min(rt), max(rt)+0.05),
       main="", xlab="", ylab="Fraction",
       pch=21,
       col="white",
       xaxt="n")
  for(i in 1:ncol(rt)){
    if(sd(rt[1:conNum,i])==0){
      rt[1,i]=0.00001
    }
    if(sd(rt[(conNum+1):(conNum+treatNum),i])==0){
      rt[(conNum+1),i]=0.00001
    }
    conData=rt[1:conNum,i]
    treatData=rt[(conNum+1):(conNum+treatNum),i]
    vioplot(conData,at=3*(i-1),lty=1,add = T,col = '#2874C5')
    vioplot(treatData,at=3*(i-1)+1,lty=1,add = T,col = '#ef5767')
    wilcoxTest=wilcox.test(conData,treatData)
    p=wilcoxTest$p.value
    if(p<0.05){
      cellPvalue=cbind(Cell=colnames(rt)[i],pvalue=p)
      outTab=rbind(outTab,cellPvalue)
    }
    mx=max(c(conData,treatData))
    lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
    #??С????ͼ???Ϸ????ϲ?????pvalue
    text(x=3*(i-1)+0.5, y=mx+0.02, labels=ifelse(p<0.001, paste0("p<0.001"), paste0("p=",sprintf("%.03f",p))), cex = 0.8)
  }
  legend("topright",
         c("Control", "Treat"),
         lwd=3,bty="n",cex=1,
         col=c("#2874C5","#ef5767"))
  text(seq(1,(3*ncol(rt)-2),3), -0.03, xpd=NA, labels=colnames(rt), cex=1, srt=45, pos=2)
  # dev.off()
  print("finished")
}







#' Title
#'
#' @param rt
#' @param group_list
#' @param Control
#' @param Treat
#'
#' @return
#' @export
#'
#' @examples
heatmap_immune <- function(rt,group_list,
                     Control='CT',
                     Treat='PC') {
  con <- group_list==Control
  treat <- group_list==Treat
  conData=rt[con,]
  treatData=rt[treat,]
  conNum=nrow(conData)
  treatNum=nrow(treatData)
  data=t(rbind(conData,treatData))
  data=data[apply(data,1,sd)>0,]
  Type=c(rep("Control",conNum), rep("Treat",treatNum))
  names(Type)=colnames(data)
  Type=as.data.frame(Type)
  p = pheatmap::pheatmap(data,
           annotation_col=Type,
           color=colorRampPalette(c(rep("#2874C5",3), "white", rep("#B3242A",3)))(50),
           cluster_cols=F,
           show_colnames=F,
           scale="row",
           fontsize = 7,
           fontsize_row=7,
           fontsize_col=7)
  return(p)
}

#' Title
#'
#' @param data
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
box_plot <- function(data,X,Y) {

  p <-   ggboxplot(data,
                   x=X, y=Y,
                   color = "group",
                   # fill = "group",
                   legend='none',ggtheme = theme_bw(),
                   # color = "black",
                   # palette = "lancet",
                   # palette =c("#2874C5","#B3242A"),
                   add = "jitter")+
    #添加p-valuep
    stat_compare_means(method = "t.test")+
    scale_color_manual(values = c("#2874C5","#B3242A"))
  # scale_fill_manual(values = c("#b7e1e9","#D5A7A9"))
  # stat_compare_means(method = "t.test")+
  # stat_compare_means(aes(label=paste0(..method..,"\n", "p=",..p.format..)), label.x = 1.5, label.y = 10)

  return(p)
}




#' Title
#'
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
ROC_plot <- function(X,Y) {

  p <-   plot.roc(Y,as.numeric(X),
                  print.auc=T,
                  legacy.axes=T,
                  print.thres=TRUE,
                  col="#B3242A")
  return(p)
}



# define function
#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
GO_plot <- function(object) {
  p <- dotplot(object, split="ONTOLOGY",showCategory = 5,label_format=50)+
    facet_grid(ONTOLOGY~.,scale="free")+
    theme(panel.grid = element_blank())+#修改主题
    theme(axis.title =element_text(size = 12, color = 'black'),
          axis.text.y =element_text(size = 12),
          legend.title=element_text(size=12))+
    # scale_color_gradient(high='#B3242A',low="#DFA0A6")+
    scale_fill_gradient(high='#B3242A',low="#f1707d")
  return(p)
}


# define function
#' Title
#'
#' @param object
#' @param showNum
#'
#' @return
#' @export
#'
#' @examples
KEGG_plot <- function(object,showNum=15) {
  p <- dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130)+
    theme(panel.grid = element_blank())+#修改主题
    theme(axis.title =element_text(size = 12, color = 'black'),
          axis.text.y =element_text(size = 12),
          legend.title=element_text(size=12))+
    scale_fill_gradient(high='#B3242A',low="#f1707d")
  return(p)
}
