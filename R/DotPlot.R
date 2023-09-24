


# define function
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
KEGG_plot <- function(object,showNum=15) {
  p <- dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130)+
    theme(panel.grid = element_blank())+#修改主题
    theme(axis.title =element_text(size = 12, color = 'black'),
          axis.text.y =element_text(size = 12),
          legend.title=element_text(size=12))+
    scale_fill_gradient(high='#B3242A',low="#f1707d")
  return(p)
}
