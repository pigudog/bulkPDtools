# define function
sankeyplot <- function(corLodes) {

  mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),5)
  mycol <- rep(c("#00AFBB","#FC4E07",'#99c355',"#e6b707","#2874C5","#f87669","#92C5DE","#F4A582","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494"),5)
  p = ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +
    #??aes.flow??????????ɫ??forward˵????ɫ??ǰ??һ?£?backward˵???ͺ???һ?¡?
    geom_flow(width = 1/8,aes.flow = "forward") +
    geom_stratum(alpha = .9,width = 0.2) +
    scale_fill_manual(values = mycol) +
    #size = 2.4???????????ֵĴ?С
    geom_text(stat = "stratum", size = 2.4,color="black") +
    xlab("") + ylab("") + theme_bw() +
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
    theme(panel.grid =element_blank()) +
    theme(panel.border = element_blank()) +
    ggtitle("") + guides(fill = FALSE)
  return(p)
}

