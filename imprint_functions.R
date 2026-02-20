# Functions for IMPRINT analysis

# Boxplot function with statistical analysis
myboxplot = function(data,xcolumn,ycolumn,palette,mytitle,myylabel){
  
  myboxplot = ggplot(data, aes(x = xcolumn, y = ycolumn, group = (xcolumn), color = xcolumn)) +
    geom_boxplot(width=0.8,outlier.shape=NA) + 
    geom_jitter(alpha = 0.25, size = 2) + 
    scale_color_manual(xcolumn, breaks=c("1", "2", "3"),values=c("1"=palette[1],"2"=palette[2],"3"=palette[3])) +
    labs(x='Cluster',y=myylabel,title=mytitle) +
    theme(legend.position="none",axis.title=element_text(size=12),title=element_text(size=12))
  
  kw = kruskal.test(ycolumn ~ xcolumn, data = data) 
  kw.pval = kw$p.value
  
  # if (kw.pval<0.05){
  my_comparisons <- list( c("1", "2") )
  
  myboxplot +
    stat_compare_means(comparisons = my_comparisons,label='p.signif',size=3) #+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = max(ycolumn)*1.1)
  # stat_compare_means(comparisons = my_comparisons,size=3,method="wilcox.test") #+ # Add pairwise comparisons p-value
  # }
}

# Stacked bar plot function with Fisher's exact test
dem.bar.plot = function(data,xcolumn,ycolumn,mytitle,myylabel){
  
  data = data[!is.na(ycolumn),]
  xcolumn = xcolumn[!is.na(ycolumn)]
  ycolumn = ycolumn[!is.na(ycolumn)]
  
  stat.summary = summary(table(ycolumn,xcolumn))
  pval = stat.summary$p.value
  
  mycolors = colorRampPalette(brewer.pal(8,'Accent'))(length(levels(ycolumn)))
  # mycolors = brewer.pal(length(levels(ycolumn)),'Accent')
  
  mybarplot = ggplot(data, aes(x = xcolumn, fill=ycolumn)) +
    geom_bar(position='fill') + 
    annotate('text',x=1.5,y=1.1,label=paste0('Chi-Squared p=',as.character(signif(pval,2))),size=3) +
    labs(x='Cluster',y=myylabel,title=mytitle) + 
    scale_fill_manual(values=mycolors) +
    theme(legend.text = element_text(size=10),axis.title=element_text(size=8),title=element_text(size=8),
          legend.key.size = unit(0.5, 'cm'))
  
  mybarplot
}


mypca = function(data,firstPC,secondPC,colors,colorBy,PCx,PCy){
  pca = prcomp((data),scale=F)
  prop.var = summary(pca)$importance[2,]
  summary(pca)
  x_scores = pca$x[,1:4]
  featnames = colnames(X)
  
  sc=ggplot(x_scores, aes(.data[[firstPC]], .data[[secondPC]], color = colorBy, col=colors)) +
    geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
    geom_point(size=3,alpha=0.75)  + scale_color_manual(colorBy, breaks=c("1","2","3"),values=c("1"=colors[1],"2"=colors[2],"3"=colors[3])) + 
    labs(x=paste0(firstPC,' (',round(prop.var[PCx]*100,1),'%)'),y=paste0(secondPC,' (',round(prop.var[PCy]*100,1),'%)')) + 
    theme(text=element_text(size=20)) + theme(legend.position = 'none' , legend.title=element_blank()) +
    stat_ellipse(level=0.95)
  
  prin_comp = prcomp((data), rank = 4)
  components = prin_comp[["x"]]
  components = data.frame(components)
  
  explained_variance = summary(prin_comp)[["sdev"]]
  explained_variance = explained_variance[1:4]
  comp = prin_comp[["rotation"]]
  
  loadings = comp
  for (i in seq(explained_variance)){
    loadings[,i] <- comp[,i] * explained_variance[i]
  }
  loadings = as.data.frame(loadings)
  loadings$Feat = featnames
  
  norm_vec <- function(x) sqrt(sum(x^2))
  mag = apply(loadings[,1:2],1,norm_vec)
  loadings$mag = mag
  loadings =  loadings[order(loadings$mag,decreasing=TRUE),]
  ld = ggplot(loadings, aes(.data[[firstPC]], .data[[secondPC]], label=rownames(loadings))) +
    geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
    geom_segment(aes(x = 0, y = 0, xend = .data[[firstPC]], yend = .data[[secondPC]]), color = "black") + geom_point(size=3, color = 'black') + 
    geom_text(hjust=0.1, vjust=-1, size=4) +
    labs(x=paste0(firstPC,' (',round(prop.var[PCx]*100,1),'%)'),y=paste0(secondPC,' (',round(prop.var[PCy]*100,1),'%)')) +
    theme(text=element_text(size=20)) + xlim(-10,10) + xlim(min(loadings$PC1)*1.1, max(loadings$PC1)*1.1) + ylim(min(loadings$PC2)*1.1, max(loadings$PC2)*1.1)
  plot_grid(sc,ld)
}

