################################################################################
# FIGURE 3A ####################################################################
################################################################################
my.pvals = matrix(NA, nrow= 2, ncol = length(feature.range))
p.adj = matrix(NA, nrow= 2, ncol = length(feature.range))
cohens.d = matrix(NA, nrow= 2, ncol = length(feature.range))

for (i in 1:2){
  for (j in feature.range){
    group.1 = as.matrix(bm.sort[(bm.sort$Vax_timing_bin == 3 | bm.sort$Vax_timing_bin == 4 | bm.sort$Vax_timing_bin == 5) &  bm.sort$Time_point==i,j])
    group.2 = as.matrix(bm.sort[bm.sort$Vax_timing_bin == 1 & bm.sort$Time_point==i,j])
    
    stat.result = wilcox.test(group.1,group.2)
    my.pvals[i,j-5] = stat.result$p.value
    
    effect_size = cohen.d(log10(group.1+1),log10(group.2+1), hedges.correction = T)
    cohens.d[i,j-5] = effect_size[["estimate"]]
    
  }
}

rownames(cohens.d)=c('2wk','6wk')
colnames(cohens.d)=colnames(bm.sort)[feature.range]
rownames(my.pvals)=c('2wk','6wk')
colnames(my.pvals)=colnames(bm.sort)[feature.range]
cohens.d = cohens.d %>% as.data.frame()
cohens.d[is.na(cohens.d)] = 0

# set up heat map parameters
paletteLength=100
myColor = colorRampPalette(c("blue4","white","firebrick"))(paletteLength)

myBreaks = c(seq(min(cohens.d),-1,length.out= 1),seq(-1+1/paletteLength, 0, length.out=ceiling(paletteLength/2)-2), seq(0+1/paletteLength, 1, length.out=floor(paletteLength/2)-1), seq(1+1/paletteLength, max(cohens.d), length.out=2))

# plot heatmap
Fig_3A = pheatmap(cohens.d,color=myColor,cluster_rows=F,cluster_cols=F,gaps_col=c(seq(from=21,to=ncol(cohens.d),by=21),ncol(cohens.d)-3),breaks=myBreaks,fontsize_col=4,  number_color='black',fontsize_number=14)

################################################################################
# FIGURE 3B-E ##################################################################
################################################################################

# Uncomment the antigen for which you want to visualize box plots.

myAntigen = 'B_Colorado'

# myAntigen = 'A_Victoria_2570_2019' #:D
myAntigen = 'H7'
# H3N2s
myAntigen = 'HongKong_4801'
myAntigen = 'Singapore'
myAntigen = 'Kansas_14'

myAntigen = 'Cambodia'
myAntigen = 'Darwin'
myAntigen = 'California'

myAntigen = 'Victoria'
myAntigen = 'HongKong_45'
myAntigen = 'B_Phuket_3073_2013'
myAntigen = 'B_Austria_1359417_2021'

rm(vax.t1,vax.t2,novax.t1,novax.t2)

vax.t1 = bm[(bm$Vax_timing_bin == 3 | bm$Vax_timing_bin == 4 | bm$Vax_timing_bin == 5) & bm$Time_point==1,] %>% dplyr::select(contains(myAntigen))
vax.t1 = log10(vax.t1+1)
vax.t1 = pivot_longer(vax.t1,cols=colnames(vax.t1),names_to='Ag',values_to='MFI')
vax.t1$Fc = sub("\\_.*", "", vax.t1$Ag)
vax.t1$group = 'Vax_2wk'

vax.t2 = bm[(bm$Vax_timing_bin == 3 | bm$Vax_timing_bin == 4 | bm$Vax_timing_bin == 5) & bm$Time_point==2,] %>% dplyr::select(contains(myAntigen))
vax.t2 = log10(vax.t2+1)
vax.t2 = pivot_longer(vax.t2,cols=colnames(vax.t2),names_to='Ag',values_to='MFI')
vax.t2$Fc = sub("\\_.*", "", vax.t2$Ag)
vax.t2$group = 'Vax_6wk'

novax.t1 = bm[(bm$Vax_timing_bin == 1) & bm$Time_point==1,] %>% dplyr::select(contains(myAntigen))
novax.t1 = log10(novax.t1+1)
novax.t1 = pivot_longer(novax.t1,cols=colnames(novax.t1),names_to='Ag',values_to='MFI')
novax.t1$Fc = sub("\\_.*", "", novax.t1$Ag)
novax.t1$group = 'No_Vax_2wk'

novax.t2 = bm[(bm$Vax_timing_bin == 1) & bm$Time_point==2,] %>% dplyr::select(contains(myAntigen))
novax.t2 = log10(novax.t2+1)
novax.t2 = pivot_longer(novax.t2,cols=colnames(novax.t2),names_to='Ag',values_to='MFI')
novax.t2$Fc = sub("\\_.*", "", novax.t2$Ag)
novax.t2$group = 'No_Vax_6wk'

plot.data = rbind(vax.t1,vax.t2,novax.t1,novax.t2)
plot.data = plot.data %>% filter(plot.data$Fc != 'ADCP' & plot.data$Fc !='ADNP')
my.comparisons = list( c("No_Vax_2wk","Vax_2wk"), c("No_Vax_6wk","Vax_6wk"))
plot.data$group = factor(plot.data$group, levels=c('Vax_2wk','No_Vax_2wk','Vax_6wk','No_Vax_6wk'))

Fig_3BE = ggplot(plot.data, aes(x = group, y =MFI, color = as.factor(group), group = as.factor(group))) + facet_wrap(. ~ Fc, nrow = 2, scales = "free") + 
  geom_boxplot(color='black',size = 0.5, outlier.shape=NA) + 
  geom_jitter(alpha=0.25,width=0.1) +
  stat_compare_means(comparisons = my.comparisons,label='p.signif',size=3,method='wilcox.test') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.75), legend.position = "none") + 
  labs(x = "",y = "MFI (log10)")

Fig_3BE
