################################################################################
# FIGURE 1A ####################################################################
################################################################################
# set up heat map input data
plot.data= log10(bm.sort[,feature.range]) %>% as.data.frame()
plot.data = scale(plot.data,center = TRUE, scale = TRUE)
plot.data[is.na(plot.data)]=0

# set up row annotations
row.labels = bm.sort$Vax_timing_bin
row.labels=recode(row.labels, '1'='No Vaccine',
                  '2'='Before Pregnancy',
                  '3'='First Trimester',
                  '4'='Second Trimester',
                  '5'='Third Trimester',
                  '6'='While Breastfeeding')
row.labels = row.labels %>% as.data.frame()
rownames(row.labels) = paste0("row_", seq(nrow(row.labels)))
rownames(plot.data)=rownames(row.labels)
colnames(row.labels)=c('Preg_vax')

# set up heat map parameters
ann_colors =brewer.pal(6, "BrBG")
ann_colors = list(Preg_vax = c('No Vaccine'=ann_colors[1],
                               'Before Pregnancy'=ann_colors[2],
                               'First Trimester'=ann_colors[3],
                               'Second Trimester'=ann_colors[4],
                               'Third Trimester'=ann_colors[5],
                               'While Breastfeeding'=ann_colors[6]))

paletteLength=100
myBreaks = c(seq(min(plot.data),-2,length.out= 1),seq(-2, 0, length.out=ceiling(paletteLength/2) + 1), 
             seq(1/paletteLength, 2-1/paletteLength, length.out=floor(paletteLength/2)),seq(2,max(plot.data),length.out= 2))
myColor = colorRampPalette(c("blue4", "white", "firebrick"))(paletteLength+2)

#plot heat map
Fig_1A = pheatmap(plot.data,color=myColor,cluster_rows=F,cluster_cols=F,
         breaks=myBreaks,gaps_col=c(seq(from=21,to=ncol(plot.data),by=21),
                                    ncol(plot.data)-3),annotation_row=row.labels, annotation_colors = ann_colors, 
         gaps_row = c(sum(bm$`Time_point`==1),(sum(bm$`Time_point`==1) + sum(bm$`Time_point`==2))), 
         fontsize_col=4, show_rownames = FALSE)

################################################################################
# FIGURE 1B ####################################################################
################################################################################


plot.data = transfer.ratios
plot.data[is.na(plot.data)]=0
plot.data=plot.data[, order(colnames(plot.data))]
plot.data = plot.data %>% select(-contains('EBV'))
plot.data = plot.data %>% select(-contains('RSV'))
plot.data = plot.data %>% as.matrix() %>% apply(2, mean) %>% t() %>% as.data.frame()

plot.data = plot.data %>% pivot_longer(.,cols=colnames(plot.data),names_to='Ag',values_to='MFI')
plot.data$Fc = sub("\\_.*", "", plot.data$Ag)
plot.data$strain = sub("\\_.*", "", plot.data$Ag)
plot.data$strain = sapply(strsplit(plot.data$Ag, "_"), function(x) x[2])    
plot.data$strain[plot.data$strain=='HongKong'|plot.data$strain=='Michigan']='NA'
plot.data = plot.data %>% arrange(match(Fc, c("IgM", "IgG1", "IgG2","IgG3",'IgG4','IgA1','IgA2','FcgR2A','FcgR2B','FcgR3A','FcgR3B','ADCD')))
brewer.pal(3,"Set2")
# myPalette= c(rep('#FC8D62',21),rep('#8DA0CB',12),rep('#66C2A5',3))
myPalette= c(rep('#FC8D59',21),rep('#FFFFBF',12),rep('#91BFDB',3))
myPalette= c(rep('#DEEBF7',21),rep('#9ECAE1',12),rep('#3182BD',3))
myPalette= c(rep('#E0ECF4',21),rep('#9EBCDA',12),rep('#8856A7',3))

Fig_1B = plot.data %>% mutate(Fc = fct_relevel(Fc, "IgM", "IgG1", "IgG2","IgG3",'IgG4','IgA1','IgA2','FcgR2A','FcgR2B','FcgR3A','FcgR3B','ADCD')) %>%
  ggplot(aes(x=Fc,y=MFI,fill = Fc,shape=strain)) + geom_boxplot(width=0.8,outlier.shape=NA,fill=myPalette) +
  geom_point(alpha=0.5,size=5,position = position_jitterdodge(dodge.width=0.8))  + theme_classic() +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20),
        legend.text = element_text(size=20),legend.title = element_text(size=20), panel.background = element_blank()) + 
  labs(x='Fc',y='MFI (log10)')

Fig_1B
