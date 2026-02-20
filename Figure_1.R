################################################################################
# FIGURE 1B ####################################################################
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
pheatmap(plot.data,color=myColor,cluster_rows=F,cluster_cols=F,
         breaks=myBreaks,gaps_col=c(seq(from=21,to=ncol(plot.data),by=21),
                                    ncol(plot.data)-3),annotation_row=row.labels, annotation_colors = ann_colors, 
         gaps_row = c(sum(bm$`Delta_weeks_bin`==1),(sum(bm$`Delta_weeks_bin`==1) + sum(bm$`Delta_weeks_bin`==2))), 
         fontsize_col=4, show_rownames = FALSE)

################################################################################
# FIGURE 1C ####################################################################
################################################################################

#populate!