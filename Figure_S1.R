################################################################################
# FIGURE S1 ####################################################################
################################################################################

X = bm.sort %>% filter(bm.sort$Time_point==1)
X = distinct(X, `SubjectID`, .keep_all = TRUE)

shared.rows = intersect(X$`SubjectID`,pca.variables$SubjectID)
X = X %>% filter(X$`SubjectID` %in% shared.rows)
pca.variables = pca.variables %>% filter(pca.variables$SubjectID %in% shared.rows)
pca.variables = pca.variables[order(pca.variables$SubjectID),]

X = X[order(X$`SubjectID`),-c(1:5)] 

plot.data = log10(X+1) %>% na.omit() %>% t() 
plot.data = plot.data[rowSums(plot.data)>0,]

colors = c("#619CFF","#F8766D","#00BA38",'black','white')

pca = prcomp(t(plot.data),scale=T)
prop.var = summary(pca)$importance[2,]
summary(pca)
x_scores = pca$x[,1:4]

sc1=ggplot(x_scores, aes(PC1, PC2, color = as.factor(pca.variables$child_sex))) +
  geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
  geom_point(size=3,alpha=0.75)  + 
  labs(x=paste0('PC1',' (',round(prop.var[1]*100,1),'%)'),y=paste0('PC2',' (',round(prop.var[2]*100,1),'%)')) + 
  theme(text=element_text(size=12)) + labs(color='Sex')

sc2=ggplot(x_scores, aes(PC1, PC2, color = as.factor(pca.variables$birth.season))) +
  geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
  geom_point(size=3,alpha=0.75)  + 
  labs(x=paste0('PC1',' (',round(prop.var[1]*100,1),'%)'),y=paste0('PC2',' (',round(prop.var[2]*100,1),'%)')) + 
  theme(text=element_text(size=12)) + labs(color='Season')
sc3=ggplot(x_scores, aes(PC1, PC2, color = as.factor(pca.variables$birth.year))) +
  geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
  geom_point(size=3,alpha=0.75)  +
  labs(x=paste0('PC1',' (',round(prop.var[1]*100,1),'%)'),y=paste0('PC2',' (',round(prop.var[2]*100,1),'%)')) + 
  theme(text=element_text(size=12)) + labs(color='Year')
sc4=ggplot(x_scores, aes(PC1, PC2, color = pca.variables$prepreg_BMI)) +
  geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
  geom_point(size=3,alpha=0.75)  + scale_color_viridis_c(na.value='#FDD700',direction=1,option='magma') + 
  labs(x=paste0('PC1',' (',round(prop.var[1]*100,1),'%)'),y=paste0('PC2',' (',round(prop.var[2]*100,1),'%)')) + 
  theme(text=element_text(size=12)) + labs(color='BMI')
sc5=ggplot(x_scores, aes(PC1, PC2, color = as.factor(pca.variables$race_white))) +
  geom_hline(yintercept=0, color = "grey") + geom_vline(xintercept=0, color = "grey") + 
  geom_point(size=3,alpha=0.75)  +
  labs(x=paste0('PC1',' (',round(prop.var[1]*100,1),'%)'),y=paste0('PC2',' (',round(prop.var[2]*100,1),'%)')) + 
  theme(text=element_text(size=12)) + labs(color='Race')

Fig_S1 = plot_grid(sc1,sc3,sc4,sc5,sc2,ncol=2)
