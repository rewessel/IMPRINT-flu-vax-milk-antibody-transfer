## SET UP CODE FOR IMPRINT ANALYSIS

## LOAD REQUIRED PACKAGES
library(ggplot2)
library(pheatmap)
library(readxl)
library(caret)
library(tidyverse)
library(writexl)
library(RColorBrewer)
library(cowplot)
library(mltools)
library(data.table)
library(Hmisc)
library(ggpubr)

# SET WORKING DIRECTORY
getwd()
setwd('Z:/users/wesselr/data/IMPRINT') #fix this typo later
getwd()

## IMPORT DATA
bm = read_xlsx('IMPRINT_breastmilk_combined_features_01142026.xlsx')

# adjust the variable names automatically here so they show up in all plots
# All A-strains should have a HhNn label attached to them
mynames = colnames(bm)
mynames[grep("A_Brisbane", mynames)] = paste0(mynames[grep("A_Brisbane", mynames)],'_H1N1')
mynames[grep("A_California_04_2009", mynames)] = paste0(mynames[grep("A_California_04_2009", mynames)],'_H1N1')
mynames[grep("A_Cambodia_e0826360", mynames)] = paste0(mynames[grep("A_Cambodia_e0826360", mynames)],'_H3N2')
mynames[grep("A_Darwin", mynames)] = paste0(mynames[grep("A_Darwin", mynames)],'_H3N2')
mynames[grep("A_Hawaii", mynames)] = paste0(mynames[grep("A_Hawaii", mynames)],'_H1N1')
mynames[grep("A_HongKong_4801", mynames)] = paste0(mynames[grep("A_HongKong_4801", mynames)],'_H3N2')
mynames[grep("A_HongKong_45", mynames)] = paste0(mynames[grep("A_HongKong_45", mynames)],'_H3N2')
mynames[grep("A_Indonesia_05", mynames)] = paste0(mynames[grep("A_Indonesia_05", mynames)],'_H5N1')
mynames[grep("A_Kansas_14_2017", mynames)] = paste0(mynames[grep("A_Kansas_14_2017", mynames)],'_H3N2')
mynames[grep("A_Michigan_45", mynames)] = paste0(mynames[grep("A_Michigan_45", mynames)],'_H1N1')
mynames[grep("A_Shanghai_1_2013", mynames)] = paste0(mynames[grep("A_Shanghai_1_2013", mynames)],'_H7N9')
mynames[grep("A_Singapore_INF", mynames)] = paste0(mynames[grep("A_Singapore_INF", mynames)],'_H3N2')
mynames[grep("A_Victoria_2570", mynames)] = paste0(mynames[grep("A_Victoria_2570", mynames)],'_H1N1')
mynames[grep("A_Wisconsin_588", mynames)] = paste0(mynames[grep("A_Wisconsin_588", mynames)],'_H1N1')
mynames[grep("A_Wisconsin_67", mynames)] = paste0(mynames[grep("A_Wisconsin_67", mynames)],'_H1N1')
colnames(bm) = mynames

# Rearrange the columns of the breast milk data set per Tal's suggestions:
# A strains (H1N1, H3N2, H5N1, H7N9) --> B strains --> NA antigens
# Then within those categories, organize by year... how to do this? extract the year and put them in sequential order

fc.types = sub("^(([^_]*_){1}).*", "\\1", mynames[feature.range = 4:282]) %>% unique()
  
mynames.full = matrix(data=NA)
# mynames.fc = list()
for (i in 1:length(fc.types)){ #length(fc.types)
  my.names.current = mynames[grep(fc.types[i],mynames)]

  mynames.A = my.names.current[grep('_A_',my.names.current)]
  mynames.A = mynames.A[order(as.numeric(str_match(mynames.A, "20\\s*(.*?)\\s*_")[, 2]))]
  mynames.A = c(mynames.A[grep('H1N1',mynames.A)],mynames.A[grep('H3N2',mynames.A)],mynames.A[grep('H5N1',mynames.A)],mynames.A[grep('H7N9',mynames.A)])
  
  mynames.B = my.names.current[grep('_B_',my.names.current)]
  mynames.B = mynames.B[order(as.numeric(str_match(mynames.B, "20\\s*(.*?)\\s*_")[, 2]))]
  
  mynames.NA = my.names.current[grep('NA',my.names.current)]
  
  mynames.fc = c(mynames.A,mynames.B,mynames.NA)
  
  mynames.full = c(mynames.full,mynames.fc)
}

mynames.full = mynames.full[-1]


# bm.full = read_xlsx('IMPRINT_breastmilk_combined_features_11102025.xlsx')
# bm = read_xlsx('IMPRINT_breastmilk_combined_features_08272025.xlsx')

# bm.merge=merge(bm[,c(2,283)],bm.full,by='Barcode')
# write_xlsx(bm.merge,'IMPRINT_breastmilk_combined_features_01142026.xlsx')

matvax = read_xlsx('G_ALTER_Matvax_20250326_REW.xlsx')

# THIS IS THE COLUMN RANGE OF SYSTEMS SEROLOGY FEATURES
feature.range = 4:282

bm[,feature.range] = bm[,feature.range] %>% lapply(as.numeric)
bm = cbind(bm[,1:3],bm[,match(mynames.full,colnames(bm))],bm[,283:290])

# ADD A COLUMN THAT BINS VACCINE TIMINGS INTO CATEGORIES
#NEW VERSION
vax.timing = data.frame(matrix(NA,nrow = nrow(bm)))
colnames(vax.timing) = 'time'
# vax.timing$time[bm$trimester_at_vax_No_Vaccine==1 | (bm$trimester_at_vax_After_Pregnancy>0 &
#                                                        bm$vaccinated_while_bf==0 &
#                                                        bm$trimester_at_vax_First_Trimester==0 & 
#                                                        bm$trimester_at_vax_Second_Trimester==0 & 
#                                                        bm$trimester_at_vax_Third_Trimester==0 &
#                                                        bm$trimester_at_vax_Before_Pregnancy==0)] = 1
# 
# vax.timing$time[bm$trimester_at_vax_Before_Pregnancy>0 &
#                   bm$trimester_at_vax_First_Trimester==0 &
#                   bm$trimester_at_vax_Second_Trimester==0 &
#                   bm$trimester_at_vax_Third_Trimester==0 &
#                   bm$vaccinated_while_bf==0] = 2 # before pregnancy
# 
# 
# vax.timing$time[bm$trimester_at_vax_First_Trimester==1] = 3
# vax.timing$time[bm$trimester_at_vax_Second_Trimester==1] = 4
# vax.timing$time[bm$trimester_at_vax_Third_Trimester==1] = 5
# 
# vax.timing$time[bm$vaccinated_while_bf>0 &
#                   bm$trimester_at_vax_First_Trimester==0 &
#                   bm$trimester_at_vax_Second_Trimester==0 &
#                   bm$trimester_at_vax_Third_Trimester==0] = 6
vax.timing$time[bm$trimester_at_vax_No_Vaccine==1] = 1

vax.timing$time[bm$trimester_at_vax_Before_Pregnancy>0 &
                  bm$trimester_at_vax_First_Trimester==0 &
                  bm$trimester_at_vax_Second_Trimester==0 &
                  bm$trimester_at_vax_Third_Trimester==0] = 2 # before pregnancy


vax.timing$time[bm$trimester_at_vax_First_Trimester==1] = 3
vax.timing$time[bm$trimester_at_vax_Second_Trimester==1] = 4
vax.timing$time[bm$trimester_at_vax_Third_Trimester==1] = 5

vax.timing$time[bm$trimester_at_vax_After_Pregnancy>0 &
                  bm$trimester_at_vax_First_Trimester==0 &
                  bm$trimester_at_vax_Second_Trimester==0 &
                  bm$trimester_at_vax_Third_Trimester==0 &
                  bm$trimester_at_vax_Before_Pregnancy==0] = 6


bm$Vax_timing_bin = vax.timing$time


# ADD A COLUMN THAT BINS VACCINE TIMINGS INTO CATEGORIES
# OLD VERSION
# vax.timing = data.frame(matrix(NA,nrow = nrow(bm)))
# colnames(vax.timing) = 'time'
# vax.timing$time[bm$trimester_at_vax_No_Vaccine==1] = 1
# vax.timing$time[bm$trimester_at_vax_First_Trimester==1] = 2
# vax.timing$time[bm$trimester_at_vax_Second_Trimester==1] = 3
# vax.timing$time[bm$trimester_at_vax_Third_Trimester==1] = 4
# vax.timing$time[(bm$trimester_at_vax_After_Pregnancy>=1 | bm$trimester_at_vax_Before_Pregnancy>=1) & 
#                   bm$trimester_at_vax_First_Trimester!=1 & bm$trimester_at_vax_Second_Trimester!=1 & 
#                   bm$trimester_at_vax_Third_Trimester!=1] = 5 # group before or after pregnancy into a single category




# CREATE AN ORGANIZED VERSION OF THE BM DATA SET, ARRANGED BY VACCINE TIMING AND TIME POINT
bm.sort = bm[order(bm$Delta_weeks_bin,bm$Vax_timing_bin),]

# SET GGPLOT THEME
theme_set(theme_bw(16))
