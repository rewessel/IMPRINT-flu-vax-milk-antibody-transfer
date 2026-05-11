# SET UP CODE FOR IMPRINT BREAST MILK ANALYSIS

# LOAD REQUIRED PACKAGES
library(ggplot2)
library(pheatmap)
library(readxl)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(Hmisc)
library(ggpubr)
library(effsize)

# IMPORT DATA
# bm.new = read_xlsx('AB updated IMPRINT_BM_LuminexADCDdata_zeroed.xlsx')
# bm.new = bm.new[-177,]
# bm.new[bm.new<0] = 0

bm = read_xlsx('Milk_serology_v3.xlsx')
matvax = read_xlsx('Maternal_flu_vaccine_records.xlsx')
pca.variables = read_xlsx('PCA_variables.xlsx')
transfer.ratios = read_xlsx('Milk_serum_transfer_ratios_v2.xlsx')

# COMBINE PBS ZEROED DATA WITH RELEVANT COLUMNS FROM bm
# bm = merge(bm[,2:5],bm.new,by='Barcode')
bm = bm %>% select(!contains('EBV') & !contains('RSV') & !contains('Ebola'))

# THIS IS THE COLUMN RANGE OF SYSTEMS SEROLOGY FEATURES
feature.range = 6:284
# feature.range = 6:323

bm[,feature.range] = bm[,feature.range] %>% lapply(as.numeric)

# CREATE AN ORGANIZED VERSION OF THE BM DATA SET, ARRANGED BY VACCINE TIMING AND TIME POINT
bm.sort = bm[order(bm$Time_point,bm$Vax_timing_bin),]

# SET GGPLOT THEME
theme_set(theme_bw(16))

# ADD vax_timing_bin column to transfer.ratios DATA
# transfer.ratios = merge(transfer.ratios,bm[,c(1,5)],by='SubjectID')
# transfer.ratios = unique(transfer.ratios)
transfer.ratios.sort = transfer.ratios[order(transfer.ratios$Time_point,transfer.ratios$Vax_timing_bin),]

