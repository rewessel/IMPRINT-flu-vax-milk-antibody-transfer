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
bm = read_xlsx('Milk_serology.xlsx')
matvax = read_xlsx('Maternal_flu_vaccine_records.xlsx')
pca.variables = read_xlsx('PCA_variables.xlsx')
transfer.ratios = read_xlsx('Milk_serum_transfer_ratios.xlsx')

# THIS IS THE COLUMN RANGE OF SYSTEMS SEROLOGY FEATURES
feature.range = 6:284

bm[,feature.range] = bm[,feature.range] %>% lapply(as.numeric)

# CREATE AN ORGANIZED VERSION OF THE BM DATA SET, ARRANGED BY VACCINE TIMING AND TIME POINT
bm.sort = bm[order(bm$Time_point,bm$Vax_timing_bin),]

# SET GGPLOT THEME
theme_set(theme_bw(16))
