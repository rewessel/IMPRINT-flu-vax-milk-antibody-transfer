# Flu-vaccination-breastmilk-transfer
Last updated by Remziye Wessel (wesselr@mit.edu) on 20 Feb 2026.

This repository contains R and MATLAB scripts to perform analysis in the manuscript "Maternal influenza vaccination boosts durable breast milk antibody-mediated protection in infants" (Wessel and Tong et al, 2026).


Files:

IMPRINT_breastmilk_analysis.R

IMPRINT_breastmilk_analysis.m

Figure_1.R - Plots Fig1A heatmap and Fig1B boxplots (requires Milk_serum_transfer_ratios.xlsx). 

Figure_2.m - Performs O-PLSDA modeling shown in Figure 2.

Figure_3.m - Plots line graphs in Fig3F-I and Figure S6 (manually tell the script which antigen you want to plot).

Figure_3.R - Calculates cohen's d effect size and plots heatmap in Fig3A, plots boxplots and performs statistical analysis in Fig 3B-E and Figure S5 (manually tell the script which antigen you want to plot).

Figure_4.m - Performs O-PLSDA modeling shown in Figure 4.

Figure_5.m - Performs O-PLSDA modeling shown in Figure 5 and Figure S7-S8 (requires Fig6_combined_serum_milk.xlsx, which is subsetted data to only infants with matched cord serum and milk samples and whose first 6 months of life coincided with a flu season).

Figure_S1.R - Performs PCA and plots scores plots colored by a select few metadata variables stored in PCA_variables.xlsx.

IMPRINT_breastmilk_setup.R - Loads all required packages and reads dataframes needed to run anaylses in R.



Software:

R version 4.5.1

MATLAB version R2024a


R Packages:

ggplot2 (4.0.0)

pheatmap (1.0.3)

readxl (1.4.5)

tidyverse (2.0.0)

writexl (1.5.4)

RColorBrewer (1.1.3)

ggpubr (0.6.2)

effsize (0.8.1)

cowplot (1.2.0)

Other dependencies:

You will need to install the MATLAB PLSR package from rewessel/PLSR-DA_MATLAB to reproduce Figures 2, 4, 5, S7, and S8.


All code was executed on a Dell machine running Windows 11 Pro (v25H2).

