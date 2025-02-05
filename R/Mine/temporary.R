# Repeat initial proteomics analysis from Alexis in my own way and with non-batch corrected data
# TEMPORARY, TO THEN BE INCORPORATED TO THE MAIN CODE


# Libraries
library(dplyr)
library(data.table)
library(limma)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tibble)
library(stringr)
library(gprofiler2)
library(xlsx)
library(stringr)
library(gridExtra)
library(grid)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(factoextra)
library(glue)
library(dendextend)
library(ggrepel)
library(ggvenn)
library(tidyr)



# Load data
proteomics_raw <- as.data.frame(fread("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/11_Preprocessing_Boris/AlternativeAAUsage-tRNA/AlternativeAAUsage-tRNA_peptidecentric_ProteinMaxLFQ_PCAoutlier_removed_complete.tsv"))
sample_layout <- as.data.frame(fread("S:/AG/AG-CF-HTMS/AG-Ralser-Share/30-0092_AndreaLehmann-AlternativeAAUsage-tRNA/05_DataAnalysis/11_Preprocessing_Boris/AlternativeAAUsage-tRNA/AlternativeAAUsage-tRNA_peptidecentric_PrecursorQuantity_filename_annotations.tsv"))
master_dataset <- as.data.frame(fread("C:/MyStuff/tRNAs/Data/GtRNAdb/master_tRNA_dataset.csv"))

































