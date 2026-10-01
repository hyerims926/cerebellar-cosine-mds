### Set working directory ###
setwd("~/data/dev_0813")

### Set Required Libraries ###
library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(SingleCellExperiment)
library(harmony)
library(RSpectra)
library(slingshot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(MASS)
library(plotly)
library(tradeSeq)
library(tibble)
library(pheatmap)
library(emdist)
library(tidytext)
library(readr)
library(msigdbr)
library(FNN)
install.packages("msigdbr")

use_python("/home/rstudio/.cache/R/reticulate/uv/cache/archive-v0/jI1FTbpQ61w1gysg0o7YX/bin/python", required = TRUE)
library(zellkonverter)
library(biomaRt)

### color ###

cols <- c(
  "progenitor" = "#6abd12",
  "MBO" = "#bce590",
  
  "VZ_neuroblast_1" = "#facc15",
  "VZ_neuroblast_2" = "#fde047", 
  "VZ_neuroblast_3" = "#fef08a",
  
  "NTZ_neuroblast_1" = "#7d5207",
  "NTZ_neuroblast_2" = "#b39054",
  "NTZ_neuroblast_3" = "#ecd09f",
  
  "GCP" = "#7f1d1d",
  "GC_diff_1" = "#f76161",
  "GC_diff_2" = "#ff9f9f",
  "GC_defined" = "#ffd5d5",
  
  "GCP/UBCP" = "#f27500",
  "GC/UBC_diff" = "#fb9f4a",
  "UBC_diff" = "#f7b272",
  "UBC_defined" = "#fdc999",
  
  "Purkinje_diff" = "#0a229b",
  "Purkinje_defined" = "#90a2f2",
  "Purkinje_maturing" = "#cbd4f9",
  
  "interneuron_diff" = "#069372",
  "interneuron_defined" = "#7cd9c3",
  "GABA_DN_defined" = "#c2ece2",
  
  "glut_DN_defined" = "#ae3970",
  "isth_N_diff" = "#fdb8e9",
  "isth_N_defined" = "#ffe1f7",
  
  "glioblast" = "#550584",
  "oligo_progenitor" = "#8c2ebf",
  "oligodendrocyte" = "#c38ce3",
  "astrocyte" = "#d7c8f0",
  
  "meningeal" = "grey10",
  "parabrachial" = "grey30")

### Set seed ###
seed <- 1234
set.seed(seed)

### Read h5ad ###

path_to_h5ad <- "~"

sce <- readH5AD(path_to_h5ad)
# sce <- readH5AD("~/data/115k_mouse_develop.h5ad")
assayNames(sce)
obj_115k <- as.Seurat(
  sce,
  counts = "X",
  data = "X")

obj_115k[["RNA"]] <- obj_115k[["originalexp"]]
DefaultAssay(obj_115k) <- "RNA"
obj_115k[["originalexp"]] <- NULL

table(obj_115k$author_stage)

tapply(obj_115k$nFeature_originalexp, obj_115k$author_stage, 
       quantile, probs = c(0.01, 0.05, 0.5, 0.95, 0.995), na.rm = TRUE)

obj_115k <- subset(
  obj_115k,
  subset =
    nFeature_originalexp >
    ave(nFeature_originalexp, author_stage,
        FUN = function(x) quantile(x, 0.02, na.rm = TRUE)) &
    nFeature_originalexp <
    ave(nFeature_originalexp, author_stage,
        FUN = function(x) quantile(x, 0.995, na.rm = TRUE)))

table(obj_115k$precisest_label)
table(obj_115k$author_stage)
table(obj_115k$dev_state)

levels_115k <- c(
  "progenitor",
  "neural_crest_progenitor",
  "MBO",
  
  "VZ_neuroblast_1",
  "VZ_neuroblast_2",
  "VZ_neuroblast_3",
  
  "NTZ_neuroblast_1",
  "NTZ_neuroblast_2",
  "NTZ_neuroblast_3",
  
  "GCP",
  "GCP/UBCP",
  "GC/UBC_diff",
  "GC_diff_1",
  "GC_diff_2",
  "GC_defined",
  "UBC_diff",
  "UBC_defined",
  
  "Purkinje_diff",
  "Purkinje_maturing",
  "Purkinje_defined",
  
  "interneuron_diff",
  "interneuron_defined",
  "GABA_DN_defined",
  
  "glut_DN_defined",
  "isth_N_diff",
  "isth_N_defined",
  "parabrachial",
  "motorneuron",
  "noradrenergic",
  
  "glioblast",
  "astrocyte",
  "oligo_progenitor",
  "oligodendrocyte",
  "ependymal",
  "mural/endoth",
  "meningeal",
  "immune",
  "erythroid")

obj_115k$anno <- factor(obj_115k$dev_state, levels = levels_115k)
table(obj_115k$anno)

saveRDS(obj_115k, "Rdata/obj_115k.RDS")
obj_115k <- readRDS("Rdata/obj_115k.RDS")

obj_E10.5 <- subset(obj_115k, subset = author_stage %in% c("E10.5"))
saveRDS(obj_E10.5, "Rdata/obj_E10.5.RDS")

obj_E11.5 <- subset(obj_115k, subset = author_stage %in% c("E11.5"))
saveRDS(obj_E11.5, "Rdata/obj_E11.5.RDS")

obj_E12.5 <- subset(obj_115k, subset = author_stage %in% c("E12.5"))
saveRDS(obj_E12.5, "Rdata/obj_E12.5.RDS")

obj_E13.5 <- subset(obj_115k, subset = author_stage %in% c("E13.5"))
saveRDS(obj_E13.5, "Rdata/obj_E13.5.RDS")

obj_E14.5 <- subset(obj_115k, subset = author_stage %in% c("E14.5"))
saveRDS(obj_E14.5, "Rdata/obj_E14.5.RDS")

obj_E15.5 <- subset(obj_115k, subset = author_stage %in% c("E15.5"))
saveRDS(obj_E15.5, "Rdata/obj_E15.5.RDS")

obj_E17.5 <- subset(obj_115k, subset = author_stage %in% c("E17.5"))
saveRDS(obj_E17.5, "Rdata/obj_E17.5.RDS")

obj_P0 <- subset(obj_115k, subset = author_stage %in% c("P0"))
saveRDS(obj_P0, "Rdata/obj_P0.RDS")

obj_P4 <- subset(obj_115k, subset = author_stage %in% c("P4"))
saveRDS(obj_P4, "Rdata/obj_P4.RDS")

obj_P7 <- subset(obj_115k, subset = author_stage %in% c("P7"))
saveRDS(obj_P7, "Rdata/obj_P7.RDS")

obj_P14 <- subset(obj_115k, subset = author_stage %in% c("P14"))
saveRDS(obj_P14, "Rdata/obj_P14.RDS")

obj_adult <- subset(obj_115k, subset = author_stage %in% c("adult (9 weeks)"))
saveRDS(obj_adult, "Rdata/obj_adult.RDS")

