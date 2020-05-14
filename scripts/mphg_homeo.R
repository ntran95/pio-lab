library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)
library(ggrepel)
library(readxl)
library(writexl)
options(future.globals.maxSize = 5000 * 1024^2)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "mphg_homeo"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
}

# =========================================================== Read in integrated macrophage seurat object =================================================

mphg_integrated <- readRDS("../data/mphg_homeo/SeurObj_anchored_update3_mphg_regen_14k_cells_seurat3_v1.0_.RDS")

meta_mphg_integrated <- mphg_integrated@meta.data

#check how many cells pertain to which treatment in the original mphg integrated object
table(mphg_integrated$data.set)
#homeo - 13937

length(rownames(mphg_integrated@meta.data))
# =========================================================== Subset for Homeo Treatment =================================================
#subset for only homeo treatments, specified by the "data.set" column in the meta.data
mphg_homeo <- subset(mphg_integrated, subset  = data.set == "homeo")

#check amt of cells selected for homeo trt is correctly subseted
table(mphg_homeo$data.set)
#homeo-13937
length(rownames(mphg_homeo@meta.data))

meta_mphg_homeo <- mphg_homeo@meta.data

# =========================================================== Recluster filtered obj after subsetting  ===========================================================

if (DefaultAssay(mphg_homeo) == "integrated"){
  print("switching assay to 'RNA'")
  DefaultAssay(mphg_homeo) <- "RNA"
}else {
  print("assay already defaulted to RNA")
  #check if normalized matrix is mo
  if(all(is.na(mphg_homeo[["RNA"]]@data)) == "FALSE"){
    print("current object already normalized")
  }else{
    mphg_homeo <- NormalizeData(mphg_homeo, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
  }
  #check if the variable features is already populated by checking if obj[["RNA]]@var.features is empty
  if (all(is.na(mphg_homeo[["RNA"]]@var.features)) == "TRUE") {
    print("populating variable features")
    mphg_homeo <- FindVariableFeatures(mphg_homeo, selection.method = "vst", nfeatures = 2000)
  }
  #check if scaled data is empty, if TRUE, scale data
  if(all(is.na(mphg_homeo[["RNA"]]@scale.data)) == "TRUE"){
    print("current object not scaled, scaling now")
    mphg_homeo <- ScaleData(mphg_homeo)
  }
}

