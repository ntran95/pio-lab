library(Seurat)
# library(plotly)
#install.packages("future")
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
# library(gridExtra)
# library(ggrepel)
library(ggrepel)
library(stringr)
library(patchwork)
options(future.globals.maxSize = 5000 * 1024^2)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "integ_DGE"
  
  figurePath <- function(filename, format){paste0("/pio-lab/scripts/", script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

# =========================================================== Load filtered integ smartseq and homeo object

obj_integrated_filtered <- readRDS("../data/filtered_adj_fpkm_1828_smartseq_integ.RDS")

# =========================================================== FindMarkers by treatments and seq method

obj_integrated_filtered$cell.type.treatment.method <- paste(Idents(obj_integrated_filtered), obj_integrated_filtered$treatment, obj_integrated_filtered$seq.method, sep = "_")

obj_integrated_filtered$cell.type.treatment.method

Idents(obj_integrated_filtered) <- "cell.type.treatment.method"

x <- FindMarkers(obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2", verbose = TRUE)


# =========================================================== Create Function for DGE table
DGEtable <- function(seurat_obj, ident.1, ident.2) {
  gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
  #if there isn't a column in the metadata specifying cluster ident by treatments, make one
  if (length(seurat_obj$cell.type.treatment.method) == 0) {
    print("creating new column in metadata")
    seurat_obj$cell.type.treatment.method <- paste(Idents(seurat_obj), seurat_obj$treatment, seurat_obj$seq.method, sep = "_")
  }
  print("switching to identity based on treatment and method")
  Idents(seurat_obj) <- "cell.type.treatment.method"
  df <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, verbose = TRUE)
  df$Gene.name.uniq <- rownames(df)
  #inspect new col name change
  names(df)
  #merge gene table with all markers generated based on uniq gene symbol
  df <- merge(df, gene_table, by = "Gene.name.uniq") 
  #reorder df in ascending order based on cluster (col #7) and then avg_logFC (col#3)
  df <- df[order( df[,3], decreasing = TRUE),]
  return(df)
}

y <-DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2")

meta <- obj_integrated_filtered@meta.data

cluster.ident <- strsplit(unique(as.character(na.omit(meta$cell.type.ident))), split = "[][']|,\\s*")

cluster.ident.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)){
  print(paste0("working on cluster: ",cluster.ident[[i]]))
  cluster.ident.list[[i]] <- DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2"))
  assign(paste0(cluster.ident[[i]],".DGE.1hrvHomeo"), DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2")))
}

