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
library(ggrepel)
#library(patchwork)
options(future.globals.maxSize = 5000 * 1024^2)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "integ_DGE"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

# =========================================================== Load filtered integ smartseq and homeo object =================================================

obj_integrated_filtered <- readRDS("../data/filtered_adj_fpkm_1828_smartseq_integ.RDS")

# =========================================================== FindMarkers by treatments and seq method

obj_integrated_filtered$cell.type.treatment.method <- paste(Idents(obj_integrated_filtered), obj_integrated_filtered$treatment, obj_integrated_filtered$seq.method, sep = "_")

obj_integrated_filtered$cell.type.treatment.method

Idents(obj_integrated_filtered) <- "cell.type.treatment.method"

x <- FindMarkers(obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2", verbose = TRUE)


# =========================================================== Create Function for DGE table =================================================
DGEtable <- function(seurat_obj, ident.1, ident.2) {
  gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
  #ensure that the default assay is switched to "RNA" prior to DGE analysis
  if (DefaultAssay(obj_integrated_filtered) == "integrated"){
    print("switching assay to RNA")
    DefaultAssay(seurat_obj) <- "RNA"
  }
  #if there isn't a column in the metadata specifying cluster ident by treatments, make one
  if(!"cell.type.treatment.method" %in% colnames(seurat_obj@meta.data)){
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
  #reorder df in ascending order based on avg_logFC (col#3)
  df <- df[order( df[,3], decreasing = TRUE),]
  return(df)
}

y <-DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2")

meta <- obj_integrated_filtered@meta.data

#sort alphabetically, extract unique cluster idents, omit NA's, turn to characters, split each idents into elements of list
cluster.ident <- strsplit(sort(unique(as.character(na.omit(meta$cell.type.ident)))), split = "[][']|,\\s*")

cluster.ident.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)){
  print(paste0("working on cluster: ",cluster.ident[[i]]))
  cluster.ident.list[[i]] <- DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2"))
  assign(paste0(cluster.ident[[i]],".smrtseq.1hrvHomeo"), DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2")))
}

names(cluster.ident.list) <- cluster.ident

# =========================================================== Read in Sungmin's Regen App data  =================================================s
setwd("../data/sungminregenappDGE/")

file_list <- list.files(path="../sungminregenappDGE/")

sungmin.regen.cluster.ident.list <- vector("list", length = length(file_list))

for (i in 1:length(file_list)){
  sungmin.regen.cluster.ident.list[[i]] <- read.table(file = file_list[[i]], header = TRUE, sep = "\t")
  assign(paste0(gsub(".tsv", "" ,file_list[[i]]),".sungmin.1hrvHomeo"), read.table(file = file_list[[i]], header = TRUE, sep = "\t"))
}

names(sungmin.regen.cluster.ident.list) <- gsub(".tsv", "", file_list)

# =========================================================== Merge smartseq and sungmin's regen DGE  =================================================s

#empty list to store merged DGE tables
merged.data <- vector("list", length = length(cluster.ident))

conserved.smrtseq <- vector("list", length = length(cluster.ident))

scatterplot.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)) {
  #merge the smrtseq data with sungmin's data, store within a list, merged.data
  merged.data[[i]] <- merge(x = cluster.ident.list[[i]], y = sungmin.regen.cluster.ident.list[[i]], suffixes = c(".smrtseq2", ".regen.app"),
                            by = c("Gene.name.uniq", "Gene.stable.ID", "Gene.name", "ZFIN.ID", "Gene.description"), sort = FALSE, all = TRUE)
  merged.data[[i]]$diff_avg_logFC <- merged.data[[i]]$avg_logFC.smrtseq2 - merged.data[[i]]$avg_logFC.regen.app
  #generate scatterplots
  png(figurePath(paste0(cluster.ident[[i]], ".png")), width = 11,
      height = 9, units = "in", res = 200)
  scatterplot.list[[i]] <- ggplot(merged.data[[i]], aes(avg_logFC.regen.app, avg_logFC.smrtseq2, label = Gene.name.uniq)) + 
    geom_point() + 
    ggtitle(paste0(cluster.ident[[i]])) + 
    geom_text(aes(label=Gene.name.uniq),hjust=0, vjust=0) + 
    coord_fixed(xlim = c(-1, 3), ylim = c(-1,3)) + 
    scale_x_continuous(name ="Sungmin Regen App - avg LogFC - 1hrVShomeo") +
    scale_y_continuous(name = "smartseq2 - avg LogFC - 1hrVShomeo")
  print(scatterplot.list[[i]])
  dev.off()
  #look through each row of the merged dataframe for each cluster ident
  for (x in 1:nrow(merged.data[[i]])){
    #search for NA's in the avg logFC for sungmin's dataset, indicating genes that are only conserved in the smrtseq data
    if (is.na(merged.data[[i]]$avg_logFC.regen.app[x]) == "TRUE"){
      conserved.smrtseq[[i]] <- merged.data[[i]] %>% select(Gene.name.uniq, avg_logFC.smrtseq2, avg_logFC.regen.app)
    }
  }
  
}
names(merged.data) <- cluster.ident
names(conserved.smrtseq) <- cluster.ident
names(scatterplot.list) <- cluster.ident


merge <- merge(x = `AP-cells.smrtseq.1hrvHomeo`, y = `AP-cells.sungmin.1hrvHomeo`, suffixes = c(".smrtseq2", ".regen.app"),
               by = c("Gene.name.uniq", "Gene.stable.ID", "Gene.name", "ZFIN.ID", "Gene.description"), sort = FALSE, all = TRUE)

amp <-merged.data$`amp-SCs`

scatterplot.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)) {
  scatterplot.list[[i]] <- ggplot(merged.data[[i]], aes(avg_logFC.regen.app, avg_logFC.smrtseq2, label = Gene.name.uniq))  + coord_fixed(xlim = c(-1, 3), ylim = c(-1,3)) + scale_x_continuous(name ="Sungmin Regen App - avg LogFC") +
    geom_label_repel(aes(label = Gene.name.uniq),
                     box.padding   = 0.35, 
                     point.padding = 0.5,
                     segment.color = 'grey50') +  theme_classic() +
    geom_point(color = "blue", size = 3)
  print(scatterplot.list[[i]])
  }
names(scatterplot.list) <- cluster.ident

scatterplot.list$`early-HCs`
