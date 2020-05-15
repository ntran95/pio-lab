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

mphg_homeo <- RunPCA(mphg_homeo, npcs = 100, verbose = TRUE, features = NULL)
mphg_homeo <- FindNeighbors(mphg_homeo, dims = 1:15, verbose = TRUE)
mphg_homeo <- FindClusters(mphg_homeo, resolution = .8, verbose = TRUE)
mphg_homeo <- RunUMAP(mphg_homeo, reduction = "pca", dims = 1:15, verbose = TRUE)

DimPlot(mphg_homeo)

# =========================================================== annotation ===========================================================
meta <- mphg_homeo@meta.data

#split the cell cluster idents extracted from metadata into a uniquely named list
cluster.ident <- strsplit(unique(as.character(meta$cell.type.ident)), split = "[][']|,\\s*")

#create empty list to store counts of cells to cluster idents
cluster.ident.list <- vector("list", length = length(cluster.ident))

#loop through the cell type ident, match it with the ident from the metadata, tally up the numberof cells expressed in that cluster
for (i in 1:length(cluster.ident)){
  cluster.ident.list[[i]] <- meta %>% filter(cell.type.ident == cluster.ident[[i]]) %>% count(seurat_clusters, sort =TRUE)
  print(cluster.ident[[i]])
}

names(cluster.ident.list) <- cluster.ident

cluster.ident.list


#desired unlabelled clusters passed through 
pos_list <- c(14, 15)
pos_clusters <- vector(mode = "list", length = length(pos_list))
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)

for (x in 1:length(pos_list)){
  print(paste0("finding markers for cluster: ", pos_list[x]))
  pos_clusters[[x]] <- FindMarkers(mphg_homeo, ident.1 = pos_list[x], only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10, verbose = TRUE)
  print(pos_clusters[[x]])
  pos_clusters[[x]]["Gene.name.uniq"] <- row.names(pos_clusters[[x]])
  pos_clusters[[x]] <- merge(pos_clusters[[x]], gene_table, by = "Gene.name.uniq")
}

names(pos_clusters) <- paste("cluster_", pos_list)
pos_clusters$`cluster_ 10`


# =========================================================== Rename filtered integ clusters  ===========================================================

colnames(meta)
cells <- list("mcamb" = 2, "rbp4" = 11, "fabp3" = 3, "f13a1b" = 0, "tspan10" = 1, "eomesa" = 12, "pcna" = c(6,9), "cldnh" =8, 
              "irg1" = 7,"runx3" = 4,"spock3" = 13, "mpx" = 5, "gata3" = 12, "stat1b" = 10, "14" = 14, "15" = 15)
meta$homeo.cell.type.ident <- factor(rep("", nrow(meta)),
                                       levels = names(cells), ordered = TRUE)
for (i in 1:length(cells)) {
  meta$homeo.cell.type.ident[meta$seurat_clusters %in% cells[[i]]] <- names(cells)[i]
}
mphg_homeo@meta.data <- meta
Idents(mphg_homeo) <- mphg_homeo@meta.data$homeo.cell.type.ident

png(figurePath(paste0("annotated.umap.png"))
    ,width = 11, height = 9, units = "in", res = 300)
DimPlot(mphg_homeo, label = TRUE) + NoLegend()
dev.off()

all.markers <- FindAllMarkers(mphg_homeo, logfc.threshold = 0.10)

# =========================================================== Rename filtered integ clusters  ===========================================================

png(figurePath(paste0("annotated.umap.png"))
    ,width = 11, height = 9, units = "in", res = 300)
DimPlot(mphg_homeo, label = TRUE, group.by = c("cell.type.ident", "homeo.cell.type.ident")) + NoLegend()
dev.off()

library(RColorBrewer)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(integer_with_Nclusters)


png(figurePath(paste0("old.cluster.ident.png"))
    ,width = 11, height = 9, units = "in", res = 300)
DimPlot(mphg_homeo, reduction = "umap", pt.size = 0.20,
         label = TRUE, label.size = 4, group.by   = "cell.type.ident", cols = gg_color_hue(16)) 
dev.off()