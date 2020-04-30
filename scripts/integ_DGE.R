setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts/")

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


if (FALSE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "integ_DGE"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

readRDS("../../obj_intregrated.RDS")

central.diff.exp <- read.table("../../diff_exp_results.tsv", sep = "\t", header = T)

top10.diff.genes <- central.diff.exp %>% top_n(10, wt = avg_logFC) %>% pull(Gene.name.uniq)


levels(Idents(obj_integrated))

#subset Central cells 

central.sub <- subset(obj_integrated, idents = "central-cells")

central.sub

central.sub <- RunPCA(central.sub, npcs = 100, verbose = TRUE, features = NULL)
central.sub <- FindNeighbors(central.sub, dims = 1:10)
central.sub <- FindClusters(central.sub, resolution = .6)
central.sub <- RunUMAP(central.sub, reduction = "pca", dims = 1:10)

d1 <- DimPlot(central.sub, reduction = "umap", group.by  = "treatment", cols = c("red", "green3"))

d2 <- DimPlot(central.sub, reduction = "umap")

d1 + d2

DimPlot(central.sub, reduction = "umap", group.by  = "data.set")

FeaturePlot(central.sub, features = c("fosab", "zgc:158463", "btg2", "si:ch73-335l21.4", "sox3", "si:ch211-243g18.2", "si:dkey-222f2.1", "her15.1*1", "si:dkey-33m11.8", "her4.2*1"), split.by = "treatment", 
            cols = c("grey", "red"))

VlnPlot(central.sub, features = top10.diff.genes, pt.size = 0.4, group.by = "treatment")

# ===========================================================  Find markers for cluster/treatment of interest

obj_integrated$celltype.treatment <- paste(Idents(obj_integrated), obj_integrated$treatment, sep = "_")
#obj_integrated$cell.type.ident <- Idents(obj_integrated)
Idents(obj_integrated) <- "celltype.treatment"
central.dge.df <- FindMarkers(obj_integrated, ident.1 = "central-cells_1hr", ident.2 = "central-cells_homeo", verbose = FALSE)
head(central.dge.df, n = 15)

#change the column name formerally called "gene" generated from FindAllMarkers() to  "Gene.name.uniq" so we can merge with gene_table
#must have same column name in order to match genes
central.dge.df$Gene.name.uniq <- rownames(central.dge.df)
#inspect new col name change
names(central.dge.df)

#merge gene table with all markers generated based on uniq gene symbol
central.dge.df <- merge(central.dge.df, gene_table, by = "Gene.name.uniq") 

#reorder df in ascending order based on cluster (col #7) and then avg_logFC (col#3)
central.dge.df <- central.dge.df[order( central.dge.df[,3], decreasing = TRUE),]
write.table(central.dge.df, file = '../data/central.dge.df.tsv', sep = '\t', row.names = F)

#extract top 10 genes DE
top10.central.dge <- central.dge.df$Gene.name.uniq[1:10]

obj_integrated@meta.data$celltype.treatment <- factor(obj_integrated@meta.data$celltype.treatment)

obj_integrated@meta.data$treatment <- factor(obj_integrated@meta.data$treatment)


xy<- FeaturePlot(obj_integrated, features = top10.central.dge, split.by  = "treatment", reduction = "umap", pt.size = 0.25)
for (i in 1:length(xy)) {
  xy[[i]] <- xy[[i]] 
}
#png(figurePath("common_features.png"), width = 40,
#   height = 80, units = "in", res = 200)
print(cowplot::plot_grid(plotlist = xy, ncol = 4))
#dev.off()

# ===========================================================  Plots 

Idents(obj_integrated) <- "cell.type.ident"

DimPlot(obj_integrated, reduction = "umap", group.by  = "treatment", cols = c("red", "green3"))

DimPlot(obj_integrated, reduction = "umap", cells = top10.central.dge)

vlnplt <- VlnPlot(obj_integrated, features = top10.central.dge, group.by = "treatment", pt.size = .0, combine = FALSE, cols = c("red", "green3"))
wrap_plots(plots = vlnplt, ncol = 4)

png(figurePath("central.dge.vln.png"),
    width = 11, height = 9, units = "in", res = 300)
print(wrap_plots(plots = vlnplt, ncol = 4))
dev.off()

ridgeplt <- RidgePlot(obj_integrated, features = top10.central.dge, group.by = "treatment", cols = c("red", "green3"))

png(figurePath("central.dge.ridge.png"),
    width = 11, height = 9, units = "in", res = 300)
print(ridgeplt)
dev.off()

dp <- DotPlot(obj_integrated, split.by  = "treatment", features = top10.central.dge, cols = c("lightblue", "orangered"), col.min = -2.0, col.max = 2.0) + RotatedAxis()

png(figurePath("central.dge.dotplt.png"),
    width = 11, height = 9, units = "in", res = 300)
print(dp)
dev.off()

Idents(obj_integrated) <- "celltype.treatment"
DoHeatmap(obj_integrated, features = top10.central.dge, group.by = "data.set")

str(obj_integrated@meta.data)

top10.central.dge

# =========================================================== Create Function for DGE table
DGEtable <- function(seurat_obj, ident.1, ident.2) {
  gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
#if there isn't a column in the metadata specifying cluster ident by treatments, make one
  if (length(seurat_obj$celltype.treatment) == 0) {
    print("hi")
    seurat_obj$celltype.treatment <- paste(Idents(seurat_obj), seurat_obj$treatment, sep = "_")
  }
  print("else")
  Idents(seurat_obj) <- "celltype.treatment"
  df <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, verbose = FALSE)
  df$Gene.name.uniq <- rownames(df)
#inspect new col name change
  names(df)
#merge gene table with all markers generated based on uniq gene symbol
  df <- merge(df, gene_table, by = "Gene.name.uniq") 
#reorder df in ascending order based on cluster (col #7) and then avg_logFC (col#3)
  df <- df[order( df[,3], decreasing = TRUE),]
  return(df)
}

central.1hr.v.homeo <- DGEtable(obj_integrated, ident.1 = "central-cells_1hr", ident.2 = "central-cells_homeo")

central.1hr.v.homeo %>% select(Gene.name.uniq,avg_logFC) %>% filter(avg_logFC > 0.0) %>% pull(Gene.name.uniq)

DGEtable(obj_integrated, ident.1 = "central-cells_1hr", ident.2 = c("central-cells_homeo", "mantle-cells_homeo"))

# =========================================================== Measure UMI Abundance
options(scipen = 100, digits = 4)
png(figurePath("UMIAbundance.png"), width = 10,
    height = 6, units = "in", res = 400)
boxplot(homeo.isl1_sib_10X@meta.data$nCount_RNA, fpkm_matrix_1828_smartseq2@meta.data$nCount_RNA, 
        notch = TRUE, log = "y", names = c("10X", "Smartseq2"), 
        main = "UMI Abundance Comparison B/T 10X and SmartSeq2", col = c("lightblue", "purple"), 
        las=1, font.axis=1)
abline(h=mean(homeo.isl1_sib_10X@meta.data$nCount_RNA), lty = 2, col = "lightblue",)
abline(h=mean(fpkm_matrix_1828_smartseq2@meta.data$nCount_RNA), lty = 2, col = "purple")
axis(side=2, at=1861, labels = TRUE, las = 1, font.axis=1)
axis(side=2, at=2470, labels = TRUE, las = 1, font.axis=1)
axis(side=2, at=904524, labels = TRUE, las = 1, font.axis=1)
dev.off()

png(figurePath("hist-10X.png"), width = 10,
    height = 6, units = "in", res = 400)
hist(homeo.isl1_sib_10X@meta.data$nCount_RNA, main = "UMI Abundance-10X",
     xlab = "10X", col = "lightblue")
dev.off()

png(figurePath("hist-smartseq2.png"), width = 10,
    height = 6, units = "in", res = 400)
hist(fpkm_matrix_1828_smartseq2@meta.data$nCount_RNA,
     main = "UMI Abundance-SmartSeq2",
     xlab = "SmartSeq2",
     col = "purple")
dev.off()
