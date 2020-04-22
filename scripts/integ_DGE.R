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

#extract top 10 genes DE
top10.central.dge <- central.dge.df$Gene.name.uniq[1:10]

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

vlnplt <- VlnPlot(obj_integrated, features = top10.central.dge, group.by = "treatment", pt.size = .0, combine = FALSE)
wrap_plots(plots = vlnplt, ncol = 4)

RidgePlot(obj_integrated, features = top10.central.dge, group.by = "treatment")

DotPlot(obj_integrated, group.by = "treatment", features = top10.central.dge) + RotatedAxis() + scale_x_reverse()

Idents(obj_integrated) <- "celltype.treatment"
DoHeatmap(obj_integrated, features = top10.central.dge, group.by = "data.set")

str(obj_integrated@meta.data)
