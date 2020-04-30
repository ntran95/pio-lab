####Set up environment####
#install.packages("Seurat")
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
options(future.globals.maxSize = 5000 * 1024^2)

setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")

# Daniel stuff
if (FALSE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "adj_isl1_sib_counts_10X"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

####Read in Data####

#load in isl1_sib_counts_10X dir containing count matrix
#files in this folder should be in .gz format
isl1_sib_10X.data <- Read10X(data.dir = "../data/isl1_sib_counts_10X/")

adj_homeo.isl1_sib_10X <- CreateSeuratObject(counts = isl1_sib_10X.data, project = "homeo.isl1.sib.10X", min.cells = 1, min.features = 1)

####Standard pre-processing workflow####
#calculate mitochondrial contamination
#do not add suffix to identify sample in "percent.mt" for metadata, will cause convergance issues in downstream integration with multiple datasets
adj_homeo.isl1_sib_10X[["percent.mt"]] <- PercentageFeatureSet(adj_homeo.isl1_sib_10X, pattern = "^mt-")

####Subsetting####
#subsetting parameters: omit cells with genes less than 3000 and greater than 200. omit mitochondrial contamination greater than 10%. omit cells with molecules greater than 10,000
adj_homeo.isl1_sib_10X <- subset(adj_homeo.isl1_sib_10X, subset = nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA <20000) 

adj_homeo.isl1_sib_10X

plan("multiprocess")
adj_homeo.isl1_sib_10X <- NormalizeData(adj_homeo.isl1_sib_10X, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
#mutliprocess for NormalizationData using 1G is executable, not same case for ScaleData()

####FindVariableFeatures####
#Identification of highly variable genes
adj_homeo.isl1_sib_10X <- FindVariableFeatures(adj_homeo.isl1_sib_10X, selection.method = "vst", nfeatures = 2000)

adj_homeo.isl1_sib_10X <- ScaleData(adj_homeo.isl1_sib_10X, features = NULL, vars.to.regress = "nCount_RNA", verbose = TRUE)

####Perform PCA####
#default pc is set to 100
adj_homeo.isl1_sib_10X <- RunPCA(adj_homeo.isl1_sib_10X, npcs = 100, verbose = TRUE, features = NULL)

VizDimLoadings(adj_homeo.isl1_sib_10X, dims = 1:2, reduction = "pca")

DimHeatmap(adj_homeo.isl1_sib_10X, dims = 1:15, cells = 500, balanced = TRUE)

DimPlot(adj_homeo.isl1_sib_10X, reduction = "pca")

####Choosing PC####

####Elbow Plot####
ElbowPlot(adj_homeo.isl1_sib_10X, ndims = 50)

ElbowPlot(adj_homeo.isl1_sib_10X, ndims = 35)

ElbowPlot(adj_homeo.isl1_sib_10X, ndims = 20)

####Clustering/UMAP####
adj_homeo.isl1_sib_10X <- FindNeighbors(adj_homeo.isl1_sib_10X, dims = 1:10, features = NULL)

adj_homeo.isl1_sib_10X <- FindClusters(adj_homeo.isl1_sib_10X, resolution = 1.2)

adj_homeo.isl1_sib_10X <- RunUMAP(adj_homeo.isl1_sib_10X, dims = 1:10, reduction = "pca")


umap.unlabeled <- DimPlot(adj_homeo.isl1_sib_10X, reduction = "umap",
                          label = TRUE, pt.size= 0.4)
umap.unlabeled

####Common Features- FeaturePlot####
#common_features <- scan(paste0("/Volumes/easystore/SIMR_2019/pio-lab/data/gene-lists/common_neuromast_features.txt"), what = "character")
e <- FeaturePlot(adj_homeo.isl1_sib_10X, features = c("wnt2", "sost", "	znf185", "si:ch211-229d2.5", "si:ch73-261i21.5","pcna", "mki67"),
                 reduction = "umap", pt.size = 0.25, combine = FALSE, label = TRUE)
for (i in 1:length(e)) {
  e[[i]] <- e[[i]] + NoLegend() + NoAxes()
}
png(figurePath("DV.AP.AMP_featuresPC10.png"), width = 40,
    height = 25, units = "in", res = 400)
print(cowplot::plot_grid(plotlist = e, ncol = 3))
dev.off()

