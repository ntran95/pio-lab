library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
options(future.globals.maxSize = 8000 * 1024^2)
library(RColorBrewer)
library(stringr)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "homeo_samples_integ"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}

}

# =========================================================== Read in Exp Matrix in Loop
sample.ls <- c("L29314", "L34727", "L34728")
#create empty list for exp matrix
samples.data.list <- vector("list", length = length(sample.ls))

for (i in 1:length(sample.ls)){
  print(paste0("reading in sample: ", sample.ls[i]))
  samples.data.list[[i]] <- Read10X(data.dir = paste0("../data/homeo_sample_", sample.ls[i], "/filtered_feature_bc_matrix/zipped/"))
  print(paste0("working dir: ", getwd()))
  print(paste0("Pulling matrix from: ../data/homeo_sample_", sample.ls[i], "/filtered_feature_bc_matrix/zipped/"))
  print(paste0("finished reading: ", sample.ls[i]))
}
names(samples.data.list) <- sample.ls

# =========================================================== Create Seurat Object
#the dims() of the count mtx in the object matches the original in samples.data.list() when using lapply
#no features or cells subsetted
obj_list <- lapply(samples.data.list, CreateSeuratObject, project = "homeo.integ")
names(obj_list) <- sample.ls

for (i in 1:length(sample.ls)){
  print(obj_list[i])
  #specify dataset
  obj_list[[i]]@meta.data$data.set <- paste0("homeo-",sample.ls[i])
  #specify treatment
  obj_list[[i]]@meta.data$treatment <- "homeo"
  #calculate percent mito
  obj_list[[i]]@meta.data[["percent.mt"]] <- PercentageFeatureSet(obj_list[[i]], pattern = "^mt-")
  obj_list[[i]]@meta.data$percent.mt <- obj_list[[i]]@meta.data$percent.mt$nCount_RNA
  #QC violin
  #png(figurePath(paste0(sample.ls[i],"QC-vln.png")),
   #   width = 11, height = 9, units = "in", res = 300)
  #print(VlnPlot(obj_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  #dev.off()
  print(summary(obj_list[[i]]$nFeature_RNA))
  print(summary(obj_list[[i]]$percent.mt))
  #print("The original object summary: ", obj_list[i])
  obj_list[[i]] <- subset(obj_list[[i]], subset = nCount_RNA < 15000 & percent.mt < 10)
  print(obj_list[[i]])
}


#obj_list <- lapply(obj_list, subset, subset = nCount_RNA < 15000 & percent.mt < 10)

#rename seur objects in obj_list
L29314_homeo <- obj_list$L29314
L34727_homeo <- obj_list$L34727
L34728_homeo <- obj_list$L34728

# =========================================================== Normalize/Find Variable Features

for (i in 1:length(obj_list)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]], normalization.method = "LogNormalize", 
                                 scale.factor = 10000, verbose = TRUE)
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = TRUE)
  
}
#check that normalization works 
obj_list$L29314[["RNA"]]@data

# =========================================================== Add the isl1 homeo data
adj_homeo_isl1_10X <- readRDS("../data/adj_homeo.isl1_sib_10X.RDS")

adj_homeo_isl1_10X@meta.data$data.set <- "homeo-10X-isl1"

sample.ls <- c("isl1_10X","L29314", "L34727", "L34728")

obj_list <- c(adj_homeo_isl1_10X,obj_list)

names(obj_list) <- sample.ls

  
# =========================================================== Integration
all_shared_genes <- lapply(obj_list, row.names) %>% Reduce(intersect, .) 

all_shared_genes

dims = 1:15
SC_transform <- FALSE

if (SC_transform) {
  for (i in 1:length(obj_list)) {
    print("hi")
    obj_list[[i]] <- SCTransform(
      obj_list[[i]], verbose = FALSE)
  }
  seurat_obj_features <- SelectIntegrationFeatures(
    object.list = obj_list, nfeatures = 2000)
  obj_list <- PrepSCTIntegration(
    object.list = obj_list, anchor.features = all_shared_genes)
  # reference = 1 is the homeostatic dataset in obj list
  obj_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                        anchor.features = seurat_obj_features, normalization.method = "SCT",
                                        dims = dims, reference = 1) 
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims,
                                  normalization.method = "SCT")
} else {
  for (i in 1:length(obj_list)) {
    obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = FALSE)
    obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
  }
  obj_anchors <- FindIntegrationAnchors(object.list = obj_list,
                                        dims = dims, reference = 1) # 1 is the homeostatic dataset in obj list
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims, features.to.integrate = all_shared_genes)
  DefaultAssay(obj_integrated) <- "integrated"
}

obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE)


# =========================================================== Adjust object integrated metadata ===========================================================

obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE)

#rename homeo samples L##### to LIMS order number
meta_obj_integ <- obj_integrated@meta.data

data.set <- as.vector(obj_integrated@meta.data$data.set)

meta_obj_integ <- meta_obj_integ%>% mutate(data.set =case_when(str_detect(data.set, "homeo-L29314") ~ "homeo-2047", 
                                                               str_detect(data.set, "homeo-L34727") ~ "homeo-2410-7",
                                                               str_detect(data.set, "homeo-L34728") ~ "homeo-2410-8",
                                                               TRUE ~ as.vector(obj_integrated@meta.data$data.set)))

obj_integrated@meta.data$data.set <- meta_obj_integ$data.set


# =========================================================== UMAP/Clustering
dim_list <- c(10,15,20,25,30)

common_features <- scan(paste0("../data/gene-lists/common_neuromast_features.txt"), what = "character")

if (!SC_transform) {
  plan("multiprocess")
  obj_integrated <- ScaleData(obj_integrated, verbose = TRUE, vars.to.regress = "nCount_RNA")
  obj_integrated <- RunPCA(obj_integrated, npcs = 100, verbose = TRUE, features = NULL)
}
for (pc in dim_list){
  obj_integrated <- FindNeighbors(obj_integrated, dims = 1:pc, verbose = TRUE)
  obj_integrated <- FindClusters(obj_integrated, resolution = 1.0, verbose = TRUE)
  obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:pc, verbose = TRUE)
  png(figurePath(paste0("umap.by.dataset.PC",pc,".png"))
      ,width = 11, height = 9, units = "in", res = 300)
  print(DimPlot(obj_integrated, group.by = "data.set", cols = brewer.pal(n = 4, name = "RdBu")) + DarkTheme()) 
  dev.off()
  png(figurePath(paste0("umap.unlabelled.PC",pc,".png"))
      ,width = 11, height = 9, units = "in", res = 300)
  print(DimPlot(obj_integrated))
  dev.off()
  DefaultAssay(obj_integrated) <- "RNA"
  
  integrated_featplt <- FeaturePlot(obj_integrated, common_features,
                                    reduction = "umap", pt.size = 0.25, combine = FALSE, label = TRUE)
  for (i in 1:length(integrated_featplt)) {
    integrated_featplt[[i]] <- integrated_featplt[[i]] + NoLegend() + NoAxes()
  }
  png(figurePath(paste0("common_features.PC",pc,".png")), width = 40,
      height = 80, units = "in", res = 200)
  print(cowplot::plot_grid(plotlist = integrated_featplt, ncol = 4))
  dev.off()
  DefaultAssay(obj_integrated) <- "integrated"
  
  
}

# =========================================================== PC 25 Annotating
obj_integrated <- FindNeighbors(obj_integrated, dims = 1:25, verbose = TRUE)
obj_integrated <- FindClusters(obj_integrated, resolution = 1.0, verbose = TRUE)
obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:25, verbose = TRUE)
DimPlot(obj_integrated, group.by = "data.set")
DimPlot(obj_integrated)

saveRDS(object = obj_integrated, file = "../data/homeo_samples_integ.RDS")

obj_integrated <- readRDS("../data/homeo_samples_integ.RDS")

meta_common_features <- read.table(file = "../data/gene-lists/meta_common_features.tsv", sep = "", header = T)

#plan("multiprocess")
#all.markers.integ <- FindAllMarkers(obj_integrated, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.01, return.thresh = 0.001, verbose = TRUE)

gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)

# =========================================================== Annotate -pos clusters ===========================================================
#desired unlabelled clusters passed through 
pos_list <- c(12,20,21,22,23,24,25)
pos_clusters <- vector(mode = "list", length = length(pos_list))

 for (x in pos_list){
   print(paste0("finding markers for cluster: ", x))
   #plan("multiprocess")
   assign(paste0("cluster_",x), FindMarkers(obj_integrated, ident.1 = x, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10, verbose = TRUE))
   
   }


for (x in 1:length(pos_list)){
  print(paste0("finding markers for cluster: ", pos_list[x]))
  plan("multiprocess")
  pos_clusters[[x]] <- FindMarkers(obj_integrated, ident.1 = pos_list[x], only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10, verbose = TRUE)
  print(pos_clusters[[x]])
  pos_clusters[[x]]["Gene.name.uniq"] <- row.names(pos_clusters[[x]])
  pos_clusters[[x]] <- merge(pos_clusters[[x]], gene_table, by = "Gene.name.uniq")
  }

names(pos_clusters) <- paste("cluster_", pos_list)

# =========================================================== Rename Clusters ===========================================================

meta <- obj_integrated@meta.data
colnames(meta)
cells <- list("mature-HCs" = 0, "early-HCs" = 14,  "HC-prog" = 16,
              "central-cells" = c(1,2,3,9,19), "DV-cells" = c(10), "AP-cells" = c(8,18),
              "amp-SCs" = c(11,17), "mantle-cells" = c(5,6), "Inm" = c(4,7,13,15),"col1a1b-pos" = 12,
              "clec14a-pos" = 21, "mfap4-pos" = 23, "c1qtnf5-pos" = 24, "apoa1b-pos" = 25, "krt91-pos" = 20, "hbbe1-pos" = 22)
meta$cell.type.ident <- factor(rep("", nrow(meta)),
                               levels = names(cells), ordered = TRUE)
for (i in 1:length(cells)) {
  meta$cell.type.ident[meta$seurat_clusters %in% cells[[i]]] <- names(cells)[i]
}
obj_integrated@meta.data <- meta
Idents(obj_integrated) <- obj_integrated@meta.data$cell.type.ident

umap.labeled <- DimPlot(obj_integrated, reduction = "umap", label = TRUE, pt.size= 0.4) + NoLegend()

png(figurePath("annotated.clusters.png"), width = 11,
    height = 9, units = "in", res = 200)
print(umap.labeled)
dev.off()

saveRDS(object = obj_integrated, file = "../data/homeo_samples_integ.RDS")

# =========================================================== Remove unspecific -pos clusters ===========================================================
Ident.list <- levels(Idents(obj_integrated))
neuromast.list <- vector(mode = "list")
for (i in 1:length(Ident.list)) {
  if (grepl("-pos", Ident.list[i]) == "FALSE") {
    neuromast.list <- append(neuromast.list, Ident.list[i])
   }
}

obj_integrated_filtered <- subset(obj_integrated, idents = neuromast.list)

levels(Idents(obj_integrated_filtered))

DimPlot(obj_integrated_filtered)

# =========================================================== Recluster after subsetting  ===========================================================
#obj_integrated_filtered <- ScaleData(obj_integrated_filtered, verbose = TRUE, vars.to.regress = "nCount_RNA")
obj_integrated_filtered <- RunPCA(obj_integrated_filtered, npcs = 100, verbose = TRUE, features = NULL)
obj_integrated_filtered <- FindNeighbors(obj_integrated_filtered, dims = 1:25, verbose = TRUE)
obj_integrated_filtered <- FindClusters(obj_integrated_filtered, resolution = 1.0, verbose = TRUE)
obj_integrated_filtered <- RunUMAP(obj_integrated_filtered, reduction = "pca", dims = 1:25, verbose = TRUE)
DimPlot(obj_integrated_filtered, group.by = "data.set")
DimPlot(obj_integrated_filtered)
