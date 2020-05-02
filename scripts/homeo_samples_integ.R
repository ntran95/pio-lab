library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

if (TRUE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "homeo_samples_integ"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}

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

# =========================================================== UMAP/Clustering
if (!SC_transform) {
  obj_integrated <- ScaleData(obj_integrated, verbose = TRUE, vars.to.regress = "nCount_RNA")
  obj_integrated <- RunPCA(obj_integrated, npcs = 100, verbose = TRUE, features = NULL)
  obj_integrated <- FindNeighbors(obj_integrated, dims = 1:15)
  obj_integrated <- FindClusters(obj_integrated, resolution = 1.2)
  obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:15)
  
}

saveRDS(object = obj_integrated, file = "../data/homeo_samples_integ.RDS")

readRDS("../data/homeo_samples_integ.RDS")
# =========================================================== Sample L29314


# =========================================================== Sample L34727


# =========================================================== Sample L34728





