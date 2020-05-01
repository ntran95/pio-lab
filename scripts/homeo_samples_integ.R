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
  png(figurePath(paste0(sample.ls[i],"QC-vln.png")),
      width = 11, height = 9, units = "in", res = 300)
  print(VlnPlot(obj_list[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()
  print(summary(obj_list[[i]]$nFeature_RNA))
  print(summary(obj_list[[i]]$percent.mt))
  #print("The original object summary: ", obj_list[i])
  #obj_list[[i]] <- subset(obj_list[[i]], subset = nCount_RNA < 15000 & percent.mt < 10)
  #print(obj_list[[i]])
}

obj_list <- lapply(obj_list, subset, subset = subset = nCount_RNA < 15000 & percent.mt < 10)

obj_list$L29314

#rename seur objects in obj_list
L29314_homeo <- obj_list$L29314
L34727_homeo <- obj_list$L34727
L34728_homeo <- obj_list$L34728

#calculate percent.mt for each obj in list
percent.mt <-lapply(obj_list, PercentageFeatureSet, pattern = "^mt-")
obj_list$L29314[["percent.mt"]] <- percent.mt$L29314
obj_list$L34727[["percent.mt"]] <- percent.mt$L34727
obj_list$L34728[["percent.mt"]] <- percent.mt$L34728

# =========================================================== Normalize/Find Variable Features

for (i in 1:length(obj_list)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]], verbose = TRUE)
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = TRUE)
}

# =========================================================== Sample L29314

L29314_homeo@meta.data
VlnPlot(L29314_homeo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(L29314_homeo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(L29314_homeo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# =========================================================== Sample L34727


# =========================================================== Sample L34728



