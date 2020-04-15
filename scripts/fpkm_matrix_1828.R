setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts/")

library(Seurat)
# library(plotly)
#install.packages("future")
# library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
# library(gridExtra)
# library(ggrepel)
library(ggrepel)
library(stringr)
script_name <- "fpkm_matrix_1828"

readRDS("../data/fpkm_matrix_1828.RDS")

#reading in .rds, contains the expression matrix 
#variable fpkm_matrix_1828 contains 32489 (features) x 649 (cells) sparse Matrix of class "dgCMatrix"
#similarly to Read10X(data.dir = "../data/etc/")
#keep original matrix/object for reference
#fpkm_matrix_1828.data <- readRDS("../data/fpkm_matrix_1828.RDS")

#fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = fpkm_matrix_1828.data, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

#fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 


###########convert rownames of fpkm expression matrix ()###########
#load gene table
fpkm_expression_mtx <- readRDS("../data/fpkm_matrix_1828.RDS")
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
as.matrix(gene_table)
row.names(gene_table) <- gene_table$Gene.stable.ID
#merge
fpkm_expression_mtx <- merge(fpkm_expression_mtx,gene_table["Gene.name.uniq"],by="row.names",all.x=TRUE, sort = FALSE)

#check to see if gene symbol [column 651] matches ENSGARG ensembl id (row 1). gene symbol will be in last column
#32489x 651
fpkm_expression_mtx[1,c(1,651)]

fpkm_expression_mtx[3, c(1,651)]

#turn Gene.name.uniq into rownames, automatically delete the colum
fpkm_expression_mtx <- data.frame(fpkm_expression_mtx, row.names = "Gene.name.uniq")

#delete ENSDARG column
fpkm_expression_mtx$Row.names <- NULL
#turn back to matrix
fpkm_expression_mtx <- as.matrix(fpkm_expression_mtx)
dim(fpkm_expression_mtx)
#dim 32489   649

###############Create Seurat Object###############

fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = fpkm_expression_mtx, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 

#############Add percent.mt#####################
fpkm_matrix_1828_smartseq2[["percent.mt"]] <- PercentageFeatureSet(fpkm_matrix_1828_smartseq2, pattern = "^mt-")


############add treatments to metadata (1hr and homeo)###########
meta_smartseq2 <- fpkm_matrix_1828_smartseq2@meta.data

barcode <- rownames(meta_smartseq2)
meta_smartseq2 <- meta_smartseq2%>% mutate(treatment =case_when(str_detect(barcode, "homeo") ~ "homeo", 
                                                                str_detect(barcode, "1hr") ~ "1hr ",
                                                                TRUE ~ barcode))
#get barcode back as rownames
rownames(meta_smartseq2) <- rownames(fpkm_matrix_1828_smartseq2@meta.data)

#add to seurat metadata
fpkm_matrix_1828_smartseq2@meta.data$treatment<- meta_smartseq2$treatment

##########Load isl1_sib_10X_homeo Data#############
load("../data/workspace_homeo_isl1_sib_10X.RData")

##########Seurat Object List#################
#isl1_sib_10X_homeo is in first index
seurat_obj_list <- c(homeo.isl1_sib_10X, fpkm_matrix_1828_smartseq2)
dims = 1:15

##############Integration Loop###############
##############Integration Loop###############
SC_transform <- FALSE
if (SC_transform) {
  for (i in 1:length(seurat_obj_list)) {
    print("hi")
    seurat_obj_list[[i]] <- SCTransform(
      seurat_obj_list[[i]], verbose = FALSE)
  }
  seurat_obj_features <- SelectIntegrationFeatures(
    object.list = seurat_obj_list, nfeatures = 2000)
  seurat_obj_list <- PrepSCTIntegration(
    object.list = seurat_obj_list, anchor.features = seurat_obj_features)
  # reference = 1 is the homeostatic dataset in obj list
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
                                        anchor.features = seurat_obj_features, normalization.method = "SCT",
                                        dims = dims, reference = 1) 
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims,
                                  normalization.method = "SCT")
} else {
  for (i in 1:length(seurat_obj_list)) {
    #seurat_obj_list[[i]] <- NormalizeData(seurat_obj_list[[i]], verbose = FALSE)
    seurat_obj_list[[i]] <- FindVariableFeatures(seurat_obj_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
  }
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
                                        dims = dims, reference = 1) # 1 is the homeostatic dataset in obj list
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims)
  DefaultAssay(obj_integrated) <- "integrated"
}
#obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE, levels = ids)
if (FALSE) {
  saveRDS(obj_integrated, dataPath(
    paste0("SeurObj_before_clust", "_", script_name,"_.RDS")))
  obj_integrated <- readRDS(
    dataPath(paste0("SeurObj_before_clust", "_", script_name,"_.RDS")))
}
# =========================================================== UMAP/Clustering
if (!SC_transform) {
  obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
}
