####Set up environment####
#install.packages("Seurat")
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
setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts/")

if (FALSE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "fpkm_matrix_1828"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}
#reading in .rds, contains the expression matrix 
#variable fpkm_matrix_1828 contains 32489 (features) x 649 (cells) sparse Matrix of class "dgCMatrix"
#similarly to Read10X(data.dir = "../data/etc/")
#keep original matrix for reference
#keep original matrix/object for reference
fpkm_matrix_1828.data <- readRDS("../data/fpkm_matrix_1828.RDS")

fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = fpkm_matrix_1828.data, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 


###########convert rownames of fpkm expression matrix ()###########
#load gene table
fpkm_expression_mtx <- readRDS("../data/fpkm_matrix_1828/fpkm_matrix_1828.RDS")
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
as.matrix(gene_table)
row.names(gene_table) <- gene_table$Gene.stable.ID
#merge
fpkm_expression_mtx <- merge(fpkm_expression_mtx,gene_table["Gene.name.uniq"],by="row.names",all.x=TRUE, sort = FALSE)

#check to see if gene symbol [column 651] matches ENSGARG ensembl id (row 1). gene symbol will be in last column
#32489x 651
fpkm_expression_mtx[1,c(1,651)]
fpkm_expression_mtx[3, c(1,651)]
#turn Gene.name.uniq into rownames
rownames(fpkm_expression_mtx) <- fpkm_expression_mtx$Gene.name.uniq
#delete ENSDARG column
fpkm_expression_mtx$Row.names <- NULL
fpkm_expression_mtx$Gene.name.uniq <- NULL
#turn back to matrix
fpkm_expression_mtx <- as.matrix(fpkm_expression_mtx)
dim(fpkm_expression_mtx)
#dim 32489   649

rownames(fpkm_expression_mtx) <- gsub("[_]", "-", rownames(fpkm_expression_mtx))
colnames(fpkm_expression_mtx) <- gsub("[.]", "-", colnames(fpkm_expression_mtx))


#check if rownames is empty
which(rownames(x = fpkm_expression_mtx) == '')
#[1] 11478
print(rownames(fpkm_expression_mtx)[11478])
#[1] ""
#delete
new_fpkm_expression_mtx <- fpkm_expression_mtx[-11478,]
#turn to sparse matrix
fpkm_expression_mtx <- as.sparse(fpkm_expression_mtx)

#check
dim(new_fpkm_expression_mtx)
#[1] 32488   649

which(rownames(x = new_fpkm_expression_mtx) == '')
#integer(0)

#####################Create Seurat Object############################
fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = new_fpkm_expression_mtx, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

#############Add percent.mt#####################
fpkm_matrix_1828_smartseq2[["percent.mt"]] <- PercentageFeatureSet(fpkm_matrix_1828_smartseq2, pattern = "^mt-")

############add treatments to metadata (1hr and homeo)###########
meta_smartseq2 <- fpkm_matrix_1828_smartseq2@meta.data

barcode <- rownames(meta_smartseq2)
meta_smartseq2 <- meta_smartseq2%>% mutate(treatment =case_when(str_detect(barcode, "homeo") ~ "homeo-smrtseq", 
                                                                str_detect(barcode, "1hr") ~ "1hr-smrtseq",
                                                                TRUE ~ barcode))
#get barcode back as rownames
rownames(meta_smartseq2) <- rownames(fpkm_matrix_1828_smartseq2@meta.data)

#add to seurat metadata
fpkm_matrix_1828_smartseq2@meta.data$treatment<- meta_smartseq2$treatment


fpkm_matrix_1828_smartseq2@meta.data


##############Find Variable Features###########################
fpkm_matrix_1828_smartseq2 <- FindVariableFeatures(fpkm_matrix_1828_smartseq2, selection.method = "vst",nfeatures = 2000, verbose = FALSE)

####Clustering/UMAP####
fpkm_matrix_1828_smartseq2 <- ScaleData(fpkm_matrix_1828_smartseq2, features = NULL, vars.to.regress = NULL, verbose = TRUE)

fpkm_matrix_1828_smartseq2 <- RunPCA(fpkm_matrix_1828_smartseq2, npcs = 100, verbose = TRUE, features = NULL)

fpkm_matrix_1828_smartseq2 <- FindNeighbors(fpkm_matrix_1828_smartseq2, dims = 1:15, features = NULL)

fpkm_matrix_1828_smartseq2 <- FindClusters(fpkm_matrix_1828_smartseq2, resolution = 1.0)

fpkm_matrix_1828_smartseq2 <- RunUMAP(fpkm_matrix_1828_smartseq2, dims = 1:15, reduction = "pca")

fpkm_smart_seq2_unlabel_cluster <- DimPlot(fpkm_matrix_1828_smartseq2, reduction = "umap")

head(fpkm_matrix_1828_smartseq2@meta.data)

png(figurePath("fpkm_unlabeled_umap.png"),
    width = 11, height = 9, units = "in", res = 300)
print(fpkm_smart_seq2_unlabel_cluster)
dev.off()
######Load########
load("../data/workspace_homeo_isl1_sib_10X.RData")

#######change treatment column of 10X Homeo#########
homeo.isl1_sib_10X@meta.data$treatment <- "homeo-10X"
#check
head(homeo.isl1_sib_10X@meta.data)


##############Integration#####################

#seurat_obj_list <- c(fpkm_matrix_1828_smartseq2, homeo.isl1_sib_10X)
split_fpkm_smartseq <- SplitObject(fpkm_matrix_1828_smartseq2, split.by = "treatment")
seurat_obj_list <- c(homeo.isl1_sib_10X, split_fpkm_smartseq[1], split_fpkm_smartseq[2])
dims = 1:15
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
  ("hi")
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, dims = dims, reference = 1, verbose = TRUE) # 1 is the homeostatic dataset in obj list
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims, features.to.integrate	= rownames(homeo.isl1_sib_10X))
  DefaultAssay(obj_integrated) <- "integrated"
  print("hi")
}


obj_integrated@meta.data$treatment <- factor(obj_integrated@meta.data$treatment, ordered = TRUE)
if (FALSE) {
  saveRDS(obj_integrated, "../data/fpkm_matrix_1828/SeurObj_before_clust_smartseq_10X.RDS")
  #obj_integrated <- readRDS(dataPath(paste0("SeurObj_before_clust", "_", script_name,"_.RDS")))
}

# =========================================================== UMAP/Clustering
if (!SC_transform) {
  obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
  obj_integrated <- RunPCA(obj_integrated, npcs = 100, verbose = TRUE, features = NULL)
  obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:15)
  obj_integrated <- FindNeighbors(obj_integrated, dims = 1:25)
  obj_integrated <- FindClusters(obj_integrated, resolution = 1.2)
  
}

####################Elbow######################
ElbowPlot(obj_integrated, ndims = 40)

####################JackStraw##################
obj_integrated <- JackStraw(obj_integrated, num.replicate = 100, dims = 35)
obj_integrated <- ScoreJackStraw(obj_integrated, dims = 1:35)
JackStrawPlot(obj_integrated, reduction = "pca", dims = 1:35)


###################DImplit#######################
obj_integrated_by_trt_umap <- DimPlot(obj_integrated, group.by  = "treatment", label = TRUE)

png(figurePath("integrated_by_trt.png"),
    width = 11, height = 9, units = "in", res = 300)
print(obj_integrated_by_trt_umap)
dev.off()

DimPlot(obj_integrated, label = TRUE)

save.image("../data/post-integrated-SeurObj.RData")

