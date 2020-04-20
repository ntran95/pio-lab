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
fpkm_matrix_1828.data <- readRDS("../data/fpkm_matrix_1828/fpkm_matrix_1828.RDS")

fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = fpkm_matrix_1828.data, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 


###########convert rownames of fpkm expression matrix ()###########
#load gene table
fpkm_expression_mtx <- readRDS("../data/fpkm_matrix_1828/fpkm_matrix_1828.RDS")
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
as.matrix(gene_table)
row.names(gene_table) <- gene_table$Gene.stable.ID
#add gene name column with "_" instead of "-"
#merge
new_fpkm_expression_mtx <- merge(fpkm_expression_mtx,gene_table["Gene.name.uniq"],by="row.names",all.x=TRUE, sort = FALSE)

#check to see if gene symbol [column 651] matches ENSGARG ensembl id (row 1). gene symbol will be in last column
#32489x 651
new_fpkm_expression_mtx[1,c(1,651)]

new_fpkm_expression_mtx[3, c(1,651)]

#turn Gene.name.uniq into rownames, automatically delete the colum
rownames(new_fpkm_expression_mtx) <- new_fpkm_expression_mtx$Gene.name.uniq

#delete ENSDARG column
new_fpkm_expression_mtx$Row.names <- NULL
new_fpkm_expression_mtx$Gene.name.uniq <- NULL
#turn back to matrix
new_fpkm_expression_mtx <- as.matrix(new_fpkm_expression_mtx)
dim(new_fpkm_expression_mtx)
#dim 32489   649

#replace "_" with "-" for gene names
#based on Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#colnames(new_fpkm_expression_mtx) <- gsub("[_]", "-", colnames(new_fpkm_expression_mtx))
rownames(new_fpkm_expression_mtx) <- gsub("[_]", "-", rownames(new_fpkm_expression_mtx))


fpkm_expression_mtx

#confirm that "_" replaced with "-"
grep("_", rownames(fpkm_expression_mtx))

###############Create Seurat Object###############

smartseq_fpkm <- CreateSeuratObject(counts = new_fpkm_expression_mtx, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

fpkm_matrix_1828_smartseq2corrected_gene
#Error: Feature names of counts matrix cannot be empty

################convert ensembl id to gene symbol after CreateSeuratObject################
counts.mtx <- merge(as.data.frame(fpkm_matrix_1828_smartseq2[["RNA"]]@counts),as.data.frame(gene_table["Gene.name.uniq"]),by="row.names",all.x=TRUE, sort = FALSE)
counts.mtx <- data.frame(counts.mtx, row.names = "Gene.name.uniq")
#delete ENSDARG column
counts.mtx$Row.names <- NULL
#turn back to matrix
counts.mtx <- as.matrix(counts.mtx)
dim(counts.mtx)
#26320   649
#fpkm_matrix_1828_smartseq2 <- SetAssayData(object = fpkm_matrix_1828_smartseq2, slot = "counts", new.data = counts.mtx)
dim(fpkm_matrix_1828_smartseq2[["RNA"]]@counts)

#replace "_" with "-" for gene names
#based on Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
corrected_counts_mtx <- gsub("[_]", "-",rownames(counts.mtx))

rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@counts) <- corrected_counts_mtx


####################Convert object[["RNA]]@data (normalized data) to gene symbol################
rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@data) <- corrected_counts_mtx

#check
fpkm_matrix_1828_smartseq2[["RNA"]]@data

#############################Convert object[["RNA]]@meta.features to gene symbol##############
rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@meta.features) <- corrected_counts_mtx


########Error: Feature names of counts matrix cannot be empty#######
#rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@counts) <- gsub("[-]", "_",rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@counts))
#rownames(fpkm_matrix_1828_smartseq2[["RNA"]]@counts)

#rownames(homeo.isl1_sib_10X[["RNA"]]@counts) <- gsub("[-]", "_",rownames(homeo.isl1_sib_10X[["RNA"]]@counts))
#rownames(homeo.isl1_sib_10X[["RNA"]]@counts)

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
SC_transform <- TRUE
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
    seurat_obj_list[[i]] <- NormalizeData(seurat_obj_list[[i]], verbose = FALSE)
    seurat_obj_list[[i]] <- FindVariableFeatures(seurat_obj_list[[i]], selection.method = "vst",nfeatures = length(rownames(seurat_obj_list[[i]])), verbose = FALSE)
    print("hi")
  }
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, dims = dims, reference = 1, verbose = TRUE) # 1 is the homeostatic dataset in obj list
  #obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims, features.to.integrate	= rownames(homeo.isl1_sib_10X))
  #DefaultAssay(obj_integrated) <- "integrated"
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

save.image("../data/pre_integrated_SeurObj.RData")

###########################assign ensdarg id to homeo data##############
homeo_counts_mtx <- as.matrix(homeo.isl1_sib_10X[["RNA"]]@counts)
new_homeo_counts_mtx <- merge(homeo_counts_mtx,gene_table["Gene.stable.ID"],by="row.names",all.x=TRUE, sort = FALSE)


