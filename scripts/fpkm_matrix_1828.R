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
options(future.globals.maxSize = 5000 * 1024^2)


if (FALSE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "fpkm_matrix_1828"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}
  
  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

# =========================================================== Load smartseq2 dataset

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

# =========================================================== Create Seurat Object############################
fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = new_fpkm_expression_mtx, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

# =========================================================== Add percent.mt#####################
fpkm_matrix_1828_smartseq2[["percent.mt"]] <- PercentageFeatureSet(fpkm_matrix_1828_smartseq2, pattern = "^mt-")


# =========================================================== add treatments to metadata (1hr and homeo)###########
meta_smartseq2 <- fpkm_matrix_1828_smartseq2@meta.data

barcode <- rownames(meta_smartseq2)
#specify treatment and sequencing technique (data.set)
meta_smartseq2 <- meta_smartseq2%>% mutate(data.set =case_when(str_detect(barcode, "homeo") ~ "homeo-smrtseq", 
                                                                str_detect(barcode, "1hr") ~ "1hr-smrtseq",
                                                                TRUE ~ barcode)) %>% mutate(treatment =case_when(str_detect(barcode, "homeo") ~ "homeo", 
                                                                                                                str_detect(barcode, "1hr") ~ "1hr",
                                                                                                                TRUE ~ barcode))
#get barcode back as rownames
rownames(meta_smartseq2) <- rownames(fpkm_matrix_1828_smartseq2@meta.data)

#add treatment to seurat metadata
fpkm_matrix_1828_smartseq2@meta.data$treatment<- meta_smartseq2$treatment
#add sequencing technique to seurat metadata
fpkm_matrix_1828_smartseq2@meta.data$data.set<- meta_smartseq2$data.set
#check
fpkm_matrix_1828_smartseq2@meta.data
# =========================================================== Normalize smartseq2 =========================
#note: normalization step for smartseq2 object is skipped, only specifically for this step
# =========================================================== Add 10X Homeo data ==========================

#isl1_sib_10X.data <- Read10X(data.dir = "../data/isl1_sib_counts_10X")

#homeo.isl1_sib_10X <- CreateSeuratObject(counts = isl1_sib_10X.data, project = "homeo.isl1.sib.10X", min.cells = 1, min.features = 1)

#homeo.isl1_sib_10X[["percent.mt"]] <- PercentageFeatureSet(homeo.isl1_sib_10X, pattern = "^mt-")

#homeo.isl1_sib_10X <- subset(homeo.isl1_sib_10X, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA <10000) 

#homeo.isl1_sib_10X <- NormalizeData(homeo.isl1_sib_10X, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)

#homeo.isl1_sib_10X

load("../data/workspace_homeo_isl1_sib_10X.RData")

# =========================================================== Modify 10X metadata ===========================

homeo.isl1_sib_10X@meta.data$data.set <- "homeo-10X"
homeo.isl1_sib_10X@meta.data$treatment <- "homeo"

#check
head(homeo.isl1_sib_10X@meta.data)

# =========================================================== Seurat Object List ================================

split_fpkm_smartseq <- SplitObject(fpkm_matrix_1828_smartseq2, split.by = "data.set")
seurat_obj_list <- c(homeo.isl1_sib_10X, split_fpkm_smartseq[1], split_fpkm_smartseq[2])
#find common/shared genes among the three seurat objects (for IntegrateData())
all_shared_genes <- lapply(seurat_obj_list, row.names) %>% Reduce(intersect, .) 

all_shared_genes

# =========================================================== Integration Loop ==================================

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
    object.list = seurat_obj_list, anchor.features = all_shared_genes)
  # reference = 1 is the homeostatic dataset in obj list
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
                                        anchor.features = seurat_obj_features, normalization.method = "SCT",
                                        dims = dims, reference = 1) 
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims,
                                  normalization.method = "SCT")
} else {
  for (i in 1:length(seurat_obj_list)) {
    seurat_obj_list[[i]] <- NormalizeData(seurat_obj_list[[i]], verbose = FALSE)
    seurat_obj_list[[i]] <- FindVariableFeatures(seurat_obj_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = FALSE)
  }
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
                                        dims = dims, reference = 1) # 1 is the homeostatic dataset in obj list
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims, features.to.integrate = all_shared_genes)
  DefaultAssay(obj_integrated) <- "integrated"
}

obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE)

# =========================================================== UMAP/Clustering
if (!SC_transform) {
  obj_integrated <- ScaleData(obj_integrated, verbose = FALSE)
  obj_integrated <- RunPCA(obj_integrated, npcs = 100, verbose = TRUE, features = NULL)
  obj_integrated <- FindNeighbors(obj_integrated, dims = 1:15)
  obj_integrated <- FindClusters(obj_integrated, resolution = 1.2)
  obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:15)
  
}

obj_integrated_by_trt_umap <- DimPlot(obj_integrated, group.by  = "treatment", pt.size = 0.4) + DarkTheme() 

png(figurePath("integrated_by_trt.png"),
    width = 11, height = 9, units = "in", res = 300)
print(obj_integrated_by_trt_umap)
dev.off()

obj_integrated_by_dataset_umap <- DimPlot(obj_integrated, group.by  = "data.set", pt.size = 0.4) + DarkTheme() 

png(figurePath("integrated_by_dataset.png"),
    width = 11, height = 9, units = "in", res = 300)
print(obj_integrated_by_dataset_umap)
dev.off()

obj_integrated_unlabeled <- DimPlot(obj_integrated, label = TRUE, pt.size = 0.4) + NoLegend()

png(figurePath("integrated_umap_unlabelled.png"),
    width = 11, height = 9, units = "in", res = 300)
print(obj_integrated_unlabeled)
dev.off()


# ===========================================================  Common Features # =========================================================== 
common_features <- scan(paste0("/Volumes/easystore/SIMR_2019/pio-lab/data/gene-lists/common_neuromast_features.txt"), what = "character")

marker.ident <- c("early-HC", "early-HC", "early-HC", "mature-HC", "early-HC",
                  "prog-HC", "prog-HC", "prog-HC", "prog-HC", "early-HC",
                  "prog-HC", "prog-HC", "prog-HC", "prog-HC", "prog-HC",
                  "prog-HC", "early-HC", "D/V Pole", "D/V Pole", "Mantle",
                  "Mantle", "Mantle", "Amplifying Support", "Amplifying Support", "Central",
                  "Central", "Central", "Central", "Central", "Central",
                  "AP Pole", "AP Pole", "AP Pole", "Common neuromast cell types", "Interneuromast",
                  "Interneuromast/Mantle")
common_features_meta <- data.frame(common_features, marker.ident)

integrated_featplt <- FeaturePlot(obj_integrated, common_features,
                                  reduction = "umap", pt.size = 0.25, combine = FALSE, label = TRUE)
for (i in 1:length(integrated_featplt)) {
  integrated_featplt[[i]] <- integrated_featplt[[i]] + NoLegend() + NoAxes()
}
png(figurePath("common_features.png"), width = 40,
    height = 80, units = "in", res = 200)
print(cowplot::plot_grid(plotlist = integrated_featplt, ncol = 4))
dev.off()

################################VlnPlot#########################################
integrated_vln <- VlnPlot(obj_integrated, common_features,pt.size = 0.0, combine = FALSE, log = TRUE)
for (i in 1:length(integrated_vln)) {
  integrated_vln[[i]] <- integrated_vln[[i]] + NoLegend() + NoAxes()
}
png(figurePath("vln_common_features.png"), width = 40,
    height = 80, units = "in", res = 200)
print(cowplot::plot_grid(plotlist = integrated_vln, ncol = 4))
dev.off()

# ===========================================================  FindAllMarkers ===========================================
#plan("multiprocess")
all.markers.integrated <- FindAllMarkers(obj_integrated, only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.01, return.thresh = 0.001, verbose = TRUE)

#change the column name formerally called "gene" generated from FindAllMarkers() to  "Gene.name.uniq" so we can merge with gene_table
#must have same column name in order to match genes
colnames(all.markers.integrated)[7] <- "Gene.name.uniq"

#inspect new col name change
names(all.markers.integrated)

#merge gene table with all markers generated based on uniq gene symbol
all.markers.integrated <- merge(all.markers.integrated, gene_table, by = "Gene.name.uniq") 
####Rename Clusters####
meta_integrated <- obj_integrated@meta.data
colnames(meta_integrated)
cells <- list("mature-HCs" = c(9,13), "early-HCs" = c(8,12),  "HC-prog" = 14,
              "central-cells" = c(7,4,1), "DV/AP Pole" = 3,
              "amplfying support" = c(10,16), "mantle-cells" = c(0,5), "col1a1b-pos" = 11,
              "c1qtnf5-pos" = 17, "interneuromast" = c(2,6), "apoa1b-pos" = 20, "mfap4-pos" = 19, "krt91-pos" = 18, "clec14a" =15)
meta_integrated$cell.type.ident <- factor(rep("", nrow(meta_integrated)),
                                          levels = names(cells), ordered = TRUE)
for (i in 1:length(cells)) {
  meta_integrated$cell.type.ident[meta_integrated$seurat_clusters %in% cells[[i]]] <- names(cells)[i]
}
obj_integrated@meta.data <- meta_integrated
Idents(obj_integrated) <- obj_integrated@meta.data$cell.type.ident

integrated.umap.labeled <- DimPlot(obj_integrated, reduction = "umap", label = TRUE, pt.size= 0.4) + NoLegend()

####print labeled plot####
png(figurePath("annotated_umap_clusters.png"),
    width = 11, height = 9, units = "in", res = 300)
print(integrated.umap.labeled)
dev.off()