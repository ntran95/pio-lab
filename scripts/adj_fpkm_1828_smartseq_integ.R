library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)
options(future.globals.maxSize = 8000 * 1024^2)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "adj_fpkm_1828_smartseq_integ"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
  
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

# =========================================================== add homeo 10X samples ===========================

homeo_samples_integ <- readRDS("../data/homeo_samples_integ.RDS")

#search for NAs
grep("TRUE",is.na(colnames(homeo_samples_integ)))

# =========================================================== Seurat Object List ================================
#split homeo data by samples (c("isl1_10X","L29314", "L34727", "L34728"))
split_homeo_samples_integ <- SplitObject(homeo_samples_integ, split.by = "data.set")
#split smartseq by treatment (c(1hr-smrtseq, homeo-smrtseq))
split_fpkm_smartseq <- SplitObject(fpkm_matrix_1828_smartseq2, split.by = "data.set")

seurat_obj_list <- c(split_homeo_samples_integ, split_fpkm_smartseq)

all_shared_genes <- lapply(seurat_obj_list, row.names) %>% Reduce(intersect, .) 

length(all_shared_genes)

length(seurat_obj_list)

for (i in 1:length(seurat_obj_list)){
  #get current assay
  #print(DefaultAssay(object = seurat_obj_list[[i]]))
  current_assay <-DefaultAssay(object = seurat_obj_list[[i]])
  #print(current_assay)
  #change "integrated" assay back to "RNA" 
  if (current_assay == "integrated"){
    #print(current_assay)
    DefaultAssay(seurat_obj_list[[i]]) <- "RNA"
  }
  print(DefaultAssay(object = seurat_obj_list[[i]]))
}

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
    seurat_obj_list[[i]] <- NormalizeData(seurat_obj_list[[i]], verbose = TRUE)
    seurat_obj_list[[i]] <- FindVariableFeatures(seurat_obj_list[[i]], selection.method = "vst",nfeatures = 2000, verbose = TRUE)
  }
  obj_anchors <- FindIntegrationAnchors(object.list = seurat_obj_list,
                                        dims = dims, reference = 1) # 1 is the homeostatic dataset in obj list
  obj_integrated <- IntegrateData(anchorset = obj_anchors, dims = dims, features.to.integrate = all_shared_genes)
  DefaultAssay(obj_integrated) <- "integrated"
  
}

obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE)

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
  print(DimPlot(obj_integrated, group.by  = "data.set") + DarkTheme())
  dev.off()
  png(figurePath(paste0("umap.unlabelled.PC",pc,".png"))
      ,width = 11, height = 9, units = "in", res = 300)
  print(DimPlot(obj_integrated))
  dev.off()
  integrated_featplt <- FeaturePlot(obj_integrated, common_features,
                                    reduction = "umap", pt.size = 0.25, combine = FALSE, label = TRUE)
  for (i in 1:length(integrated_featplt)) {
    integrated_featplt[[i]] <- integrated_featplt[[i]] + NoLegend() + NoAxes()
  }
  png(figurePath(paste0("common_features.PC",pc,".png")), width = 40,
      height = 80, units = "in", res = 200)
  print(cowplot::plot_grid(plotlist = integrated_featplt, ncol = 4))
  dev.off()
  
}
