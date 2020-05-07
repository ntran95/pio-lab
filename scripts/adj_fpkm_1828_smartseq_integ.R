library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(stringr)
options(future.globals.maxSize = 8000 * 1024^2)
library(RColorBrewer)


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

# =========================================================== Change Default Assay  ================================
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

# =========================================================== Adjust object integrated metadata

obj_integrated@meta.data$data.set <- factor(obj_integrated@meta.data$data.set, ordered = TRUE)

#rename homeo samples L##### to LIMS order number
meta_obj_integ <- obj_integrated@meta.data

data.set <- as.vector(obj_integrated@meta.data$data.set)

meta_obj_integ <- meta_obj_integ%>% mutate(data.set =case_when(str_detect(data.set, "homeo-L29314") ~ "homeo-2047", 
                                                               str_detect(data.set, "homeo-L34727") ~ "homeo-2410-7",
                                                               str_detect(data.set, "homeo-L34728") ~ "homeo-2410-8",
                                                               TRUE ~ as.vector(obj_integrated@meta.data$data.set)))

obj_integrated@meta.data$data.set <- meta_obj_integ$data.set

# =========================================================== UMAP/Clustering ===========================================================

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
  print(DimPlot(obj_integrated, group.by = "data.set", cols = brewer.pal(n = 6, name = "RdBu")) + DarkTheme()) 
  dev.off()
  png(figurePath(paste0("umap.unlabelled.PC",pc,".png"))
      ,width = 11, height = 9, units = "in", res = 300)
  print(DimPlot(obj_integrated))
  dev.off()
  png(figurePath(paste0("umap.by.treatment.PC", pc,".png"))
      ,width = 11, height = 9, units = "in", res = 300)
  print(DimPlot(obj_integrated, group.by = "treatment") + DarkTheme())
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

# =========================================================== PC 25 Annotating ===========================================================
obj_integrated <- FindNeighbors(obj_integrated, dims = 1:25, verbose = TRUE)
obj_integrated <- FindClusters(obj_integrated, resolution = 1.0, verbose = TRUE)
obj_integrated <- RunUMAP(obj_integrated, reduction = "pca", dims = 1:25, verbose = TRUE)
DimPlot(obj_integrated, group.by = "data.set")
print(DimPlot(obj_integrated, label = TRUE))

#saveRDS(object = obj_integrated, file = paste0("../data/", script_name,".RDS"))

meta_common_features <- read.table(file = "../data/gene-lists/meta_common_features.tsv", sep = "", header = T)

gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)

meta_common_features <- read.table(file = "../data/gene-lists/meta_common_features.tsv", sep = "", header = T)


# =========================================================== Annotate -pos clusters ===========================================================
#desired unlabelled clusters passed through 
pos_list <- c(13,20,21,22,23,24,25)
pos_clusters <- vector(mode = "list", length = length(pos_list))

for (x in 1:length(pos_list)){
  print(paste0("finding markers for cluster: ", pos_list[x]))
  pos_clusters[[x]] <- FindMarkers(obj_integrated, ident.1 = pos_list[x], only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10, verbose = TRUE)
  print(pos_clusters[[x]])
  pos_clusters[[x]]["Gene.name.uniq"] <- row.names(pos_clusters[[x]])
  pos_clusters[[x]] <- merge(pos_clusters[[x]], gene_table, by = "Gene.name.uniq")
}

names(pos_clusters) <- paste("cluster_", pos_list)

# =========================================================== Annotate neuromast clusters
meta <- obj_integrated@meta.data

cluster.ident <- strsplit(unique(na.omit(meta$cell.type.ident)), split = "[][']|,\\s*")

cluster.ident.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)){
  cluster.ident.list[[i]] <- meta %>% filter(cell.type.ident == cluster.ident[[i]]) %>% count(seurat_clusters, sort =TRUE)
  print(cluster.ident[[i]])
}

names(cluster.ident.list) <- cluster.ident



colnames(meta)
cells <- list("mature-HCs" = c(6,12), "early-HCs" = c(17,19),  "HC-prog" = 16,
              "central-cells" = c(0,1,4), "D/V-cells" = c(9,8), "A/P-cells" = 7,
              "amplfying support" = c(11,18), "mantle-cells" = c(3,5), "interneuromast" = c(2,10,14,15),"col1a1b-pos" = 13,
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


saveRDS(object = obj_integrated, file = paste0("../data/", script_name,".RDS"))
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


