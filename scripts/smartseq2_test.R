setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts/")
script_name <- "fpkm_matrix_1828"
readRDS("../data/fpkm_matrix_1828/fpkm_matrix_1828.RDS")
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
fpkm_expression_mtx <- readRDS("../fpkm_matrix_1828/data/fpkm_matrix_1828.RDS")
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

#rownames(fpkm_expression_mtx) <- gsub("[:]", "_", rownames(fpkm_expression_mtx))
#rownames(fpkm_expression_mtx) <- gsub("[*]", "-", rownames(fpkm_expression_mtx))
#rownames(fpkm_expression_mtx) <- gsub("[_]", "", rownames(fpkm_expression_mtx))


gene_vector <- rownames(fpkm_expression_mtx)
is.null(rownames(fpkm_expression_mtx))

fpkm_expression_mtx <- as.sparse(fpkm_expression_mtx)

#export
write.table(fpkm_expression_mtx, file="../data/fpkm_matrix_1828/corr_gene_fpkm_matrix_1228.tsv", col.names=TRUE, row.names = TRUE, sep = "\t")

counts.mtx <- read.table("../data/fpkm_matrix_1828/corr_gene_fpkm_matrix_1228.tsv", sep = "\t", header = TRUE)

saveRDS(fpkm_expression_mtx, file = "../data/fpkm_matrix_1828/corr_gene_fpkm_matrix_1228.RDS")
x <- read.table(file = system.file("../data/fpkm_matrix_1828/corr_gene_fpkm_matrix_1228.tsv", package = 'Seurat'),as.is = TRUE)

###############Create Seurat Object###############
fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = as.sparse(fpkm_expression_mtx), project = "fpkm_smartseq2", min.cells = 1, min.features = 1)
fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 



###########################assign ensdarg id to homeo data##############
isl1_sib_10X_count_matrix <- as.matrix(readRDS("../data/isl1_sib_counts_10X/isl1_sib_counts_10_matrix.RDS"))
head(rownames(isl1_sib_10X_count_matrix))
dim(isl1_sib_10X_count_matrix)


#isl1_sib_10X_count_matrix<- merge(isl1_sib_10X_count_matrix,gene_table["Gene.stable.ID"],by="row.names",all.x=TRUE, sort = FALSE)

head(rownames(isl1_sib_10X_count_matrix))

head(isl1_sib_10X_count_matrix$Gene.stable.ID)

row.names(isl1_sib_10X_count_matrix) <- isl1_sib_10X_count_matrix$Gene.stable.ID

isl1_sib_10X_count_matrix$Row.names <- NULL
isl1_sib_10X_count_matrix$Gene.name.uniq <- NULL

dim(isl1_sib_10X_count_matrix)


######
load("../data/workspace_homeo_isl1_sib_10X.RData")

isl1_info <- gene_table[gene_table$Gene.name.uniq %in% rownames(isl1_sib_10X_count_matrix),]
new_isl1_sib_10X_matrix <- isl1_sib_10X_count_matrix[match(isl1_info$Gene.stable.ID, rownames(isl1_sib_10X_count_matrix)),]
rownames(new_isl1_sib_10X_matrix) <- isl1_info$Gene.stable.ID

saveRDS(new_isl1_sib_10X_matrix, file = "../data/isl1_sib_counts_10_matrix.RDS")

isl1_test <- CreateSeuratObject(counts = new_isl1_sib_10X_matrix, project = "isl1_sib_homeo", min.cells = 0, min.features = 0)

