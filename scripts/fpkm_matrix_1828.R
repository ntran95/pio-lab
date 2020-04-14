setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts/")
script_name <- "fpkm_matrix_1828"

readRDS("../data/fpkm_matrix_1828.RDS")

#reading in .rds, contains the expression matrix 
#variable fpkm_matrix_1828 contains 32489 (features) x 649 (cells) sparse Matrix of class "dgCMatrix"
#similarly to Read10X(data.dir = "../data/etc/")
fpkm_matrix_1828.data <- readRDS("../data/fpkm_matrix_1828.RDS")

fpkm_matrix_1828_smartseq2 <- CreateSeuratObject(counts = fpkm_matrix_1828.data, project = "fpkm_smartseq2", min.cells = 1, min.features = 1)

fpkm_matrix_1828_smartseq2
#26320 features across 649 samples within 1 assay 

# convert rownames of fpkm expression matrix ()
#use test variable first
fpkm_expression_mtx <- readRDS("../data/fpkm_matrix_1828.RDS")
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
as.matrix(gene_table.m)
row.names(gene_table) <- gene_table$Gene.stable.ID
#load gene table
#fpkm_expression_mtx.merge <- merge(fpkm_expression_mtx, gene_table, by.x = rownames(fpkm_expression_mtx))

#matched <- cbind(fpkm_expression_mtx, gene_table[, "Gene.name.uniq"][match(rownames(fpkm_expression_mtx), rownames(gene_table))])

new.mtx <- merge(fpkm_expression_mtx,gene_table["Gene.name.uniq"],by="row.names",all.x=TRUE)

#check to see if gene symbol [column 651] matches ENSGARG ensembl id (row 1). gene symbol will be in last column
#32489x 651
new.mtx[1,c(1,651)]

new.mtx[3, c(1,651)]

#turn Gene.name.uniq into rownames
new.mtx <- data.frame(new.mtx, row.names = "Gene.name.uniq")

#delete ENSDARG column
new.mtx$Row.names <- NULL



#add treatments to metadata (1hr and homeo)