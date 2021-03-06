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
library(ggrepel)
library(readxl)
library(writexl)
options(future.globals.maxSize = 5000 * 1024^2)


if (TRUE) {
  #setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  setwd("/home/ntran2/bgmp/pio-lab/scripts")
  
  script_name <- "integ_DGE"
  
  figurePath <- function(filename, format){paste0(script_name, "_figures/", filename)}
}

# =========================================================== Load filtered integ smartseq and homeo object =================================================

obj_integrated_filtered <- readRDS("../data/filtered_adj_fpkm_1828_smartseq_integ.RDS")

# =========================================================== Test FindMarkers by treatments and seq method  =================================================
DefaultAssay(obj_integrated_filtered) <- "RNA"
obj_integrated_filtered$cell.type.treatment.method <- paste(Idents(obj_integrated_filtered), obj_integrated_filtered$treatment, obj_integrated_filtered$seq.method, sep = "_")
obj_integrated_filtered$cell.type.treatment.method
Idents(obj_integrated_filtered) <- "cell.type.treatment.method"
x <- FindMarkers(obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2", verbose = TRUE)
x$Gene.name.uniq <- rownames(x)
#merge gene table with all markers generated based on uniq gene symbol
x <- merge(x, gene_table, by = "Gene.name.uniq") 
#reorder df in ascending order based on avg_logFC (col#3)
x <- x[order( x[,3], decreasing = TRUE),]


# =========================================================== Create Function for DGE table =================================================
#this function automates the findmarkers() by treatment and seqmethod for each cluster identity (like the section above) and returns the 
#df in a list. The seurat object is modified in function however does not affect the object outside the environment
DGEtable <- function(seurat_obj, ident.1, ident.2) {
  gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)
  #ensure that the default assay is switched to "RNA" prior to DGE analysis
  if (DefaultAssay(obj_integrated_filtered) == "integrated"){
    print("switching assay to RNA")
    DefaultAssay(seurat_obj) <- "RNA"
  }
  #if there isn't a column in the metadata specifying cluster ident by treatments, make one
  if(!"cell.type.treatment.method" %in% colnames(seurat_obj@meta.data)){
    print("creating new column in metadata")
    seurat_obj$cell.type.treatment.method <- paste(Idents(seurat_obj), seurat_obj$treatment, seurat_obj$seq.method, sep = "_")
  }
  print("switching to identity based on treatment and method")
  Idents(seurat_obj) <- "cell.type.treatment.method"
  df <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2, verbose = TRUE)
  df$Gene.name.uniq <- rownames(df)
  #inspect new col name change
  names(df)
  #merge gene table with all markers generated based on uniq gene symbol
  df <- merge(df, gene_table, by = "Gene.name.uniq") 
  #reorder df in ascending order based on avg_logFC (col#3)
  df <- df[order( df[,3], decreasing = TRUE),]
  return(df)
}

y <-DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = "central-cells_1hr_smartseq2", ident.2 = "central-cells_homeo_smartseq2")

meta <- obj_integrated_filtered@meta.data

#sort alphabetically, extract unique cluster idents, omit NA's, turn to characters, split each idents into elements of list
cluster.ident <- strsplit(sort(unique(as.character(na.omit(meta$cell.type.ident)))), split = "[][']|,\\s*")

cluster.ident.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)){
  print(paste0("working on cluster: ",cluster.ident[[i]]))
  cluster.ident.list[[i]] <- DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2"))
  assign(paste0(cluster.ident[[i]],".smrtseq.1hrvHomeo"), DGEtable(seurat_obj = obj_integrated_filtered, ident.1 = paste0(cluster.ident[[i]],"_1hr_smartseq2"), ident.2 = paste0(cluster.ident[[i]],"_homeo_smartseq2")))
}

names(cluster.ident.list) <- cluster.ident

# =========================================================== Read in Sungmin's Regen App data  =================================================s
setwd("../data/sungminregenappDGE/")

file_list <- list.files(path="../original_sungminregenappDGE/")

sungmin.regen.cluster.ident.list <- vector("list", length = length(file_list))

for (i in 1:length(file_list)){
  sungmin.regen.cluster.ident.list[[i]] <- read.table(file = file_list[[i]], header = TRUE, sep = "\t")
  assign(paste0(gsub(".tsv", "" ,file_list[[i]]),".sungmin.1hrvHomeo"), read.table(file = file_list[[i]], header = TRUE, sep = "\t"))
}

names(sungmin.regen.cluster.ident.list) <- gsub(".tsv", "", file_list)

setwd("/home/ntran2/bgmp/pio-lab/scripts")


# =========================================================== Merge smartseq and sungmin's regen DGE  =================================================s

#empty list to store merged DGE tables
merged.data <- vector("list", length = length(cluster.ident))

conserved.smrtseq <- vector("list", length = length(cluster.ident))

scatterplot.list <- vector("list", length = length(cluster.ident))

for (i in 1:length(cluster.ident)) {
  #merge the smrtseq data with sungmin's data, store within a list, merged.data
  merged.data[[i]] <- merge(x = cluster.ident.list[[i]], y = sungmin.regen.cluster.ident.list[[i]], suffixes = c(".smrtseq2", ".regen.app"),
                            by = c("Gene.name.uniq", "Gene.stable.ID", "Gene.name", "ZFIN.ID", "Gene.description"), sort = FALSE, all = TRUE)
  merged.data[[i]]$diff_avg_logFC <- merged.data[[i]]$avg_logFC.smrtseq2 - merged.data[[i]]$avg_logFC.regen.app
  #generate scatterplots
  png(figurePath(paste0(cluster.ident[[i]], ".png")), width = 11,
      height = 9, units = "in", res = 200)
  scatterplot.list[[i]] <- ggplot(merged.data[[i]], aes(avg_logFC.regen.app, avg_logFC.smrtseq2, label = Gene.name.uniq)) + 
    geom_point() + 
    ggtitle(paste0(cluster.ident[[i]])) + 
    geom_text(aes(label=Gene.name.uniq),hjust=0, vjust=0) + 
    coord_fixed(xlim = c(-1, 3), ylim = c(-1,3)) + 
    scale_x_continuous(name ="Sungmin Regen App - avg LogFC - 1hrVShomeo") +
    scale_y_continuous(name = "smartseq2 - avg LogFC - 1hrVShomeo")
  print(scatterplot.list[[i]])
  dev.off()
  #filter avg_logFC in sungmin's regen dataset for NA's, rearrange to show avg_logFC between smartseq and sungmins
  conserved.smrtseq[[i]] <- filter(merged.data[[i]], is.na(avg_logFC.regen.app) == "TRUE") %>% select(Gene.name.uniq, avg_logFC.smrtseq2, avg_logFC.regen.app)
  
  
}
names(merged.data) <- cluster.ident
names(conserved.smrtseq) <- cluster.ident
names(scatterplot.list) <- cluster.ident


merge <- merge(x = `AP-cells.smrtseq.1hrvHomeo`, y = `AP-cells.sungmin.1hrvHomeo`, suffixes = c(".smrtseq2", ".regen.app"),
               by = c("Gene.name.uniq", "Gene.stable.ID", "Gene.name", "ZFIN.ID", "Gene.description"), sort = FALSE, all = TRUE)

amp <-merged.data$`amp-SCs`


for (i in 1:length(cluster.ident)) {
  conserved.smrtseq[[i]] <- filter(merged.data[[i]], is.na(avg_logFC.regen.app) == "TRUE") %>% select(Gene.name.uniq, avg_logFC.smrtseq2, avg_logFC.regen.app)
  }
names(conserved.smrtseq) <- cluster.ident

# =========================================================== Volcano Plot  =================================================s

for (i in 1:length(cluster.ident)){
  #generate volcano plots for the smartseq2 data
  #label genes that are upregulated
  upreg.genes.to.label <- cluster.ident.list[[i]]$Gene.name.uniq[1:20]
  #label genes that are downregulated
  downreg.genes.to.label <- tail(cluster.ident.list[[i]]$Gene.name.uniq, 20)
  
  png(figurePath(paste0("volcanoplot.",cluster.ident[[i]],".smartseq2.png")), width = 11,
      height = 9, units = "in", res = 200)
  with(cluster.ident.list[[i]], plot(avg_logFC, -log10(p_val), pch=20, main=paste0("Analysis of Cluster: ", cluster.ident[[i]]), col = "grey"))
  grid(NULL,NULL, lty = 6, col = "lightgrey") 
  mtext("Volcano Plot - Smartseq2 Data -1hrVShomeo", line = 0)
  with(subset(cluster.ident.list[[i]] , p_val<.05), points(avg_logFC, -log10(p_val), pch=20,col="black")) #color significant genes
  with(subset(cluster.ident.list[[i]], (avg_logFC)>1), points(avg_logFC, -log10(p_val), pch=20, col="green")) #color upregulated genes
  with(subset(cluster.ident.list[[i]],(avg_logFC)<(1*-1)), points(avg_logFC, -log10(p_val), pch=20, col="red")) #color downregulated genes
  #specify points to label - upreg
  text(cluster.ident.list[[i]]$avg_logFC[1:20], -log10(cluster.ident.list[[i]]$p_val[1:20]),labels=upreg.genes.to.label)
  #specify points to label - downreg
  text(tail(cluster.ident.list[[i]]$avg_logFC, 20), -log10(tail(cluster.ident.list[[i]]$p_val, 20)), labels = downreg.genes.to.label)
  dev.off()
  
  #generate volcano plots for the Sungmin's data
  #label genes that are upregulated
  upreg.genes.to.label <- sungmin.regen.cluster.ident.list[[i]]$Gene.name.uniq[1:20]
  #label genes that are downregulated
  downreg.genes.to.label <- tail(sungmin.regen.cluster.ident.list[[i]]$Gene.name.uniq, 20)
  
  png(figurePath(paste0("volcanoplot.",cluster.ident[[i]], ".sungmin.regen.png")), width = 11,
      height = 9, units = "in", res = 200)
  with(sungmin.regen.cluster.ident.list[[i]], plot(avg_logFC, -log10(p_val), pch=20, main=paste0("Analysis of Cluster: ", cluster.ident[[i]]), col = "grey"))
  grid(NULL,NULL, lty = 6, col = "lightgrey") 
  mtext("Volcano Plot - Sungmin's Regeneration Data - 1hrVShomeo", line = 0)
  with(subset(sungmin.regen.cluster.ident.list[[i]] , p_val<.05), points(avg_logFC, -log10(p_val), pch=20,col="black")) #color significant genes
  with(subset(sungmin.regen.cluster.ident.list[[i]], (avg_logFC)>1), points(avg_logFC, -log10(p_val), pch=20, col="green")) #color upregulated genes
  with(subset(sungmin.regen.cluster.ident.list[[i]],(avg_logFC)<(1*-1)), points(avg_logFC, -log10(p_val), pch=20, col="red")) #color downregulated genes
  #specify points to label - upreg
  text(sungmin.regen.cluster.ident.list[[i]]$avg_logFC[1:20], -log10(sungmin.regen.cluster.ident.list[[i]]$p_val[1:20]),labels=upreg.genes.to.label)
  #specify points to label - downreg
  text(tail(sungmin.regen.cluster.ident.list[[i]]$avg_logFC, 20), -log10(tail(sungmin.regen.cluster.ident.list[[i]]$p_val, 20)), labels = downreg.genes.to.label)
  dev.off()
  
}

#label genes that are upregulated
upreg.genes.to.label <- amp$Gene.name.uniq[1:20]
#label genes that are downregulated
downreg.genes.to.label <- tail(amp$Gene.name.uniq, 20)

with(amp, plot(avg_logFC, -log10(p_val), pch=20, main="Analysis of Cluster: ", col = "grey"))
grid(NULL,NULL, lty = 6, col = "lightgrey") 
mtext("Volcano Plot - Smartseq2 Data: ", line = 0)
with(subset(amp , p_val<.05), points(avg_logFC, -log10(p_val), pch=20,col="black"))
with(subset(amp, (avg_logFC)>1), points(avg_logFC, -log10(p_val), pch=20, col="green"))
with(subset(amp,(avg_logFC)<(1*-1)), points(avg_logFC, -log10(p_val), pch=20, col="red"))
#specify points to label - upreg
text(amp$avg_logFC[1:20], -log10(amp$p_val[1:20]),labels=upreg.genes.to.label)
#specify points to label - downreg
text(tail(amp$avg_logFC, 20), -log10(tail(amp$p_val, 20)), labels = downreg.genes.to.label)
legend("bottomleft", legend = c("upreg", "downreg", "not sig"),
       col = c("green","red", "grey"), pch = 19,
       bty = "n", 
       pt.cex = 2, 
       cex = .5, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))


with(sungmin.regen.cluster.ident.list$`amp-SCs`, plot(avg_logFC, -log10(p_val), pch=20, main=paste0("Analysis of Cluster: amp-SCs"), col = "grey"))
grid(NULL,NULL, lty = 6, col = "lightgrey") 
mtext("Volcano Plot - Sungmin's Regeneration Data", line = 0)
with(subset(sungmin.regen.cluster.ident.list$`amp-SCs` , p_val<.05), points(avg_logFC, -log10(p_val), pch=20,col="black")) #color significant genes
with(subset(sungmin.regen.cluster.ident.list$`amp-SCs`, (avg_logFC)>1), points(avg_logFC, -log10(p_val), pch=20, col="green")) #color upregulated genes
with(subset(sungmin.regen.cluster.ident.list$`amp-SCs`,(avg_logFC)<(1*-1)), points(avg_logFC, -log10(p_val), pch=20, col="red")) #color downregulated genes
#specify points to label - upreg
text(sungmin.regen.cluster.ident.list[[i]]$avg_logFC[1:20], -log10(sungmin.regen.cluster.ident.list[[i]]$p_val[1:20]),labels=upreg.genes.to.label)
#specify points to label - downreg
text(tail(sungmin.regen.cluster.ident.list$`amp-SCs`$avg_logFC, 20), -log10(tail(sungmin.regen.cluster.ident.list$`amp-SCs`$p_val, 20)), labels = downreg.genes.to.label)
legend("bottomleft", legend = c("upreg", "downreg", "not sig"),
       col = c("green","red", "grey"), pch = 19,
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = F , 
       inset = c(0.1, 0.1))

save.image("../data/integ_DGE.RData")
save(merged.data, file = "../data/integ_DGE.RDS")


# =========================================================== export  =================================================s
#creates a single excel file with multiple sheets/tabs pertaining to DGE analysis based on cluster IDs
write_xlsx(list("amp-SCs" = merged.data$`amp-SCs`,
                "AP-cells" = merged.data$`AP-cells`,
                "central-cells" = merged.data$`central-cells`,
                "DV-cells" = merged.data$`DV-cells`,
                "early-HCs" = merged.data$`early-HCs`,
                "HC-prog" = merged.data$`HC-prog`,
                "Inm" = merged.data$Inm,
                "mantle-cells" = merged.data$`mantle-cells`,
                "mature-HCs" = merged.data$`mature-HCs`), path = "../data/integ_DGE/exported_files/IDGE-smrtseqVsungminregen-1hrVhomeo.xlsx")

