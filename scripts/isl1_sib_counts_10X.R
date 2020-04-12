###############################scRNA-seq Analysis of 10X Homeostatic isl1 Sibling#######################

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
options(future.globals.maxSize = 5000 * 1024^2)

setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")

# Daniel stuff
if (FALSE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")

  script_name <- "isl1_sib_counts_10X"

  figurePath <- function(filename){paste0("/n/projects/ddiaz/Analysis",
  "/Scripts/Nicole/pio-lab/scripts/", script_name, "_figures/", filename)}

  #devtools::load_all("/n/projects/ddiaz/Analysis/Scripts/SeuratExtensions")
}

####Read in Data####

#load in isl1_sib_counts_10X dir containing count matrix
#files in this folder should be in .gz format
isl1_sib_10X.data <- Read10X(data.dir = "../data/isl1_sib_counts_10X/")

homeo.isl1_sib_10X <- CreateSeuratObject(counts = isl1_sib_10X.data, project = "homeo.isl1.sib.10X", min.cells = 1, min.features = 1)

####Standard pre-processing workflow####
#calculate mitochondrial contamination
#do not add suffix to identify sample in "percent.mt" for metadata, will cause convergance issues in downstream integration with multiple datasets
homeo.isl1_sib_10X[["percent.mt"]] <- PercentageFeatureSet(homeo.isl1_sib_10X, pattern = "^mt-")

####QC metrics-VlnPlot####
#vlnplots of number of genes, molecule count, and mitochondrial % in cell
#summary tables used to determine minimal lower threshold nCounts_RNA for subsetting

#nFeature_RNA
#will not be delimitting threshold for nFeature_RNA, we want all possible genes
nFeature.vln.homeo.isl1.sib.10X <-VlnPlot(object = homeo.isl1_sib_10X, features = "nFeature_RNA") + geom_hline(yintercept = 3000, color = 'red') + geom_hline(yintercept = 200, color = 'red') +  scale_y_continuous(breaks=c(200, 3000))
summary(homeo.isl1_sib_10X$nFeature_RNA)

nFeature.vln.homeo.isl1.sib.10X

#nCount_RNA
nCount.vln.homeo.isl1.sib.10X <- VlnPlot(object = homeo.isl1_sib_10X, features = "nCount_RNA") + ylim(0, 15000) + geom_hline(yintercept = 10000, color = 'red') 
nCount.vln.homeo.isl1.sib.10X
summary(homeo.isl1_sib_10X$nCount_RNA)

#percent.mito %
pct.mito.vln.homeo.isl1.sib.10X <- VlnPlot(object = homeo.isl1_sib_10X, features = "percent.mt") + geom_hline(yintercept = 10, color = 'red') + scale_y_continuous(breaks=10)
summary(homeo.isl1_sib_10X$percent.mt)

#CombinePlots(plots = list(nFeature.vln.homeo.isl1.sib.10X + nCount.vln.homeo.isl1.sib.10X + pct.mito.vln.homeo.isl1.sib.10X), ncol = 3)

#boxplot to visulize lower bound outliers in nFeature, set threshold to ~200
boxplot(homeo.isl1_sib_10X$nFeature_RNA, lwd= 2, ylim = c(0,500))
stripchart(homeo.isl1_sib_10X$nFeature_RNA, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue', bg= 'black')

#head(homeo.isl1_sib_10X$nCount_RNA)

####QC metrics-FeatureScatter####
featurescatter.nCountXpercent.mt.homeo.isl1.sib.10X <- FeatureScatter(homeo.isl1_sib_10X, feature1 = "nCount_RNA", feature2 = "percent.mt") + xlim(0, 15000)
featurescatter.nCountXnFeature.homeo.isl1.sib.10X <- FeatureScatter(homeo.isl1_sib_10X, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + xlim(0, 15000)
featurescatter.nCountXpercent.mt.homeo.isl1.sib.10X + featurescatter.nCountXnFeature.homeo.isl1.sib.10X

####Subsetting####
#subsetting parameters: omit cells with genes less than 3000 and greater than 200. omit mitochondrial contamination greater than 10%. omit cells with molecules greater than 10,000
homeo.isl1_sib_10X <- subset(homeo.isl1_sib_10X, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10 & nCount_RNA <10000) 

homeo.isl1_sib_10X
#16607 cells preserved from 17713 in original data

####Normalize####
#call the plan() function to initiate the multiprocess/parallelization for NormalizeData()
plan("multiprocess")
homeo.isl1_sib_10X <- NormalizeData(homeo.isl1_sib_10X, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)
#mutliprocess for NormalizationData using 1G is executable, not same case for ScaleData()

####FindVariableFeatures####
#Identification of highly variable genes
homeo.isl1_sib_10X <- FindVariableFeatures(homeo.isl1_sib_10X, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(homeo.isl1_sib_10X), 10)

#plot 2000 variable features with labels on top 10 most variable genes
VariableFeatureplot1.homeo.isl1_sib_10X <- VariableFeaturePlot(homeo.isl1_sib_10X) 
VariableFeature.labedlplot2.homeo.isl1_sib_10X <- LabelPoints(plot = VariableFeatureplot1.homeo.isl1_sib_10X, points = top10, repel = TRUE) + ggtitle("Top 10 Variable Genes Using VariableFeaturePlot")

VariableFeature.labedlplot2.homeo.isl1_sib_10X

####Scale####
all.genes <- rownames(homeo.isl1_sib_10X)
#25622 genes total
#another method of accessing unique, nonduplciated gene symbols in counts data.
#returns 0 if there are no duplicates
test<- all.genes %in% all.genes[duplicated(all.genes)] 
sum(grepl(pattern = "TRUE", x = test))

plan("sequential")
#options(future.globals.maxSize = 3500 * 1024^2)
#plan("multiprocess")
homeo.isl1_sib_10X <- ScaleData(homeo.isl1_sib_10X, features = NULL, vars.to.regress = NULL, verbose = TRUE)

####Save Environment####
#default save RData env
save.image("../data/workspace_homeo_isl1_sib_10X.RData")
#3.1G

####Perform PCA####
#default pc is set to 100
homeo.isl1_sib_10X <- RunPCA(homeo.isl1_sib_10X, npcs = 100, verbose = TRUE, features = NULL)

VizDimLoadings(homeo.isl1_sib_10X, dims = 1:2, reduction = "pca")

DimHeatmap(homeo.isl1_sib_10X, dims = 1:15, cells = 500, balanced = TRUE)

DimPlot(homeo.isl1_sib_10X, reduction = "pca")

####Choosing PC####

####Elbow Plot####
ElbowPlot(homeo.isl1_sib_10X, ndims = 50)

ElbowPlot(homeo.isl1_sib_10X, ndims = 35)

ElbowPlot(homeo.isl1_sib_10X, ndims = 20)

####JackStraw####
#dims is defaulted to 20
homeo.isl1_sib_10X <- JackStraw(homeo.isl1_sib_10X, num.replicate = 100, dims = 35)
homeo.isl1_sib_10X <- ScoreJackStraw(homeo.isl1_sib_10X, dims = 1:35)
JackStrawPlot.PC20.homeo.isl1_sib_10X<- JackStrawPlot(homeo.isl1_sib_10X, dims = 1:35)
JackStrawPlot.PC20.homeo.isl1_sib_10X

####Clustering/UMAP####
homeo.isl1_sib_10X <- FindNeighbors(homeo.isl1_sib_10X, dims = 1:15, features = NULL)

homeo.isl1_sib_10X <- FindClusters(homeo.isl1_sib_10X, resolution = 1.2)

homeo.isl1_sib_10X <- RunUMAP(homeo.isl1_sib_10X, dims = 1:15, reduction = "pca")


umap.homeo.isl1.sib.10X <- DimPlot(homeo.isl1_sib_10X, reduction = "umap",
  label = TRUE, pt.size= 0.4)

umap.homeo.isl1.sib.10X
png(figurePath("umap_clusters.png"),
  width = 11, height = 9, units = "in", res = 300)
  print(umap.homeo.isl1.sib.10X)
dev.off()

####Annotate Clusters####

#import gene info list pulled from ensembl
gene_table <- read.table("../data/Danio_Features_unique_Ens98_v1.tsv", sep = "\t", header = TRUE)

####Findallmarkers####

#mirrored as closely in parameters to Lush, 2019. as possible (see "Quality control, dimensional reduction, and cell classification" under Methods)
#fold change greater than 0.10, or less than 􏰀0.10
#only return p-value <0.001 (see return.thres())
all.markers.homeo.isl1_sib_10X <- FindAllMarkers(homeo.isl1_sib_10X, only.pos = FALSE, min.pct = 0.01, logfc.threshold = 0.01, return.thresh = 0.001, verbose = TRUE)

#change the column name formerally called "gene" generated from FindAllMarkers() to  "Gene.name.uniq" so we can merge with gene_table
#must have same column name in order to match genes
colnames(all.markers.homeo.isl1_sib_10X)[7] <- "Gene.name.uniq"

#inspect new col name change
names(all.markers.homeo.isl1_sib_10X)

#merge gene table with all markers generated based on uniq gene symbol
all.markers.homeo.isl1_sib_10X <- merge(all.markers.homeo.isl1_sib_10X, gene_table, by = "Gene.name.uniq") 

#reorder df in ascending order based on cluster (col #7) and then avg_logFC (col#3)
all.markers.homeo.isl1_sib_10X <- all.markers.homeo.isl1_sib_10X[order( all.markers.homeo.isl1_sib_10X[,7], all.markers.homeo.isl1_sib_10X[,3] ),]

#top 50
top_50.all.marker.homeo.isl1 <- all.markers.homeo.isl1_sib_10X %>% group_by(cluster) %>% top_n(n = 50, wt = (avg_logFC))

#top 50 Heat Map
HeatMap.Top50Markers.Homeo.isl1<- DoHeatmap(homeo.isl1_sib_10X, features = top_50.all.marker.homeo.isl1$Gene.name.uniq)

#export in data/ dir
#dir.create('../data/results_homeo.isl1_sib_10X/')
write.table(all.markers.homeo.isl1_sib_10X, file = '../data/results_homeo.isl1_sib_10X/findallmarkers.tsv', sep = '\t', row.names = F)

####FindMarkers####
for (x in levels(Idents(homeo.isl1_sib_10X))) {
  cluster_name <- paste('cluster', x, sep="_")
  print(paste("finding conserved markers and exporting genes from cluster", x, sep=" "))
  temp_cluster <- assign(paste0("cluster", x, sep = "_"),FindMarkers(homeo.isl1_sib_10X, ident.1 = x, only.pos = FALSE, min.pct = 0.10, logfc.threshold = 0.10, verbose = TRUE))
  #cluster_name$Gene.name.uniq <- rownames(cluster_name)
  #cluster_name <- merge(cluster_name, gene_table, by = "Gene.name.uniq")
  #write.table(temp_cluster, paste(cluster_name, ".tsv", sep=""), sep="\t", col.names=TRUE)
}

####Identify Cluster Identity####
#canonical marker lists
polar.cells.markers <- c("sost", "adcyap1b", "six2b", "fsta", "wnt2", "srrt")
hair.cells.markers <- c("atoh1a", "dld", "her4.1", "tekt3")
ap.cells.markers <- c("cx44.2", "fap", "fgf10a", "hmx2", "hmx3a", "mcf2lb", "nfasca", "pdgfaa", "prom1b", "ompb", "pip5k1bb", "si:ch73-261i21.5", "tnfsf10l3", "wdpcp", "wnt11r", "znf185")
mantle.cells.markers <- c("cldne", "crb3b", "crip1", "cts12", "fat1b", "fgfr2", "map7d3", "mcama", "ovgp1", "pkhd1l1", "pknox2", "ponzr6", "rgs4", "sfrp1a", "stat3",
                          "sulf1",
                          "tmem119b",
                          "tnfb",
                          "tnfsf10",
                          "tspan1")
ring.cell.markers <- c("alpl", "ckma",  "cmah" , "fndc7rs4","hopx", "nr4a3", "tcf7l2", "tfap2a")
central.cell.markers <- c("ccl20a.3","crabp2a","cxcl20", "dhrs3a", "dkk2", "fabp11a", "fabp7a",
"fgf3", "fhdc2", "gata2a", "hbgefb", "hey2", "igfbp5a"  ,"isl1", "itgb4", "kitlga", "krt92", "lima1a", "mb", "mhc1lda", "ptprz1a", "rdh10a", "si:ch73-380n15.2", "sox21b","sult3st2",  "sult3st3", "abcg4b", "AL954146.1",
"abcg4b","bmp4", "CABZ01067232.1l", "ebf3a", "si:dkey-4p15.5 fgf22", "gata2b", "glula", "her15.1", "her6",
"hey2", "HSPA8", "krt18", "lfng", "mb", "mfng", "palm1a", "proca1","prox1a", "slc1a3a/glasta", "soul5", "tppp3", "vegfab", "zgc:165423")


#si:dkey-4p15.5 = ENSDARG00000086272

#find top 5 avg_logFC for ap cell markers, top avg_logFC arranged as 1st row
filter(all.markers.homeo.isl1_sib_10X, Gene.name.uniq == ap.cells.markers) %>% top_n(n = 5, wt = (avg_logFC)) %>% arrange(desc(avg_logFC))


####Polar Cells####
polar.cell.featureplot <-FeaturePlot(homeo.isl1_sib_10X, features = polar.cells.markers, label = TRUE) + ggtitle("Expression Levels of Polar Cell Markers")
polar.cell.vlnplot <- VlnPlot(homeo.isl1_sib_10X, features = polar.cells.markers, pt.size  = 0)+ ggtitle("Expression Levels of Polar Cell Markers")
#ars2 does not show up in all.marker.homeo, perhaps signal is too low?
#0,8,11,6?

####Hair Cell Lineage####
hair.cell.featureplot <- FeaturePlot(homeo.isl1_sib_10X, features = hair.cells.markers, label = TRUE) + ggtitle("Expression Levels of Hair Cell Markers")
hair.cell.vlnplt <- VlnPlot(homeo.isl1_sib_10X, features = hair.cells.markers, pt.size  = 0)+ ggtitle("Expression Levels of Hair Cell Markers")

####A-P Cells####
ap.cells.featureplt <- FeaturePlot(homeo.isl1_sib_10X, features = ap.cells.markers, label = TRUE) + ggtitle("Expression Levels of A-P Cell Markers")
ap.cells.vlnplt <- VlnPlot(homeo.isl1_sib_10X, features = ap.cells.markers, pt.size = 0) + ggtitle("Expression Levels of A-P Cell Markers")

####Mantle Cells####
mantle.cells.featureplt <- FeaturePlot(homeo.isl1_sib_10X, features = mantle.cells.markers, label = TRUE)+ ggtitle("Expression Levels of Mantle Cell Markers")
mantle.cells.vlnplt <-  VlnPlot(homeo.isl1_sib_10X, features = mantle.cells.markers, pt.size = 0) + ggtitle("Expression Levels of Mantle Cell Markers")
#ERROR: The following requested variables were not found: lye, si:dkey−103g5.3


####Ring#### 
ring.cell.featureplt <- FeaturePlot(homeo.isl1_sib_10X, features = ring.cell.markers, label = TRUE)+ ggtitle("Expression Levels of Ring Cell Markers")
ring.cell.vlnplt <- VlnPlot(homeo.isl1_sib_10X, features = ring.cell.markers, pt.size = 0)

####Central Cells####
#Part 1
central.cells.featureplt <- FeaturePlot(homeo.isl1_sib_10X, features = central.cell.markers, label = TRUE)+ ggtitle("Expression Levels of Central Cell Markers")
central.cells.vlnplt <- VlnPlot(homeo.isl1_sib_10X, features = central.cell.markers, pt.size = 0)+ ggtitle("Expression Levels of Central Cell Markers")
#Error:The following requested variables were not found: hbgefb


#####Top10GeneCellIdentityFunction#####
Top10GeneCellIdentity <- function(marker.list, df){
  #This function will return a adjusted dataframe for the genes with the 10 (or less) highest avg_logFC
  #and its corresponding cluster
  #it will also generate a FeaturePlot and Violin Plot corresponding
  #this function may be helpful as most canonical markers have low signal, this will discern what signal is highest to help identify clusters
  top10.df <-filter(df, Gene.name.uniq == marker.list) %>% top_n(n = 10, wt = (avg_logFC)) %>% arrange(desc(avg_logFC))
  gene.list <- unique(top10.df$Gene.name.uniq)
  feature.plot <- FeaturePlot(homeo.isl1_sib_10X, features = gene.list, label = TRUE)
  vln.plt <- VlnPlot(homeo.isl1_sib_10X, features = gene.list, pt.size = 0)
  return(list(top10.df, feature.plot, vln.plt))
  }
Top10GeneCellIdentity(marker.list = ap.cells.markers, df = all.markers.homeo.isl1_sib_10X)

FeaturePlotToPdf <- function(marker.list) {
  for (x in marker.list){
    to.pdf <- FeaturePlot(homeo.isl1_sib_10X, features = x, label = TRUE)
    #png(figurePath("umap_clusters.png"),
    #width = 11, height = 9, units = "in", res = 300)
    #print(t)
    #dev.off()
    print(to.pdf)
  }
}
FeaturePlotToPdf(ap.cells.markers)


####Save Figures to PDF####
pdf("./isl1_sib_counts_10X_figures/pct.mito.vln.homeo.isl1.sib.10X.pdf")
pct.mito.vln.homeo.isl1.sib.10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/nFeature.vln.homeo.isl1.sib.10X.pdf")
nFeature.vln.homeo.isl1.sib.10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/nCount.vln.homeo.isl1.sib.10X.pdf")
nCount.vln.homeo.isl1.sib.10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/featurescatterQC.pdf", width=10, height=5)
featurescatter.nCountXpercent.mt.homeo.isl1.sib.10X + featurescatter.nCountXnFeature.homeo.isl1.sib.10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/VariableFeature.labedlplot2.homeo.isl1_sib_10X.pdf")
VariableFeature.labedlplot2.homeo.isl1_sib_10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/HeatMapPC15.pdf")
DimHeatmap(homeo.isl1_sib_10X, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf("./isl1_sib_counts_10X_figures/ElblowPlot.pdf")
ElbowPlot(homeo.isl1_sib_10X, ndims = 20)
dev.off()

pdf("./isl1_sib_counts_10X_figures/JackStrawPlot.PC20.homeo.isl1_sib_10X.pdf")
JackStrawPlot.PC20.homeo.isl1_sib_10X
dev.off()

pdf("./isl1_sib_counts_10X_figures/UMAP.homeo.isl1_sib_10X.pdf", width=10, height=5)
DimPlot(homeo.isl1_sib_10X, reduction = "umap", label = TRUE) 
dev.off()

pdf("./isl1_sib_counts_10X_figures/HeatMap.Top50Markers.Homeo.isl1.pdf")
HeatMap.Top50Markers.Homeo.isl1
dev.off()
