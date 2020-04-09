###############################scRNA-seq Analysis of 10X Homeostatic isl1 Sibling#######################

####Set up environment####
#install.packages("Seurat")
library(Seurat)
library(plotly)
#install.packages("future")
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggrepel)
options(future.globals.maxSize = 5000 * 1024^2)

setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")

####Read in Data####

#load in isl1_sib_counts_10X dir containing count matrix
#files in this folder should be in .gz format
isl1_sib_10X.data <- Read10X(data.dir = "../data/isl1_sib_counts_10X/")

homeo.isl1_sib_10X <- CreateSeuratObject(counts = isl1_sib_10X.data, project = "homeo.isl1.sib.10X", min.cells = 1, min.features = 1)

homeo.isl1_sib_10X