library(Seurat)
library(future)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggrepel)

if (TRUE) {
  setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")
  
  script_name <- "homeo_samples_integ"
  
  figurePath <- function(filename, format){paste0("/Volumes/easystore/SIMR_2019",
                                                  "/pio-lab/scripts/", script_name, "_figures/", filename)}

}

sample.ls <- c("L29314", "L34727", "L34728")

for (file in 1:length(sample.ls)){
  
  
  
}





