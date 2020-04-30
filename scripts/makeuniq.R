# =========================================================== R script that will make the gene names uniq on the cell ranger output
setwd("/Volumes/easystore/SIMR_2019/pio-lab/scripts")

sample.ls <- c("L29314", "L34727", "L34728")

#loop thru each
for (file in 1:length(sample.ls)){
  setwd(paste0("/Volumes/easystore/SIMR_2019/pio-lab/data/homeo_sample_", sample.ls[file], sep = ""))
  setwd("./filtered_feature_bc_matrix/")
  print(getwd())
  # Start R in directory with features.tsv.gz
  if(file.exists("./features.tsv.gz")) {
    system("mkdir backup zipped unzipped ; cp -n *.gz backup/ ; gunzip *.gz")
  }
  features <- read.delim("./features.tsv", header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)
  # Check for repeated genes
  grep(" \\(1 of many\\)", features[,2], value = TRUE)
  grep(" ", features[,2], value = TRUE)
  # Remove "(1 of many)" from common gene names
  features[,2] <- gsub(" \\(1 of many\\)", "", features[,2])
  # find number of duplicates
  table(duplicated(features[,2]))
  # make repeats unique
  features[,2] <- make.unique(features[,2], sep = "*")
  if("Gene Expression" %in% features[,3]) {
    features <- features[,-3]
  }
  write.table(features, "./features.tsv",
              quote = FALSE, sep = "\t", row.names = FALSE)
  system("cp -n * unzipped/ ; gzip * ; mv *gz zipped")

}

#output will save the original filtered mtx in "backups"
#will create the features.tsv file with uniq gene names in "unzipped" and "zipped" files



