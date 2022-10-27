#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

# reading in new lineage reports
lineage_reports <- read.csv(list.files(path = ".", pattern = "all_lineage_reports*")[1])
lineage_reports$RBD_level <- NA

# assemble a list of variant tables
variant_tables <- list.files(path = ".", pattern = "*vcf")

# looping through each variant table and counting RBD mutations
for (i in 1:length(variant_tables)){
  
  strain_name <- basename(variant_tables[i]) %>%
    str_remove(".vcf") %>%
    str_replace_all("_", "/")
  
  row <- as.numeric(rownames(lineage_reports[lineage_reports$taxon==strain_name,]))
  
  mutations <- read.delim(variant_tables[i], skip = 55)
  # mutations <- mutations[mutations$REF!="-",]
  mutations <- mutations[mutations$ALT!="N",]
  mutations <- mutations[!grepl("-", mutations$ALT),]
  
  rbd <- mutations[mutations$POS > (21554+957) &
                     mutations$POS < (21554+1623),]
  
  if (nrow(rbd)>0){
    mutations <- nrow(rbd)
  } else {
    mutations <- 0
  }
  
  lineage_reports$RBD_level[row] <- mutations
  
}

# constructing new filename
new_filename <- str_replace(list.files(path = ".", pattern = "all_lineage_reports*")[1], 
                            "all_",
                            "rbd_classified_")

# Exporting the final lineage report
write.csv(lineage_reports, new_filename, quote = F, row.names = F)
