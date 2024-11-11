
#=============== CRC Metagenomics Project ===============#

# Aim
  # To merge per-sample relative abundance data into a single abundance matrix


options(
  
  java.parameters = "-Xmx64g",
  stringsAsFactors = F
  
) 

setwd(dir = "D:/2-Research/2-CRC metagenomics/") 


# install.packages("tidyverse")
# install.packages("ggpubr")
# install.packages("magrittr")
# install.packages("broom")
# install.packages("remotes")
# remotes::install_github("jbisanz/qiime2R")

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")


library(phyloseq)
library(tidyverse)
library(ggpubr)
library(magrittr)


# 1. Species

file_paths <- list.files(path = "species/", pattern = "\\.tsv$", 
                         full.names = TRUE)

data_list <- lapply(file_paths, read_tsv)

merged_data <- Reduce(function(x, y) merge(x, y,
                                           by = "clade_name", 
                                           all = TRUE), 
                      data_list); merged_data[is.na(merged_data)] <- 0

write.csv(merged_data, file = "species.csv", row.names = FALSE)
  

# 2. Strain
  
file_paths <- list.files(path = "strain/", pattern = "\\.tsv$", 
                         full.names = TRUE)

data_list <- lapply(file_paths, read_tsv)

merged_data <- Reduce(function(x, y) merge(x, y,
                                           by = "clade_name", 
                                           all = TRUE), 
                      data_list); merged_data[is.na(merged_data)] <- 0

write.csv(merged_data, file = "strain.csv", row.names = FALSE)



# 3. Genus

file_paths <- list.files(path = "genus/", pattern = "\\.tsv$", 
                         full.names = TRUE)

data_list <- lapply(file_paths, read_tsv)

merged_data <- Reduce(function(x, y) merge(x, y,
                                           by = "clade_name", 
                                           all = TRUE), 
                      data_list); merged_data[is.na(merged_data)] <- 0

write.csv(merged_data, file = "genus.csv", row.names = FALSE)