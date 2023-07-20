library(tidyverse)
library(here)



Ida2021_counts <- readRDS(here('Ida2021/res_mOTUs3_genus_l75g3.rds')) %>% 
  as.data.frame() 
Ida2021_counts <- Ida2021_counts[-1, ]
write.table(Ida2021_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Ida2021_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

Qin2014_counts <- readRDS(here('Qin2014//res_mOTUs3_genus_l75g3.rds')) %>% 
  as.data.frame() 
Qin2014_counts <- Qin2014_counts[-1, ]
write.table(Qin2014_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Qin2014_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
