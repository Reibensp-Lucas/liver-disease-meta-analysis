####  taking all RDS files of 16S count matrices on genus level and convert 
# them in the right format while moving to the right folder  

library(tidyverse)
library(here)

setwd("/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/16S")

Cortez2021_counts <- readRDS(here('Cortez_2021/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame() 
Cortez2021_counts <- Cortez2021_counts[-1, ]
write.table(Cortez2021_counts, 
              "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Cortez2021_counts.txt",
              sep = "\t", row.names = TRUE, col.names = TRUE)
  
Duan2019_counts <- readRDS(here('Duan_2019/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
    as.data.frame() 
Duan2019_counts <- Duan2019_counts[-1, ]
write.table(Duan2019_counts, 
              "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Duan2019_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
  

Furukawa2020_counts <- readRDS(here('Furukawa_2020/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame() 
Furukawa2020_counts <-  Furukawa2020_counts[-1, ]
write.table(Furukawa2020_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Furukawa2020_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

Iwasawa2017_counts <- readRDS(here('Iwasawa_2017/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame() 
Iwasawa2017_counts <-  Iwasawa2017_counts[-1, ]
write.table(Iwasawa2017_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Iwasawa2017_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)


Liu2019_counts <- readRDS(here('Liu_2019/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame() 
Liu2019_counts <- Liu2019_counts[-1, ]
write.table(Liu2019_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Liu2019_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

Wang2017_counts <- readRDS(here('Wang_2017/tax_level_profiles/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame()
Wang2017_counts <- Wang2017_counts[-1, ]  
write.table(Wang2017_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Wang2017_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)

Lin2023_counts <- readRDS(here('Lin_2023/QC_tax_level/res_mapseq_genus_standardAM.rds')) %>% 
  as.data.frame()
Lin2023_counts <-  Lin2023_counts[-1, ]
write.table(Lin2023_counts, 
            "/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/raw_profiles/Lin2023_counts.txt",
            sep = "\t", row.names = TRUE, col.names = TRUE)
