library(tidyverse)
library(here)
library(vegan)
library(gridExtra)
library(patchwork)
library(SIAMCAT)
# install.packages("ggpubr")
# install.packages("ggrepel")
# BiocManager::install("SIAMCAT")

TODO: #Load the 16S samples with and without QC and check differences in total reads and relAbs for each QC/noQC - pair. 
  # Think about visualization
  # check whtehter PCoA of sample pairs cluster or whether clustering is happening via method
  # differntial abuncance test of single genera with or without QC SIAMCAT cancer samples 

TODO: # load the WGS counts and check whether the   

# 
# PRJNA518071 <- readRDS("/Users/zellergroup/Google Drive/My Drive/Documents/HCC_gut_microbiome/data/PRJNA517801/res_mOTUs3_mOTU_l75g3.rds")
# View(PRJNA518071)
# 
# PRJNA518071_species <- readRDS("/Users/zellergroup/Google Drive/My Drive/Documents/HCC_gut_microbiome/data/PRJNA517801/res_mOTUs3_species_l75g3.rds")
# 
# View(PRJNA518071_species)
# Qin_2014_species <- readRDS("/Users/zellergroup/Google Drive/My Drive/Documents/HCC_gut_microbiome/data/WGS/Qin_2014/res_mOTUs3_species_l75g3.rds")
# Ida_2021_species <- readRDS("/Users/zellergroup/Google Drive/My Drive/Documents/HCC_gut_microbiome/data/WGS/Ida_2021/res_mOTUs3_species_l75g3.rds")
# Lin2023 <- readRDS("/Users/zellergroup/Google Drive/My Drive/Documents/HCC_gut_microbiome/data/16S/Lin2023/res_mapseq_genus_standardAM.rds")
# colnames(Lin2023) 
# 
# 
# 
# 
# Lin2023
# write.table(Lin2023, file = "Lin2023_counts.txt", sep = "\t")

qc <- readRDS("/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/16S/Lin2023/QC_tax_level/res_mapseq_genus_standardAM.rds")
no_qc <- readRDS("/Users/zellergroup/Library/CloudStorage/GoogleDrive-lucas.reibenspies@embl.de/My Drive/Documents/HCC_gut_microbiome/data/16S/Lin2023/noQC_tax_level/res_mapseq_genus_standardAM.rds")

#1: Compare total Bacterial reads
bac_reads_no_qc <- no_qc["Bacteria",] %>% 
  enframe() %>% 
  rename(sample_id = name) %>% 
  mutate(condition = "no QC") %>% 
  pivot_wider(names_from = sample_id, values_from = value)


bac_reads_wide <- qc["Bacteria",] %>% 
  enframe() %>% 
  rename(sample_id = name) %>% 
  mutate(condition = "QC") %>% 
  pivot_wider(names_from = sample_id, values_from = value) %>% 
  rbind(. , bac_reads_no_qc) 
  
bac_reads_long <- pivot_longer(bac_reads_wide, !condition ,names_to = "sample_id", values_to = "value")

plot_bac_reads <- bac_reads_long %>% 
  ggplot(aes(x = condition, y = value, fill = condition)) +
  geom_boxplot() +
  scale_fill_brewer(type = "seq", palette = "Greens") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  theme_classic() +
  ggtitle("Bacterial reads in QC vs. no QC") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.position = "none") +
  geom_jitter(aes(shape = condition, color = condition)) +
  scale_shape() +
  scale_color_manual(values = c("no QC" = "cyan4", "QC" = "brown")) +
  labs(x = "Method", y = "total bacteria reads per sample", shape = "single values of method") +
  # theme(legend.position = "right") +
    NULL
          
ggsave(plot_bac_reads, filename = here("../boxplot_bac_reads.pdf"), width = 5, height = 5)

##### load metadata ###
Lin_meta <- read_tsv("../data/metadata/Lin2023_meta_tidy.tsv") 
Lin_meta %>% 
  group_by(condition) %>% 
  summarise(count = n())  
  
Lin_meta <- Lin_meta %>% 
  rename(sample_id = Run) %>% 
  rename(disease = condition)
#### Scatterplot ####

scatterplot_bac_reads <- bac_reads_long %>%
  relocate(sample_id, .before = condition) %>% 
  full_join(., Lin_meta, by = "sample_id") %>% 
#  separate(condition, c("QC", "no QC"))
  pivot_wider(names_from = condition, values_from = value) %>% 
  # rename()
  ggplot(aes(y = log10(`no QC` +1), x = log10(QC +1), color = disease, shape = disease, group = sample_id)) +
  geom_point() + 
  labs(x = "log10 of total bacteria reads with QC", y = "log10 of total bacteria reads without QC", title = "Scatterplot of total bacterial reads by method") +
  geom_line() +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  tune::coord_obs_pred()+
  
  NULL


ggsave(scatterplot_bac_reads, filename = here("../Scatter_reads.pdf"), width = 5, height = 5)

#2: Ordination -- PCOA
ord_qc <- qc %>% 
  as_tibble(rownames = "bac") %>% 
  pivot_longer(-bac,names_to = "sample_id",values_to = "counts") %>% 
  filter(bac != 'Bacteria') %>% 
  mutate(sample_id = paste0(sample_id,'_QC')) %>% 
  group_by(sample_id) %>% 
  mutate(rel = counts/sum(counts)) %>% 
  ungroup()
ord_no_qc <- no_qc %>% 
  as_tibble(rownames = "bac") %>% 
  pivot_longer(-bac,names_to = "sample_id",values_to = "counts") %>% 
  filter(bac != 'Bacteria') %>% 
  mutate(sample_id = paste0(sample_id,'_no_QC')) %>% 
  group_by(sample_id) %>% 
  mutate(rel = counts/sum(counts)) %>% 
  ungroup()

ordination <-  bind_rows(ord_qc, ord_no_qc) %>% 
  select(-counts) %>% 
  pivot_wider(names_from = 'bac', values_from = 'rel', values_fill = 0) %>% 
  as.data.frame() %>% 
  column_to_rownames("sample_id") %>% 
  as.matrix()
class(ordination)
  
  # column_to_rownames("sample_id")


# bac_reads_long %>% 
#   mutate(sample_id_new = paste(sample_id,condition,sep = '_')) %>% 
#   select(-sample_id,-condition) %>% 
#   pivot_wider(names_from = sample_id_new,values_from = value)
qc_distances <- vegdist(x = log10(ordination + 1E-6), method = "manhattan") # maybe insert log10 formula 
# qc_distances <- vegdist(x = (ordination), method = "bray") # maybe insert log10 formula 
pcoa_object <- cmdscale(d = qc_distances, k = 2,eig = T) %>% 
  as.data.frame()

eigenvalues <- pcoa_object$eig
var <- eigenvalues/(sum(eigenvalues))*100


pcoa_object1 <-  pcoa_object %>% 
  rownames_to_column("sampleID") %>% 
  as_tibble()%>% 
  separate(sampleID, into = c("sample_id", "suffix"), sep = "_") %>% 
  rename(method = suffix) %>% 
  mutate(method = gsub("no", "no_QC", method))
pcoa_object2 <-  pcoa_object %>% 
   rownames_to_column("sampleID") %>% 
   as_tibble() %>% 
   separate(sampleID, into = c("sample_id", "suffix"), sep = "_") %>% 
   rename(method = suffix) %>% 
   mutate(method = gsub("no", "no_QC", method)) %>% 
   pivot_wider(names_from = method, values_from = c(V1, V2)) %>% 
  distinct(sampleID, .keep_all = TRUE)


plot_object_meta <- inner_join(pcoa_object1, select(Lin_meta, sample_id, disease), by = "sample_id") %>% 
  rename(v1 = V1) %>% 
  rename(v2 = V2)

pcoa_plot <- plot_object_meta %>% 
  ggplot(aes(x = v1, y = v2, col = method)) +
  geom_point(aes(shape = method)) +
  # geom_line(aes(group = sample_id)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on condition and QC status")+ 
  xlab(paste0("PC1 (", round(var[1], 2), "% Variance explained)")) +
  ylab(paste0("PC1 (", round(var[2], 2), "% Variance explained)")) +

    NULL
pcoa_plot

ggsave(pcoa_plot, filename = here("../pcoa_plot_method.pdf"), width = 5, height = 5)



 
#3 Scatter: genus by genus comparison 
thresh_for_prev <- 1e-6
thresh_prev <- 0.05
rrafection <- 1000

scatter_qc <-  qc %>% 
  as_tibble(rownames = 'bac') %>% 
  filter(bac != "Bacteria") %>% 
  pivot_longer(-bac,names_to = 'Sample_id',values_to = "count") %>% 
  mutate(rel = count/sum(count)) %>% 
  rename(count_QC = count) %>% 
  rename(rel_QC = rel) %>% 
  group_by(bac) %>% 
  mutate(N_samples = n(),
         prev_QC = sum(rel_no_QC>thresh_for_prev)/n()) %>% 
  filter(prev_QC > thresh_prev) %>% 
  ungroup()

hist(log10(no_qc),breaks = 100)

scatter_no_qc <- no_qc %>% 
  as_tibble(rownames = 'bac') %>% 
  filter(bac != "Bacteria") %>% 
  pivot_longer(-bac,names_to = 'Sample_id',values_to = "count_no_QC") %>%
  mutate(rel_no_QC = count_no_QC/sum(count_no_QC)) %>% 
  group_by(bac) %>% 
  mutate(N_samples = n(),
         prev_no_QC = sum(rel_no_QC>thresh_for_prev)/n()) %>% 
  filter(prev_no_QC > thresh_prev) %>% 
  ungroup()

 
tmp <- scatter_no_qc %>% select(bac,prev) %>% distinct()
summary(tmp$prev>0.05)


scatter_object <-  full_join(scatter_no_qc, scatter_qc, by = c("Sample_id", "bac"))
  
scatterplot_count  <- scatter_object %>% 
 ggplot(aes(x = log10(count_no_QC + 1), y = log10(count_QC + 1))) +
#  ggplot(aes(x = count_no_QC, y = count_QC)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(title = "total reads of Bacterial Genera per sample", 
       x = "log10(counts without QC)", 
       y = "log10(counts with QC)") +
  tune::coord_obs_pred() +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  ggpubr::stat_cor(method = "spearman") +
  NULL

ggsave(scatterplot_count, filename = here("../Scatte_count.pdf"), width = 5, height = 5)

hist(log10(scatter_object$rel_no_QC +1E-7), breaks = 100)
scatterplot_rel  <- scatter_object %>% 
  ggplot(aes(x = log10(rel_no_QC +1e-6), y = log10(rel_QC +1e-6))) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(title = "Relative Abundance of Bacterial Genera per sample", 
       x = "log10(relative abundance without QC)", 
       y = "log10(relative abundance with QC)") +
  tune::coord_obs_pred() +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  ggpubr::stat_cor(method = "spearman") +
  # ggpubr::stat_cor(method = "pearson") +
  # ggrepel::geom_text_repel()
  NULL
                                                      
scatterplot_rel
scatterplot_count + scatterplot_rel
ggsave(scatterplot_rel, filename = here("../Scatter_rel.pdf"), width = 5, height = 5)

#compare with rrarefied samples: 
rrscatter_no_qc <- ordination %>% 
  rrarefy(1000)
