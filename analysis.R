#################################
### 1. Load packages and define params
#################################

# install.packages("pheatmap")
# install.packages("patchwork")
# install.packages("gridExtra")
library(RColorBrewer)
library(tidyverse)
library(here)
library(vegan)
library(gridExtra)
library(patchwork)
library(SIAMCAT)
library(pheatmap)
# library(wesanderson)

# notice that combine() is masked from 'package:dplyr' now. 

depth1 <- 2000
depth2 <- 10000
depth3 <- 50000


#############################################
### 2. Load and double-check profiles + meta
#############################################

# define basepaths
pathToProfiles <- here('..', 'data', 'raw_profiles')
pathToMetadata <- here('..', 'data', 'metadata', 'tidy_data')
pathToData <- here('..', 'data')

####### Instructions ##########
# Step 1: Get list/vector of files in basePath. Hint: use list.files() function
# Step 2: Read in each profile as a dataframe into a list. Give respective entry of list the name of the profile.
### list.files(basePath)
## You can use read_tsv() function to read in a data frame
## Figure out how to read in data frame that has an additional column
# Step 3: Implement a few sanity checks:
## print number of samples (expected to be in col)
## print a few examples of samples 
## print number of taxa (expected to be in rows)
## print a few examples of taxa
## print largest and smallest entry in table

##### read count profiles as tibble into a list #########

#TODO: see that new 16S profiles(genus_level) are converted to the rright format and include them, too

readProfiles <- function(pathToProfiles) {
  profiles_file_names <- list.files(pathToProfiles, full.names = TRUE) # Step 1
  profiles_list <- list()
  for (i in seq_along(profiles_file_names)) {  # Step 2
    file_name <- basename(profiles_file_names[i])
    df_name <- gsub("_counts.txt", "", file_name)
    df <- read.table(profiles_file_names[i], sep = "\t", header = TRUE)
    df$Taxa <- rownames(df)
    rownames(df) <- NULL 
    df <- as_tibble(df)
    assign(df_name, df)
    profiles_list[[df_name]] <- df
  }
  count_profiles <- profiles_list
  if (length(count_profiles) == length(profiles_file_names)) {
    cat("All files from", pathToProfiles, "were added to the list 'count_profiles'", "\n")
    cat("\n", "Here is a summary of each data.frame:", "\n")
    for (i in seq_along(count_profiles)) {
      df <- count_profiles[[i]]
      cat("\n ########### Summary for", names(count_profiles)[i], "###########", "\n") #title of summary
      cat("Number of samples:", ncol(df), "\n") #number of samples (columns)
      cat("Sample examples:", paste(colnames(df)[1:3], collapse=", "), "\n") #examples for samples
      cat("Number of taxa:", nrow(df), "\n") #number of taxa (rows)
      cat("Taxa examples:", paste(df$Taxa[1:3], collapse=", "), "\n") # examples for taxw
      cat("Highest abundancy:", df %>% 
            select_if(is.numeric) %>% 
            max(na.rm = TRUE), paste("(",df$Taxa[which(df == max(df %>% select_if(is.numeric), na.rm = TRUE), arr.ind = T)[1]], ")", "\n"))
      cat("Lowest Abundancy:",  df %>% 
            select_if(is.numeric) %>% 
            min(na.rm = TRUE), paste("(",nrow(which(df == min(df %>% select_if(is.numeric), na.rm = TRUE), arr.ind = T)), "Taxa)"), "\n")
    }
  } else {
    stop("Not all files in directory were added to 'count_profiles'")
  }
  
  return(count_profiles)
}

count_profiles <- readProfiles(pathToProfiles)

## Inspect each dataframe #####

summariseProfiles <- function(profiles_list) {
  for (i in seq_along(profiles_list)) {
    df <- profiles_list[[i]]
    cat("\n ########### Summary for", names(profiles_list)[i], "###########", "\n") #title of summary
    cat("Number of samples:", ncol(df), "\n") #number of samples (columns)
    cat("Sample examples:", paste(colnames(df)[1:3], collapse=", "), "\n") #examples for samples
    cat("Number of taxa:", nrow(df), "\n") #number of taxa (rows)
    cat("Taxa examples:", paste(df$Taxa[1:3], collapse=", "), "\n") # examples for taxa
    cat("Highest abundancy:", df %>% 
          select_if(is.numeric) %>% 
          max(na.rm = TRUE), paste("(",df$Taxa[which(df == max(df %>% select_if(is.numeric), na.rm = TRUE), arr.ind = T)[1]], ")", "\n"))
    cat("Lowest Abundancy:",  df %>% 
          select_if(is.numeric) %>% 
          min(na.rm = TRUE), paste("(",nrow(which(df == min(df %>% select_if(is.numeric), na.rm = TRUE), arr.ind = T)), "Taxa)"), "\n")
  }
}
#)
# delete the ".singles" where it occurs

remove.singles <- function(df_list) {
  for (i in seq_along(count_profiles)) {
    df <- count_profiles[[i]]
    colnames(df) <- gsub(".singles", "", colnames(df))
    df_list[[i]] <- df
  }
  cat("Run names were modified. New Names do not contain '.singles' anymore.", "\n")
  return(df_list)
}

count_profiles <- remove.singles(count_profiles)

####### Read metadata #############

readMetadata <- function(pathToMetadata) {
    meta_file_names <- list.files(pathToMetadata, full.names = TRUE) # Step 1
   meta_list <- list()
    for (i in seq_along(meta_file_names)) {  # Step 2
      file_name <- basename(meta_file_names[i])
      #df_name <- gsub(".tsv", "", file_name)
      df_name <- gsub("_meta_tidy.tsv", "", file_name)
      df <- read_tsv(meta_file_names[i])
      assign(df_name, df)
      meta_list[[df_name]] <- df
    }
    meta_data <- meta_list 
    if (length(meta_list) == length(meta_file_names)) {
      cat("\n","All files from", pathToMetadata, "were added to the list 'meta_data'", "\n")
      cat("\n", "Here is a summary of each file:", "\n")
      for (i in seq_along(meta_data)) {
        df <- meta_data[[i]]
        cat("\n ########### Summary for", names(meta_data)[i], "###########", "\n") #title of summary
        cat("Number of Runs:", nrow(df), "\n") #number of samples (rows)
        cat("Run examples:", paste(df[1:3, "Run"], collapse=", "), "\n") #examples for samples
        cat("AssayType:", paste(df$AssayType)[2], "\n") #AssayType
        cat("Average spot length:", mean(df$AvgSpotLen), "\n")
        cat("Library Layout:", paste(df$LibraryLayout[2]), "\n")
        cat("Metadata parameters:", paste(colnames(df), collapse=", "), "\n") # all parameters
        cat("Conditions: ")
        cat(paste(df %>% 
          group_by(condition) %>% 
          summarise(count = n())), "\n") # examples for Condition  
        ##### Sometimes there is a second sample names column
        if ("SampleName_2" %in% colnames(df)) {
          df$SampleName_2 <- gsub("[0-9]+$", "", df$SampleName_2)
          cat("Conditions_2: ")
          cat(paste(df %>% 
                      group_by(SampleName_2) %>% 
                      summarise(count = n())), "\n")# examples for Condition 
        } 
      }
    } else {
      stop("Not all files in directory were added to 'meta_data'")
    }
    
    return(meta_data)
}

metadata <- readMetadata(pathToMetadata)
#### check if sample names are consistent along raw counts and metadata #####

## Manually throwing out data or a  profile.
UC_samples <- tibble(metadata = metadata) %>% 
  unnest(cols = c(metadata)) %>% 
  select(Run, condition) %>% 
  filter(str_detect(condition, "UC")) %>% 
  pull(Run)
match_UC <- sapply(metadata, function(df) any(df$Run %in% UC_samples))


young_samples <- tibble(metadata = metadata) %>% 
  unnest(cols = c(metadata)) %>% 
  select(Run, Age) %>% 
  filter(Age < 6) %>% 
  pull(Run)
match_young <- sapply(metadata, function(df) any(df$Run %in% young_samples))
metadata <- lapply(metadata, function(df) df[!df$Run %in% c(young_samples,UC_samples), ])

all_data <- tibble(counts = count_profiles,
                   dataset_name = names(count_profiles)) %>%
  inner_join(tibble(metadata = metadata,
                    dataset_name = names(metadata)), by = 'dataset_name') %>%
  relocate(counts, metadata)



# Loop through the dataframes in both lists and check if their names match
align_samples <- function(meta_list, count_list) { for (i in seq_along(meta_list)) {
  # Get the dataframe name from meta_list
  name_meta <- names(meta_list)[i]
    # Get the corresponding dataframe name from count_list
  name_count <- gsub("_meta", "", name_meta)
    # Check if the Run column values from meta_list match the colnames from count_list
  run_col <- meta_list[[i]][["Run"]]
  print(str_c("Number of runs in metadata: ", length(run_col)))
  print(str_c("Number of samples in profile: ", dim(count_list[[i]])[[2]] - 1))
  # TODO
  ## Replace print statements either with warnings (or break the code completely).
  if (all(run_col %in% colnames(count_list[[i]]))) {
    print(paste("Sample names match for", name_meta, "and", name_count))
  } else {
    print(paste("Sample names do not match for", name_meta, "and", name_count))
    print(str_c("Overlap size between meta and profiles: ", sum(run_col %in% colnames(count_list[[i]]))))
  }
}
}  

align_samples(all_data$metadata, all_data$counts)

conditions <- all_data %>% 
  select(metadata, dataset_name) %>% 
  unnest(cols = c(metadata)) %>% 
  relocate(condition, .after = dataset_name) %>% 
  select(dataset_name, condition) %>% 
  table() %>% 
  as_tibble() %>% 
  filter(n != 0) 
  
summary_table <- conditions %>%
  group_by(condition) %>% 
  summarize(dataset_names = paste(dataset_name, collapse = ", "),
            n_values = paste(n, collapse = ", "),
            sum_n = sum(n)) 
summary_table <- summary_table %>% 
  add_row(condition = "all",
          dataset_names = "all",
          n_values = "all",
          sum_n = sum(summary_table$sum_n))


condition_colors <- c("Control" = "grey90", "Resolved HBV" = "grey69", "CHB" = "cyan", "CIR" = "cyan3", "HCC" = "dodgerblue4", "CIR+HCC" = "cyan4",
                      "AUD" = "yellow", "ALD" = "tan1", "AH" = "tan3", "HBVC" = "dodgerblue3", "NHBVC" = "dodgerblue", "NAFLD" = "darkgreen", 
                      "PBC" = "orangered", "PSC" = "red3")
all_conditions <- summary_table$condition[1:15]
all_conditions <- c(all_conditions[7], all_conditions[-7])
case_conditions <- all_conditions[-1]
###############
### 3. analysis
###############

all_data <- all_data %>% 
  #select(counts, dataset_name) %>% 
  mutate(counts_long = map(counts, function(x) return(x %>%
                                                            pivot_longer(-Taxa) %>%
                                                            rename(sampleID = name, count = value))))

#### count histograms of all samples ####

depth_histogram_data <- all_data %>%
  #select(-profiles) %>%
  select(dataset_name, counts_long) %>%
  unnest(cols = counts_long) %>%
  group_by(sampleID, dataset_name) %>%
  summarize(depth = sum(count))

samplesFiltered1 <- depth_histogram_data %>%
  filter(depth <= depth1) %>%
  group_by(dataset_name) %>%
  tally() %>%
  mutate(x = 1.25E5, y = 0)
#  mutate(x = 2E5, y = 25)
samplesFiltered2 <- depth_histogram_data %>%
  filter(depth <= depth2) %>%
  group_by(dataset_name) %>%
  tally() %>%
  mutate(x = 1.25E5, y = 0)
#  mutate(x = 2E5, y = 25)
samplesFiltered3 <- depth_histogram_data %>%
  filter(depth <= depth3) %>%
  group_by(dataset_name) %>%
  tally() %>%
  mutate(x = 1.25E5, y = 0)
#  mutate(x = 2E5, y = 25)

p1 <- ggplot() +
  #geom_histogram(data = depth_histogram_data, fill = "#852020",aes(x = depth)) +
   geom_histogram(data = depth_histogram_data, fill = "red",aes(x = depth)) +
  facet_grid(dataset_name~., scales = "free_y") +
  labs(x = "sequencing depth", y = "number of runs",
       title = ggtitle("Histogram of sequencing depth")) +
  xlim(-8000, 150000) +
  theme_light() +
  geom_vline(xintercept = depth1) +
  geom_vline(xintercept = depth2) +
  geom_vline(xintercept = depth3) +
  # geom_text(data = data.frame(depth1 = depth1), 
  #           aes(x = depth1 + 15000, y = 28, label = depth1),
  #           stat = "unique") + 
  annotate(geom = "text", x = depth1 - 2000, y = 100, label = depth1) +
  annotate(geom = "text", x = depth2 + 2000, y = 100, label = depth2) +
  annotate(geom = "text", x = depth3 + 2000, y = 100, label = depth3 ) +
  # geom_text(data = samplesFiltered1,
  #        aes(x = x, y = y, label = n)) +
  # geom_text(data = samplesFiltered2,
  #          aes(x = x, y = y - 3, label = n)) +
  geom_text(data = samplesFiltered1, nudge_y = 75,
            aes(x= x, y = y, label = str_c("Samples removed at depth1: ",
                                         n))) +
  geom_text(data = samplesFiltered2, nudge_y = 50,
            aes(x= x, y = y, label = str_c("Samples removed at depth2: ",
                                           n))) +
  geom_text(data = samplesFiltered3, nudge_y = 25,
              aes(x= x, y = y, label = str_c("Samples removed at depth3: ",
                                             n)))


p1


#### filter by sequencing depth ####
filtered_samples <- depth_histogram_data %>% 
  filter(depth >= depth1) %>% 
  select(sampleID) %>% 
  pull(sampleID)
all_data <- all_data %>% 
  mutate(counts_long_high_depth = map(counts_long, function(x) return(x %>% 
                                        filter(sampleID %in% filtered_samples))))

#### Subsampling + relative abundance####
all_data <- all_data %>%
  mutate(profiles_subsampled_hd = map(counts_long_high_depth, function(x) return(x %>%
                                        pivot_wider(id_cols = sampleID, 
                                                        names_from = Taxa, 
                                                       values_from = count) %>%
                                        column_to_rownames('sampleID') %>%
                                        as.matrix() %>%
                                        rrarefy(depth1) %>%
                                        as.data.frame() %>%
                                        rownames_to_column('sampleID') %>%
                                        pivot_longer(-sampleID) %>%
                                        rename(Taxa = name,count = value) %>%
                                        relocate(sampleID, Taxa) %>% 
                                        group_by(sampleID) %>% 
                                        mutate(relAb = proportions(count)))))

# join with metadata 
all_data <- all_data %>% 
  mutate(profiles_long_with_meta = map2(counts_long_high_depth, metadata, function(pro, met) {
    return(inner_join(pro, met, by = c('sampleID' = "Run")) %>%
             group_by(sampleID) %>%
             mutate(relAb = proportions(count)))
  })) %>% mutate(profiles_subsampled_with_meta = map2(profiles_subsampled_hd, metadata, function(pro, met) {
    return(inner_join(pro, met, by = c('sampleID' = "Run")))}))
  
# filter unresolved taxa 
 all_data <- all_data %>% 
 mutate(profiles_long_meta_resolved = map(profiles_long_with_meta, function (x) return(x %>% 
                                                          filter(Taxa != "not_resolved")))) %>% 
   mutate(profiles_subsampled_meta_resolved = map(profiles_subsampled_with_meta, function (x) return(x %>% 
                                                                                           filter(Taxa != "not_resolved"))))
# make df with metadata that we want to use for PCoA 
 meta_all <- all_data %>%
   select(profiles_long_meta_resolved, dataset_name) %>% 
   unnest() %>% select(sampleID, condition, dataset_name, Continent, Age, BMI, Sex, comorbidity) %>%  #add other meta-information if necessary
   distinct()
    

### get wide profile of all_data in one matrix
wide_matrix <- all_data %>%
  select(profiles_long_meta_resolved, dataset_name) %>% 
  unnest() %>% 
  pivot_wider(id_cols = Taxa, names_from = sampleID, values_from = relAb, values_fill = 0) %>%
  as.data.frame() %>% 
  mutate(Taxa = str_replace(Taxa, "-", '_')) %>%
  mutate(Taxa = str_replace(Taxa, "/", '_')) %>%
  mutate(Taxa = ifelse(Taxa == "beta proteobacterium Mena 23_3_3c[genus]", "beta proteobacterium", Taxa)) %>%
  column_to_rownames('Taxa')

### distance matrix calculation and turning into matrix
pairwise_distances <- vegdist(x = t(log10(wide_matrix + 1E-5)), method = "manhattan") # maybe insert log10 formula 
pcoa_object <- cmdscale(d = pairwise_distances, k = 2)
pcoa_object <- pcoa_object %>% 
  as.data.frame() %>% 
  rownames_to_column("sampleID") %>% 
  as_tibble()
head(pcoa_object)
# TODO: make pc <-  formula to round down to the 2nd decimal (0.08 becomes 0.01)

# Join metadata back (disease type, dataset name, ...)
plot_object_meta <- inner_join(pcoa_object, meta_all, by = "sampleID")
plot_object_meta <- plot_object_meta %>% 
  rename(data = dataset_name)#%>%
  # rename(v1 = V1) %>% 
  # rename(v2 = V2)
  age_breaks <- seq(0, max(plot_object_meta$Age, na.rm = TRUE) + 10, by = 10)
  
plot_object_meta <- plot_object_meta %>%  
  mutate(age_group = cut(Age, breaks = age_breaks, labels = 
                           paste(age_breaks[-length(age_breaks)], "-", age_breaks[-1])))

# Plot PCoA based on distance matrix with different color codings side-by-side
p2_condition <- plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, col = condition)) +
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette = "Set3") +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on condition")

p3_data <- plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, col = data)) +
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette = "Set3") +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on dataset")

# p4_LD <- plot_object_meta %>% 
#   ggplot(aes(x = v1, y = v2, color =  ifelse(condition == "Control", "Control",
#                                         ifelse(condition %in% c("Resolved HBV", "CHB"), "Liver disease", 
#                                           ifelse(condition %in% c("HCC", "CIR", "CIR_HCC"), "advanced liver disease", condition))))) +
#   geom_point() +
#   theme_classic() +
#   theme(legend.position = "bottom") +
#   # scale_color_manual(values = c("Control" = "#009999", "liver disease" = "#FF3300"))
#   scale_color_discrete(name = "condition", labels = c("Control", "Liver disease", "Advanced liver disease"))

p4_LD <- plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, color = case_when(
    condition %in% c("Control") ~ "Control",
    condition %in% c("Resolved HBV", "CHB", "ALD", "AUD", "ALD", "NAFLD", "PBC", "AH", "HBV") ~ "Liver disease",
    condition %in% c("HCC", "CIR", "CIR+HCC", "PSC", "HBVC", "NHBVC") ~ "Advanced liver disease",
    TRUE ~ as.character(condition)
  ))) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(color = "liver disease type" ) +
  # scale_color_discrete(name = "condition", labels = all_conditions) +
  # scale_color_manual(
  #   name = "condition",
  #   values = c("Control" = "cyan3",
  #              "Liver disease" = "coral1",
  #              "Advanced liver disease" = "forestgreen"),
  #   labels = c("Control", "Liver disease", "Advanced liver disease")
  # )
NULL

condition_groups <-  c(Control = "gray90", Cholangitis = "cyan2", Alcoholic = "chocolate2", Hepatitis = "deeppink2", Carcinoma = "darkolivegreen3")

p_condition_group <- plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, color = case_when(
    condition %in% c("Control") ~ "Control",
    condition %in% c("Resolved HBV", "CHB", "AH") ~ "Hepatitis", 
    condition %in% c("ALD", "AH", "AUD") ~ "Alcoholic",
    condition %in% c("PSC", "PBC") ~ "Cholangitis",
    condition %in% c("HBVC", "HCC", "NHBVC", "CIR+HCC") ~ "Carcinoma",
    condition %in% c("CIR", "CIR+HCC") ~ "Cirrhosis", 
    condition %in% c("NAFLD") ~ "NAFLD",
    TRUE ~ as.character(condition)
  ))) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on condition group") +
  scale_color_manual(values = condition_groups) +
  labs(color = "condition group" ) +
  NULL

p_etiology <- plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, color = case_when(
    condition %in% c("Control") ~ "Control",
    condition %in% c("Resolved HBV", "CHB", "HBVC") ~ "Viral", 
    condition %in% c("ALD", "AH", "AUD") ~ "Alcoholic",
    condition %in% c("PSC", "PBC") ~ "Autoimmune",
    condition %in% c("Cirrhosis", "HCC", "NAFLD", "NHBVC", "CIR+HCC") ~ "unknown",
    TRUE ~ as.character(condition)
  ))) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on etiology") +
  labs(color = "etiology" ) +
  # scale_color_discrete(name = "condition", labels = all_conditions) +
  # scale_color_manual(
  #   name = "condition",
  #   values = c("Control" = "cyan3",
  #              "Liver disease" = "coral1",
  #              "Advanced liver disease" = "forestgreen"),
  #   labels = c("Control", "Liver disease", "Advanced liver disease")
  # )
  NULL



p_Continent <- plot_object_meta %>% 
  ggplot(aes(x = 'V1', y = 'V2', col = Continent, rm.na = TRUE)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on geographic location")

# grid.arrange(p2_condition, p3_data, p4_LD, p_Continent, ncol = 2)
grid.arrange(p2_condition, p3_data, p_Continent, ncol = 2)

p_BMI <-  plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, col = BMI)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on BMI")

p_sex <-  plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, col = Sex)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on sex")

p_Age <-  plot_object_meta %>% 
  ggplot(aes(x = V1, y = V2, col = age_group)) +
  geom_point() + #aes(color = is.na(age_group))) +
#  scale_color_manual(values = c("TRUE" = "gray95")) +
  theme_classic() +
  theme(legend.position = "bottom") +
  ggtitle("distance analysis based on Age")


(p_Continent + p_BMI) / (p_sex + p_Age)
############# Richness & Diversity Analysis ##############

# use the vegan::diversity function
# important to check format
# TODO: rrafection is lost and should definitely be recovered. 
subsampled_matrix <- all_data %>%
  select(profiles_subsampled_meta_resolved, dataset_name) %>% 
  unnest() %>% 
  pivot_wider(id_cols = Taxa, names_from = sampleID, values_from = relAb, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames('Taxa') 
  

#### Richness #####
richness <- specnumber(t(subsampled_matrix))
richness <-  data.frame(sampleID = names(richness), value = as.numeric(richness))
richness <- richness %>% 
  as_tibble() %>% 
  rename(richness = value)
plot_object_meta <- inner_join(plot_object_meta, richness, by = "sampleID")

summary_richness <- plot_object_meta %>% 
  select(condition, richness) %>% 
  group_by(condition)

wilcox_richness <-  pairwise.wilcox.test(summary_richness$richness, 
                                 summary_richness$condition, p.adjust.method = "none") %>% 
  broom::tidy() %>%
  filter(group1 == "Control" | group2 == "Control") %>%
  mutate(condition = ifelse(group1 == "Control", group2, group1)) %>%
  select(-group1, -group2) %>%
  relocate(condition)

p7_richness <- plot_object_meta %>% 
  ggplot(aes(x = condition, y = richness, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors, breaks = all_conditions) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  theme_classic() +
  labs(x = "Condition", y = "Genus richness per sample") +
  scale_y_continuous(limits = c(30 ,90), breaks = seq(0, 180, by = 20)) +
  ggtitle("Genus richness by condition") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.position = "none") +
  scale_x_discrete(limits = all_conditions) +
  geom_text(data = wilcox_richness, aes(x = condition, y = 75, label = str_c("p = ", round(p.value, digits = 4)))) +
  NULL


#### shannon diversity ####
diversity_matrix <- t(subsampled_matrix)
shannon <-  diversity(diversity_matrix, "shannon")
shannon <- data.frame(sampleID = names(shannon), value = as.numeric(shannon))
shannon <- shannon %>% 
  as_tibble() %>% 
  rename(shannon = value)
  
  
# first, add shannon values to plot_object_meta(?)
plot_object_meta <- inner_join(plot_object_meta, shannon, by = "sampleID")


### simpson diversity #####
diversity_matrix <- t(subsampled_matrix)
simpson <-  diversity(diversity_matrix, "simpson")
simpson <- data.frame(sampleID = names(simpson), value = as.numeric(simpson))
simpson <- simpson %>% 
  as_tibble() %>% 
  rename(simpson = value)


# first, add simpson values to plot_object_meta(?)
plot_object_meta <- inner_join(plot_object_meta, simpson, by = "sampleID")

#### wilcoxon test of Controls vs. each Condition ####

summary_shannon <- plot_object_meta %>% 
  select(condition, shannon) %>% 
  group_by(condition)

wilcox1 <-  pairwise.wilcox.test(summary_shannon$shannon, 
                                 summary_shannon$condition, p.adjust.method = "none") %>% 
  broom::tidy() %>%
  filter(group1 == "Control" | group2 == "Control") %>%
  mutate(condition = ifelse(group1 == "Control", group2, group1)) %>%
  select(-group1, -group2) %>%
  relocate(condition)





#  map(plot_object_meta, function(x) return(x %>% 
#                                             select(condition, shannon) %>%
#                                             )) 


#### plot diversity ####
p6_simpson <- plot_object_meta %>% 
  ggplot(aes(x = condition, y = simpson, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors, breaks = all_conditions) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  theme_classic() +
  labs(x = "Condition", y = "Simpson index") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  ggtitle("Simpson diversity by condition") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.position = "none") +
  scale_x_discrete(limits = all_conditions) +
  NULL

p5_shannon <- plot_object_meta %>% 
  ggplot(aes(x = condition, y = shannon, fill = condition)) +
  geom_boxplot() +
  scale_fill_manual(values = condition_colors, breaks = all_conditions) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  theme_classic() +
  labs(x = "Condition", y = "Shannon index") +
  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, by = 0.5)) +
  ggtitle("Shannon diversity by condition") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5), 
        legend.position = "none") +
  scale_x_discrete(limits = all_conditions) +
  geom_text(data = wilcox1, aes(x = condition, y = 3.8, label = str_c("p = ", round(p.value, digits = 4)))) +
  NULL

 

 p7_richness / p5_shannon 


##############
## 4. SIAMCAT
#############
#  filter(dataset_name == c("Lin2023", "Yun2019", "Ni2019"))
meta_all_sc <- column_to_rownames(meta_all, var = 'sampleID')
# meta_all_sc <- column_to_rownames(meta_all, var = 'sampleID')
#TODO: create label and use wide/matrix to create siamcat object
sc_label <- create.label(meta = meta_all_sc, label = 'condition',
                         case = case_conditions, control = 'Control')
sc.obj <- siamcat(feat = wide_matrix, label = sc_label, meta = meta_all_sc)
show(sc.obj)
sc.obj <- filter.features(sc.obj, filter.method = 'abundance', cutoff = 0.001)
sc.obj <- check.associations(sc.obj, log.n0 = 1e-06, alpha = 0.05)
association.plot(sc.obj, sort.by = 'fc', 
                 panels = c('fc', 'prevalence', 'auroc'), max.show = 40)
# sc_label_hcc <- create.label(meta = meta_first_sc, label = 'condition', 
#                              case = 'HCC', control = 'Control')
# sc.obj_hcc <- siamcat(feat = wide_matrix, label = sc_label_hcc, meta = meta_first_sc) 
# sc.obj_hcc <- filter.features(sc.obj_hcc, filter.method = 'abundance', cutoff = 0.001)
# sc.obj_hcc <- check.associations(sc.obj_hcc, log.n0 = 1e-06, alpha = 0.05)
#association.plot(sc.obj_hcc, sort.by = 'fc', 
#                 panels = c('fc', 'prevalence', 'auroc'))
#TODO: ggsave for plots to use for retreat presentaiton

# ggsave(filename = here("../SIAMCAT_retreat.png"), width = 10, height = 5)

# wide_matrix %>% 
#   as_tibble(rownames = 'bac') %>% 
#   pivot_longer(-bac,names_to = 'Sample_ID',values_to = 'relAb') %>% 
#   filter(str_detect(bac,pattern = 'upriav|epidi|ordonia|esorhiz|zcobacter')) %>% 
#   ggplot(aes(x=bac,y=log10(relAb+1e-5)))+
#   geom_point()




# SIAMCAT analysis of single datasets
# Schritt 1: SIAMCAT für jeden Datensatz einzeln (map() function kann hier wieder behilflich sein)
  # Schritt 1a wide_matrix fuer jeden Daatensatz einzeln erstellen und in tibble einspeichern
# Schritt 2: Assoziierte Features (test results) extrahieren associations_hcc <- associations(sc.obj_hcc) %>% 
# + rownames_to_column('Taxa') %>% 
#   + select(Taxa, fc)
# > associations_hcc
# Schritt 3: jede Matrix ins long format bringen, Genus, Datensatz, p-value, fold-change
# Scrhitt 4: diese alle mit rbind() zusammenfuegen
# ______________
# Schritt 5 : heatmap dann machen
all_data <- all_data %>%
  #select(profiles_long_meta_resolved, dataset_name) %>%
  mutate(wide_matrices = map(profiles_long_meta_resolved, function(x) return(x %>%
                                                                               pivot_wider(id_cols = Taxa, names_from = sampleID, values_from = relAb, values_fill = 0) %>%
                                                                               as.data.frame() %>%
                                                                               column_to_rownames('Taxa')))) 
# split_tibbles <- split(as.data.frame(meta_all_sc), meta_all_sc$dataset_name)
# SIAMCAT_data <- SIAMCAT_data %>% 
#   mutate(meta_data = split_tibbles) 

all_data <- all_data %>%
  filter(dataset_name != "Furukawa2020") %>% 
  mutate(sc_label = map(metadata, function(x)  {
    x <- x %>%
      mutate(condition = case_when(condition == "Control" ~ "Control",
                                   condition != "Control" ~ "Case")) %>%
      select(Run, condition) %>%
      as.data.frame() %>% 
      column_to_rownames('Run') 
    
    #print(table(x$condition))
    #return(x)
    # label <- x$condition
    # case <- setdiff(unique(label), "Control")
    # control <- "Control"
    # label <- setNames(label, label)
    # case <- setNames(case, case)
    # control <- setNames(control, control)
    # create.label(meta = x$metadata, label = label, case = case, control = control)
    return(create.label(meta = x, label = 'condition', case = 'Case', control = "Control"))
  }))

all_data <- all_data %>%
  mutate(metadata = map(metadata, function(x) return(x %>% as.data.frame() %>% select(Run, condition) %>% column_to_rownames('Run')))) %>%
  #mutate(siamcat = map(.x = wide_matrices, .f = ~ siamcat(feat = .x$wide_matrices, label= sc_label, meta = .x$metadata)))
  mutate(siamcat = pmap(list(wide_matrices, sc_label, metadata), function(mat, lab, meta) return(siamcat(feat = mat, meta = meta, label = lab))))

label <- 'condition'

all_data <- all_data %>%
  mutate(siamcat = map(siamcat, function(x) {
    x <- x %>%
      filter.features() %>% 
      check.associations(log.n0 = 1e-06, alpha = 0.05)
    return(x)
  })) %>%
  mutate(associations = map2(siamcat, dataset_name, function(s, n) {
    associations(s) %>% 
      rownames_to_column('Taxa') %>% 
      select(Taxa, fc, p.adj) %>% 
      filter(Taxa != 'Bacteria') %>% 
      mutate(dataset_name = n) 
  }))

p_vals_selected <- all_data %>%
  select(dataset_name, associations) %>%
  unnest() %>%
  pivot_wider(id_cols = Taxa, names_from = dataset_name, values_from = p.adj)
p_vals_selected <- p_vals_selected[apply(p_vals_selected[, 2:dim(p_vals_selected)[2]], 1, function(x) {
  #x <- x[2:length(x)]
  #print(x)
  #adsads
  if(all(is.na(x))) {
    return(FALSE)
  }
  x <- x[!is.na(x)]
  return(any(x < 0.1))
}), ] 

fcs_selected <- all_data %>%
  select(dataset_name, associations) %>%
  unnest() %>%
  pivot_wider(id_cols = Taxa, names_from = dataset_name, values_from = fc) %>%
  inner_join(p_vals_selected %>% select(Taxa))


p_vals_selected <- p_vals_selected %>%
  column_to_rownames('Taxa') %>%
  as.matrix()

fcs_selected <- fcs_selected %>%
  column_to_rownames('Taxa') %>%
  as.matrix()

p_vals_selected_asterisk <- p_vals_selected
p_vals_selected_asterisk <- ifelse(is.na(p_vals_selected_asterisk), 100000, p_vals_selected_asterisk)
p_vals_selected_asterisk <- ifelse(p_vals_selected_asterisk < 0.1, "*", "")


make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}





heatmap_groups = data.frame(Groups = as.factor(c("Cholangitis", "Alcoholic", "Cholangitis", "Hepatitis", "Carcinoma", "Carcinoma", "Hepatitis", "Carcinoma")))
ann_colors = list(
  Groups = c(Cholangitis = "cyan2", Alcoholic = "chocolate2", Hepatitis = "deeppink2", Carcinoma = "darkolivegreen3"))
rownames(heatmap_groups) <- colnames(fcs_selected)
heatmap_col_label <- c( "Cortez (n = 27)", "Duan (n = 64)", "Iwasawa (n = 39)", "Lin (n = 149)", "Liu (n = 90)", "Ni (n = 109)", "Wang (n = 263)","Yun2019 (n = 37)")
p8_heatmap <- pheatmap(fcs_selected, cluster_rows = F, cluster_cols = F, annotation_col = heatmap_groups, annotation_colors = ann_colors, 
                       main = "differentially abundant taxa in LD", cellwidth = 20, na_col = "lightgrey", labels_col = heatmap_col_label,
                       angle_col = 45, labels_row = make_bold_names(fcs_selected, rownames, c("Veillonella", "Fusobacterium", "Enterococcus", "Prevotella", "Parabacteroides",
                                                                                              "Streptococcus", "Lactobacillus", "Eubacterium", "Akkermansia")), display_numbers = p_vals_selected_asterisk)

boxplot_data <- all_data %>%
  select(dataset_name, profiles_long_meta_resolved) %>%
  unnest() %>%
  select(dataset_name, Taxa, relAb, sampleID, condition) %>%
  mutate(condition = case_when(condition == "Control" ~ "Control",
                               condition != "Control" ~ "Case")) %>%
  filter(Taxa != "Bacteria") 
  
plot_stuff_boxplot <- function(boxplot_data, genus = NULL) {
  tmp <- boxplot_data %>% filter(Taxa == 'Fusobacterium')
  p <- ggplot() +
    geom_boxplot(data = tmp, aes(x = dataset_name, y = relAb + 1E-5, color = condition), outlier.color = NA, position = 'dodge', alpha = 0.3) +
    geom_point(data = tmp, aes(x = dataset_name, y = relAb + 1E-5, color = condition, group = condition), position = position_jitterdodge()) +
    theme_classic() +
    scale_y_log10()
  return(p)
}


get_quantile_plot <- function(inputData, axisColumn, labelColumn, valueColumn, expectedNumLevels = 3, xlab = 'Genera', ylab = "Relative Abundances (log10)") {
  
  colnames(inputData)[colnames(inputData) == axisColumn] <- "genus"
  colnames(inputData)[colnames(inputData) == labelColumn] <- "label"
  colnames(inputData)[colnames(inputData) == valueColumn] <- "relAb"
  
  
  print(head(inputData))
  #groupColors <- c("red", "blue")
  #quantiles <- c(0.5, 0.7, 0.9, 0.95)
  quantileData <- list()
  #SSS <- inputData
  # inputData <- inputData %>%
  #  filter(type == "reference (CV)")
  # medians <- inputData %>%
  # group_by(genus, label, plotGroup) %>%
  # summarize(median = median(relAb))
  for (quantile in quantiles) {
    tmp <- inputData %>% 
      #select(cyl, wt) %>% 
      group_by(genus, label, plotGroup) %>%
      #mutate(cyl) %>%
      summarize(value_max = quantile(relAb, probs = quantile),
                value_min = quantile(relAb, probs = 1-quantile)) %>%
      mutate(quantile = quantile)
    quantileData[[length(quantileData) + 1]] <- tmp
  }
  
  quantileData <- do.call('rbind', quantileData) %>%
    #left_join(medians %>% select(-plotGroup), by = c("genus", "label")) %>%
    mutate(genus = as.factor(genus)) %>%
    arrange(desc(quantile))  %>%
    #mutate(Quantiles = map2_chr(quantile, label, function(x, y) return(str_c((1-x) * 100, '% - ', x * 100, "% - ", y))))
    mutate(Quantiles = map2_chr(quantile, label, function(x, y) {
      return(str_c((x * 100), "% - ", y))
    }))
  
  # print(head(quantileData))
  
  
  # print(head(quantileData))
  # return(quantileData)
  l <- quantileData %>%
    ungroup() %>%
    select(Quantiles, quantile, label) %>%
    distinct() %>%
    group_by(label) %>%
    nest() %>%
    mutate(data = map(data, function(x) return(x %>% arrange(quantile)))) %>%
    unnest() %>%
    # arrange(quantile) %>%
    ungroup() %>%
    pull(Quantiles)
  
  names(colorVec) <- l
  
  quantileData <- quantileData %>%
    mutate(Quantiles = factor(Quantiles, levels = (l))) %>%
    arrange(desc(quantile))
  
  
  labelMap <- levels(quantileData$genus)  
  names(labelMap) <- 1:length(labelMap)
  
  #print(levels(quantileData$Quantiles))
  print(length(unique(quantileData$label)))
  print(expectedNumLevels)
  if (length(unique(quantileData$label)) != expectedNumLevels) {
    asdaddads
  }
  
  groupLevels <- unique(quantileData$label)
  
  colorVec <- c(colorRampPalette(c(groupColors[1], "white"))(length(quantiles)),
                colorRampPalette(c(groupColors[2], "white"))(length(quantiles)))
  names(colorVec) <- l
  # print(head(quantileData))
  
  p <- ggplot(data = quantileData %>%
                mutate(Quantiles = factor(Quantiles, levels = l))) +
    geom_rect(data = quantileData %>%
                filter(label == groupLevels[1]) %>%
                mutate(Quantiles = factor(Quantiles, levels = l)),
              aes(xmin = as.integer(genus) -0.1 - 0.15, xmax = as.integer(genus) + 0.1 - 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
    geom_rect(data = quantileData %>%
                filter(label == groupLevels[2]) %>%
                mutate(Quantiles = factor(Quantiles, levels = l)),
              aes(xmin = as.integer(genus) -0.1 + 0.15, xmax = as.integer(genus) + 0.1 + 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
    geom_point(data = quantileData %>%
                 filter(label == groupLevels[1]) %>%
                 filter(str_detect(Quantiles, "50"))%>%
                 mutate(Quantiles = factor(Quantiles, levels = l)),
               aes(x = as.integer(genus)-0.15, y = value_max), fill = 'darkgreen', size = 3, pch = 23, color = 'white') +
    geom_point(data = quantileData %>%
                 filter(label == groupLevels[2]) %>%
                 filter(str_detect(Quantiles, "50"))%>%
                 mutate(Quantiles = factor(Quantiles, levels = l)),
               aes(x = as.integer(genus)+0.15, y = value_max), fill = 'darkgreen', size = 3, pch = 23, color = 'white') +            
    scale_x_continuous(breaks = as.integer(names(labelMap)),  labels = labelMap) +
    #scale_fill_manual(breaks = c("Dark","DarkLight","Medium","LightDark","Light"),
    #                values=c("red", "orange","yellow","cadetblue2","dodgerblue"))
    scale_fill_manual(values = colorVec[!str_detect(names(colorVec), '50')], breaks = names(colorVec[!str_detect(names(colorVec), '50')])) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    #coord_flip() +
    facet_wrap(~plotGroup, nrow = 2, scales = "free") +
    xlab(xlab) +
    ylab(ylab) +
    NULL 
  return(p)
}

quantiles <- c(0.5, 0.6, 0.7, 0.8, 0.9)
#groupColors <- c("#852020", '#21520e', "#090938")
groupColors <- c("#852020", "#090938")
colorVec <- c(colorRampPalette(c(groupColors[1], "white"))(length(quantiles)),
              colorRampPalette(c(groupColors[3], "white"))(length(quantiles)),
              colorRampPalette(c(groupColors[2], "white"))(length(quantiles)))
colorVecGreen <- colorRampPalette(c(groupColors[2], "white"))(length(quantiles))

p_fuso <- get_quantile_plot(boxplot_data %>%
                    filter(Taxa == "Fusobacterium") %>%
                    mutate(plotGroup = 'Fusobacterium') %>%
                    mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_veillonella <- get_quantile_plot(boxplot_data %>%
                    filter(Taxa == "Veillonella") %>%
                    mutate(plotGroup = 'Veillonella') %>%
                    mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)


p_enterococcus <- get_quantile_plot(boxplot_data %>%
                                     filter(Taxa == "Enterococcus") %>%
                                     mutate(plotGroup = 'Enterococcus') %>%
                                     mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)
p_prevotella <-  get_quantile_plot(boxplot_data %>%
                                     filter(Taxa == "Prevotella") %>%
                                     mutate(plotGroup = 'Prevotella') %>%
                                     mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_eubac <-  get_quantile_plot(boxplot_data %>%
                                     filter(Taxa == "Eubacterium") %>%
                                     mutate(plotGroup = 'Eubacterium') %>%
                                     mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_akkermansia <-  get_quantile_plot(boxplot_data %>%
                                     filter(Taxa == "Akkermansia") %>%
                                     mutate(plotGroup = 'Akkermansia') %>%
                                     mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_lactobac <-  get_quantile_plot(boxplot_data %>%
                                   filter(Taxa == "Lactobacillus") %>%
                                   mutate(plotGroup = 'Lactobacillus') %>%
                                   mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_streptoc <-  get_quantile_plot(boxplot_data %>%
                                   filter(Taxa == "Streptococcus") %>%
                                   mutate(plotGroup = 'Streptococcus') %>%
                                   mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)

p_parabac <-  get_quantile_plot(boxplot_data %>%
                                   filter(Taxa == "Parabacteroides") %>%
                                   mutate(plotGroup = 'Parabacteroides') %>%
                                   mutate(relAb = log10(relAb + 1E-5)), 'dataset_name', 'condition', 'relAb', expectedNumLevels = 2)


(p_veillonella + p_fuso + p_eubac) / (p_prevotella + p_enterococcus + p_akkermansia) / (p_streptoc + p_lactobac + p_parabac) + 
  plot_layout(guides = 'collect')
