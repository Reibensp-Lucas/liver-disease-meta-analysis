library(tidyverse)
library(here)

##### Ni2019 #######
Ni2019_meta <- read_tsv(here('..', 'data', 'metadata', "Ni2019_meta.tsv"))
colnames(Ni2019_meta)
Ni2019_meta_tidy <- Ni2019_meta %>% 
  select("Run", "AssayType", "AvgSpotLen", "Bases", "BioProject", "BioSample", "CenterName", "GeoLocNameCountry", "GeoLocNameCountryContinent", "Host", "Instrument", "LibraryName", "LibraryLayout", "Organism", "Platform", "SampleName", "SraStudy") %>%
  rename(condition = LibraryName) %>% 
  mutate(condition = gsub("\\d+", "", condition)) %>%
  relocate(SraStudy, .before = CenterName) %>% 
  mutate(comorbidity = case_when(
    condition == "Hypertension" ~ "Hypertension",
    condition == "Diabetes" ~ "Diabetes",
    grepl("HCC", condition) ~ case_when(
      grepl("hyper", condition, ignore.case = TRUE) ~ "Hypertension",
      grepl("diab", condition, ignore.case = TRUE) ~ "Diabetes",
      grepl("hyperdiab", condition, ignore.case = TRUE) ~ "Hypertension, Diabetes",
      TRUE ~ NA_character_
    ),
    TRUE ~ NA_character_
  )) %>% 
  mutate(condition = gsub("HCC.*", "HCC", condition)) %>% 
  mutate(condition = gsub("Hypertension", "Control", condition)) %>% 
  mutate(condition = gsub("Diabetes", "Control", condition)) %>% 
  mutate(condition = gsub("Health", "Control", condition)) %>% 
  relocate(comorbidity, .before = LibraryLayout) %>% 
  rename(Country = GeoLocNameCountry) %>% 
  rename(Continent = GeoLocNameCountryContinent)
#made the condition "Hypertension" and "Diabetes" to controls and saved it into the new "comorbidity" column

write_tsv(Ni2019_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Ni2019_meta_tidy.tsv"))



##### Yun2019 ####
Yun2019_meta <- read_tsv(here('..', 'data', 'metadata', "Yun2019_meta.tsv"))
colnames(Yun2019_meta)
Yun2019_meta_tidy <- Yun2019_meta %>% 
  select("Run", "AssayType", "AvgSpotLen", "Bases", "BioProject", "BioSample", "CenterName", "GeoLocNameCountry", "GeoLocNameCountryContinent", "Instrument", "LibraryLayout", "Organism", "Platform", "SampleName", "SampleName_2","SraStudy") %>%
  rename(condition = SampleName_2) %>% 
  mutate(condition = gsub("\\d+", "", condition)) %>%
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(condition, .before = LibraryLayout) %>% 
  mutate(comorbidity = case_when(
    condition == "Hypertension" ~ "Hypertension",
    condition == "Diabetes" ~ "Diabetes",
    grepl("HCC", condition) ~ case_when(
      grepl("hyper", condition, ignore.case = TRUE) ~ "Hypertension",
      grepl("diab", condition, ignore.case = TRUE) ~ "Diabetes",
      grepl("hyperdiab", condition, ignore.case = TRUE) ~ "Hypertension, Diabetes",
      TRUE ~ NA_character_
    ),
    TRUE ~ NA_character_
  )) %>% 
  relocate(comorbidity, .before = LibraryLayout) %>% 
  mutate(condition = gsub("CN", "Control", condition)) %>% 
  rename(Country = GeoLocNameCountry) %>% 
  rename(Continent = GeoLocNameCountryContinent)

write_tsv(Yun2019_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Yun2019_meta_tidy.tsv"))


##### Ren2017 ########

#Ren2017_meta <- read_tsv(here('..', 'data', "Ren2017_meta.tsv"))
## link to the Condition is not possible from this metadata, unfortunately. Will have to omit this study.

##### Lin2023 ############
Lin_2023_meta <- read.csv(here('..', 'data', 'metadata',"Lin_2023_SraRunTable.txt"))
colnames(Lin_2023_meta)                          
Lin_2023_meta_tidy <- Lin_2023_meta %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "geo_loc_name_country", "geo_loc_name_country_continent", "Instrument", "LibraryLayout", "Organism", "Platform", "Sample.Name","SRA.Study") %>% 
  rename(condition = Sample.Name) %>% 
  mutate(condition = gsub("Healthy.*", "Control", condition)) %>% 
  mutate(condition = gsub("Chron.*", "CHB", condition)) %>% 
  mutate(condition = gsub("Advan.*", "CIR+HCC", condition)) %>% 
  mutate(condition = gsub("Resolv.*", "Resolved HBV", condition)) %>% 
  rename(SraStudy = SRA.Study) %>% 
  rename(CenterName = Center.Name) %>% 
  rename(Country = geo_loc_name_country) %>% 
  rename(Continent = geo_loc_name_country_continent) %>% 
  rename(AssayType = Assay.Type) %>% 
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(condition, .before = LibraryLayout)
write_tsv(Lin_2023_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Lin2023_meta_tidy.tsv"))

######### Duan2019 ######
Duan_2019_meta  <- read.csv(here('..', 'data', 'metadata', "Duan_2019_meta.txt"), na.strings = c("", "missing"))
colnames(Duan_2019_meta)
Duan_2019_meta_tidy <- Duan_2019_meta %>% 
  filter(HOST == 'Homo sapiens') %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name",
         "Instrument", "LibraryLayout", "Organism", "Platform","SRA.Study",
         "host_sex", "Host_age", "Host_disease", "host_BMI", "sample_source") %>% 
  rename(condition = Host_disease) %>% 
  mutate(Country = case_when(
    sample_source == NA ~ sample_source,
    sample_source %in% c("Spain", "Belgium", "France", "Mexico", "United Kingdom", "USA") ~ sample_source,
    grepl("University", sample_source, ignore.case = TRUE) ~ "USA")) %>% 
   mutate(Continent = case_when(
     Country  %in% c("Spain", "Belgium", "France", "United Kingdom") ~ "Europe",
     Country  %in% c("Mexico", "USA") ~ "North America")) %>% 
  rename(SraStudy = SRA.Study) %>% 
  rename(CenterName = Center.Name) %>%  
  rename(AssayType = Assay.Type) %>% 
  rename(Sex = host_sex) %>% 
  rename(BMI = host_BMI) %>% 
  rename(Age = Host_age) %>% 
  mutate(condition = case_when(
    condition == "alcoholic cirrhosis" ~"ALD",
    condition == "alcoholic liver disease" ~ "ALD",
    condition == "alcohol use disorder" ~ "AUD",
    condition == "alcoholic hepatitis" ~ "AH",
    condition == "healthy liver control" ~ "Control")) %>% 
#  mutate(BMI = case_when(BMI == c("", "missing") ~ NA)) %>% 
  mutate(Age = gsub("*.years", "", Age)) %>% 
#  mutate(Sex = case_when(Sex == "missing" ~ NA)) %>%
#  mutate(Age = case_when(Age == "missing" ~ NA)) %>%
  mutate_at(c('BMI', 'Age'), as.numeric) %>% 
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(condition, .before = LibraryLayout) %>% 
  relocate(Continent, .before = Instrument) %>% 
  relocate(Country, .before = Continent)
write_tsv(Duan_2019_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Duan2019_meta_tidy.tsv"))

run_ids_Duan <- Duan_2019_meta_tidy$Run

####### Iwasawa_2017 ##### 
Iwasawa2017 <-  read.csv(here('..', 'data', 'metadata', "Iwasawa_2017_meta.txt"), na.strings = c("", "missing"))
colnames(Iwasawa2017)
Iwasawa2017_tidy <-  Iwasawa2017 %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "geo_loc_name_country",
         "geo_loc_name_country_continent", "Host_age", "host_disease_stat", "host_sex",  "Instrument", "LibraryLayout", 
         "Organism", "Platform", "Sample.Name","SRA.Study", "chem_administration") %>% 
  rename(condition = host_disease_stat, SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, 
         Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Age = Host_age, 
         Sex = host_sex, SampleName = Sample.Name, Intervention = chem_administration) %>%
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(condition, .before = LibraryLayout) %>% 
  relocate(Continent, .before = Instrument) %>% 
  relocate(Country, .before = Continent) %>% 
  relocate(Age, .before = SampleName) %>% 
  relocate(Sex, .before = Age) %>% 
  relocate(Intervention, .before = SampleName) %>% 
  mutate(condition = case_when(
    condition == "Primary Sclesrosing Cholangitis" ~ "PSC",
    condition == "Ulcerative Colitis" ~ "UC", 
    condition == "healthy control" ~ "Control")) 
write_tsv(Iwasawa2017_tidy, here('..', 'data', 'metadata', 'tidy_data', "Iwasawa2017_meta_tidy.tsv"))


####### Cortez_2021 #####
Cortez2021 <-  read.csv(here('..', 'data', 'metadata', "Cortez_2021_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
colnames(Cortez2021)
Cortez2021_tidy <-  Cortez2021 %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "geo_loc_name_country",
         "geo_loc_name_country_continent", "Host_age", "Host_disease", "host_sex",  "Instrument", "LibraryLayout", 
         "Organism", "Platform", "Sample.Name","SRA.Study") %>%
  rename(condition = Host_disease, SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, 
         Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Age = Host_age, 
         Sex = host_sex, SampleName = Sample.Name) %>%
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(condition, .before = LibraryLayout) %>% 
  relocate(Continent, .before = Instrument) %>% 
  relocate(Country, .before = Continent) %>% 
  relocate(Age, .before = SampleName) %>% 
  relocate(Sex, .before = Age) %>%
  mutate(condition = case_when(
    condition == "None" ~ "Control", 
    condition == "Primary Sclerosing Cholangitis" ~ "PSC", 
    condition == "Ulcerative Colitis" ~ "UC", 
    condition == "Primary Sclerosing Cholangitis + Ulcerative Colitis" ~ "PSC+UC")) 
write_tsv(Cortez2021_tidy, here('..', 'data', 'metadata', 'tidy_data', "Cortez2021_meta_tidy.tsv"))


###### Furukawa_2020 #######
Furukawa2020 <-  read.csv(here('..', 'data', 'metadata', "Furukawa_2020_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
colnames(Furukawa2020)
Furukawa2020_tidy <- Furukawa2020 %>% 
  select("Run", "Assay.Type", "body_mass_index", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "geo_loc_name_country",
         "geo_loc_name_country_continent", "Age", "host_disease_stat", "sex",  "Instrument", "LibraryLayout", 
         "Organism", "Platform", "Sample.Name","SRA.Study", "chem_administration", "gastrointest_disord", "host_height") %>% 
  rename(comorbidity = host_disease_stat, SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, 
         Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, Sex = sex, SampleName = Sample.Name, 
         Height = host_height, Intervention = chem_administration, BMI = body_mass_index, Surgery = gastrointest_disord) %>% 
  mutate(condition = str_extract(comorbidity, "^[^,]+"))%>% 
  mutate(condition = gsub("primary biliary cirrhosis", "PBC", condition)) %>% 
  mutate(comorbidity = gsub("primary biliary cirrhosis", "", comorbidity)) %>%
  mutate(comorbidity = na_if(trimws(comorbidity), "")) %>%
  mutate(comorbidity = ifelse(!is.na(comorbidity), str_sub(comorbidity, 2), comorbidity)) %>% 
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(comorbidity, .before = LibraryLayout) %>% 
  relocate(condition, .before = comorbidity) %>% 
  relocate(Continent, .before = Instrument) %>% 
  relocate(Country, .before = Continent) %>% 
  relocate(Height, .before = Surgery) %>% 
  relocate(BMI, .before = Height) %>% 
  relocate(Age, .before = BMI) %>% 
  relocate(Sex, .before = Age) %>% 
  relocate(SampleName, .after = Height) %>% 
  relocate(Intervention, .before = SampleName) %>%
  relocate(Surgery, .before = Intervention)
write_tsv(Furukawa2020_tidy, here('..', 'data', 'metadata', 'tidy_data', "Furukawa2020_meta_tidy.tsv"))

##### Wang_2017 ##### 
Wang2017_meta <-  read.csv(here('..', 'data', 'metadata', "Wang_2017_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
Wang_pat_meta <-  readxl::read_excel(here('..', 'data', 'metadata', "Wang_2017_patient_metadata.xls"), col_names = T)
colnames(Wang2017_meta)
Wang_pat_meta_new <-  Wang_pat_meta
colnames(Wang_pat_meta_new) <- Wang_pat_meta[1, ]  
Wang_pat_meta_new <-  Wang_pat_meta_new[-1, ]
Wang_pat_meta_new <- Wang_pat_meta_new %>% 
  select(-1) %>% 
  rename(condition = Disease, ChildPughScore = Child.Pugh.score, ChildPughLevel = Child.Pugh.level, 'Run' = 'SRA accession number', 'Gender' = 'Gender (1=m)') %>% 
  relocate(Run, .before = condition) %>% 
  mutate(Sex = case_when(
    Gender == 1 ~ "male", 
    Gender == 2 ~ "female")) %>% 
  select(-Gender) %>% 
  mutate(condition = case_when(
    condition == "healthy" ~ "Control",
    condition == "HBV" ~ "HBV")) %>% 
  mutate_at(c('Age', 'PT', 'ALT', 'AST', 'GGT', 'TBIL', 'DBIL', 'IDBIL', 'ALB', 'ChildPughScore', 'BMI'), as.numeric)

Wang2017_meta_tidy <- Wang2017_meta %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "geo_loc_name_country",
         "geo_loc_name_country_continent", "Host_disease", "Instrument", "LibraryLayout", 
         "Organism", "Platform", "Sample.Name","SRA.Study") %>% 
  rename(condition = Host_disease, SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, 
         Country = geo_loc_name_country, Continent = geo_loc_name_country_continent, SampleName = Sample.Name) %>% 
  mutate(condition = case_when(
    condition == "healthy" ~ "Control",
    condition == "HBV" ~ "HBV",
    condition == "NAFLD" ~ "NAFLD"))
Wang2017_meta_tidy <- left_join(Wang2017_meta_tidy, Wang_pat_meta_new, by = "Run") %>% 
  select(-condition.y) %>% 
  rename(condition = condition.x)  %>% 
  relocate(SraStudy, .before = CenterName) %>% 
  relocate(Instrument, .before = condition) %>%
  relocate(Sex, .before = Age) %>% 
  relocate(BMI, .after = Age) %>% 
  relocate(SampleName, .after = ChildPughLevel)
write_tsv(Wang2017_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Wang2017_meta_tidy.tsv"))  

##### Liu2019 #########
Liu2019_meta <-  read.csv(here('..', 'data', 'metadata', "Liu2019_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
colnames(Liu2019_meta)
Liu2019_meta_tidy <-  Liu2019_meta %>% 
  select(Run, Assay.Type, AvgSpotLen, Bases, BioProject, BioSample, Center.Name, disease_state, Instrument, Organism, Platform, Sample.Name, SRA.Study, LibraryLayout) %>% 
  rename(condition = disease_state, SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, SampleName = Sample.Name) %>% 
  mutate(condition = case_when(
    condition == "healthy control" ~ "Control",
    condition == "HBVC" ~ "HBVC",
    condition == "NHBVC" ~ "NHBVC")) %>% 
  mutate(Country = "China") %>% 
  mutate(Continent = "Asia") %>% 
  mutate(CenterName = "Nanjing University Hospital") %>% 
  relocate(SraStudy, .after = BioSample) %>% 
  relocate(condition, .before = LibraryLayout) %>% 
  relocate(Continent, .before = Instrument) %>% 
  relocate(Country, .before = Continent) %>% 
  relocate(Platform, .after = LibraryLayout) %>% 
  relocate(Organism, .before = Platform) %>% 
  relocate(SampleName, .after = Platform)
write_tsv(Liu2019_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Liu2019_meta_tidy.tsv")) 


###### Qin2014 ######### 
Qin2014_meta  <- read.csv(here('..', 'data', 'metadata', "Qin2014_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
colnames(Quin2014_meta)
Qin2014_clinical <-  readxl::read_excel(here('..', 'data', 'metadata', "Quin2014_clinical_meta.xlsx"), na = c("", "missing", "NA", "not_applicable","None", "-"))
colnames(Quin2014_clinical)                   
Qin2014_clinical_new <- Qin2014_clinical %>% 
  select(-Stage, -...22, -...23) %>% 
  rename(SampleName = `Sample ID`, BMI = `BMI (kg/m2)`, Alb = `Alb (g/L)`, TB = `TB (umol/L)`, PT = `PT (S)`, 
         Cirrhosis = `Cirrhotic(Y or N)`, HBV = `HBV related  (Y or N)`, Alcoholic = `Alcohol related (Y or N)`,
         Other = `Other causes related`) 

map_other_to_abbr <- function(value) {
  if (value == "Hepatitis E Virus related") {
    return("HEV")
  } else if (value == "Hepatitis C Virus related") {
    return("HCV")
  } else if (value == "Hepatitis D Virus related") {
    return("HDV")
  } else if (value == "autoimmune related") {
    return("Autoimmune")
  } else if (value == "shistosoma related") {
    return("Shistosoma")
  } else if (value == "schistosoma，Hepatitis E virus related") {
    return("Shistosoma+HEV")
  } else if (value == "hepatolenticular degeneration related") {
    return("WD")
  } else if (value == "primary biliary cirrhosis & autoimmune related") {
    return("PBC+Autoimmune")
  } else if (value == "schistosoma，Hepatitis E virus and autoimmune related") {
    return("Shistosoma+HEV+Autoimmune")  
  } else {
    return("")
  }
}

Qin2014_clinical_new <- Qin2014_clinical_new %>% 
  rowwise() %>%
  mutate(
    condition = ifelse(Cirrhosis == "Y", "Cirrhosis", "Control")) %>%
  mutate(etiology = 
      paste(ifelse(HBV == "Y", "HBV", ""), collapse = "+") %>%
      paste(ifelse(Alcoholic == "Y", "Alcoholic", ""), collapse = "+") %>%
      paste(ifelse(Other == "autoimmune related", "Autoimmune", ""), collapse = "+") %>%
      paste(map_other_to_abbr(Other), collapse = "+")
  )
Qin2014_clinical_new <- Qin2014_clinical_new %>% 
  select(-Cirrhosis, -HBV, -Alcoholic, -Other) %>% 
  as_tibble()


Qin2014_meta_tidy <-  Qin2014_meta %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "Instrument", "LibraryLayout", 
          "Organism", "Platform", "Sample.Name", "Sample_Name","SRA.Study") %>% 
  rename( SraStudy = SRA.Study, AssayType = Assay.Type, CenterName = Center.Name, SampleName = Sample.Name) %>% 
  mutate(SampleName = gsub("*._Run*", "", SampleName)) %>% 
  mutate(SampleName = gsub("([A-Z]+)([0-9]+)", "\\1-\\2", SampleName)) # %>% 
  # filter(!duplicated(SraStudy) | duplicated(SraStudy, fromLast = TRUE))
  # distinct(SampleName, .keep_all = TRUE)
Qin2014_meta_tidy <- left_join(Qin2014_meta_tidy, Qin2014_clinical_new, by = "SampleName") %>% 
  filter(!is.na(condition))
write_tsv(Qin2014_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Qin2014_meta_tidy.tsv"))


########## Ida2021 ###########
Ida2021_meta  <- read.csv(here('..', 'data', 'metadata', "Ida2021_meta.txt"), na.strings = c("", "missing", "NA", "not_applicable"))
colnames(Ida2021_meta)
Ida2021_meta_tidy <- Ida2021_meta %>% 
  select("Run", "Assay.Type", "AvgSpotLen", "Bases", "BioProject", "BioSample", "Center.Name", "Instrument", "LibraryLayout", 
         "Organism", "Platform", "Source_material_ID", "geo_loc_name_country", "geo_loc_name_country_continent")  %>% 
  rename(condition = Source_material_ID, AssayType = Assay.Type, CenterName = Center.Name, Country = geo_loc_name_country, Continent = geo_loc_name_country_continent) %>% 
  mutate(condition = gsub("\\d+", "", condition)) %>% 
  mutate(condition = gsub("HD", "Control", condition))
write_tsv(Ida2021_meta_tidy, here('..', 'data', 'metadata', 'tidy_data', "Ida2021_meta_tidy.tsv")) 
