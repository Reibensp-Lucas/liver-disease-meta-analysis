# Liver disease meta analysis
Scripts and other files that were used during my internship at Zeller team

This is a description of the project and where to find the files and scripts that I worked with during my internship at Zeller group.
wherever possible, I also added the directory path for `ssh username@nile.emble.de server` , where most of the files can be found in.

The project consits of one main part - the liver disease meta analysis - and one minor part - QC impact on microbial profiles in low quality samples . The latter one was a small exploratio, together with Fabian, that came up when profiling the 16S data for the meta anlysis. In short: 

## QC difference testing
Fabian noticed that some of the datasets were of bad quality and where thus subjected to heavy filtering when performing quality control as part of the profiling. This resulted in lost reads. We ware asking whether this affected the outcome of subsequent analyses and investigated this shortly. 

The Script **'QC_difference_test.R'** can be found in `g/scb/zeller/reibensp/LD_meta_analysis/scripts/QC_difference_test.R` and the **qc_diffence_test** commit in this repository. 

Data can be found in `g/scb/zeller/reibensp/LD_meta_analysis/data/16S/Lin_2023/` with the subfolders `noQC_tax_level/` and `QC_tax_level/`. 

## Liver disease meta analysis ###### 
### dataset colleciton and literature search 
dataset collection that resulted from the literature search can be found in [Microbiome-disease-association spreadsheet](https://docs.google.com/spreadsheets/d/1lqbHJrT2GXUAKYTmvlJhKC2kbS035wCi6HPpivI93AU/edit) 
In the _lines 106 - 140_, Green BioProject IDs represent data that was downloaded and profiled. 
Green metatadata columns represent files where metadata that allows the link to the case or control group of the study cohort is available (usually SRA Accession metadata, if not other specified) pubmed search script  and results can be found under commit **literature_search**. 

### raw data after profiling 
ALl data can be found in `g/scb/zeller/reibensp/LD_meta_analysis/data/` There are subfolders here:
`../16S/` contains folders with the names of the datasets, each containing the profiling results in different tax levels, and the **create_count_dfs.R** snippet to manually convert the tax level profiles into the needed count matrix file for the analysis, plus the default 16S.Rproj file.
In the same way `../WGS/` contains the tax level profiles in the dataset_name_folders, as well as the **create_count_matrices.R** script to convert them to the rght format, as well as the WGS.Rproj
`../raw_profiles/` contains these raw_profiles that can be read in from there to the meta anlaysis script (**analysis.R**)


### cleaning metadata 
original metadata can be found in`g/scb/zeller/reibensp/LD_meta_analysis/data/metadata/` as .txt or .tsv files. Data cleaning was done with R - **cleaning_metadta.R** can be accessed via `g/scb/zeller/reibensp/LD_meta_analysis/scripts/cleaning_metadat.R` or found in this repo. 
Resulting tidy versions of the metadata were saved to `g/scb/zeller/reibensp/LD_meta_analysis/data/metadata/tidy_data/` as .tsv files.
