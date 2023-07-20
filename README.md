# Liver disease meta analysis
Scripts and other files that were used during my internship at Zeller team

This is a description of the project and where to find the files and scripts that I worked with during my internship at Zeller group.
wherever possible, I also added the directory path for `ssh username@nile.emble.de server` , where most of the files can be found in.

The project consits of one main part - the liver disease meta analysis - and one minor part - QC impact on microbial profiles in low quality samples . The latter one was a small exploratio, together with Fabian, that came up when profiling the 16S data for the meta anlysis. In short: 

## QC difference testing
Fabian noticed that some of the datasets were of bad quality and where thus subjected to heavy filtering when performing quality control as part of the profiling. This resulted in lost reads. We ware asking whether this affected the outcome of subsequent analyses and investigated this shortly. 

The Script **'QC_difference_test.R'** can be found in `g/scb/zeller/reibensp/LD_meta_analysis/scripts/QC_difference_test.R` and the **scripts folder** in this repository. 

Files can be found in `g/scb/zeller/reibensp/LD_meta_analysis/data/16S/Lin_2023/` with the subfolders `noQC_tax_level` and `QC_tax_level`. and the respective **folders** in this repo.

## Liver disease meta analysis ###### 
### dataset colleciton and literature search 
dataset collection that resulted from the literature search can be found in [Microbiome-disease-association spreadsheet](https://docs.google.com/spreadsheets/d/1lqbHJrT2GXUAKYTmvlJhKC2kbS035wCi6HPpivI93AU/edit) 
In the _lines 106 - 140_, Green BioProject IDs represent data that was downloaded and profiled. 
Green metatadata columns represent files where metadata that allows the link to the case or control group of the study cohort is available (usually SRA Accession metadata, if not other specified) pubmed search script  and results can be found in folder **literature_search**. 
