# liver-disease-deta-analysis-
Scripts and other files that were used during my internship at Zeller team

This is a description of the project and where to find the files and scripts that I worked with during my internship at Zeller group.
wherever possible, I also added the directory path for @nile.emble.de server, where most of the files can be found in.

The project consits of one main part - the liver disease meta analysis - and one minor part - QC impact on microbial profiles in low quality samples . The latter one was a small exploratio, together with Fabian, that came up when profiling the 16S data for the meta anlysis. In short: 

Fabian noticed that some of the datasets were of bad quality and where thus subjected to heavy filtering when performing quality control as part of the profiling. This resulted in lost reads. We ware asking whether this affected the outcome of subsequent analyses and investigated this shortly. The Script  
'QC_difference_test.R' can be found in 'g/scb/zeller/reibensp/LD_meta_analysis/scripts/QC_difference_test.R' and the scripts folder in this repo. 
Files can be found in 'g/scb/zeller/reibensp/LD_meta_analysis/data/16S/Lin_2023/' with the subfolders 'noQC_tax_level' and 'QC_tax_level'. and the respective folders in this repo.

#### Liver disease meta analysis ###### 

dataset collection that resulted from the literature search can be found in 'https://docs.google.com/spreadsheets/d/1lqbHJrT2GXUAKYTmvlJhKC2kbS035wCi6HPpivI93AU/edit' , lines 
