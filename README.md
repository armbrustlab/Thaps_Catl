# Thaps_Catl

GitHub page contiaing code for transcrtiptome processing described in paper 

"Flavobacterial exudates disrupt cell cycle progression and metabolism of the diatom Thalassiosira pseudonana" 

# Code files:

- WGCNA.R: Script that runs WGCNA clustering.  
- GLM_DE.R: Script that does differential expression analysis using EdgeR. 
- enrichment.R: Script that performs functional enrichemnt on genes in WGCNA modules and differentially expressed genes for the three time comparisons.

R version used: R version 4.1.0 (2021-05-18)

# Instructions for running scripts:

1. Download all data .csv files from Zenodo (10.5281/zenodo.6672614). There are 5 data files needed for the code to run: all_abundance_TPM.csv, DE.csv, Experimental_parameters.csv, thaps_GO.txt, WGCNA_network.csv
3. Create a folder on your desktop called Thaps_Catl, and save all files and code into that folder (~/Desktop/Thaps_Catl). 
4. Open code in R and run scripts. 
