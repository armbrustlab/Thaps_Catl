# Zinka 2020_06_16
# Functional enrichemnt analysis 
#
# Inputs: 
# - T. pseudonana GO annotations file (thaps_GO.txt)
# - Table of WGCNA network summary (WGCNA_network.csv)
# - Table of differentially expressed genes at each time point (DE.csv)
#
# Outputs: 
# - Summary of enriched terms for the module or DE time comparison being analyzed 
#   (for topGO enrichment: GO_enrich_summary, for KEGG enrichment: KEGG_enrich_summary)

library(topGO)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)

################################################################
# Importing files:
################################################################

##### Importing T. pseudonana GO annotations from a txt file:
thaps_GO <- readMappings(file = '~/Desktop/Thaps_Catl/thaps_GO.txt')
geneNames <-  names(thaps_GO) #Setting the gene universe genes

##### Importing table of WGCNA network:
WGCNA <- read.csv('~/Desktop/Thaps_Catl/WGCNA_network.csv')
WGCNA = as.data.frame(WGCNA)

##### Importing table of differentially expresed genes:
DE <- read.csv('~/Desktop/Thaps_Catl/DE.csv')

################################################################
# Selecting the module or time comparison for enrichment analysis 
################################################################
#Specify the WGCNA module or DE time comparison to do enrichemnt analysis on. 
#Un-comment the DE section to do enrichemnt on differentally expressed genes. 

##### For enrichment of WGCNA modules:
#Set module to look at (all possible modules: 1,5,6,7,9,11,17,21,23,33,34,35)
module = "9" 
#Define list of genes that belongs to that module:
MyInterestingGenes = WGCNA$Gene_ID [ WGCNA$ModuleLabel == module]

##### For enrichment of DE genes at each time point:
#Set the timepoint to look at (all possible timepoints: t24, t72, t120)
module = "t24"
#Define list of genes that belongs to DE timepoint:
if (module == "t24") {
  MyInterestingGenes <- DE$t24_geneID
}
if (module == "t72") {
  MyInterestingGenes <- DE$t72_geneID
}
if (module == "t120") {
  MyInterestingGenes <- DE$t120_geneID
}

################################################################
# Functional enrichemnt with topGO
################################################################

#Creating a geneList object that is a named factor indicatig which genes are interesting and which are not. 
geneList <- factor(as.integer(geneNames %in% MyInterestingGenes))
names(geneList) <- geneNames
str(geneList)

#Building the topGOdata object: 
ontology = "BP" #Select what type of GO enrichemnt I want to do (BP, MF, CC)
GOdata <- new("topGOdata", 
              ontology = ontology, 
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = thaps_GO,
              nodeSize = 5)
GOdata

#Classic Fisher's exact test:
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)
resultFisher

#Summary table:
GO_enrich_summary <- GenTable(GOdata,
                   classicFisher = resultFisher,
                   #weightFisher = resultWeight,
                   orderBy = "classicFisher", 
                   ranksOf = "classicFisher", 
                   topNodes = length(resultFisher@score),
                   numChar = 1000)

################################################################
# Functional enrichemnt with KEGG
################################################################

KEGG_enrich <- enrichKEGG(gene = MyInterestingGenes, 
                              organism = 'tps',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.1)

cnetplot(KEGG_enrich, categorySize = "pvalue") #Visualize enriched term and associated genes. 

# Create summary:
KEGG_enrich_summary = subset(KEGG_enrich@result , KEGG_enrich@result$qvalue <= 0.1) 
