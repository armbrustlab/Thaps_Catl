# Zinka 2020_06_16
# Script for analyzing the expression data using network correlation analysis with WGCNA 
#
#Inputs:
# - Table with TPM values for all genes in all samples (all_abundnace_TPM.csv)
# - Table with measured experimental paratmeters (Experimental_parameters.csv)
#
#Outputs:
# - Plots of module contruction and relation to expermental traits, plots of module and gene selection. 
# - Data table with module membership and gene significnace for filtrate treatment for all analyzed smaples and genes. 

library(WGCNA)
options(stringsAsFactors = FALSE)
library(dplyr)
library(tidyverse)

################################################################
# Importing files:
################################################################

#Importing the raw TPM expression table:
ExprTable_raw = read.csv('~/Desktop/Thaps_Catl/all_abundnace_TPM.csv')

#Importing the raw experimental traits table:
ExpTraits = read.csv('~/Desktop/Thaps_Catl/Experimental_parameters.csv')

################################################################
# Pre-processing the data:
################################################################
ExprTable_raw <- ExprTable_raw %>% remove_rownames %>% column_to_rownames(var = "Gene.ID") %>% as.data.frame()

##### Log2 transforming the raw data ##### 
ExprTable_log = log2(ExprTable_raw + 1)
ExprTable_log <- as.matrix(ExprTable_log)

###### Filter out low level genes ###### 
# Calculate median TPM across all genes and samples, keep genes for which greater than 25% of the samples have expression above the cutoff. 
med = median(ExprTable_log)
threshold = 0.25
ExprTable_log_filt = ExprTable_log[which(sapply(1:nrow(ExprTable_log), function(x) { sum(ExprTable_log[x,]>=med)/ncol(ExprTable_log) } ) >= threshold),]

#Boxplots for before and after removing low expression genes:  
par(mfrow=c(1,1), mar = c(16,8,3,3) +0.1)
boxplot(ExprTable_log, main="All samples Log2 TPM data including low-level genes", cex.axis = 1, las = 3, title.cex = 1, ylab = "log2 (TPM)", cex.lab = 1)
boxplot(ExprTable_log_filt, main="All samples Log2 TPM data after deletion low-level genes", cex.axis = 1, las = 3, title.cex = 1, ylab = "log2 (TPM)", cex.lab = 1)

#Renaming final expression table to ExprTable.
ExprTable = as.data.frame(t(ExprTable_log_filt))

################################################################
# Network construction and module identification:
################################################################

########## Check data for missing values ###########
gsg = goodSamplesGenes(ExprTable, verbose = 3);
gsg$allOK #if this returns TRUE, all genes have passed the cut. 

########### Sample dendrogram to check for outliers ###########
#Cluster the samples to check for outliers:
sampleTree = hclust(dist(ExprTable), method = "average");
sizeGrWindow(12,9)
par(cex = 1.5);
par(mar = c(2,6,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1,
     cex.axis = 1, cex.main = 1, cex = 1)

########### Removing outlier samples ###########
ExprTable <- as.data.frame(t(ExprTable))
ExprTable$t24_cont_3 <- NULL
ExprTable$t72_cont_3 <- NULL
ExprTable = as.data.frame(t(ExprTable))

#Re-plot sample dendrogram:
sampleTree = hclust(dist(ExprTable), method = "average");
sizeGrWindow(12,9)
par(cex = 1.5);
par(mar = c(2,6,2,0))
plot(sampleTree, main = "Sample clustering after deleting outliers", sub="", xlab="", cex.lab = 1,
     cex.axis = 1, cex.main = 1, cex = 1)

########### Set adjacency and topological overlap matrix ###########
softPower = 18;
adjacency = adjacency(ExprTable, power = softPower, type = "signed");
#Toplogical overlap matrix:
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Calculate scale free topology and soft conncectivity:
k = softConnectivity(datE = ExprTable, power = 18)
sizeGrWindow(10,5)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k, main = "Scale free topology")

########### Cluster to get initial modules ###########
#Clustering using TOM:
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

#Setting minimum module size:
minModuleSize = 30;

#Optimizing DeepSplit value
mColorh=NULL
for (ds in 0:4){
  dynamicMods_optimize = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                       deepSplit = ds, pamRespectsDendro = FALSE,
                                       minClusterSize = minModuleSize, pamStage = TRUE);
  mColorh=cbind(mColorh,(dynamicMods_optimize));
}
sizeGrWindow(10,10)
par(mfrow=c(1,1))
plotDendroAndColors(geneTree, mColorh, paste("dpSplt =",0:4), main = "pamStage = TRUE, pamRespectsDendro = FALSE", dendroLabels=FALSE);
#DeepSplit=1 was used. 

#Module identification using dynamic tree cut:
deepSplit = 1
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = deepSplit, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = dynamicMods
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, paste("Dynamic, ds =", deepSplit, sep = " "),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


########### Merge modules with similar co-expression ###########
#Calculate eigengenes
MEList = moduleEigengenes(ExprTable, colors = dynamicColors)
MEs = MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
#Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
#Plot
sizeGrWindow(7, 6)
par(cex = 1);
par(mar = c(2,6,2,0))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red", cex = 2)

#Merge modules
merge = mergeCloseModules(ExprTable, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;

#Calculate and plot new dendogram of eigengenes:
mergedMEDiss = 1-cor(mergedMEs)
mergedMETree = hclust(as.dist(mergedMEDiss), method = "average")
sizeGrWindow(7, 6)
par(mar = c(0,5,2,1))
par(cex = 1);
plot(mergedMETree, main = "Clustering of merged module eigengenes",
     xlab = "", sub = "", cex = 1, cex.main = 1, cex.lab = 1, cex.axis = 1)

#Plot new dendrogram:
sizeGrWindow(50, 50)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.lab =2, cex.axis = 2, cex.colorLabels = 1.5, cex.main = 2)

# Rename module colors and labels
moduleColors = mergedColors
moduleLabels = moduleColors
MEs = mergedMEs;

################################################################
# Relating modules to other experimental parameters:
################################################################

########### Manipulating experimental data ###########
ExpTraits <- ExpTraits %>% remove_rownames %>% column_to_rownames(var = "Sample") %>% as.data.frame()
ExpTraits <- as.data.frame(t(ExpTraits))
# Remove any samples that were removed from the sample dendrogram analysis:
ExpTraits$t24_cont_3 <- NULL
ExpTraits$t72_cont_3 <- NULL
# Re-formating the table:
ExpTraits <- as.data.frame(t(ExpTraits))

########### Visualize how samples correlate with experimental data ###########
# Re-cluster samples
sampleTree2 = hclust(dist(ExprTable), method = "average")
# Convert experimental measurnments to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(ExpTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(7, 6)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(ExpTraits),
                    main = "Sample dendrogram and trait heatmap",
                    cex.dendroLabels = 1.8, cex.rowText = 1, cex.main = 1, cex.lab = 1, cex.axis = 1)

########### Qunatifying module-trait relationships ############ 
nGenes = ncol(ExprTable);
nSamples = nrow(ExprTable);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(ExprTable, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, ExpTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(12,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 8.5, 3, 6));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(ExpTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               cex.axis = 1, cex.lab = 1, cex.main = 1)


########### Gene Significance and Module Membership calculations ############ 
#Select the trait I want to examine
trait = 'Filtrate treatment'
weight = as.data.frame(ExpTraits$SA60.treatment);
names(weight) = "Treatment"

#Calculate module membership and gene significance for the trait. 
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(ExprTable, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(ExprTable, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

#Calculate gene significance for single trait. 
trait= ExpTraits$SA60.treatment
GS1 = as.numeric(cor (ExprTable, trait, use = "p"))
GeneSignificance = abs(GS1) 
#Calcualte module significance (it is the average gene significance):
ModuleSignificance = tapply(GeneSignificance, moduleColors, mean, na.rm = T)
sizeGrWindow(8,7)
par(mfrow = c(1,1), mar = c(5,5,5,5))
plotModuleSignificance(GeneSignificance,moduleColors, cex.lab=1, cex.main=1, cex.axis=1, cex =1,
                       ylab = "Module significance",
                       main = "Module significance for filtrate treatment")

########### Identifying genes with high Gene Significance and Module Membership ############ 
#Select module and trait to look at:
module = "9"
trait = 'Filtrate treatment'
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1), mar = c(5,6,5,15));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in","module", module),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = "black")
abline(h = 0.8, col = "red", lwd = 2)
abline(v = 0.8, col = "red", lwd = 2)

################################################################
# Comparing intramodular connectivity to module membership and gene significance 
################################################################

# Calculating intramodular connectivity for each gene:
IntraMod = intramodularConnectivity(adjacency, moduleColors)

# Plotting relationshjip between gene significance and intramodular connectivity:
colorlevels=unique(moduleColors) 
sizeGrWindow(30,30)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,4.5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(IntraMod$kWithin[restrict1],
                     GeneSignificance[restrict1], col="grey",
                     main=paste("Mod",whichmodule, sep = " "),
                     xlab = "kWithin Connectivity", ylab = "Gene Significance", abline = TRUE)
}

#Plotting relationsip between module membership and intramodular connectivity.
sizeGrWindow(8,6)
par(mfrow=c(1,3))
which.color="9";
restrictGenes=moduleColors==which.color
verboseScatterplot(IntraMod$kWithin[ restrictGenes],
                   (geneModuleMembership[restrictGenes, paste("MM", which.color, sep="")])^6,
                   col="black",
                   main = paste("Mod", which.color, sep = " "),
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color="5";
restrictGenes=moduleColors==which.color
verboseScatterplot(IntraMod$kWithin[ restrictGenes],
                   (geneModuleMembership[restrictGenes, paste("MM", which.color, sep="")])^6,
                   col="black",
                   main = paste("Mod", which.color, sep = " "),
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color="21";
restrictGenes=moduleColors==which.color
verboseScatterplot(IntraMod$kWithin[ restrictGenes],
                   (geneModuleMembership[restrictGenes, paste("MM", which.color, sep="")])^6,
                   col="black",
                   main = paste("Mod", which.color, sep = " "),
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")


################################################################
# Saving results
################################################################

# Create final sumary table to export:
genes = names(ExprTable) 
GS_All = corAndPvalue(ExprTable, ExpTraits, use = "p"); #Gene significance table
MM_All = corAndPvalue(ExprTable, MEs, use = "p"); #Module membership table

# Form matrices holding the GS and MM (called GSmat and MMmat):
GSmat = rbind(GS_All$cor, GS_All$p)
nTraits = ncol(ExpTraits)
traitNames = colnames(ExpTraits)
dim(GSmat) = c(nGenes, 2*nTraits)
probes = genes
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
  c("GS.", "p.GS."),
  rep(traitNames, rep(2, nTraits)))

MMmat = rbind(MM_All$cor, MM_All$p)
nME = ncol(MEs)
MENames = colnames(MEs)
dim(MMmat) = c(nGenes, 2*nME)
probes = genes
rownames(MMmat) = probes;
colnames(MMmat) = spaste(
  c("MM.", "p.MM."),
  rep(MENames, rep(2, nME)))

#Save into a csv file:
summary_info = data.frame(Gene_ID = genes,
                  ModuleLabel = moduleColors,
                  ModuleColor = labels2colors(moduleColors),
                  GSmat,
                  MMmat);

# Order the genes first by module color, then by geneTraitSignificance for the Treatment trait. 
GeneOrder = order(summary_info$ModuleLabel, -abs(summary_info$GS.SA60.treatment))
summary_info = summary_info[GeneOrder,]
head(summary_info)

#Output table:
write.csv(summary_info, file = "~/Desktop/Thaps_Catl/Network_summary_info.csv", row.names = FALSE, quote = FALSE)
