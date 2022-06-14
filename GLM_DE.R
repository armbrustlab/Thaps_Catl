# Zinka 2020_06_16
# Differential expression analysis using EdgeR. 
#
# Inputs: 
# - Table with TPM values for all genes in all samples (all_abundnace_TPM.csv)
#
# Outputs: 
# - PLots: Mean square deviation plot, dispertion plots, differential expression plots. 
# - Data table with logFC and p-values for differential expression of trated vs. control T. pseudonana cultures at each time point. 

library(edgeR)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)

################################################################
# Importing files:
################################################################

#Reading in count file containing all experimental groups with only raw TPM counts, no annotations! 
rawdata<- read.csv(file = '~/Desktop/Thaps_Catl/all_abundnace_TPM.csv',
                   check.names=FALSE, 
                   stringsAsFactors = FALSE, 
                   header = TRUE)
head(rawdata)

################################################################
# Pre-processing data:
################################################################

#Format
counts <- select(rawdata, -"Gene ID") %>% 
  as.matrix()
genes <- select(rawdata, "Gene ID")

#Create metadata table for sample characteristics (time, group, replicate)
samples <- as.array(colnames(counts)) %>% 
  strsplit2( "[_]") 
colnames(samples) <- c("time", "treatment", "rep")

#Put data into DGEList format:
y <- DGEList(counts=counts, samples = samples, genes=genes)

#TMM normalization:
y <- calcNormFactors(y) #Estimates normalization factors (used in Anders 2013)

# Mean Square deviation (MSD) plot:
pdf(file = "~/Desktop/Thaps_Catl/MSD.pdf", width = 10, height = 8)
plotMDS(y)
dev.off()

################################################################
# GLM model for differential expression analysis
################################################################

#Making the design matrix:
samples <- data.frame(samples)
Group <- factor(paste(samples$time, samples$treatment, sep="."))
samples <- cbind(samples,Group=Group)
design <- model.matrix(~0+Group)
colnames(design) <- levels(Group) 

#Estimating dispersion:
disp = estimateDisp(y, design)
pdf(file = "~/Desktop/Thaps_Catl/Estimated_dispersion.pdf", width = 10, height = 8)
plotBCV(disp)
dev.off()

#Fitting GLMs to my data with the design matrix superimposed:
fit <- glmQLFit(disp, design) 
pdf(file ="~/Desktop/Thaps_Catl/Fitted_dispersion.pdf", width = 10, height = 8)
plotQLDisp(fit)
dev.off()

################################################################
# Setting contrasts
################################################################

#Contrasts for pairwise comparisons between treated and control transcripts:
#Un-commemnt the one to be analyzed. 

# Compare t24.filtrate to t24.cont: (-1*t24.cont, +1*t24.filtrate)
contrast = c(0,0,-1,1,0,0)
name = 't24.filtrate_vs_t24.cont'

# # Compare t72.filtrate to t72.cont: (-1*t72.cont, +1*t72.filtrate)
# contrast = c(0,0,0,0,-1,1)
# name = 't72.filtrate_vs_t72.cont'

# # Compare t120.filtrate to t120.cont: (-1*t120.cont, +1*t120.filtrate)
# contrast = c(-1,1,0,0,0,0)
# name = 't120.filtrate_vs_t120.cont'

################################################################
# Differential expression analsyis
################################################################

#DE with Quasi-Likelihood F-tests:
de = glmQLFTest(fit, contrast = contrast)
topTags(de)

#Calculate the total number of differentially expressed genes at 5% FDR
DEnumbers <- summary(decideTests(de))
DEnumbers = data.frame(DEnumbers)
pdf(file = "~/Desktop/Thaps_Catl/Num_DE_genes.pdf", width = 10, height = 8)
ylim <- c(0,1.1*max(DEnumbers$Freq))
xx <- barplot(DEnumbers$Freq, 
              main = DEnumbers$Var2[1], 
              names.arg = DEnumbers$Var1, 
              ylab = "Number of DE genes",
              ylim = ylim)
DEnumbers$Freq <- as.numeric(DEnumbers$Freq)
text(x = xx, y = DEnumbers$Freq, label = DEnumbers$Freq, pos = 3)
dev.off()

#Plot log-fold change against log-counts per million, with DE genes highlighed:
pdf(file = "~/Desktop/Thaps_Catl/plotMD.pdf", width = 10, height = 8)
plotMD(de)
abline(h=c(-1,1), col="blue") #adding a blue line that indicates a 2-fold change. 
dev.off()

#Export csv table of DE genes for the selected contrast:
file = paste("~/Desktop/Thaps_Catl/", name, ".csv", sep = "")
write.csv(top_genes$table, file=file)

