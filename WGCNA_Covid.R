library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(gridExtra)
library(corrplot)
allowWGCNAThreads()         
data <- read.delim("C:/Users/minot/Downloads/GSE152418_p20047_Study1_RawCounts/GSE152418_p20047_Study1_RawCounts.txt", header = T)
#get metadata
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,46:50)]

# prepare data ------------------------------------------------------------
data[1:10,1:10]
data <- data %>% 
  gather(key = 'samples',value = 'counts', -ENSEMBLID) %>% 
  mutate(samples = gsub('\\.', '-', samples)) %>% 
  inner_join(., phenoData, by = c('samples' = 'title')) %>% 
  select(1,3,4) %>% 
  spread(key = 'geo_accession', value = 'counts') %>% 
  column_to_rownames(var = 'ENSEMBLID')

# Quality control --------------------------------------------------------
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are outliers ------------------------------------------
data <- data[gsg$goodGenes == TRUE,]
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pcp - method 2  ---------------------------------------------------------
pca <- prcomp(t(data))
pca.dat <- pca$x
# variance calculation ----------------------------------------------------
pca.var <- pca$sdev^2
pca.dat <- as.data.frame(pca.dat)
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
while (!is.null(dev.list())) dev.off() 
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# Exclude outlier data ----------------------------------------------------

samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4614995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# normalize the data ------------------------------------------------------
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))
names(colData)

# row names identical to column names -------------------------------------
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# create dds dataset ------------------------------------------------------

sum(is.na(data.subset))
data.subset[is.na(data.subset)] <- 0
str(data.subset)
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1)

# Variance stabilization  -------------------------------------------------
install.packages("tidyverse")
library(tidyverse)

dds_norm <- vst(dds75)
norm.counts <- assay(dds_norm) %>% 
  t()

# Network construction  ---------------------------------------------------
library(WGCNA)
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
sft.data <- sft$fitIndices

# install package ggplot2  ----------------------------------------------
install.packages("ggplot2")

# install package gridExtra -----------------------------------------------
install.packages("gridExtra")
# visualize to select power  ----------------------------------------------
library(ggplot2)
library(gridExtra)
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()
a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric  ----------------------------------------------
library(WGCNA)
soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

# model eigen genes -------------------------------------------------------
cor <- temp_cor
module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)

# get number of genes in each module  -------------------------------------
table(bwnet$colors)

# plot dendrogram and the module colors before and after merging u --------
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# install dplyr -----------------------------------------------------------

install.packages("dplyr")

# create traits file - binarize categorical variables  --------------------
library(dplyr)
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>% 
  select(8)

# binarize categorical columns  -------------------------------------------
library(WGCNA)
colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))
severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)
traits <- cbind(traits, severity.out)

# Define number of samples and genes  -------------------------------------

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# install tibble package  -------------------------------------------------
install.packages("tibble")

# visualize module-trait association as a heatmap -------------------------
library(tibble)
install.packages("corrplot")
library(corrplot)
heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

cor_matrix <- cor(heatmap.data, use = "pairwise.complete.obs")

corrplot(cor_matrix, 
         method = "color", 
         col = c("blue1", "skyblue", "white", "pink", "red"), 
         tl.col = "black", 
         tl.srt = 45)

# dataframe for segregating module turquoise ---------------------------------------------------------------
install.packages("dplyr")   # if not already installed
library(dplyr)
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(bwnet.colors == 'turquoise') %>% 
  rownames()
view(module.gene.mapping)

# module membership / intramolecular connectivity  ------------------------

module.membership.measures <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
module.membership.measure.pvals[1:10,1:10]

# calculating the gene significance and associated p-values  --------------

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)

# new ---------------------------------------------------------------------
METurquoise <- module_eigengenes[, "MEturquoise"]
MEGreen <- module_eigengenes[, "MEgreen"]
# new1 --------------------------------------------------------------------
write.csv(METurquoise, "turquoise_genes.csv", row.names = FALSE)
write.csv(MEGreen, "green_genes.csv", row.names = FALSE)

# new2 --------------------------------------------------------------------

write.csv(traits, "traits_data.csv", row.names = TRUE)

# module eigengene and metadata extarction  -------------------------------
library(GEOquery)
library(WGCNA)
library(dplyr)

# Define GEO dataset ID
geo_id <- "GSE152418"

# Download and extract metadata
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))

# Select relevant columns
phenoData <- phenoData[, c(1,2,46:50)]

# Save metadata to CSV
write.csv(phenoData, "GSE152418_metadata.csv", row.names = FALSE)

# Extract module eigengenes
module_eigengenes <- bwnet$MEs

# Ensure target labels are correctly assigned
target_labels <- phenoData$disease_state  # Adjust if column name differs
if (is.null(target_labels)) {
  stop("Error: Target labels column not found. Check column name in phenoData.")
}

# Add target labels to module eigengenes
module_eigengenes <- cbind(module_eigengenes, target_labels)

# Extract turquoise and green module eigengenes
turquoise_genes <- module_eigengenes %>% select(starts_with("MEturquoise"), target_labels)
green_genes <- module_eigengenes %>% select(starts_with("MEgreen"), target_labels)

# Save filtered data to CSV
write.csv(turquoise_genes, "turquoise_genes_with_targets.csv", row.names = TRUE)
write.csv(green_genes, "green_genes_with_targets.csv", row.names = TRUE)

# Print first few rows to verify
head(turquoise_genes)
head(green_genes)

# check -------------------------------------------------------------------
library(GEOquery)
library(WGCNA)

# Define GEO dataset ID
geo_id <- "GSE152418"

# Download and extract metadata
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(gse[[1]])

# Select relevant columns
phenoData <- phenoData[, c("geo_accession", "disease state:ch1")]

# Extract module eigengenes
module_eigengenes <- bwnet$MEs

# Ensure row names match between module_eigengenes and phenoData
common_samples <- intersect(rownames(module_eigengenes), phenoData$geo_accession)
module_eigengenes <- module_eigengenes[common_samples, ]
target_labels <- phenoData[match(common_samples, phenoData$geo_accession), "disease state:ch1"]

# Extract only green and turquoise genes
green_genes <- module_eigengenes$MEgreen
turquoise_genes <- module_eigengenes$MEturquoise

# Create data frames with aligned target values
green_data <- data.frame(GreenGenes = green_genes, Target = target_labels)
turquoise_data <- data.frame(TurquoiseGenes = turquoise_genes, Target = target_labels)

# Save as CSV files
write.csv(green_data, "green_genes_with_targets.csv", row.names = FALSE)
write.csv(turquoise_data, "turquoise_genes_with_targets.csv", row.names = FALSE)

# Print success message
cat("Feature extraction and saving completed successfully!")

# noq ---------------------------------------------------------------------
mapping_data <- read.table("module.gene.mapping", header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Extract IDs corresponding to turquoise color
turquoise_ids <- mapping_data$ENSG[mapping_data$bwnet$colors == "turquoise"]
write.table(turquoise_ids, file="turquoise_ids.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Extract IDs corresponding to green color
green_ids <- mapping_data$ENSG[mapping_data$bwnet$colors == "green"]
write.table(green_ids, file="green_ids.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Print the number of genes in each module
cat("Number of genes in turquoise module:", length(turquoise_ids), "\n")
cat("Number of genes in green module:", length(green_ids), "\n")