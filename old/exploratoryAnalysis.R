## Unsupervised data analysis
## 1. MA plots (or scatterplots) between samples
## 2. Principal component analysis for samples
## 3. Hierarchical clustering of highly variant peptides/proteins
## 4. Entries can be selected by CV (coefficient of variation) or MAD (mean absolute deviation)
##    e.g. all peptides/proteins, 5% of highly variable peptides/proteins, etc.
rm(list = ls())
library(limma)
library(readxl)
library(ggplot2)

level = "protein"
inFile = "id_uni_prot_quan.xlsx"
if (length(grep("pep", inFile))) {
  tbl = read_excel(inFile, skip = 4) # Peptide publication table
} else {
  tbl = read_excel(inFile, skip = 1) # Protein publication table
}
rawData = as.data.frame(tbl)
if (level == "peptide") {
  entry = rawData[, 1]
} else if (level == "protein") {
  entry = rawData[, 2]
}

## CV or MAD calcuation is based on log2-transformed intensities
## but output format is raw-intensity scale
colInd = grep("sig", colnames(rawData))
data = log(rawData[, colInd], 2)
cv = apply(data, 1, sd) / rowMeans(data)
mad = apply(abs(data - apply(data, 1, median)), 1, median)
threshold = 0.4 ## Threshold percentage
rowInd = NULL
rowInd = cv > quantile(cv, prob = 1 - threshold)

## Return data for the following analyses
## Column 1: entry (either protein accession or peptide sequence)
## Column 2~ : log2-transformed intensity values for reporter ions
entry = entry[rowInd]
data = cbind(entry, log(rawData[rowInd, colInd], 2))

# ## Selection of a subset of data
# cv = apply(data, 1, sd) / rowMeans(data)
# mad = apply(abs(data - apply(data, 1, median)), 1, median)
# threshold = 10/ 100 ## Threshold percentage
# rowInd = NULL
# rowInd = cv > quantile(cv, prob = 1 - threshold)
# 
# rawData[rowInd, ]











# ## PCA plot
# graphics.off()
# resPCA = prcomp(t(subData), center = TRUE, scale = TRUE)
# eigs = resPCA$sdev ^ 2
# xlabPCA = paste0("PC1 (", round((eigs[1] / sum(eigs)) * 100, 2),"%)")
# ylabPCA = paste0("PC2 (", round((eigs[2] / sum(eigs)) * 100, 2),"%)")
# resPCA = data.frame(resPCA$x[, 1:2])
# ggplot(data = resPCA[, 1:2], aes(PC1, PC2)) + 
#   geom_jitter(size = 3) + 
#   geom_text(aes(label = rownames(resPCA)), vjust = "inward", hjust = "inward", size = 4) +
#   labs(x = xlabPCA, y = ylabPCA)

# ## Hierarchical clustering
# color <- colorRampPalette(c("blue", "white", "red"))(n = 100)
# mat = as.matrix(subData)
# mat = t(scale(t(mat), center = T, scale = F)) # Only mean-centering
# matBreaks = seq(-3, 3, length.out = 100)
# pheatmap(mat = mat, color = color, breaks = matBreaks, show_rownames = F)





