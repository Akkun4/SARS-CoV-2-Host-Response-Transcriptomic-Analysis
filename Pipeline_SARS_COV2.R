# Load data into a matrix
data <- matrix(
  c(531.4478663, 894.9485, 11095.82, 937.9433, 31307.15, 630.6218, 998.3089113, 1735.505, 4294.123, 76.71832, 236.6921, 6224.914554, 1297.278146, 1099.780765, 1863.263, 9152.999855, 38.3318690),
  nrow = 17,
  ncol = 6,
  dimnames = list(
    c("CSF3", "ROR1-AS1", "CXCL3", "CXCL5", "CXCL1", "CXCL16", "ZCCHC2", "NUAK2", "MT1E", "TNFSF15", "TUBB3", "NA", "TRAF1", "BHLHE41", "NA", "NA", "PRKCG"),
    c("baseMean", "log2FoldChange", "IfcSE", "stat", "pvalue", "padj")
  )
)


heatmap(transposed_data, Rowv = NA, Colv = NA, scale = "column", margins = c(5, 10), main = "Genes vs Fold Change Heatmap")


# Install required packages if not already installed
if (!require("ComplexHeatmap")) install.packages("ComplexHeatmap")
if (!require("heatmaply")) install.packages("heatmaply")

# Generate heatmap using ComplexHeatmap
library(ComplexHeatmap)
library(heatmaply)

library(clusterProfiler)

# Define color palette
my_palette <- colorRampPalette(c("green", "yellow", "red"))

# Create heatmap
row_dendrogram <- FALSE
col_dendrogram <- FALSE
heatmap(
  data,
  cluster = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  col = my_palette,
  row_dendrogram = row_dendrogram,
  col_dendrogram = col_dendrogram,
  main = "Heatmap of data",
  cell_border = "gray"
)

# Add title and axis labels
title("Heatmap of example data")
rownames(data) <- rownames(data)
axis.title(x = "Genes", y = "log2FoldChange")

# Display heatmap using heatmaply (optional)
# heatmaply(data, main = "Heatmap with hover tooltip")

# Load necessary libraries
install.packages("clusterProfiler")
library(GOplot)  # Replace with the appropriate organism database

# Example gene list
gene_list = read.csv("results_retinopathy.csv")
Gene = gene_list$Gene.name
FoldChange = gene_list$log2FoldChange

# Create a matrix with genes and fold changes
heatmap_data <- matrix(gene_data$FoldChange, ncol = 1)

# Set row names to gene names
rownames(heatmap_data) <- gene_data$Gene

# Generate heatmap
heatmap(heatmap_data, Rowv = NA, Colv = NA, scale = "column", margins = c(5, 10), main = "Genes vs Fold Change Heatmap")