# Gene Expression Analysis Script
# Author: Taofeeq Togunwa
# Date: November 21, 2024
# Purpose: Analyze glioblastoma gene expression data, visualize results, and perform pathway enrichment.

# ==== 1. Load Required Libraries ====
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt")

install.packages("gplots", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)

# Load libraries
library(biomaRt)
library(gplots)
library(dplyr)
library(ggplot2)

# ==== 2. Load the Gene Expression Data ====
# Path to the dataset
data_path <- "Data/glioblastoma.csv" # Update the path if needed
data <- read.csv(data_path, header = TRUE, row.names = 1)

# Examine data structure
cat("Dataset Summary:\n")
print(summary(data))
cat("Number of rows (genes):", nrow(data), "\n")
cat("Number of columns (samples):", ncol(data), "\n")

# ==== 3. Generate Heatmaps ====
# Function to generate and save heatmaps
generate_heatmap <- function(data_matrix, output_path, color_palette = "Blues3", scale = "row", dendrogram = "both") {
  png(output_path, width = 1000, height = 800)
  heatmap.2(
    as.matrix(data_matrix),
    trace = "none",
    scale = scale,
    dendrogram = dendrogram,
    Colv = TRUE,
    Rowv = TRUE,
    col = hcl.colors(100, palette = color_palette),
    main = paste("Heatmap with", color_palette, "palette")
  )
  dev.off()
  cat("Heatmap saved to:", output_path, "\n")
}

# Heatmap without scaling
generate_heatmap(data, "plots/heatmap_no_scaling.png", color_palette = "Blues3", scale = "none")

# Heatmap with scaling and diverging palette
generate_heatmap(data, "plots/heatmap_diverging.png", color_palette = "green-brown", scale = "row")

# ==== 4. Fold Change and P-Value Calculation ====
# Define groups
group1 <- c(1:5) # First five columns
group2 <- c(6:10) # Last five columns

# Extract group data
group1_data <- data[, group1]
group2_data <- data[, group2]

# Calculate fold change and log2 fold change
group1_mean <- rowMeans(group1_data)
group2_mean <- rowMeans(group2_data)
fold_change <- (group2_mean - group1_mean) / group1_mean
logFC <- log2(fold_change)

# Calculate p-values
pvalues <- apply(data, 1, function(row) {
  t.test(row[group1], row[group2])$p.value
})

# Combine results into a single data frame
results <- data.frame(Gene = rownames(data), logFC = logFC, pvalues = pvalues)
write.csv(results, "results/fold_change_pvalues.csv")
cat("Fold change and p-value results saved to: results/fold_change_pvalues.csv\n")

# ==== 5. Identify Upregulated and Downregulated Genes ====
# Define thresholds
logFC_up_cutoff <- 1
logFC_down_cutoff <- -1
pvalue_cutoff <- 0.05

# Filter significant genes
upregulated_genes <- results %>% filter(logFC > logFC_up_cutoff & pvalues < pvalue_cutoff)
downregulated_genes <- results %>% filter(logFC < logFC_down_cutoff & pvalues < pvalue_cutoff)

# Save results
write.csv(upregulated_genes, "results/upregulated_genes.csv")
write.csv(downregulated_genes, "results/downregulated_genes.csv")
cat("Upregulated genes saved to: results/upregulated_genes.csv\n")
cat("Downregulated genes saved to: results/downregulated_genes.csv\n")

# ==== 6. Pathway Enrichment Analysis ====
# Example: Retrieve gene symbols using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(upregulated_genes) # Example Ensembl IDs
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Save gene conversion results
write.csv(gene_conversion, "results/gene_conversion.csv")
cat("Gene conversion results saved to: results/gene_conversion.csv\n")

# ==== 7. Pathway Visualization ====
# Read pathway enrichment data
enrichment_data <- read.csv("pathway_R/enrichment.csv")

pathway = c("Ribosomal small subunit assembly", 
            "Maturation of SSU-rRNA from tricistronic rRNA", 
            "Proteolysis", 
            "Cellular sodium ion homeostasis", 
            "Maturation of SSU-rRNA", 
            "Ribosome assembly", 
            "Regulation of defense response to virus by virus", 
            "Sodium ion homeostasis", 
            "Ribosomal small subunit biogenesis", 
            "Extracellular matrix disassembly")

nGenes = c(1, 1, 3, 1, 1, 1, 1, 1, 1, 1)
FDR = c(1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.4e-01, 1.7e-01, 1.7e-01)
fold_enrichment = c(200, 88.4, 5.6, 172.7, 66.7, 59.4, 126.6, 65.5, 45.2, 40.4)



# Create a lollipop plot for pathways
png("plots/pathway_lollipop.png", width = 1200, height = 800)
ggplot(enrichment_data, aes(y = reorder(pathway, nGenes), x = nGenes)) +
  geom_segment(aes(yend = pathway, xend = 0), color = "gray") +
  geom_point(aes(size = -log10(FDR), color = fold_enrichment), alpha = 0.7) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    x = "Number of Genes", y = "Pathway",
    title = "Pathway Enrichment Analysis",
    size = "-log10(FDR)", color = "Fold Enrichment"
  ) +
  theme_minimal()
dev.off()
cat("Pathway lollipop plot saved to: plots/pathway_lollipop.png\n")

# ==== End of Script ====
