# **Gene Expression Data Visualization and Pathway Enrichment Analysis for Glioblastoma**

## **Project Overview**

This project focuses on analyzing gene expression data from glioblastoma multiforme (GBM) samples to identify differentially expressed genes, visualize their expression patterns, and perform pathway enrichment analysis. Glioblastoma is a highly aggressive and lethal brain tumor, with poor prognosis and limited treatment options. The insights from this analysis aim to contribute to a better understanding of GBM's molecular mechanisms and potential therapeutic targets.

### **Objectives**
1. Identify upregulated and downregulated genes from gene expression data.
2. Visualize expression patterns using heatmaps with clustering.
3. Perform functional pathway enrichment to highlight biological processes and pathways associated with glioblastoma progression.
4. Relate findings to current advancements in glioblastoma research.

---

## **Findings**
### **1. Gene Expression Analysis**
- **Upregulated Genes**: We identified several genes with significant upregulation (log2 fold change > 1, p-value < 0.05), indicating potential overactivation in glioblastoma samples.
- **Downregulated Genes**: Genes with log2 fold change < -1 and p-value < 0.05 were classified as downregulated, suggesting reduced activity in tumor conditions.

### **2. Pathway Enrichment Analysis**
The functional pathway enrichment analysis identified several key pathways that are significantly enriched in glioblastoma samples. These include:

- **Ribosomal Small Subunit Assembly**:  
  A highly enriched pathway (**Fold Enrichment: ~200**), reflecting altered ribosomal biogenesis, a known hallmark of cancer.  
  - **Gene involved**: PWP2.

- **Proteolysis**:  
  Involves degradation of proteins, critical for tumor microenvironment remodeling and invasion.  
  - **Genes involved**: TMPRSS3, TPSAB1, TPSB2.

- **Cellular Sodium Ion Homeostasis**:  
  Highlighted with a high fold enrichment (**~172.7**), potentially linking to ionic transport dysregulation in glioblastoma cells.  
  - **Gene involved**: TMPRSS3.


### **3. Heatmap Visualizations**
Heatmaps revealed distinct expression patterns across samples:
- Clustering of upregulated genes showed coherent patterns linked to specific GBM phenotypes.
- Diverging color palettes effectively visualized both upregulation and downregulation trends.

---

## **Relevance in Current GBM Research**
1. **Ribosomal Pathways**: Elevated ribosome biogenesis has been identified as a driver of GBM growth, making these pathways potential therapeutic targets.
2. **Proteolysis and ECM Remodeling**: Findings reinforce the role of extracellular matrix remodeling in facilitating tumor invasion, consistent with studies on glioblastoma's invasive phenotype.
3. **Ion Transport Pathways**: Dysregulated ion homeostasis contributes to tumor cell survival under hypoxic and acidic conditions, which are characteristic of GBM.

By focusing on these pathways, our analysis underscores key molecular vulnerabilities in GBM that may inform the development of targeted therapies or biomarker discovery.

---

## **Methods**
### **1. Data Preprocessing**
- **Input**: A glioblastoma gene expression dataset containing over 500 differentially expressed genes.
- **Process**: Filtering, grouping by sample conditions, and normalization for downstream analysis.

### **2. Heatmap Generation**
- **Tools**: `heatmap.2()` from the **gplots** R package.
- **Visualizations**: Diverging (green-brown) and sequential (Blues) palettes for gene regulation patterns.
- **Clustering**: Hierarchical clustering of genes and samples to identify co-expression groups.

### **3. Fold Change and Statistical Analysis**
- Log2 fold change calculated to quantify gene expression differences between groups.
- P-values computed using t-tests for statistical significance.

### **4. Pathway Enrichment**
- **Input**: Upregulated and downregulated gene lists (Ensembl IDs).
- **Tools**: ShinyGO (v0.741) for functional enrichment.
- **Outputs**: Top pathways ranked by FDR and fold enrichment.

---

## **Files in Repository**
### **Directories**
- `data/`: Contains input datasets and files related to pathway enrichment.
- `plots/`: Heatmaps and other visualizations generated from the analysis.
- `src/`: Modular R scripts for preprocessing, visualization, and analysis.

### **Key Files**
1. `glioblastoma_gene.R`: Core script containing the entire pipeline.
2. `results/`: CSV files with fold changes and p-values.

### **Future Directions**
1. Validate identified pathways using in vitro or in vivo experiments.
2. Explore additional datasets to confirm reproducibility across different GBM samples.
3. Investigate therapeutic compounds targeting enriched pathways (e.g., ribosomal inhibitors, ion channel blockers).

---

## **How to Use**
### **Dependencies**
- **R version ≥ 4.4.1**
- Required packages: `gplots`, `dplyr`, `ggplot2`, `biomaRt`

### **Steps to Reproduce**
1. Clone the repository:
   ```bash
   git clone https://github.com/Taofeeq-T/Gene-Expression-Data-Visualization-and-Pathway-Enrichment-Analysis-for-Glioblastoma.git
