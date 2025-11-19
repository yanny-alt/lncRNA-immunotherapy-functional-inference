# Technical Methodology: Guilt-by-Association Analysis

## Analysis Pipeline

### 1. Input Data Preparation
- **Source:** Correlation matrix from team members (lncRNAs vs protein-coding genes)
- **Scope:** 3,645 lncRNAs analyzed against complete protein-coding transcriptome
- **Filtering:** Focused on positive Pearson correlation coefficients (r > 0)

### 2. Protein Partner Identification
- For each lncRNA, extracted top 50 most strongly co-expressed protein-coding genes
- Selection based on highest Pearson correlation values
- Generated 50 gene sets for functional enrichment analysis

### 3. Functional Enrichment
- **Tool:** clusterProfiler (v4.10.0) in R/Bioconductor
- **Databases:** 
  - Gene Ontology Biological Process (GO BP)
  - KEGG Pathways
- **Statistical Thresholds:**
  - FDR < 0.05
  - q-value < 0.1
  - Minimum gene set size: 5 genes

### 4. Results Processing
- Filtered for lncRNAs with significant enrichment (FDR < 0.05)
- Calculated functional diversity scores: SUM(GO_BP_terms + KEGG_pathways)
- Prioritized top 15 candidates based on diversity scores
- Generated comprehensive annotation summaries

## Technical Specifications

### Software & Packages
- **R version:** 4.3.1
- **Key packages:** clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2, tidyverse
- **Visualization:** ggplot2, pheatmap, RColorBrewer

### Statistical Methods
- Pearson correlation analysis
- Fisher's exact test for enrichment (clusterProfiler)
- Benjamini-Hochberg FDR correction
- Multiple testing correction across all lncRNAs

### Output Metrics
- Functional annotation completeness: 50/3,645 lncRNAs (1.4%)
- Average GO BP terms per lncRNA: 475.1
- Average KEGG pathways per lncRNA: 39.4
- Immune regulation bias: 39/50 lncRNAs (78%)
