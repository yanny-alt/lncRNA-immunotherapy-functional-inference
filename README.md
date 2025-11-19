# Guilt-by-Association Functional Inference of lncRNA Biomarkers in Hepatocellular Carcinoma

**Author:** Favour O. Igwezeke  
**Affiliation:** Faculty of Pharmaceutical Sciences, University of Nigeria, Nsukka  
**Context:** Capstone Project | Bioinformatics Research Internship at HackBio  
**Internship Period:** September 2025 – November 2025 (Remote)  
**Supervisor:** Dr. Adewale Ogunleye, PhD  
**Manuscript Status:** In Preparation 

---

## Project Overview

This repository documents my **primary contribution** to a collaborative research project focused on discovering long non-coding RNA (lncRNA) biomarkers for immunotherapy response in Hepatocellular Carcinoma (HCC). As part of my Bioinformatics Capstone Project with **HackBio**, I was specifically responsible for the **functional inference component**—moving beyond simple correlation to provide mechanistic insights into lncRNA biology using the guilt-by-association method.

> **My Role:** While team members handled initial RNA-seq processing (FASTQC, alignment, quantification), I led the downstream analytical workflow that gave biological meaning to our candidate lncRNAs.

---

## Scientific Background

Immunotherapy has revolutionized cancer treatment, but biomarkers for predicting patient response in HCC remain inadequate. Our project hypothesized that lncRNAs could serve as novel biomarkers. My work addressed the critical question: **"What biological functions do these correlated lncRNAs actually regulate?"**

The **guilt-by-association** method operates on the principle that functionally related genes are often co-expressed. By identifying the protein-coding partners of each lncRNA and analyzing their enriched pathways, we can infer the lncRNA's potential biological role.

---

## Technical Implementation

### 📁 Repository Structure

├── 📄 README.md # You are here\
├── 📁 scripts/\
│ └── guilt\_by\_association\_analysis.R # Main analysis pipeline\
├── 📁 results/\
│ ├── figures/ # Publication-ready visualizations\
│ └── tables/ # Analysis outputs & summaries\
└── 📁 docs/\
├── methodology.md # Technical details\
└── contribution\_statement.md # My specific role



### Analysis Pipeline
I developed and executed a comprehensive analytical workflow:

1.  **Input Processing:** Utilized correlation matrices from team members
2.  **Partner Identification:** Extracted top 50 protein-coding partners for 3,645 lncRNAs
3.  **Functional Enrichment:** Performed GO & KEGG analysis using `clusterProfiler`
4.  **Statistical Rigor:** Applied FDR correction (q-value < 0.05)
5.  **Results Synthesis:** Generated functional annotations and candidate prioritization
6.  **Visualization:** Created publication-quality figures and tables

### Quick Reproduction
```bash
# Clone this repository
git clone https://github.com/yanny-alt/lncRNA-immunotherapy-functional-inference.git

# Execute the full analysis pipeline
Rscript scripts/guilt_by_association_analysis.R
``` 

## Key Findings & Impact

### Summary Statistics

| Metric                              | Result                           |
| :---------------------------------- | :------------------------------- |
| Total lncRNAs Analyzed              | 3,645                            |
| lncRNAs with Significant Enrichment | 50                               |
| Average GO BP Terms per lncRNA      | 475.1                            |
| Average KEGG Pathways per lncRNA    | 39.4                             |
| Primary Functional Category         | Immune Response Regulation (78%) |


### Biological Insights

My analysis revealed that the majority of functionally enriched lncRNAs are involved in immune regulation, with strong connections to:

- PD-1/PD-L1 checkpoint pathway

- JAK-STAT signaling

- T-cell activation and differentiation

- Cytokine-cytokine receptor interaction


### Top Candidate Biomarkers

I identified and prioritized 15 high-confidence lncRNA biomarkers based on functional diversity scores. The top candidate, `ENSG00000242588`, showed enrichment in 823 GO terms and 472 KEGG pathways with a confidence score of 1,295.

***


## Results Gallery

| Functional Distribution | Top Candidate lncRNAs |
| :--- | :--- |
| <img src="https://github.com/yanny-alt/lncRNA-immunotherapy-functional-inference/raw/main/results/figures/functional_distribution_pie.png" width="400"> | <img src="https://github.com/yanny-alt/lncRNA-immunotherapy-functional-inference/raw/main/results/figures/summary_top_lncRNAs.png" width="400"> |
| _78% of annotated lncRNAs regulate immune response_ | _Top 15 candidates ranked by functional diversity_ |

| Co-expression Network | Pathway Enrichment |
| :--- | :--- |
| <img src="https://github.com/yanny-alt/lncRNA-immunotherapy-functional-inference/raw/main/results/figures/top_lncRNAs_correlation_heatmap.png" width="400"> | <img src="https://github.com/yanny-alt/lncRNA-immunotherapy-functional-inference/raw/main/results/figures/summary_correlation_vs_enrichment.png" width="400"> |
| _Modular organization of lncRNA-protein partnerships_ | _Correlation vs Enrichment_ |

***


## Skills Demonstrated

Bioinformatics & Computational Biology

- RNA-seq downstream analysis

- Functional enrichment (GO, KEGG)

- Correlation network analysis

- Statistical analysis & multiple testing correction

Programming & Data Science

- R programming (tidyverse, Bioconductor)

- Data visualization (ggplot2, pheatmap)

- Reproducible research practices

- Version control with Git/GitHub

Scientific Research

- Hypothesis-driven analysis

- Biological interpretation

- Research documentation

- Collaborative science

***


## Citation & Attribution

This work represents my independent contribution to a collaborative project conducted during my Bioinformatics Research Internship at **HackBio**. When our manuscript is published, the citation will be updated here.

> Acknowledgments: I thank my research team members for their contributions to data preprocessing and Dr. Adewale Ogunleye at HackBio for supervision and guidance throughout this capstone project.


***


## Contact

Favour O. Igwezeke\
Faculty of Pharmaceutical Sciences\
University of Nigeria, Nsukka\
\beingfave@gmail.com | \[[LinkedIn Profile]](https://www.linkedin.com/in/favourokechukwu/details/experience/) 

***

_This project was completed as part of the Bioinformatics Research Internship Capstone Project at HackBio (Sept 2025 – Nov 2025)._
