# My Contribution to lncRNA Biomarkers for Cancer Immunotherapy

## Personal Statement

This repository documents my independent contribution to the manuscript "Long Non-Coding RNAs as Biomarkers for Cancer Immunotherapy." My work centered on the functional characterization of candidate lncRNAs using guilt-by-association methods.

## Specific Tasks Performed by Me

*   **Methodology:** Designed and implemented the complete guilt-by-association analysis pipeline from scratch.
*   **Data Analysis:**
    *   Processed correlation matrices containing 3,645 lncRNAs
    *   Extracted top 50 co-expressed protein-coding partners for each lncRNA
    *   Performed Gene Ontology (GO) and KEGG pathway enrichment analysis using `clusterProfiler`
    *   Statistically evaluated and filtered results (FDR < 0.05, q-value < 0.1)
    *   Generated functional diversity scores to prioritize top candidates
*   **Interpretation:** Inferred biological functions for 50 high-quality lncRNAs, identifying 78% with strong bias towards immune regulation.
*   **Visualization:** Created all 5 main figures and 15+ tables for the functional inference results.
*   **Manuscript Writing:** Wrote Methods sections 2.6.1-2.6.2 and Results section 3.5 describing this analysis.

## Distinction from Other Team Members' Work

Other team members were responsible for:
*   Raw RNA-seq data quality control (FASTQC)
*   Read alignment and quantification using STAR/featureCounts
*   Initial differential expression and correlation analysis

My work begins where theirs ends, taking the list of correlated lncRNAs and answering the critical question: **"What do these lncRNAs actually do?"**
