




if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
install.packages(c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer"))

install.packages(c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer"))
unlink(tempdir(), recursive = TRUE)
unlink(tempdir(), recursive = TRUE)

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("textshaping", type = "source")
install.packages("ragg", type = "source")
install.packages("ggplot2")
install.packages("tidyverse", type = "source")
install.packages("pheatmap")

options(repos = c(CRAN = "http://cloud.r-project.org"))
install.packages("RColorBrewer")


install.packages("RColorBrewer")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

# Load required libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


setwd("/Users/favourigwezeke/Personal_System/Research/HACKBIO_2025/IncRNAs/data")


getwd()


dir.create("results/guilt_by_association", recursive = TRUE, showWarnings = FALSE)
dir.create("results/guilt_by_association/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/guilt_by_association/tables", recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("Functional Inference Analysis - Guilt by Association\n")
cat("============================================================\n\n")

# ============================================================================
# MODIFIED STEP 1: REBUILD CORRELATION MATRIX WITH lncRNAs + PROTEIN-CODING
# ============================================================================

cat("STEP 1: Building correlation matrix with lncRNAs AND protein-coding genes...\n")

# Load the original count data
cat("Loading count data from featureCounts...\n")
count_data <- read.table(
  "featureCounts.out.txt", 
  header = TRUE, 
  comment.char = "#", 
  row.names = 1,
  stringsAsFactors = FALSE
)

# Identify count columns (adjust indices based on your actual data structure)
count_columns <- count_data[, 6:ncol(count_data)]  
cat("Raw count data dimensions:", dim(count_columns), "\n")

# Debug: Check what the row names actually look like
cat("First few row names from count data:\n")
print(head(rownames(count_columns)))

# Load gene annotations 
cat("Loading gene annotations for filtering...\n")
feature_info <- read.table(
  "featureCounts.out.txt",
  header = TRUE,
  comment.char = "#",
  stringsAsFactors = FALSE
)

# Create proper gene lookup
gene_lookup <- feature_info %>%
  dplyr::select(Geneid) %>%
  distinct() %>%
  mutate(Geneid_clean = sub("\\..*$", "", Geneid))

cat("First few Geneid_clean from lookup:\n")
print(head(gene_lookup$Geneid_clean))

# Check if there's a mismatch in ID formats
cat("Checking ID overlap...\n")
count_rownames_clean <- sub("\\..*$", "", rownames(count_columns))
overlap <- intersect(count_rownames_clean, gene_lookup$Geneid_clean)
cat("Genes overlapping between count data and lookup:", length(overlap), "\n")

# Get biomaRt annotations
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
annotations <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id", 
  values = gene_lookup$Geneid_clean,
  mart = mart
)

gene_lookup <- left_join(gene_lookup, annotations, by = c("Geneid_clean" = "ensembl_gene_id"))
gene_lookup_unique <- gene_lookup %>% 
  dplyr::select(Geneid_clean, hgnc_symbol, gene_biotype) %>% 
  distinct(Geneid_clean, .keep_all = TRUE)

# Filter genes to include only lncRNAs and protein-coding
genes_to_include <- gene_lookup_unique %>%
  filter(gene_biotype %in% c("lncRNA", "protein_coding")) %>%
  mutate(match_id = case_when(
    !is.na(hgnc_symbol) & hgnc_symbol != "" ~ hgnc_symbol,
    TRUE ~ Geneid_clean
  ))

cat("Genes selected for correlation matrix:\n")
cat("  - lncRNAs:", sum(genes_to_include$gene_biotype == "lncRNA"), "\n")
cat("  - Protein-coding:", sum(genes_to_include$gene_biotype == "protein_coding"), "\n")

# FIXED: Use cleaned row names for matching
count_rownames_clean <- sub("\\..*$", "", rownames(count_columns))
counts_filtered <- count_columns[count_rownames_clean %in% genes_to_include$Geneid_clean, ]

# Update row names to cleaned version for proper mapping
rownames(counts_filtered) <- count_rownames_clean[count_rownames_clean %in% genes_to_include$Geneid_clean]

cat("Filtered count matrix dimensions:", dim(counts_filtered), "\n")

if(nrow(counts_filtered) == 0) {
  cat("ERROR: No genes matched after filtering. Checking why...\n")
  cat("Sample count row names (cleaned):", head(count_rownames_clean), "\n")
  cat("Sample genes_to_include Geneid_clean:", head(genes_to_include$Geneid_clean), "\n")
  stop("No overlap between count data and gene annotations. Check ID formats.")
}

# Remove genes with zero variance (constant expression)
nonzero_var_genes <- apply(counts_filtered, 1, var) > 0
counts_filtered <- counts_filtered[nonzero_var_genes, ]
cat("After removing zero-variance genes:", dim(counts_filtered), "\n")

# Convert to matrix and log2 transform
count_matrix <- as.matrix(counts_filtered)
log2_counts <- log2(count_matrix + 1)  # log2(x+1) transformation

cat("STRATEGIC APPROACH: Multi-tier correlation analysis\n")

# TIER 1: Quick discovery with highly variable genes
cat("TIER 1: Highly variable genes for initial discovery...\n")
top_genes <- 10000  # Balanced approach
gene_variance <- apply(log2_counts, 1, var)
high_var_genes <- names(sort(gene_variance, decreasing = TRUE))[1:top_genes]

# TIER 2: Ensure your lncRNAs of interest are included
lncRNA_list <- read.csv("lncRNA_extracted.csv", stringsAsFactors = FALSE) %>%
  mutate(ensembl_gene_id_clean = sub("\\..*$", "", ensembl_gene_id))

target_lncRNAs <- lncRNA_list$ensembl_gene_id_clean

# Combine strategies: high variance + your lncRNAs
final_genes <- unique(c(
  high_var_genes,
  target_lncRNAs[target_lncRNAs %in% rownames(log2_counts)]
))

cat("Final gene set:", length(final_genes), "genes\n")
cat("  - High variance genes:", sum(final_genes %in% high_var_genes), "\n")
cat("  - Your lncRNAs:", sum(final_genes %in% target_lncRNAs), "\n")

# Build matrix
log2_final <- log2_counts[final_genes, ]
correlation_matrix <- cor(t(log2_final))

# Set names
id_mapping <- setNames(genes_to_include$match_id, genes_to_include$Geneid_clean)
common_genes <- rownames(correlation_matrix)[rownames(correlation_matrix) %in% names(id_mapping)]
rownames(correlation_matrix) <- id_mapping[common_genes]
colnames(correlation_matrix) <- id_mapping[common_genes]

cat("✓ STRATEGIC matrix built:", dim(correlation_matrix), "\n")









# ============================================================================
# STEP 2: Loading gene annotations (already done above, but keeping for clarity)
# ============================================================================

cat("STEP 2: Gene annotations already loaded during matrix construction\n")
cat("  Total annotated genes:", nrow(gene_lookup), "\n")
cat("  Unique genes in lookup:", nrow(gene_lookup_unique), "\n\n")

# ============================================================================
# STEP 3: Identifying gene types (UPDATED for new matrix)
# ============================================================================

cat("STEP 3: Identifying gene types in NEW correlation matrix...\n")

# --- Step 3a: Clean lncRNA list Ensembl IDs ---
lncRNA_list <- read.csv("lncRNA_extracted.csv", stringsAsFactors = FALSE) %>%
  mutate(ensembl_gene_id_clean = sub("\\..*$", "", ensembl_gene_id))
cat("✓ Loaded", nrow(lncRNA_list), "candidate lncRNAs\n")

# --- Step 3b: Map lncRNA Ensembl IDs to HGNC symbol if available ---
lncRNA_symbols <- lncRNA_list %>%
  dplyr::select(ensembl_gene_id, ensembl_gene_id_clean) %>%
  left_join(
    gene_lookup_unique %>% dplyr::select(Geneid_clean, hgnc_symbol),
    by = c("ensembl_gene_id_clean" = "Geneid_clean")
  )

# --- Step 3c: Create match ID ---
lncRNA_symbols$match_id <- lncRNA_symbols$ensembl_gene_id_clean
has_hgnc <- !is.na(lncRNA_symbols$hgnc_symbol) & lncRNA_symbols$hgnc_symbol != ""
lncRNA_symbols$match_id[has_hgnc] <- lncRNA_symbols$hgnc_symbol[has_hgnc]

# --- Step 3d: Identify lncRNAs in NEW correlation matrix ---
lncRNA_ids_in_matrix <- intersect(lncRNA_symbols$match_id, rownames(correlation_matrix))
cat("✓ lncRNAs found in NEW correlation matrix:", length(lncRNA_ids_in_matrix), "\n")

# --- Step 3e: Protein-coding genes in NEW correlation matrix ---
protein_coding_in_matrix <- genes_to_include %>%
  filter(gene_biotype == "protein_coding") %>%
  pull(match_id) %>%
  intersect(rownames(correlation_matrix))

cat("✓ Protein-coding genes in NEW correlation matrix:", length(protein_coding_in_matrix), "\n")

cat("\n=== MATRIX COMPOSITION ===\n")
cat("Total genes in matrix:", nrow(correlation_matrix), "\n")
cat("lncRNAs:", length(lncRNA_ids_in_matrix), "\n") 
cat("Protein-coding:", length(protein_coding_in_matrix), "\n")
cat("Ratio:", round(length(protein_coding_in_matrix) / length(lncRNA_ids_in_matrix), 2), "protein-coding per lncRNA\n")

# Save the new matrix for future use
write.csv(correlation_matrix, "results/guilt_by_association/lncRNA_protein_coding_correlation_matrix.csv")
cat("✓ New correlation matrix saved for future use\n\n")


# ============================================================================
# STEP 4: Smart Guilt-by-Association Analysis
# ============================================================================

cat("STEP 4: Performing smart guilt-by-association analysis...\n")
cat("Processing all", length(lncRNA_ids_in_matrix), "lncRNAs with top 50 partners...\n")
cat("(Estimated time: 10-15 minutes)\n\n")

top_coexpressed <- list()
summary_stats <- data.frame()

for (i in seq_along(lncRNA_ids_in_matrix)) {
  lnc_id <- lncRNA_ids_in_matrix[i]
  
  # Progress indicator
  if (i %% 100 == 0) {
    cat("  Processed", i, "of", length(lncRNA_ids_in_matrix), "lncRNAs\n")
  }
  
  # Extract correlations with protein-coding genes
  lnc_correlations <- correlation_matrix[lnc_id, protein_coding_in_matrix]
  
  # Create results data frame
  cor_df <- data.frame(
    gene = protein_coding_in_matrix,
    correlation = as.numeric(lnc_correlations),
    stringsAsFactors = FALSE
  )
  
  # Keep only positive correlations and get top 50
  cor_df <- cor_df %>%
    filter(correlation > 0) %>%
    arrange(desc(correlation)) %>%
    head(50)
  
  # Store results
  top_coexpressed[[lnc_id]] <- cor_df
  
  # Save summary
  if (nrow(cor_df) > 0) {
    summary_stats <- rbind(summary_stats, data.frame(
      lncRNA = lnc_id,
      n_partners = nrow(cor_df),
      mean_correlation = mean(cor_df$correlation),
      max_correlation = max(cor_df$correlation),
      stringsAsFactors = FALSE
    ))
  }
}

cat("✓ Found co-expressed partners for", length(top_coexpressed), "lncRNAs\n")

# Save ALL results in one file (not individual CSVs)
saveRDS(top_coexpressed, "results/guilt_by_association/all_lncRNA_partners.rds")
write.csv(summary_stats, "results/guilt_by_association/summary_stats.csv", row.names = FALSE)

cat("✓ All results saved in single RDS file\n\n")



# ============================================================================
# STEP 5: FIXED Functional Enrichment Analysis
# ============================================================================

cat("STEP 5: Performing functional enrichment analysis...\n")

enrichment_results <- list()
successful_analyses <- 0

for (i in seq_along(top_lncRNAs)) {
  lnc_id <- top_lncRNAs[i]
  partners_df <- top_coexpressed[[lnc_id]]
  
  lnc_mean_corr <- top_lncRNAs_df$mean_correlation[top_lncRNAs_df$lncRNA == lnc_id]
  
  cat("  [", i, "/", length(top_lncRNAs), "] Analyzing ", lnc_id, 
      " (avg corr: ", round(lnc_mean_corr, 3), ")\n", sep = "")
  
  # FIX: Remove self-correlation (the lncRNA itself)
  partners_df <- partners_df %>%
    filter(gene != lnc_id)  # Remove the lncRNA from its own partners
  
  gene_symbols <- partners_df$gene
  
  # Remove any NA or empty strings
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  
  if (length(gene_symbols) < 10) {
    cat("    ⚠ Only", length(gene_symbols), "genes available - skipping\n")
    next
  }
  
  cat("    Using", length(gene_symbols), "protein-coding partners\n")
  
  tryCatch({
    # Convert gene symbols to Entrez IDs
    entrez_ids <- clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    
    if (nrow(entrez_ids) < 10) {
      cat("    ⚠ Only", nrow(entrez_ids), "Entrez IDs mapped\n")
      next
    }
    
    cat("    Mapped", nrow(entrez_ids), "genes to Entrez IDs\n")
    
    # GO Biological Process enrichment
    ego_bp <- enrichGO(
      gene = entrez_ids$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1,
      readable = TRUE
    )
    
    # KEGG Pathway enrichment
    kegg_enrich <- enrichKEGG(
      gene = entrez_ids$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
    
    # Store results
    enrichment_results[[lnc_id]] <- list(
      GO_BP = ego_bp,
      KEGG = kegg_enrich,
      input_genes = gene_symbols,
      entrez_ids = entrez_ids,
      n_genes_input = length(gene_symbols),
      n_genes_mapped = nrow(entrez_ids),
      mean_correlation = lnc_mean_corr
    )
    
    # Save individual results if significant
    if (!is.null(ego_bp) && nrow(ego_bp@result) > 0) {
      write.csv(ego_bp@result, 
                paste0("results/guilt_by_association/tables/", 
                       lnc_id, "_GO_BP.csv"),
                row.names = FALSE)
      cat("    ✓ GO BP:", nrow(ego_bp@result), "terms\n")
    } else {
      cat("    ○ GO BP: No significant terms\n")
    }
    
    if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
      write.csv(kegg_enrich@result, 
                paste0("results/guilt_by_association/tables/", 
                       lnc_id, "_KEGG.csv"),
                row.names = FALSE)
      cat("    ✓ KEGG:", nrow(kegg_enrich@result), "pathways\n")
    } else {
      cat("    ○ KEGG: No significant pathways\n")
    }
    
    successful_analyses <- successful_analyses + 1
    
  }, error = function(e) {
    cat("    ✗ Error:", e$message, "\n")
  })
  
  cat("\n")
}

# Save all enrichment results
saveRDS(enrichment_results, 
        "results/guilt_by_association/enrichment_results_top50.rds")

cat("\n✓ Enrichment analysis completed for", successful_analyses, 
    "out of", length(top_lncRNAs), "lncRNAs\n")

# Create summary of enrichment results
enrichment_summary <- data.frame()
for (lnc_id in names(enrichment_results)) {
  res <- enrichment_results[[lnc_id]]
  summary_row <- data.frame(
    lncRNA = lnc_id,
    n_genes_input = res$n_genes_input,
    n_genes_mapped = res$n_genes_mapped,
    mean_correlation = res$mean_correlation,
    GO_BP_terms = ifelse(!is.null(res$GO_BP), nrow(res$GO_BP@result), 0),
    KEGG_pathways = ifelse(!is.null(res$KEGG), nrow(res$KEGG@result), 0),
    stringsAsFactors = FALSE
  )
  enrichment_summary <- rbind(enrichment_summary, summary_row)
}

write.csv(enrichment_summary, 
          "results/guilt_by_association/enrichment_summary.csv",
          row.names = FALSE)

cat("✓ Enrichment summary saved\n\n")


# ============================================================================
# STEP 6: COMPREHENSIVE VISUALIZATION SUITE
# ============================================================================

cat("STEP 6: Generating comprehensive publication-quality visualizations...\n\n")

# Create figures directory if it doesn't exist
dir.create("results/guilt_by_association/figures", recursive = TRUE, showWarnings = FALSE)

# Load your enrichment results
enrichment_results <- readRDS("results/guilt_by_association/enrichment_results_top50.rds")
enrichment_summary <- read.csv("results/guilt_by_association/enrichment_summary.csv")

# ============================================================================
# 6A: SUMMARY LEVEL VISUALIZATIONS (Across all lncRNAs)
# ============================================================================

cat("6A: Creating summary-level visualizations across all lncRNAs...\n")

# 1. Distribution of Enrichment Results
tryCatch({
  p_summary1 <- ggplot(enrichment_summary, aes(x = GO_BP_terms)) +
    geom_histogram(fill = "steelblue", alpha = 0.8, bins = 20) +
    labs(title = "Distribution of GO BP Terms per lncRNA",
         x = "Number of GO BP Terms",
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave("results/guilt_by_association/figures/summary_GO_BP_distribution.png", 
         p_summary1, width = 10, height = 6, dpi = 300)
  cat("✓ Summary: GO BP distribution plot saved\n")
}, error = function(e) cat("⚠ Could not create GO BP distribution plot\n"))

# 2. Correlation vs Enrichment Strength
tryCatch({
  p_summary2 <- ggplot(enrichment_summary, aes(x = mean_correlation, y = GO_BP_terms)) +
    geom_point(alpha = 0.7, color = "darkred", size = 2) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    labs(title = "Correlation Strength vs Functional Enrichment",
         x = "Mean Correlation with Partners",
         y = "Number of GO BP Terms") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave("results/guilt_by_association/figures/summary_correlation_vs_enrichment.png", 
         p_summary2, width = 10, height = 6, dpi = 300)
  cat("✓ Summary: Correlation vs enrichment plot saved\n")
}, error = function(e) cat("⚠ Could not create correlation vs enrichment plot\n"))

# 3. Top lncRNAs by Enrichment Strength
tryCatch({
  top_lncRNAs <- enrichment_summary %>%
    arrange(desc(GO_BP_terms)) %>%
    head(15)
  
  p_summary3 <- ggplot(top_lncRNAs, aes(x = reorder(lncRNA, GO_BP_terms), y = GO_BP_terms)) +
    geom_col(fill = "purple", alpha = 0.8) +
    coord_flip() +
    labs(title = "Top 15 lncRNAs by GO BP Enrichment",
         x = "lncRNA",
         y = "Number of GO BP Terms") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave("results/guilt_by_association/figures/summary_top_lncRNAs.png", 
         p_summary3, width = 12, height = 8, dpi = 300)
  cat("✓ Summary: Top lncRNAs bar plot saved\n")
}, error = function(e) cat("⚠ Could not create top lncRNAs plot\n"))

# ============================================================================
# 6B: INDIVIDUAL lncRNA VISUALIZATIONS
# ============================================================================

cat("\n6B: Creating individual lncRNA visualizations...\n")

successful_plots <- 0

for (lnc_name in names(enrichment_results)) {
  cat("  Creating plots for", lnc_name, "...\n")
  
  results <- enrichment_results[[lnc_name]]
  
  # GO Biological Process Visualizations
  if (!is.null(results$GO_BP) && nrow(results$GO_BP@result) > 0) {
    
    # 1. Dotplot (Top 15 terms)
    tryCatch({
      p1 <- dotplot(results$GO_BP, 
                    showCategory = 15, 
                    title = paste(lnc_name, "- GO Biological Process"),
                    font.size = 10) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      ggsave(paste0("results/guilt_by_association/figures/", 
                    lnc_name, "_GO_BP_dotplot.png"),
             p1, width = 12, height = 8, dpi = 300)
      cat("    ✓ GO BP dotplot saved\n")
    }, error = function(e) {
      cat("    ⚠ Could not create GO BP dotplot for", lnc_name, "\n")
    })
    
    # 2. Barplot (Top 15 terms)
    tryCatch({
      p2 <- barplot(results$GO_BP, 
                    showCategory = 15, 
                    title = paste(lnc_name, "- GO Biological Process"),
                    font.size = 10) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      ggsave(paste0("results/guilt_by_association/figures/", 
                    lnc_name, "_GO_BP_barplot.png"),
             p2, width = 12, height = 8, dpi = 300)
      cat("    ✓ GO BP barplot saved\n")
    }, error = function(e) {
      cat("    ⚠ Could not create GO BP barplot for", lnc_name, "\n")
    })
    
    # 3. Network plot (if enough terms)
    if (nrow(results$GO_BP@result) >= 8) {
      tryCatch({
        p3 <- cnetplot(results$GO_BP, 
                       showCategory = 8,
                       colorEdge = TRUE,
                       circular = FALSE,
                       node_label = "all") +
          ggtitle(paste(lnc_name, "- GO BP Gene-Concept Network"))
        
        ggsave(paste0("results/guilt_by_association/figures/", 
                      lnc_name, "_GO_BP_network.png"),
               p3, width = 14, height = 10, dpi = 300)
        cat("    ✓ GO BP network plot saved\n")
      }, error = function(e) {
        cat("    ⚠ Could not create GO BP network plot for", lnc_name, "\n")
      })
    }
  }
  
  # KEGG Pathway Visualizations
  if (!is.null(results$KEGG) && nrow(results$KEGG@result) > 0) {
    
    # 1. Dotplot (Top 15 pathways)
    tryCatch({
      p4 <- dotplot(results$KEGG, 
                    showCategory = 15, 
                    title = paste(lnc_name, "- KEGG Pathways"),
                    font.size = 10) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      ggsave(paste0("results/guilt_by_association/figures/", 
                    lnc_name, "_KEGG_dotplot.png"),
             p4, width = 12, height = 8, dpi = 300)
      cat("    ✓ KEGG dotplot saved\n")
    }, error = function(e) {
      cat("    ⚠ Could not create KEGG dotplot for", lnc_name, "\n")
    })
    
    # 2. Barplot (Top 15 pathways)
    tryCatch({
      p5 <- barplot(results$KEGG, 
                    showCategory = 15, 
                    title = paste(lnc_name, "- KEGG Pathways"),
                    font.size = 10) +
        theme(plot.title = element_text(size = 12, face = "bold"))
      
      ggsave(paste0("results/guilt_by_association/figures/", 
                    lnc_name, "_KEGG_barplot.png"),
             p5, width = 12, height = 8, dpi = 300)
      cat("    ✓ KEGG barplot saved\n")
    }, error = function(e) {
      cat("    ⚠ Could not create KEGG barplot for", lnc_name, "\n")
    })
  }
  
  successful_plots <- successful_plots + 1
  cat("    ✓ Completed", lnc_name, "-", successful_plots, "/50\n\n")
}




tryCatch({
  # Load ALL required data files
  correlation_matrix <- read.csv("results/guilt_by_association/lncRNA_protein_coding_correlation_matrix.csv", 
                                 row.names = 1)
  top_coexpressed <- readRDS("results/guilt_by_association/all_lncRNA_partners.rds")
  enrichment_summary <- read.csv("results/guilt_by_association/enrichment_summary.csv")
  
  cat("✓ All data files loaded successfully\n")
  cat("  enrichment_summary dimensions:", dim(enrichment_summary), "\n")
  cat("  top_coexpressed contains", length(top_coexpressed), "lncRNAs\n")
  
  # Select top 10 lncRNAs by enrichment strength
  top_10_lncRNAs <- enrichment_summary[order(-enrichment_summary$GO_BP_terms), "lncRNA"][1:10]
  
  cat("  Selected top 10 lncRNAs:", paste(top_10_lncRNAs, collapse = ", "), "\n")
  
  # Get top 15 partners for each lncRNA
  heatmap_data <- data.frame()
  
  for (lncRNA in top_10_lncRNAs) {
    if (lncRNA %in% names(top_coexpressed)) {
      partners <- top_coexpressed[[lncRNA]][1:15, ]
      partners$lncRNA <- lncRNA
      heatmap_data <- rbind(heatmap_data, partners)
    }
  }
  
  cat("  Collected", nrow(heatmap_data), "lncRNA-partner pairs\n")
  
  if (nrow(heatmap_data) > 0) {
    # Get the most common partners across all lncRNAs
    gene_counts <- table(heatmap_data$gene)
    top_partners <- names(sort(gene_counts, decreasing = TRUE))[1:30]
    
    cat("  Selected top 30 most frequent partner genes\n")
    
    # Filter heatmap data to only include these top partners
    heatmap_filtered <- heatmap_data[heatmap_data$gene %in% top_partners, ]
    
    # Create matrix using simpler method
    heatmap_matrix <- matrix(0,  # Initialize with 0 instead of NA
                             nrow = length(top_10_lncRNAs),
                             ncol = length(top_partners),
                             dimnames = list(top_10_lncRNAs, top_partners))
    
    # Fill the matrix with correlation values
    for (i in 1:nrow(heatmap_filtered)) {
      lnc <- heatmap_filtered$lncRNA[i]
      gene <- heatmap_filtered$gene[i]
      corr <- heatmap_filtered$correlation[i]
      
      if (lnc %in% rownames(heatmap_matrix) && gene %in% colnames(heatmap_matrix)) {
        heatmap_matrix[lnc, gene] <- corr
      }
    }
    
    cat("  Created heatmap matrix:", dim(heatmap_matrix), "\n")
    cat("  Matrix summary - Min:", min(heatmap_matrix), "Max:", max(heatmap_matrix), 
        "Mean:", round(mean(heatmap_matrix), 3), "\n")
    
    if (nrow(heatmap_matrix) > 1 && ncol(heatmap_matrix) > 1) {
      # Plot heatmap
      pheatmap(
        heatmap_matrix,
        filename = "results/guilt_by_association/figures/top_lncRNAs_correlation_heatmap.png",
        main = "Top lncRNAs - Protein Coding Correlations",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        show_rownames = TRUE,
        show_colnames = TRUE,
        fontsize = 8,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        width = 14,
        height = 8
      )
      cat("✓ Correlation heatmap saved successfully!\n")
    } else {
      cat("⚠ Heatmap matrix too small after filtering\n")
    }
  } else {
    cat("⚠ No heatmap data collected\n")
  }
}, error = function(e) {
  cat("✗ Error in main heatmap:", e$message, "\n")
})


# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("VISUALIZATION COMPLETE!\n")
cat("============================================================\n")
cat("Generated for", successful_plots, "lncRNAs:\n")
cat("- 3 summary plots (across all lncRNAs)\n") 
cat("- 3-5 individual plots per lncRNA\n")
cat("- 1 comprehensive correlation heatmap\n")
cat("\nTotal files created: ~", successful_plots * 4 + 4, "visualization files\n")
cat("All saved in: results/guilt_by_association/figures/\n")
cat("============================================================\n")




# ============================================================================
# STEP 7: COMPREHENSIVE FUNCTIONAL ANNOTATION SUMMARY
# ============================================================================

cat("STEP 7: Generating comprehensive functional annotation summary...\n\n")

# Load required data
enrichment_results <- readRDS("results/guilt_by_association/enrichment_results_top50.rds")
enrichment_summary <- read.csv("results/guilt_by_association/enrichment_summary.csv")
top_coexpressed <- readRDS("results/guilt_by_association/all_lncRNA_partners.rds")

# Initialize summary table
functional_summary <- data.frame()

for (lnc_name in names(enrichment_results)) {
  results <- enrichment_results[[lnc_name]]
  
  # Extract top 3 GO BP terms with p-values
  if (!is.null(results$GO_BP) && nrow(results$GO_BP@result) > 0) {
    top_go_bp <- head(results$GO_BP@result, 3)
    top_go_bp_text <- paste(
      paste0(top_go_bp$Description, " (p=", 
             formatC(top_go_bp$pvalue, format = "e", digits = 2), ")"),
      collapse = "; "
    )
  } else {
    top_go_bp_text <- "No significant enrichment"
  }
  
  # Extract top 3 KEGG pathways with p-values
  if (!is.null(results$KEGG) && nrow(results$KEGG@result) > 0) {
    top_kegg <- head(results$KEGG@result, 3)
    top_kegg_text <- paste(
      paste0(top_kegg$Description, " (p=", 
             formatC(top_kegg$pvalue, format = "e", digits = 2), ")"),
      collapse = "; "
    )
  } else {
    top_kegg_text <- "No significant enrichment"
  }
  
  # Infer biological function based on enrichment results
  inferred_function <- "Unknown function"
  if (!is.null(results$GO_BP) && nrow(results$GO_BP@result) > 0) {
    all_descriptions <- tolower(paste(results$GO_BP@result$Description, collapse = " "))
    
    # Define functional categories
    immune_keywords <- c("immune", "t.cell", "t cell", "cytokine", "interferon", 
                         "antigen", "lymphocyte", "leukocyte", "activation",
                         "chemokine", "interleukin", "inflammation", "defense")
    epigenetic_keywords <- c("chromatin", "histone", "methylation", "acetylation",
                             "transcription", "dna binding", "rna polymerase", "rna processing")
    metabolic_keywords <- c("metabol", "biosynthetic", "catabolic", "oxidation",
                            "glycolysis", "phosphorylation", "atp", "respir")
    signaling_keywords <- c("signal", "pathway", "kinase", "receptor", "transduction")
    development_keywords <- c("development", "differentiation", "morphogenesis", "growth")
    
    # Score each category
    scores <- c(
      immune = sum(sapply(immune_keywords, function(kw) grepl(kw, all_descriptions))),
      epigenetic = sum(sapply(epigenetic_keywords, function(kw) grepl(kw, all_descriptions))),
      metabolic = sum(sapply(metabolic_keywords, function(kw) grepl(kw, all_descriptions))),
      signaling = sum(sapply(signaling_keywords, function(kw) grepl(kw, all_descriptions))),
      development = sum(sapply(development_keywords, function(kw) grepl(kw, all_descriptions)))
    )
    
    # Assign function based on highest score
    max_score <- max(scores)
    if (max_score > 0) {
      top_category <- names(scores)[which.max(scores)]
      inferred_function <- switch(top_category,
                                  immune = "Immune response regulation",
                                  epigenetic = "Epigenetic/chromatin regulation", 
                                  metabolic = "Metabolic process regulation",
                                  signaling = "Cell signaling pathway",
                                  development = "Development/differentiation"
      )
    } else {
      inferred_function <- "General cellular regulation"
    }
  }
  
  # Get top 5 correlated protein partners
  if (lnc_name %in% names(top_coexpressed)) {
    top_partners <- head(top_coexpressed[[lnc_name]], 5)
    top_partners_text <- paste(
      paste0(top_partners$gene, " (r=", round(top_partners$correlation, 3), ")"),
      collapse = ", "
    )
  } else {
    top_partners_text <- "No partner data"
  }
  
  # Calculate functional diversity score
  functional_diversity <- ifelse(!is.null(results$GO_BP), 
                                 nrow(results$GO_BP@result), 0) +
    ifelse(!is.null(results$KEGG), 
           nrow(results$KEGG@result), 0)
  
  # Create comprehensive summary row
  summary_row <- data.frame(
    lncRNA = lnc_name,
    Top_GO_BP_Terms = top_go_bp_text,
    Top_KEGG_Pathways = top_kegg_text,
    Inferred_Function = inferred_function,
    Top_Correlated_Partners = top_partners_text,
    Functional_Diversity_Score = functional_diversity,
    N_GO_BP_Terms = ifelse(!is.null(results$GO_BP), nrow(results$GO_BP@result), 0),
    N_KEGG_Pathways = ifelse(!is.null(results$KEGG), nrow(results$KEGG@result), 0),
    Mean_Correlation = results$mean_correlation,
    N_Coexpressed_Genes = results$n_genes_input,
    stringsAsFactors = FALSE
  )
  
  functional_summary <- rbind(functional_summary, summary_row)
}

# Sort by functional diversity score (most functionally rich first)
functional_summary <- functional_summary[order(-functional_summary$Functional_Diversity_Score), ]

# Save functional summary
write.csv(functional_summary, 
          "results/guilt_by_association/tables/functional_annotation_summary.csv",
          row.names = FALSE)

cat("✓ Functional annotation summary saved\n")




# ============================================================================
# STEP 7B: CREATE MASTER ANALYSIS SUMMARY
# ============================================================================

cat("\nSTEP 7B: Creating master analysis summary...\n")

# Load additional data for comprehensive summary
correlation_matrix <- read.csv("results/guilt_by_association/lncRNA_protein_coding_correlation_matrix.csv", 
                               row.names = 1)
summary_stats <- read.csv("results/guilt_by_association/summary_stats.csv")

# Create master summary
master_summary <- data.frame(
  Analysis_Date = Sys.Date(),
  Total_Genes_in_Matrix = nrow(correlation_matrix),
  Total_lncRNAs_Analyzed = length(unique(summary_stats$lncRNA)),
  Total_Protein_Coding_Genes = length(unique(unlist(lapply(top_coexpressed, function(x) x$gene)))),
  lncRNAs_with_Functional_Enrichment = nrow(functional_summary),
  Average_GO_BP_Terms_per_lncRNA = round(mean(functional_summary$N_GO_BP_Terms), 1),
  Average_KEGG_Pathways_per_lncRNA = round(mean(functional_summary$N_KEGG_Pathways), 1),
  Average_Correlation_Strength = round(mean(functional_summary$Mean_Correlation), 3),
  Top_Functional_Category = names(sort(table(functional_summary$Inferred_Function), decreasing = TRUE))[1],
  Most_Functionally_Diverse_lncRNA = functional_summary$lncRNA[1],
  Functional_Diversity_Score = functional_summary$Functional_Diversity_Score[1],
  stringsAsFactors = FALSE
)

write.csv(master_summary, 
          "results/guilt_by_association/tables/master_analysis_summary.csv",
          row.names = FALSE)

cat("✓ Master analysis summary saved\n")



# ============================================================================
# STEP 7C: CREATE TOP CANDIDATES TABLE
# ============================================================================

cat("\nSTEP 7C: Creating top candidates table...\n")

# Select top 15 lncRNAs based on functional diversity
top_candidates <- head(functional_summary, 15)

# Add additional metrics
top_candidates$Confidence_Score <- round(
  (top_candidates$N_GO_BP_Terms / max(top_candidates$N_GO_BP_Terms)) *
    (top_candidates$N_KEGG_Pathways / max(top_candidates$N_KEGG_Pathways)) *
    (top_candidates$Mean_Correlation / max(top_candidates$Mean_Correlation)) * 100,
  1
)

# Reorder columns for better presentation
top_candidates_final <- top_candidates[, c(
  "lncRNA", "Inferred_Function", "Confidence_Score", 
  "N_GO_BP_Terms", "N_KEGG_Pathways", "Mean_Correlation",
  "Top_GO_BP_Terms", "Top_KEGG_Pathways", "Top_Correlated_Partners"
)]

write.csv(top_candidates_final, 
          "results/guilt_by_association/tables/top_lncRNA_candidates.csv",
          row.names = FALSE)

cat("✓ Top candidates table saved\n")




# Load your functional summary
functional_summary <- read.csv("results/guilt_by_association/tables/functional_annotation_summary.csv")

# Count functional categories
library(ggplot2)

function_counts <- table(functional_summary$Inferred_Function)
function_df <- data.frame(
  Category = names(function_counts),
  Count = as.numeric(function_counts),
  Percentage = round(as.numeric(function_counts) / sum(function_counts) * 100, 1)
)

# Create pie chart
p_pie <- ggplot(function_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Distribution of Inferred Biological Functions\nAmong 50 Functionally Enriched lncRNAs") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

ggsave("results/guilt_by_association/figures/Figure11_functional_distribution_pie.png",
       p_pie, width = 10, height = 8, dpi = 300)

cat("✓ Figure 11 (pie chart) created!\n")











# ============================================================================
# FINAL COMPREHENSIVE SUMMARY
# ============================================================================

cat("\n============================================================\n")
cat("GUILT-BY-ASSOCIATION ANALYSIS COMPLETELY FINISHED!\n")
cat("============================================================\n")
cat("FINAL RESULTS SUMMARY:\n")
cat("- Analyzed", master_summary$Total_lncRNAs_Analyzed, "lncRNAs total\n")
cat("- Functionally annotated", master_summary$lncRNAs_with_Functional_Enrichment, "high-quality lncRNAs\n")
cat("- Average:", master_summary$Average_GO_BP_Terms_per_lncRNA, "GO terms &", 
    master_summary$Average_KEGG_Pathways_per_lncRNA, "KEGG pathways per lncRNA\n")
cat("- Top functional category:", master_summary$Top_Functional_Category, "\n")
cat("- Most promising lncRNA:", master_summary$Most_Functionally_Diverse_lncRNA, "\n")
cat("- Overall confidence score:", master_summary$Functional_Diversity_Score, "\n")
cat("\nKEY OUTPUT FILES:\n")
cat("- functional_annotation_summary.csv: Detailed functional profiles\n")
cat("- master_analysis_summary.csv: Overall analysis metrics\n")  
cat("- top_lncRNA_candidates.csv: Prioritized lncRNAs for validation\n")
cat("- All results in: results/guilt_by_association/\n")
cat("============================================================\n")
cat("ANALYSIS COMPLETE AND READY FOR PUBLICATION! 🎉\n")
cat("============================================================\n")






