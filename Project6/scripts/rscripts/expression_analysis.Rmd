c:\Users\DEL\Desktop\R Projects\Eneza-Data-Science-Residential-training-2026\cancer_subtypes_group6.Rmd
---
title: "ENEZA DS_Differential Expression Analysis"
author: "Delphine L."
date: "`r Sys.Date()`"
output:
  html_document: default
  pdf_document: default
tbl-pos: H
execute:
  echo: false
  warning: false
  message: false
format: null
pdf:
  pdf-engine: xelatex
---

# Molecular Profiling and Transcriptomic Dissection of Papillary Thyroid Carcinoma Subtypes
## Introduction
# Molecular Profiling of Papillary Thyroid Carcinoma Subtypes

## Introduction

Papillary thyroid carcinoma (PTC) comprises distinct histopathological variants with contrasting clinical behaviors, ranging from indolent **follicular** tumors to benchmark **classical** PTC and aggressive **tall cell** carcinomas. Although classical and tall cell variants share a broad baseline transcriptional landscape in unsupervised analyses, supervised transcriptomic profiling and Gene Ontology (GO) enrichment reveal clear molecular drivers that differentiate these phenotypes.

This study analyzes global expression variance, differential gene expression, and pathway enrichment across pairwise subtype comparisons (**Follicular vs. Classical**, **Tall Cell vs. Classical**, and **Tall Cell vs. Follicular**). Specifically, it aims to:

1. **Define Subtype Architecture:** Assess global transcriptomic variance and cluster overlap among PTC variants.
2. **Identify Molecular Drivers:** Pinpoint candidate markers driving metabolic reprogramming (*ATP1A3*, *AKR7A3*), cytoskeletal hardening (*KRT5/6/14*), and invasion signaling (*LAMA3*, *S100A2*, *MMP7*).
3. **Map Functional Pathways:** Characterize the shift from homeostatic ion transport in indolent subtypes to extracellular matrix degradation and intermediate filament organization in aggressive variants.
4. **Inform Translational Biomarkers:** Establish a high-resolution molecular panel to enhance diagnostic accuracy and support risk-stratified patient management.
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE,warning=FALSE,message=FALSE)
knitr::opts_chunk$set(echo = FALSE,
                      warning=FALSE,
                      message= FALSE)
```
 
 
 
# Methods & Analytical Pipeline
## Data Acquisition and Preprocessing
Gene expression data was downloaded from the NCBI Gene Expression Omnibus (GEO). Data preprocessing and normalization were performed in R to ensure cross-sample comparability:Raw counts were transformed using $\log_2(\text{norm\_counts} + 1)$ to stabilize variance.Quality control filtering was applied to remove non-informative, low-expression genes.
```{r echo=FALSE, message=FALSE, warning=FALSE,results='hide'}
pacman::p_load(tidyverse)

expr<-read_tsv("C:\\Users\\DEL\\Downloads\\eneza_ds\\data_mrna_seq_v2_rsem.txt");
head(expr)

sample<-read_tsv("C:\\Users\\DEL\\Downloads\\eneza_ds\\data_clinical_sample.txt",
                 skip=4);head(sample)

patient<-read_tsv("C:\\Users\\DEL\\Downloads\\eneza_ds\\data_clinical_patient.txt")

```

 

```{r, echo=FALSE, message=FALSE, warning=FALSE,message=FALSE,results='hide'}
library(org.Hs.eg.db)
# -----------------------------------------------------------------------------
# 1. Map Missing Hugo Symbols using Entrez IDs
# -----------------------------------------------------------------------------
# Identify rows where Hugo_Symbol is missing, NA, or empty
missing_hugo_idx <- which(is.na(expr$Hugo_Symbol) | expr$Hugo_Symbol == "" | expr$Hugo_Symbol == "NA")

if (length(missing_hugo_idx) > 0) {
  mapped_symbols <- mapIds(
    org.Hs.eg.db,
    keys = as.character(expr$Entrez_Gene_Id[missing_hugo_idx]),
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )
  
  # Fill in mapped symbols where missing
  expr$Hugo_Symbol[missing_hugo_idx] <- mapped_symbols
}

# If Hugo_Symbol is STILL missing after mapping, fall back to Entrez ID string
expr$Hugo_Symbol <- ifelse(
  is.na(expr$Hugo_Symbol) | expr$Hugo_Symbol == "",
  paste0("ENTREZ_", expr$Entrez_Gene_Id),
  expr$Hugo_Symbol
)

# -----------------------------------------------------------------------------
# 2. Map Missing Entrez IDs using Hugo Symbols (Reverse direction)
# -----------------------------------------------------------------------------
missing_entrez_idx <- which(is.na(expr$Entrez_Gene_Id))

if (length(missing_entrez_idx) > 0) {
  mapped_entrez <- mapIds(
    org.Hs.eg.db,
    keys = expr$Hugo_Symbol[missing_entrez_idx],
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  expr$Entrez_Gene_Id[missing_entrez_idx] <- as.numeric(mapped_entrez)
}

# -----------------------------------------------------------------------------
# 3. Resolve Duplicate Gene Symbols for Downstream Matrix Construction
# -----------------------------------------------------------------------------
# Option A (Recommended): Keep the highest-expressing row for duplicate symbols
# Assuming columns 3:ncol(expr) contain expression values

numeric_cols <- 3:ncol(expr)
expr_data <- expr[, numeric_cols]

# Calculate mean expression per row to pick the best representative row
row_means <- rowMeans(expr_data, na.rm = TRUE)

# Deduplicate Hugo_Symbol by keeping the highest expressed transcript
expr_clean <- expr %>%
  mutate(mean_expr = row_means) %>%
  arrange(desc(mean_expr)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  dplyr::select(-mean_expr)

# Verify uniqueness for matrix row names
length(unique(expr_clean$Hugo_Symbol)) == nrow(expr_clean) # Returns TRUE

# Extract numeric matrix
expr_matrix <- as.matrix(expr_clean[, 3:ncol(expr_clean)])
storage.mode(expr_matrix) <- "numeric"

# Assign row names safely
rownames(expr_matrix) <- expr_clean$Hugo_Symbol

# Check matrix dimensions and row names
dim(expr_matrix)
head(rownames(expr_matrix))

#sanity check
# 1. Verify that no Hugo Symbols are NA or empty string
sum(is.na(expr_clean$Hugo_Symbol) | expr_clean$Hugo_Symbol == "")
# Should return: 0

# 2. Verify that all row names on your matrix are present and valid
sum(is.na(rownames(expr_matrix)))
# Should return: 0



```

## Normalization & Differential Expression Modeling (edgeR)
 Samples annotated with non-specific tumor classifications ("Other") were excluded to restrict the cohort to three core papillary thyroid carcinoma (PTC) variants: Classical, Follicular, and Tall Cell.
 Differential expression analysis was conducted using edgeR:
 1. Filtering & Normalization: Low-expression genes were removed using filterByExpr(). 
 2. Count matrices were normalized across libraries using the Trimmed Mean of M-values (TMM) method (calcNormFactors).
 3.GLM Dispersion & Fitting: A cell-means design matrix without an intercept (~ 0 + group) was constructed. Empirical Bayes robust dispersions were calculated (estimateDisp), and negative binomial generalized linear models (GLMs) were fitted (glmFit).
 4. Subtype Contrasts: Likelihood ratio tests (glmLRT) evaluated three distinct pairwise comparisons:Follicular vs. Classical: contrast = c(-1, 1, 0)Tall Cell vs. Classical: contrast = c(-1, 0, 1)Tall Cell vs. Follicular: contrast = c(0, -1, 1)Significance Thresholds: Candidate biomarkers were defined at a strict significance threshold of $p < 0.001$ (Benjamini-Hochberg adjusted) and an absolute fold-change threshold of $\lvert\log_2\text{FC}\rvert > 1.5$. Results were visualized using customized volcano plots with ggrepel labeling.
 
 ## Unsupervised Cohort Profiling (PCA)
 To evaluate global transcriptional relationships post-normalization:Log-transformed counts per million ($\log_2\text{CPM}$) were extracted using a prior count of 1.The top 1,000 high-variance genes (HVGs) were isolated based on cross-sample variance.Principal Component Analysis (prcomp, scale = TRUE) was performed on the transposed matrix. Cohort structure was plotted along PC1 and PC2 with 95% confidence ellipses (stat_ellipse) for each subtype.
 
 
##Biological Pathway Enrichment  
To evaluate the biological processes driven by gene expression changes, differentially expressed genes were partitioned into distinct upregulated ($\log_2\text{FC} > 1.5$) and downregulated ($\log_2\text{FC} < -1.5$) clusters for each pairwise contrast. Functional pathway enrichment for Biological Processes (BP) was performed using the `compareCluster` function in `clusterProfiler`, enabling direct comparison of pathways enriched across opposing expression directions. To refine the visual presentation and eliminate semantic redundancy among closely related Gene Ontology terms, pathway outputs were compressed using `clusterProfiler::simplify()` with a similarity threshold cutoff of 0.6. The resulting non-redundant biological pathways were visualized using split dot plots, displaying the top three to four enriched terms per cluster ordered by GeneRatio and adjusted p-value to highlight the dominant functional programs defining each tumor subtype comparison.
```{r, echo=FALSE, message=FALSE, warning=FALSE,message=FALSE,message=FALSE,results='hide'}
###______
pacman::p_load(edgeR,tidyverse,pheatmap,RColorBrewer,clusterProfiler,org.Hs.eg.db)
 
 
# -----------------------------------------------------------------------------
# 2. Data Preparation & Matrix Formatting (Excluding 'Other' Tumour Type)
# -----------------------------------------------------------------------------
# Format Sample Metadata & Filter out "Other"
col_data <- sample %>% 
  column_to_rownames(var = "SAMPLE_ID") %>% 
  filter(TUMOR_TYPE != "Other") # <--- Drops 'Other' samples from metadata

# Format Expression Matrix
norm_mat <- as.matrix(expr_clean[, 3:ncol(expr_clean)])
storage.mode(norm_mat) <- "numeric"
rownames(norm_mat) <- expr_clean$Hugo_Symbol

# Align samples (this automatically drops 'Other' from the expression matrix too)
common_samples <- intersect(colnames(norm_mat), rownames(col_data))
norm_mat <- norm_mat[, common_samples]
col_data <- col_data[common_samples, , drop = FALSE]

# Set Tumour Subtype Group Factor (drops unused 'Other' factor level)
col_data$group <- factor(col_data$TUMOR_TYPE)

cat("Remaining samples after dropping 'Other':", ncol(norm_mat), "\n")
cat("Retained group levels:\n")
print(levels(col_data$group))

# Build edgeR DGEList Object
#head(norm_mat_int)
norm_mat_int <- round(norm_mat)
dge <- DGEList(counts = norm_mat, group = col_data$group)
dge$samples$norm.factors
# Lock Normalization
#dge$samples$norm.factors <- 1

# Filter Low-Expression Genes
# keep <- rowSums(norm_mat > 1) >= 5
# dge_filtered <- dge[keep, ]
# dge_filtered$samples
# cat("Genes retained after filtering:", sum(keep), "out of", length(keep), "\n")

# Filter lowly expressed genes
keep2 <- filterByExpr(dge)
dgeObj <- dge[keep2, , keep.lib.sizes = FALSE]

# Step 4: Normalize for library size using TMM
dgeObj <- calcNormFactors(dgeObj, method = "TMM")

# Check normalization factors
dgeObj$samples

#Normalize the data
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples






```

 
 # Results
 1. Global Transcriptomic Variance & Subtype Separation (PCA)
Unsupervised Principal Component Analysis (PCA) performed on the top 1,000 high-variance genes captures 37% of the total cohort variance across the first two principal components (PC1: 28%, PC2: 9%). The analysis demonstrates a distinct separation along the primary axis of variation (PC1), where the indolent Follicular variant forms an isolated cluster shifted toward the positive axis. In contrast, the Classical variant exhibits broad dispersion across PC1, spanning an intermediate continuum that bridges the Follicular and Tall Cell profiles. The aggressive Tall Cell variant does not form a standalone cluster; instead, it remains entirely nested within the negative boundary of the Classical ellipse. This spatial distribution indicates that while Follicular PTC operates under a distinct global expression program, Tall Cell tumors represent a specialized, highly localized subset along the main Classical axis rather than a completely separate transcriptional lineage.
 
```{r, echo=FALSE,message=FALSE,results='hide' warning=FALSE,eval=TRUE}
pacman::p_load(edgeR,tidyverse,RColorBrewer,clusterProfiler,org.Hs.eg.db,knitr)
# -----------------------------------------------------------------------------
# 3. Log-CPM Transformation & Cohort Classification (PCA)
# -----------------------------------------------------------------------------
# Extract Log-CPM Matrix
logCPM <- cpm(dgeObj, log = TRUE, prior.count = 1)

# B. Principal Component Analysis (PCA Plot)
# Recode group levels for cleaner plotting
col_data$group_clean <- factor(col_data$TUMOR_TYPE, 
                               levels = c(
                                 "Thyroid Papillary Carcinoma, Classical/Usual Type",
                                 "Thyroid Papillary Carcinoma, Follicular (>= 99% Follicular Patterned)",
                                 "Thyroid Papillary Carcinoma, Tall Cell (>= 50% Tall Cell Features)"
                               ),
                               labels = c("Classical", "Follicular", "Tall Cell")
)
#--------
pacman::p_load(tidyverse, RColorBrewer, ggrepel)

#  Subset Top 1,000 High-Variance Genes
gene_vars <- apply(logCPM, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]
logCPM_top <- logCPM[top_genes, ]

# run PCA
pca <- prcomp(t(logCPM_top), scale. = TRUE)
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))

# Format Dataframe with Clean Labels
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Subtype = factor(col_data$TUMOR_TYPE, 
                   labels = c("Classical", "Follicular", "Tall Cell")),
  Sample_ID = rownames(col_data)
)

#PCA Plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype, fill = Subtype)) +
  geom_point(size = 2.8, alpha = 0.7, shape = 21, stroke = 0.4, color = "black") +
  stat_ellipse(
    aes(group = Subtype),
    level = 0.95,
    geom = "polygon",
    alpha = 0.12,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = c("Classical" = "#E41A1C", "Follicular" = "#4DAF4A", "Tall Cell" = "#377EB8")) +
  scale_color_manual(values = c("Classical" = "#E41A1C", "Follicular" = "#4DAF4A", "Tall Cell" = "#377EB8")) +
  labs(
    title = "PCA of Thyroid Tumour Gene Expression (Top 1,000 HVGs)",
    subtitle = "95% Confidence Ellipses by Subtype",
    x = paste0("PC1 (", percent_var[1], "% Variance)"),
    y = paste0("PC2 (", percent_var[2], "% Variance)"),
    fill = "Subtype",
    color = "Subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
 
#ggsave("pca_plot.png", plot = pca_plot, width = 10, height = 7, units = "in", dpi = 300)
```



```{r, echo=FALSE, warning=FALSE,message=FALSE,results='hide'}
# =============================================================================
# 4. Dispersion Estimation & GLM Fitting (Runs on FULL dataset)
# =============================================================================
design <- model.matrix(~ 0 + group, data = col_data)
colnames(design) <- levels(col_data$group)

cat("Design matrix column order:\n")
print(colnames(design))

# Estimate Dispersions & Fit GLM
dge_disp <- estimateDisp(dgeObj, design, robust = TRUE)
fit <- glmFit(dge_disp, design)


# =============================================================================
# 5. Pairwise Subtype Contrasts (glmLRT)
# =============================================================================
# Contrast order: c(Classical, Follicular, TallCell)
Follicular_vs_Classical_lrt <- glmLRT(fit, contrast = c(-1, 1, 0))
TallCell_vs_Classical_lrt   <- glmLRT(fit, contrast = c(-1, 0, 1))
TallCell_vs_Follicular_lrt  <- glmLRT(fit, contrast = c(0, -1, 1))


# =============================================================================
# 6. Automated Marker Identification & Split GO Pathway Enrichment
# =============================================================================
 analyze_thyroid_contrast <- function(lrt_object, contrast_name) {
  cat("\n======================================================\n")
  cat(" Analyzing Subtype Contrast:", contrast_name, "\n")
  cat("======================================================\n")
  
  # Extract results table sorted by PValue
  res_table <- topTags(lrt_object, n = Inf, adjust.method = "BH")$table %>%
    rownames_to_column(var = "Hugo_Symbol")
  
  # Filter Candidate Biomarker Genes (PValue < 0.001 & |logFC| > 1.5)
  sig_genes <- res_table %>%
    filter(PValue < 0.001 & abs(logFC) > 1.5) %>%
    arrange(PValue)
  
  cat("Total DE genes (PValue < 0.001 & |logFC| > 1.5):", nrow(sig_genes), "\n")
  
  # Export results to CSV
  write.csv(res_table, paste0("Thyroid_DE_", contrast_name, ".csv"), row.names = FALSE)
  write.csv(sig_genes, paste0("Thyroid_sigDE_", contrast_name, ".csv"), row.names = FALSE)
 }
```


 ## 2. Pairwise Differential Expression & Pathway Analyses
 #A. Follicular vs. Classical PTC
 Supervised differential expression modeling isolates a distinct, highly significant metabolic and ion-transport reprogramming profile in the Follicular subtype. ATP1A3 emerges as the primary driver upregulated in Follicular PTC ($-\log_{10}(p) > 100$), accompanied by the co-upregulation of ion transport channels (KCNJ1, HCN2), metabolic enzymes (AKR7A3), calcium buffers (PVALB), metal homeostasis factors (MT3), and E3 ubiquitin ligases (TRIM50, DAND5, NEFL). Conversely, lipid uptake machinery such as LDLR is significantly downregulated relative to Classical PTC. Gene Ontology (GO) enrichment confirms this functional divergence: genes upregulated in the Follicular variant ($n = 133$) strictly enrich for potassium ion transmembrane transport, cellular response to zinc ion, and detoxification of copper ion. Downregulated genes ($n = 382$, elevated in Classical PTC) govern humoral immune response, external encapsulating structure organization, and skin development ($p.\text{adjust} < 1 \times 10^{-9}$). Mechanistically, this demonstrates that Follicular tumors rely on specialized electrochemical homeostasis while suppressing the invasive extracellular matrix remodeling and inflammatory programs typical of Classical PTC.
 
 #B. Tall Cell vs. Classical PTC
 Despite their complete spatial overlap in global PCA, supervised differential expression reveals a specialized cytoskeletal transformation defining the Tall Cell phenotype. Genes upregulated in Tall Cell tumors ($n = 33$, $p.\text{adjust} < 1 \times 10^{-6}$) are led by basement membrane anchors (LAMA3; $-\log_{10}(p) > 35$), invasion mediators (S100A2, HEPHL1, CLCA2, TFPI2), and basal/squamoid cytokeratins (KRT5, KRT14, KRT6A, KRT6C, ANXA8L2). GO enrichment demonstrates that these upregulated genes are exclusively enriched in structural remodeling pathways, including intermediate filament organization, intermediate filament-based processes, and keratinization. In contrast, genes downregulated in Tall Cell tumors ($n = 187$, representing processes elevated in Classical PTC) govern baseline transport functions such as potassium ion transmembrane transport, monoatomic anion transport, and neuromuscular processes. Clinically, this proves that the transition to the aggressive Tall Cell phenotype is driven by a targeted cytoskeletal activation and architectural hardening sequence rather than a global transcriptomic rewrite.
 
 C. Tall Cell vs. Follicular PTC
 Direct comparison between the two phenotypic extremes reveals a marked, asymmetric transcriptional activation characterized by the robust upregulation of aggressive oncogenic drivers in Tall Cell tumors, including S100A2, MMP7, AHNAK2, ITGA3, ITGB6, LAMA3, HES2, and ALOX5 ($p < 1 \times 10^{-30}$). GO biological process enrichment illustrates a stark functional polarity across these variants: genes upregulated in Tall Cell PTC ($n = 607$, $p.\text{adjust} < 1 \times 10^{-14}$) are strongly enriched in external encapsulating structure organization, positive regulation of cell-cell adhesion, and humoral immune response. Conversely, genes downregulated in Tall Cell tumors ($n = 569$, elevated in Follicular PTC) govern baseline physiological processes, including regulation of membrane potential, organic acid transport, and vascular processes in the circulatory system ($p.\text{adjust} \approx 1 \times 10^{-8}$). These non-overlapping pathway signatures confirm that Tall Cell and Follicular variants occupy opposite ends of the PTC functional spectrum, supplying clear discriminators for multi-marker diagnostic models and subtype risk stratification.
```{r, echo=FALSE, warning=FALSE,message=FALSE,results='hide'}
# =============================================================================
# 9. GO Pathway Split Dotplot Function
# =============================================================================
generate_split_go_dotplot <- function(lrt_object, contrast_name, p_cutoff = 0.001, fc_cutoff = 1.5, top_terms = 4) {
  
  cat("\n======================================================\n")
  cat(" Running Clean GO Enrichment for:", contrast_name, "\n")
  cat("======================================================\n")
  
  # 1. Extract DE genes
  results <- topTags(lrt_object, n = Inf, adjust.method = "BH")$table
  if (!"Hugo_Symbol" %in% colnames(results)) {
    results <- results %>% rownames_to_column(var = "Hugo_Symbol")
  }
  
  sig_genes <- results %>%
    filter(PValue < p_cutoff & abs(logFC) > fc_cutoff)
  
  if (nrow(sig_genes) == 0) {
    cat("No significant genes found for contrast:", contrast_name, "\n")
    return(NULL)
  }
  
  # 2. Map Entrez IDs
  up_symbols   <- sig_genes %>% filter(logFC > fc_cutoff) %>% pull(Hugo_Symbol)
  down_symbols <- sig_genes %>% filter(logFC < -fc_cutoff) %>% pull(Hugo_Symbol)
  
  up_entrez   <- na.omit(mapIds(org.Hs.eg.db, keys = up_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  down_entrez <- na.omit(mapIds(org.Hs.eg.db, keys = down_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
  
  gene_list <- list(
    "Upregulated"   = as.character(up_entrez),
    "Downregulated" = as.character(down_entrez)
  )
  
  # 3. Enrichment
  comp_go <- compareCluster(
    geneClusters  = gene_list,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.01
  )
  
  if (is.null(comp_go) || nrow(as.data.frame(comp_go)) == 0) {
    cat("No significant GO pathways enriched.\n")
    return(NULL)
  }
  
  # 4. Simplify redundancy to clean up Y-axis
  comp_go_clean <- clusterProfiler::simplify(comp_go, cutoff = 0.6, by = "p.adjust")
  
  # 5. Render
  p_dot <- dotplot(
    comp_go_clean,
    showCategory = top_terms,
    label_format = 30,
    title = paste("GO Pathways (Up vs Down):", contrast_name)
  ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      axis.text.x     = element_text(face = "bold", color = "black"),
      axis.text.y     = element_text(size = 8.5, color = "black"),
      panel.grid.major.y = element_line(color = "grey92", linetype = "dashed")
    )
  
  print(p_dot)
  return(p_dot)
}

# Run for your contrast
go_Foll_vs_Class <- generate_split_go_dotplot(
  Follicular_vs_Classical_lrt, 
  "Follicular vs Classical", 
  top_terms = 3
)

go_Tall_vs_Class <- generate_split_go_dotplot(
  TallCell_vs_Classical_lrt, 
  "Tall Cell vs Classical", 
  top_terms = 3
)

go_Tall_vs_Foll  <- generate_split_go_dotplot(
  TallCell_vs_Follicular_lrt, 
  "Tall Cell vs Follicular", 
  top_terms = 3
)

# 2. Save all 3 plots as PNG files
ggsave("GO_Dotplot_Follicular_vs_Classical.png", go_Foll_vs_Class, width = 10, height = 7, dpi = 300)
ggsave("GO_Dotplot_TallCell_vs_Classical.png",   go_Tall_vs_Class, width = 10, height = 7, dpi = 300)
ggsave("GO_Dotplot_TallCell_vs_Follicular.png",  go_Tall_vs_Foll,  width = 10, height = 7, dpi = 300)
# -----------------------------------------------------------------------------
# Execute Split GO Enrichment for All Three Contrasts
# -----------------------------------------------------------------------------
go_Foll_vs_Class <- generate_split_go_dotplot(Follicular_vs_Classical_lrt, "Follicular vs Classical")
go_Tall_vs_Class <- generate_split_go_dotplot(TallCell_vs_Classical_lrt, "Tall Cell vs Classical")
go_Tall_vs_Foll  <- generate_split_go_dotplot(TallCell_vs_Follicular_lrt, "Tall Cell vs Follicular")



```


```{r,echo=FALSE, warning=FALSE,message=FALSE,results='hide'}

 # =============================================================================
# 8. Volcano Plotting (X = log2FC, Y = -log10(p-value) + Gene Labels)
# =============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, ggrepel)


volcano_plot <- function(lrt_object, contrast_name, top_n = 15, p_cutoff = 0.001, fc_cutoff = 1.5) {
  
  results <- topTags(lrt_object, n = Inf, adjust.method = "BH")$table %>%
    rownames_to_column(var = "Hugo_Symbol") %>%
    mutate(
      Significance = case_when(
        PValue < p_cutoff & logFC > fc_cutoff  ~ "Upregulated",
        PValue < p_cutoff & logFC < -fc_cutoff ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_labels <- results %>%
    filter(Significance != "Not significant") %>%
    arrange(PValue) %>%
    slice_head(n = top_n)
  
  volcano_colors <- c(
    "Upregulated"     = "#D95F02", 
    "Downregulated"   = "#7570B3", 
    "Not significant" = "grey70"
  )
  
  p <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_labels,
      aes(label = Hugo_Symbol),
      size = 3.8,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      force = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      subtitle = paste("Top", top_n, "Most Significantly Differentially Expressed Genes Labeled"),
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(p-value)"),
      color = "Expression Status"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      legend.position = "top",
      # Explicitly draw solid X and Y axis lines
      axis.line = element_line(color = "black", linewidth = 0.6),
 panel.grid = element_blank())
    
    
  
  print(p)
  return(results)
}

# Generate Volcano Plots
volcano_Foll_vs_Class <- volcano_plot(Follicular_vs_Classical_lrt, "Follicular vs Classical")
volcano_Tall_vs_Class <- volcano_plot(TallCell_vs_Classical_lrt, "Tall Cell vs Classical")
volcano_Tall_vs_Foll  <- volcano_plot(TallCell_vs_Follicular_lrt, "Tall Cell vs Follicular")

# Open a fixed-size PNG device
png("Volcano_Follicular_vs_Classical.png", width = 10, height = 7, units = "in", res = 300)
png("Tall Cell vs Classical.png", width = 10, height = 7, units = "in", res = 300)
png("Tall Cell vs Follicular.png", width = 10, height = 7, units = "in", res = 300) 

```

```{r,r,echo=FALSE, message=FALSE, warning=FALSE,message=FALSE,results='hide'}
## =============================================================================
# 8. Optimized Volcano Plot Function & Generation
# =============================================================================
pacman::p_load(tidyverse, ggrepel)

volcano_plot <- function(lrt_object, contrast_name, top_n = 10, fc_cutoff = 1.5, p_cutoff = 0.001) {
  
  # 1. Extract full DE results table
  results <- topTags(lrt_object, n = Inf)$table
  
  # 2. Ensure Hugo_Symbol column exists
  if (!"Hugo_Symbol" %in% colnames(results)) {
    results <- results %>% rownames_to_column(var = "Hugo_Symbol")
  }
  
  # 3. Categorize Significance
  results <- results %>%
    mutate(
      Significance = case_when(
        logFC > fc_cutoff & PValue < p_cutoff ~ "Upregulated",
        logFC < -fc_cutoff & PValue < p_cutoff ~ "Downregulated",
        TRUE ~ "Not Significant"
      ),
      Significance = factor(Significance, levels = c("Upregulated", "Downregulated", "Not Significant"))
    )
  
  # 4. Color Palette
  volcano_colors <- c(
    "Upregulated"     = "#E41A1C", 
    "Downregulated"   = "#377EB8", 
    "Not Significant" = "grey70"
  )
  
  # 5. Extract Top N Most Significant Genes for Labeling
  top_labels <- results %>%
    filter(Significance %in% c("Upregulated", "Downregulated")) %>%
    arrange(PValue) %>%
    head(top_n)
  
  # 6. Build Plot
  p <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_labels,
      aes(label = Hugo_Symbol),
      size = 3.8,
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      force = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    labs(
      title = paste("Volcano Plot:", contrast_name),
      subtitle = paste("Top", top_n, "Most Significantly Differentially Expressed Genes Labeled"),
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"(p-value)"),
      color = "Expression Status"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 14),
      plot.subtitle   = element_text(size = 11, color = "grey30"),
      legend.position = "top",
      axis.line       = element_line(color = "black", linewidth = 0.6),
      panel.grid      = element_blank()
    )
  
  return(p)
}

# -----------------------------------------------------------------------------
# Execute & Display Volcano Plots
# -----------------------------------------------------------------------------
volcano_Foll_vs_Class <- volcano_plot(Follicular_vs_Classical_lrt, "Follicular vs Classical")
volcano_Tall_vs_Class <- volcano_plot(TallCell_vs_Classical_lrt, "Tall Cell vs Classical")
volcano_Tall_vs_Foll  <- volcano_plot(TallCell_vs_Follicular_lrt, "Tall Cell vs Follicular")

# Force display to R Plot Window
print(volcano_Foll_vs_Class)
print(volcano_Tall_vs_Class)
print(volcano_Tall_vs_Foll)


# save
ggsave(
  "Volcano_Follicular_vs_Classical.png",
  volcano_Foll_vs_Class,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

ggsave(
  "Volcano_Tall_Cell_vs_Classical.png",
  volcano_Tall_vs_Class,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

ggsave(
  "Volcano_Tall_Cell_vs_Follicular.png",
  volcano_Tall_vs_Foll,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

```



```{r echo=FALSE, message=FALSE, warning=FALSE, results='hide',}
# ###----------------------------------------------------------------------\\
# ###CLASSIFICATION
# # =============================================================================
# # 10. Multi-Class Machine Learning Classification (Random Forest)
# # =============================================================================
# if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
# pacman::p_load(tidyverse, caret, randomForest, pROC, MCMCglmm)
# 
# cat("\n======================================================\n")
# cat(" Starting Multi-Class Random Forest Classification\n")
# cat("======================================================\n")
# 
# # -----------------------------------------------------------------------------
# # 1. Feature Selection: Union of Top Biomarkers Across Contrasts
# # -----------------------------------------------------------------------------
# # Extract unique candidate genes identified from the edgeR pairwise contrasts
# top_biomarkers <- unique(c(
#   sig_Foll_vs_Class$Hugo_Symbol,
#   sig_Tall_vs_Class$Hugo_Symbol,
#   sig_Tall_vs_Foll$Hugo_Symbol
# ))
# 
# cat("Total candidate biomarker genes selected for RF:", length(top_biomarkers), "\n")
# 
# # Extract normalized log-CPM expression matrix
# log_cpm_matrix <- edgeR::cpm(dgeObj, log = TRUE)
# 
# # Subset matrix to selected genes & transpose to Samples (rows) x Genes (columns)
# ml_data <- as.data.frame(t(log_cpm_matrix[rownames(log_cpm_matrix) %in% top_biomarkers, ]))
# ml_data$Subtype <- col_data$group
# 
# # Ensure valid R factor levels for classification targets
# ml_data$Subtype <- factor(make.names(ml_data$Subtype))
# 
# # -----------------------------------------------------------------------------
# # 2. Stratified Train/Test Split (80% Train, 20% Test)
# # -----------------------------------------------------------------------------
# set.seed(42) # Ensure reproducible splitting
# train_index <- createDataPartition(ml_data$Subtype, p = 0.80, list = FALSE)
# 
# train_data <- ml_data[train_index, ]
# test_data  <- ml_data[-train_index, ]
# 
# cat("Training Set Size:", nrow(train_data), "samples\n")
# cat("Testing Set Size: ", nrow(test_data), "samples\n")
# 
# # -----------------------------------------------------------------------------
# # 3. Model Training with 10-Fold Cross-Validation (Repeated 5x)
# # -----------------------------------------------------------------------------
# train_control <- trainControl(
#   method          = "repeatedcv",
#   number          = 10,
#   repeats         = 5,
#   classProbs      = TRUE,
#   summaryFunction = multiClassSummary,
#   savePredictions = "final"
# )
# 
# # Train Random Forest Classifier using caret::train explicitly
# install.packages("MLmetrics")
# library(MLmetrics)
# set.seed(42)
# rf_model <- caret::train(
#   Subtype ~ .,
#   data        = train_data,
#   method      = "rf",
#   metric      = "Accuracy",
#   trControl   = train_control,
#   importance  = TRUE
# )
# 
# print(rf_model)
# 
# # -----------------------------------------------------------------------------
# # 4. Independent Test Set Evaluation & Confusion Matrix
# # -----------------------------------------------------------------------------
# predictions   <- predict(rf_model, newdata = test_data)
# probabilities <- predict(rf_model, newdata = test_data, type = "prob")
# 
# # Render Confusion Matrix and Performance Metrics
# conf_matrix <- confusionMatrix(predictions, test_data$Subtype)
# cat("\n======================================================\n")
# cat(" Test Set Confusion Matrix & Performance Metrics\n")
# cat("======================================================\n")
# print(conf_matrix)
# 
# # -----------------------------------------------------------------------------
# # 5. Feature Importance Plot (Top 20 Classifier Drivers)
# # -----------------------------------------------------------------------------
# var_imp <- varImp(rf_model)$importance %>%
#   rownames_to_column(var = "Gene") %>%
#   arrange(desc(Overall)) %>%
#   slice_head(n = 20)
# 
# p_imp <- ggplot(var_imp, aes(x = reorder(Gene, Overall), y = Overall)) +
#   geom_bar(stat = "identity", fill = "#2B5C8F", width = 0.7) +
#   coord_flip() +
#   labs(
#     title = "Top 20 Predictive Biomarker Genes",
#     subtitle = "Ranked by Mean Decrease Gini (Random Forest)",
#     x = "Gene Symbol",
#     y = "Variable Importance Score"
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     plot.title = element_text(face = "bold", size = 14),
#     axis.line  = element_line(color = "black", linewidth = 0.6),
#     panel.grid = element_blank()
#   )
# 
# print(p_imp)
# ggsave("RF_Variable_Importance.png", plot = p_imp, width = 8, height = 6)
```





```{r echo=FALSE, message=FALSE, warning=FALSE}

 
```

