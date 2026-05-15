# =========================================================
# DIET PROJECT — DESEQ2 MAIN EFFECT ANALYSIS
# Genotype + Diet + Timepoint
# =========================================================

BiocManager::install("DESeq2")

install.packages("pheatmap")

# -------------------------
# LOAD LIBRARIES
# -------------------------
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# -------------------------
# SET WORKING DIRECTORY
# -------------------------
setwd("~/GitHub/Caloric-Density-and-Recombination")

# -------------------------
# IMPORT COUNT MATRIX
# -------------------------
expr <- read.csv(
  "Raw Data/RNA/gene_count.csv",
  row.names = 1,
  check.names = FALSE
)

# -------------------------
# KEEP ONLY SAMPLE COLUMNS
# (first 24 columns are samples)
# -------------------------
count_matrix <- expr[, 1:24]

# Convert to matrix
count_matrix <- as.matrix(count_matrix)

# Ensure integers
storage.mode(count_matrix) <- "integer"

# -------------------------
# CREATE SAMPLE METADATA
# -------------------------
sample_info <- data.frame(
  sample = colnames(count_matrix)
)

# -------------------------
# GENOTYPE
# -------------------------
sample_info$genotype <- ifelse(
  grepl("^M42", sample_info$sample),
  "S42",
  "S399"
)

# -------------------------
# DIET
# -------------------------
sample_info$diet <- case_when(
  grepl("05x", sample_info$sample) ~ "D05",
  grepl("1x", sample_info$sample)  ~ "D1",
  grepl("2x", sample_info$sample)  ~ "D2"
)

# -------------------------
# TIMEPOINT
# -------------------------
sample_info$timepoint <- case_when(
  grepl("_2_", sample_info$sample) ~ "Day2",
  grepl("_6_", sample_info$sample) ~ "Day6"
)

# Set rownames
rownames(sample_info) <- sample_info$sample

# -------------------------
# CONVERT TO FACTORS
# -------------------------
sample_info$genotype <- factor(sample_info$genotype)
sample_info$diet <- factor(sample_info$diet)
sample_info$timepoint <- factor(sample_info$timepoint)

# -------------------------
# SET REFERENCE LEVELS
# -------------------------
sample_info$genotype <- relevel(
  sample_info$genotype,
  ref = "S42"
)

sample_info$diet <- relevel(
  sample_info$diet,
  ref = "D1"
)

sample_info$timepoint <- relevel(
  sample_info$timepoint,
  ref = "Day2"
)

# -------------------------
# CHECK SAMPLE TABLE
# -------------------------
print(sample_info)

# =========================================================
# CREATE DESEQ2 OBJECT
# =========================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
 # design = ~ timepoint + genotype + diet
  design = ~ timepoint + genotype + diet + genotype:diet
)

# -------------------------
# FILTER LOW COUNT GENES
# -------------------------
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

cat("Genes retained after filtering:", nrow(dds), "\n")

# =========================================================
# RUN DESEQ2
# =========================================================

dds <- DESeq(dds)

# -------------------------
# SHOW MODEL COEFFICIENTS
# -------------------------
resultsNames(dds)

# =========================================================
# MAIN GENOTYPE EFFECT
# S399 vs S42
# =========================================================

res_genotype <- results(
  dds,
  name = "genotype_S399_vs_S42"
)

# =========================================================
# MAIN DIET EFFECTS
# =========================================================

# -------------------------
# D05 vs D1
# -------------------------
res_D05 <- results(
  dds,
  name = "diet_D05_vs_D1"
)

# -------------------------
# D2 vs D1
# -------------------------
res_D2 <- results(
  dds,
  name = "diet_D2_vs_D1"
)

# =========================================================
# ADD GENE ANNOTATIONS
# =========================================================

annotations <- expr[, c(
  "gene_name",
  "gene_chr",
  "gene_start",
  "gene_end",
  "gene_biotype",
  "gene_description"
)]

# =========================================================
# CONVERT RESULTS TO TABLES
# =========================================================

# -------------------------
# GENOTYPE
# -------------------------
genotype_table <- as.data.frame(res_genotype) %>%
  rownames_to_column("gene_id") %>%
  left_join(
    annotations %>%
      rownames_to_column("gene_id"),
    by = "gene_id"
  )

# -------------------------
# D05
# -------------------------
diet_D05_table <- as.data.frame(res_D05) %>%
  rownames_to_column("gene_id") %>%
  left_join(
    annotations %>%
      rownames_to_column("gene_id"),
    by = "gene_id"
  )

# -------------------------
# D2
# -------------------------
diet_D2_table <- as.data.frame(res_D2) %>%
  rownames_to_column("gene_id") %>%
  left_join(
    annotations %>%
      rownames_to_column("gene_id"),
    by = "gene_id"
  )

# =========================================================
# FILTER SIGNIFICANT DEGS
# =========================================================

padj_cutoff <- 0.05
log2fc_cutoff <- 1

# -------------------------
# GENOTYPE SIGNIFICANT
# -------------------------
genotype_sig <- genotype_table %>%
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) > log2fc_cutoff
  )

# -------------------------
# D05 SIGNIFICANT
# -------------------------
diet_D05_sig <- diet_D05_table %>%
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) > log2fc_cutoff
  )

# -------------------------
# D2 SIGNIFICANT
# -------------------------
diet_D2_sig <- diet_D2_table %>%
  filter(
    !is.na(padj),
    padj < padj_cutoff,
    abs(log2FoldChange) > log2fc_cutoff
  )

# =========================================================
# SUMMARY STATISTICS
# =========================================================

cat("\n==============================\n")
cat("MAIN EFFECT RESULTS SUMMARY\n")
cat("==============================\n\n")

cat("Genotype significant DEGs:", nrow(genotype_sig), "\n")
cat("Diet D05 vs D1 significant DEGs:", nrow(diet_D05_sig), "\n")
cat("Diet D2 vs D1 significant DEGs:", nrow(diet_D2_sig), "\n")

# =========================================================
# EXPORT RESULTS
# =========================================================


# -------------------------
# FULL RESULTS
# -------------------------
write.csv(
  genotype_table,
  "Summaries and Outputs/RNA/genotype_main_effect.csv",
  row.names = FALSE
)

write.csv(
  diet_D05_table,
  "Summaries and Outputs/RNA/diet_D05_vs_D1.csv",
  row.names = FALSE
)

write.csv(
  diet_D2_table,
  "Summaries and Outputs/RNA/diet_D2_vs_D1.csv",
  row.names = FALSE
)

# -------------------------
# SIGNIFICANT ONLY
# -------------------------
write.csv(
  genotype_sig,
  "Summaries and Outputs/RNA/genotype_main_effect_significant.csv",
  row.names = FALSE
)

write.csv(
  diet_D05_sig,
  "Summaries and Outputs/RNA/diet_D05_vs_D1_significant.csv",
  row.names = FALSE
)

write.csv(
  diet_D2_sig,
  "Summaries and Outputs/RNA/diet_D2_vs_D1_significant.csv",
  row.names = FALSE
)

# =========================================================
# TOP GENES
# =========================================================

cat("\n==============================\n")
cat("TOP 10 GENOTYPE DEGs\n")
cat("==============================\n\n")

print(
  genotype_sig %>%
    arrange(padj) %>%
    select(
      gene_id,
      gene_name,
      log2FoldChange,
      padj
    ) %>%
    head(10)
)

cat("\n==============================\n")
cat("TOP 10 D05 vs D1 DEGs\n")
cat("==============================\n\n")

print(
  diet_D05_sig %>%
    arrange(padj) %>%
    select(
      gene_id,
      gene_name,
      log2FoldChange,
      padj
    ) %>%
    head(10)
)

cat("\n==============================\n")
cat("TOP 10 D2 vs D1 DEGs\n")
cat("==============================\n\n")

print(
  diet_D2_sig %>%
    arrange(padj) %>%
    select(
      gene_id,
      gene_name,
      log2FoldChange,
      padj
    ) %>%
    head(10)
)

cat("\n✓ MAIN EFFECT ANALYSIS COMPLETED SUCCESSFULLY\n")


res_interaction_D05 <- results(
  dds,
  name="genotypeS399.dietD05"
)

res_interaction_D2 <- results(
  dds,
  name="genotypeS399.dietD2"
)
