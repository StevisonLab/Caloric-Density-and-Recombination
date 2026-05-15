# =========================================================
# DIET PROJECT — FULL DESEQ2 PIPELINE (MAIN + INTERACTION)
# Genotype + Diet + Timepoint
# =========================================================

# -------------------------
# INSTALL (run once only)
# -------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#BiocManager::install("DESeq2")
#install.packages(c("pheatmap", "ggplot2", "dplyr", "tibble"))

#BiocManager::install("apeglm")

# -------------------------
# LOAD LIBRARIES
# -------------------------
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(apeglm)

# =========================================================
# 1. LOAD DATA
# =========================================================

setwd("~/GitHub/Caloric-Density-and-Recombination")

expr <- read.csv(
  "Raw Data/RNA/gene_count.csv",
  row.names = 1,
  check.names = FALSE
)

count_matrix <- as.matrix(expr[, 1:24])
storage.mode(count_matrix) <- "integer"

# =========================================================
# 2. SAMPLE METADATA
# =========================================================

sample_info <- data.frame(sample = colnames(count_matrix))
rownames(sample_info) <- sample_info$sample

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

# FACTORS
sample_info$genotype <- factor(sample_info$genotype)
sample_info$diet <- factor(sample_info$diet)
sample_info$timepoint <- factor(sample_info$timepoint)

# REFERENCE LEVELS
sample_info$genotype <- relevel(sample_info$genotype, ref = "S42")
sample_info$diet <- relevel(sample_info$diet, ref = "D1")
sample_info$timepoint <- relevel(sample_info$timepoint, ref = "Day2")

print(sample_info)

# =========================================================
# 3. DESEQ2 OBJECT
# =========================================================

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ timepoint + genotype + diet + genotype:diet
)

# FILTER LOW EXPRESSION GENES
dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("Genes retained:", nrow(dds), "\n")

# =========================================================
# 4. RUN MODEL
# =========================================================

dds <- DESeq(dds)

resultsNames(dds)

# =========================================================
# 5. MAIN EFFECTS
# =========================================================

res_genotype <- results(dds, name = "genotype_S399_vs_S42")

res_D05 <- results(dds, contrast = c("diet", "D05", "D1"))
res_D2  <- results(dds, contrast = c("diet", "D2", "D1"))

# =========================================================
# 6. INTERACTION TERMS (IMPORTANT FIX)
# =========================================================
# DO NOT GUESS NAMES — always use resultsNames(dds)

res_interaction_D05 <- results(dds, name = "genotypeS399.dietD05")
res_interaction_D2  <- results(dds, name = "genotypeS399.dietD2")

# SHRINK INTERACTIONS (paper-quality)
res_interaction_D05_shrunk <- lfcShrink(
  dds,
  coef = "genotypeS399.dietD05",
  type = "apeglm"
)

res_interaction_D2_shrunk <- lfcShrink(
  dds,
  coef = "genotypeS399.dietD2",
  type = "apeglm"
)

# =========================================================
# 7. ANNOTATIONS
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
# 8. TABLES FUNCTION
# =========================================================

make_table <- function(res, annotations) {
  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    left_join(
      annotations %>% rownames_to_column("gene_id"),
      by = "gene_id"
    )
}

genotype_table <- make_table(res_genotype, annotations)
diet_D05_table <- make_table(res_D05, annotations)
diet_D2_table  <- make_table(res_D2, annotations)

interaction_D05_table <- make_table(res_interaction_D05_shrunk, annotations)
interaction_D2_table  <- make_table(res_interaction_D2_shrunk, annotations)

# =========================================================
# 9. FILTER SIGNIFICANT GENES
# =========================================================

padj_cutoff <- 0.05
log2fc_cutoff <- 1

filter_sig <- function(df) {
  df %>%
    filter(!is.na(padj),
           padj < padj_cutoff,
           abs(log2FoldChange) > log2fc_cutoff)
}

genotype_sig <- filter_sig(genotype_table)
diet_D05_sig <- filter_sig(diet_D05_table)
diet_D2_sig  <- filter_sig(diet_D2_table)

interaction_D05_sig <- filter_sig(interaction_D05_table)
interaction_D2_sig  <- filter_sig(interaction_D2_table)

# =========================================================
# 10. SUMMARY
# =========================================================

cat("\n============================\n")
cat("SUMMARY\n")
cat("============================\n")

cat("Genotype DEGs:", nrow(genotype_sig), "\n")
cat("Diet D05 DEGs:", nrow(diet_D05_sig), "\n")
cat("Diet D2 DEGs:", nrow(diet_D2_sig), "\n")

cat("Interaction D05:", nrow(interaction_D05_sig), "\n")
cat("Interaction D2:", nrow(interaction_D2_sig), "\n")

# =========================================================
# 11. EXPORT RESULTS
# =========================================================

dir.create("Summaries and Outputs/RNA", recursive = TRUE, showWarnings = FALSE)

write.csv(genotype_table, "Summaries and Outputs/RNA/genotype.csv", row.names = FALSE)
write.csv(diet_D05_table, "Summaries and Outputs/RNA/diet_D05.csv", row.names = FALSE)
write.csv(diet_D2_table,  "Summaries and Outputs/RNA/diet_D2.csv", row.names = FALSE)

write.csv(interaction_D05_table, "Summaries and Outputs/RNA/interaction_D05.csv", row.names = FALSE)
write.csv(interaction_D2_table,  "Summaries and Outputs/RNA/interaction_D2.csv", row.names = FALSE)

write.csv(genotype_sig, "Summaries and Outputs/RNA/genotype_sig.csv", row.names = FALSE)
write.csv(interaction_D05_sig, "Summaries and Outputs/RNA/interaction_D05_sig.csv", row.names = FALSE)
write.csv(interaction_D2_sig, "Summaries and Outputs/RNA/interaction_D2_sig.csv", row.names = FALSE)

# =========================================================
# DONE
# =========================================================

cat("\n✓ FULL ANALYSIS COMPLETE\n")


vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = c("genotype", "diet"))
plotPCA(vsd, intgroup = c("timepoint", "genotype"))

sampleDists <- dist(t(assay(vsd)))

pheatmap(
  as.matrix(sampleDists),
  annotation_col = sample_info
)

genotype_table$signif <- genotype_table$padj < 0.05

ggplot(genotype_table, aes(log2FoldChange, -log10(padj), color = signif)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Genotype DEGs")

gene <- interaction_D2_sig$gene_id[1]

plotCounts(
  dds,
  gene = gene,
  intgroup = c("genotype", "diet")
)

#----------------
library(clusterProfiler)
library(org.Mm.eg.db)

ego <- enrichGO(
  gene = genotype_sig$gene_id,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)

dotplot(ego)