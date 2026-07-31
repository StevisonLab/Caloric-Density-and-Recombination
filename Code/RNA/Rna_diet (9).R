# RNA-seq ANALYSIS — CALORIC DENSITY & RECOMBINATION PROJECT
# Two Genotypes: S42 and S399 | Diet: High Calorie (HC), Low Calorie (LC), Control (C)
# Days: 2 and 6

setwd("~/GitHub/Caloric-Density-and-Recombination/")


# LIBRARIES

library(UpSetR)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(cowplot)
library(ggnewscale)
library(grid)
library(gtable)
library(scales)


# PATHS

deg_path  <- "Raw Data/RNA/DEG/"
go_path   <- "Raw Data/RNA/GO/"
kegg_path <- "Raw Data/RNA/KEGG/"

out_deg  <- "Summaries and Outputs/RNA/DEG/"
out_go   <- "Summaries and Outputs/RNA/GO/"
out_kegg <- "Summaries and Outputs/RNA/KEGG/"

# Create output folders if they don't exist
for (d in c(out_deg, out_go, out_kegg)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# HELPER: Read files
# Your .xls files are plain tab-separated text — read them directly.

read_csv_file <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "")
}


# THRESHOLDS

p_threshold      <- 0.05
log2fc_threshold <- 1      # |log2FC| > 1  →  ≥ 2-fold change


# CONDITION LABELS

cond_hc      <- "High Calorie"
cond_lc      <- "Low Calorie"
cond_ctrl    <- "Control"


# PART 1 — LOAD OVERALL (POOLED) DEG FILES--------------------

HC_vs_Ctrl  <- read_csv_file(paste0(deg_path, "HighCalorievsControl_deg.xls"))
LC_vs_Ctrl  <- read_csv_file(paste0(deg_path, "LowCalorievsControl_deg.xls"))
LC_vs_HC <- read_csv_file(paste0(deg_path, "LowCalorievsHighCalorie_deg.xls"))
# Note: S42vsS399_deg is a genotype comparison — kept separate
S42_vs_S399 <- read_csv_file(paste0(deg_path, "S42vsS399_deg.xls"))

HC_vs_Ctrl$Comparison  <- "High Calorie vs Control"
LC_vs_Ctrl$Comparison  <- "Low Calorie vs Control"
LC_vs_HC$Comparison <- "Low Calorie vs High Calorie"

S42_vs_S399$Comparison <- "S42 vs S399"


# FILTER FUNCTION

filter_deg <- function(df) {
  df %>%
    mutate(
      log2FC    = log2FoldChange,
      neglog10p = -log10(padj)
    ) %>%
    filter(
      !is.na(log2FC),
      !is.na(padj),
      !is.infinite(neglog10p),
      padj < p_threshold,
      abs(log2FC) > log2fc_threshold
    )
}

HC_vs_Ctrl_sig  <- filter_deg(HC_vs_Ctrl)
LC_vs_Ctrl_sig  <- filter_deg(LC_vs_Ctrl)
LC_vs_HC_sig <- filter_deg(LC_vs_HC)
S42_vs_S399_sig <- filter_deg(S42_vs_S399)

#all_genes_overall <- bind_rows(HC_vs_Ctrl_sig, LC_vs_Ctrl_sig )

all_genes_overall <- bind_rows(HC_vs_Ctrl_sig, LC_vs_Ctrl_sig,LC_vs_HC_sig )

cat("\n=== OVERALL DEG SUMMARY (pooled across genotypes) ===\n")
cat("Thresholds: |log2FC| >", log2fc_threshold, "& padj <", p_threshold, "\n\n")
cat("High Calorie vs Control: ", nrow(HC_vs_Ctrl_sig),  " significant DEGs\n")
cat("Low Calorie  vs Control: ", nrow(LC_vs_Ctrl_sig),  " significant DEGs\n")
cat("Low Calorie vs High Calorie: ", nrow(LC_vs_HC_sig), " significant DEGs\n")
cat("S42 vs S399:             ", nrow(S42_vs_S399_sig), " significant DEGs\n")
cat("Total unique genes (diet comparisons):", length(unique(all_genes_overall$gene_id)), "\n\n")


# UPSET PLOT — OVERALL DIET COMPARISONS--------

# --- Update your UpSet Plot lists to include LC vs HC ---
deg_list_overall <- list(
  "High Calorie vs Control" = HC_vs_Ctrl_sig$gene_id,
  "Low Calorie vs Control"  = LC_vs_Ctrl_sig$gene_id,
  "Low Calorie vs High Calorie" = LC_vs_HC_sig$gene_id
)

png(paste0(out_deg, "upset_overall_diet_DEGs1.png"), width = 1200, height = 900, res = 150)
upset(
  fromList(deg_list_overall), 
  sets = c("Low Calorie vs High Calorie" ,"Low Calorie vs Control", "High Calorie vs Control"), 
  keep.order = TRUE, 
  main.bar.color = "gray38", 
  matrix.color = "red",
  point.size = 5, 
  line.size = 2
)
dev.off()

# COMMON DEGs — OVERALL----------
common_overall <- intersect(HC_vs_Ctrl_sig$gene_id, LC_vs_Ctrl_sig$gene_id)
cat("Genes shared by HC vs Ctrl AND LC vs Ctrl:", length(common_overall), "\n")

common_overall_info <- all_genes_overall %>%
  filter(gene_id %in% common_overall) %>%
  mutate(
    Direction  = ifelse(log2FoldChange > 0, "Up", "Down"),
    gene_label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  ) %>%
  select(gene_id, gene_label, Comparison, log2FoldChange, padj, Direction) %>%
  arrange(gene_id, Comparison)

write.csv(common_overall_info, paste0(out_deg, "common_DEGs_HC_LC_vs_Ctrl.csv"), row.names = FALSE)
cat("✓ Saved: common_DEGs_HC_LC_vs_Ctrl.csv\n\n")

# Export full summary
summary_overall <- all_genes_overall %>%
  mutate(
    Direction  = ifelse(log2FoldChange > 0, "Up", "Down"),
    gene_label = ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  ) %>%
  select(gene_id, gene_label, Comparison, log2FoldChange, padj, Direction) %>%
  arrange(Comparison, padj)

write.csv(summary_overall, paste0(out_deg, "summary_overall_significant_DEGs.csv"), row.names = FALSE)
cat("✓ Saved: summary_overall_significant_DEGs.csv\n")


# PART 2 — LOAD S42 TIME-SERIES DEG FILES (Day 2 & Day 6)-----------------

S42_HC_D2_vs_C_D2 <- read_csv_file(paste0(deg_path, "S42_HC_D2vsS42_C_D2_deg.xls"))
S42_HC_D6_vs_C_D6 <- read_csv_file(paste0(deg_path, "S42_HC_D6vsS42_C_D6_deg.xls"))
S42_LC_D2_vs_C_D2 <- read_csv_file(paste0(deg_path, "S42_LC_D2vsS42_C_D2_deg.xls"))
S42_LC_D6_vs_C_D6 <- read_csv_file(paste0(deg_path, "S42_LC_D6vsS42_C_D6_deg.xls"))
S42_LC_D2_vs_HC_D2 <- read_csv_file(paste0(deg_path, "S42_LC_D2vsS42_HC_D2_deg.xls"))
S42_LC_D6_vs_HC_D6 <- read_csv_file(paste0(deg_path, "S42_LC_D6vsS42_HC_D6_deg.xls"))


# LOAD S399 TIME-SERIES DEG FILES (Day 2 & Day 6)

S399_HC_D2_vs_C_D2  <- read_csv_file(paste0(deg_path, "S399_HC_D2vsS399_C_D2_deg.xls"))
S399_HC_D6_vs_C_D6  <- read_csv_file(paste0(deg_path, "S399_HC_D6vsS399_C_D6_deg.xls"))
S399_LC_D2_vs_C_D2  <- read_csv_file(paste0(deg_path, "S399_LC_D2vsS399_C_D2_deg.xls"))
S399_LC_D6_vs_C_D6  <- read_csv_file(paste0(deg_path, "S399_LC_D6vsS399_C_D6_deg.xls"))
S399_LC_D2_vs_HC_D2 <- read_csv_file(paste0(deg_path, "S399_LC_D2vsS399_HC_D2_deg.xls"))
S399_LC_D6_vs_HC_D6 <- read_csv_file(paste0(deg_path, "S399_LC_D6vsS399_HC_D6_deg.xls"))


# PART 3 — UP/DOWN REGULATION SUMMARY TABLE (per genotype, per day)--------------

count_deg_direction <- function(df, comparison_name,
                                padj_cutoff = 0.05,
                                log2fc_thresh = 1) {
  df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      direction = case_when(
        padj < padj_cutoff & log2FoldChange >  log2fc_thresh ~ "Upregulated",
        padj < padj_cutoff & log2FoldChange < -log2fc_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    ) %>%
    count(direction) %>%
    mutate(comparison = comparison_name)
}

results_updown <- bind_rows(
  # --- S42 Day 2 ---
  count_deg_direction(S42_HC_D2_vs_C_D2,  "S42 Day2: HC vs Control",  p_threshold, log2fc_threshold),
  count_deg_direction(S42_LC_D2_vs_C_D2,  "S42 Day2: LC vs Control",  p_threshold, log2fc_threshold),
  count_deg_direction(S42_LC_D2_vs_HC_D2, "S42 Day2: LC vs HC",       p_threshold, log2fc_threshold),
  # --- S42 Day 6 ---
  count_deg_direction(S42_HC_D6_vs_C_D6,  "S42 Day6: HC vs Control",  p_threshold, log2fc_threshold),
  count_deg_direction(S42_LC_D6_vs_C_D6,  "S42 Day6: LC vs Control",  p_threshold, log2fc_threshold),
  count_deg_direction(S42_LC_D6_vs_HC_D6, "S42 Day6: LC vs HC",       p_threshold, log2fc_threshold),
  # --- S399 Day 2 ---
  count_deg_direction(S399_HC_D2_vs_C_D2,  "S399 Day2: HC vs Control", p_threshold, log2fc_threshold),
  count_deg_direction(S399_LC_D2_vs_C_D2,  "S399 Day2: LC vs Control", p_threshold, log2fc_threshold),
  count_deg_direction(S399_LC_D2_vs_HC_D2, "S399 Day2: LC vs HC",      p_threshold, log2fc_threshold),
  # --- S399 Day 6 ---
  count_deg_direction(S399_HC_D6_vs_C_D6,  "S399 Day6: HC vs Control", p_threshold, log2fc_threshold),
  count_deg_direction(S399_LC_D6_vs_C_D6,  "S399 Day6: LC vs Control", p_threshold, log2fc_threshold),
  count_deg_direction(S399_LC_D6_vs_HC_D6, "S399 Day6: LC vs HC",      p_threshold, log2fc_threshold)
)

updown_summary <- results_updown %>%
  tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0) %>%
  mutate(Total_Significant = Upregulated + Downregulated) %>%
  select(comparison, Upregulated, Downregulated, Total_Significant)

cat("\nThresholds applied: |log2FC| >", log2fc_threshold, "& padj <", p_threshold, "\n\n")
print(updown_summary)

write.csv(updown_summary, paste0(out_deg, "deg_up_down_summary_by_genotype_day.csv"), row.names = FALSE)
cat("✓ Saved: deg_up_down_summary_by_genotype_day.csv\n")

# PART 4 — UPSET PLOTS PER GENOTYPE PER DAY-------------

make_upset <- function(list_of_gene_sets, set_order, outfile,
                       matrix_col = "red") {
  
  # ---- CLEAN ----
  list_of_gene_sets <- lapply(list_of_gene_sets, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    unique(x)
  })
  
  cat("\nUPSET CHECK:\n")
  print(sapply(list_of_gene_sets, length))
  
  # ---- BUILD MATRIX FIRST (IMPORTANT FIX) ----
  upset_input <- fromList(list_of_gene_sets)
  
  png(outfile, width = 1400, height = 900, res = 160)
  
  p <- UpSetR::upset(
    upset_input,
    sets = set_order,
    keep.order = TRUE,
    order.by = "freq",
    main.bar.color = "gray38",
    matrix.color = matrix_col,
    sets.bar.color = rep("gray40", length(set_order)),
    point.size = 6,
    line.size = 2.5,
    text.scale = c(2, 2, 2, 1.5, 2, 1.5),
    mb.ratio = c(0.55, 0.45)
  )
  
  print(p)   # 🔥 CRITICAL FIX (forces rendering)
  
  dev.off()
}

# ---------- DEG FILTER + GENE EXTRACTION ----------
filter_deg_genes <- function(df) {
  
  df %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      padj < p_threshold,
      abs(log2FoldChange) > log2fc_threshold
    ) %>%
    pull(gene_id) %>%
    unique()
  
}


# S42 — DAY 2


S42_D2_list <- list(
  "HC vs Control" = filter_deg_genes(S42_HC_D2_vs_C_D2),
  "LC vs Control" = filter_deg_genes(S42_LC_D2_vs_C_D2),
  "LC vs HC"      = filter_deg_genes(S42_LC_D2_vs_HC_D2)
)

cat("\nS42 Day2 DEG counts:\n")
print(sapply(S42_D2_list, length))

make_upset(
  S42_D2_list,
  c("LC vs HC", "LC vs Control", "HC vs Control"),
  paste0(out_deg, "upset_S42_Day2.png"),
  matrix_col = "red"
)

cat("✓ Saved: upset_S42_Day2.png\n")



# S42 — DAY 6


S42_D6_list <- list(
  "HC vs Control" = filter_deg_genes(S42_HC_D6_vs_C_D6),
  "LC vs Control" = filter_deg_genes(S42_LC_D6_vs_C_D6),
  "LC vs HC"      = filter_deg_genes(S42_LC_D6_vs_HC_D6)
)

cat("\nS42 Day6 DEG counts:\n")
print(sapply(S42_D6_list, length))

make_upset(
  S42_D6_list,
  c("LC vs HC", "LC vs Control", "HC vs Control"),
  paste0(out_deg, "upset_S42_Day6.png"),
  matrix_col = "red4"
)

cat("✓ Saved: upset_S42_Day6.png\n")



# S399 — DAY 2


S399_D2_list <- list(
  "HC vs Control" = filter_deg_genes(S399_HC_D2_vs_C_D2),
  "LC vs Control" = filter_deg_genes(S399_LC_D2_vs_C_D2),
  "LC vs HC"      = filter_deg_genes(S399_LC_D2_vs_HC_D2)
)

cat("\nS399 Day2 DEG counts:\n")
print(sapply(S399_D2_list, length))

make_upset(
  S399_D2_list,
  c("LC vs HC", "LC vs Control", "HC vs Control"),
  paste0(out_deg, "upset_S399_Day2.png"),
  matrix_col = "skyblue"
)

cat("✓ Saved: upset_S399_Day2.png\n")



# S399 — DAY 6


S399_D6_list <- list(
  "HC vs Control" = filter_deg_genes(S399_HC_D6_vs_C_D6),
  "LC vs Control" = filter_deg_genes(S399_LC_D6_vs_C_D6),
  "LC vs HC"      = filter_deg_genes(S399_LC_D6_vs_HC_D6)
)

cat("\nS399 Day6 DEG counts:\n")
print(sapply(S399_D6_list, length))

make_upset(
  S399_D6_list,
  c("LC vs HC", "LC vs Control", "HC vs Control"),
  paste0(out_deg, "upset_S399_Day6.png"),
  matrix_col = "blue"
)

cat("✓ Saved: upset_S399_Day6.png\n")

# PART 4B — UPSET PLOTS PER GENOTYPE, COMBINING DAY 2 & DAY 6
# Each of the 3 comparisons appears twice (D2 + D6), colored by day
# ------------------------------------------------------------------

make_upset_byday <- function(list_of_gene_sets, set_order, day_colors, outfile) {

  list_of_gene_sets <- lapply(list_of_gene_sets, function(x) {
    x <- as.character(x)
    x <- x[!is.na(x)]
    unique(x)
  })

  cat("\nUPSET (by day) CHECK:\n")
  print(sapply(list_of_gene_sets, length))

  upset_input <- fromList(list_of_gene_sets)

  png(outfile, width = 1700, height = 1000, res = 160)

  p <- UpSetR::upset(
    upset_input,
    sets = set_order,
    keep.order = TRUE,
    order.by = "freq",
    main.bar.color = "gray38",
    matrix.color = "black",
    sets.bar.color = day_colors,   # one color per set, matched to set_order
    point.size = 6,
    line.size = 2.5,
    text.scale = c(2, 2, 2, 1.5, 2, 1.5),
    mb.ratio = c(0.55, 0.45)
  )

  print(p)
  dev.off()
}

# Shared set order and colors (2 shades per comparison: dark = D2, light = D6)
byday_set_order <- c(
  "LC vs HC (D2)",      "LC vs HC (D6)",
  "LC vs Control (D2)", "LC vs Control (D6)",
  "HC vs Control (D2)", "HC vs Control (D6)"
)

byday_colors <- c(
  "#E41A1C", "#FB9A99",   # LC vs HC        — D2 dark red   / D6 light red
  "#377EB8", "#A6CEE3",   # LC vs Control   — D2 dark blue  / D6 light blue
  "#4DAF4A", "#B2DF8A"    # HC vs Control   — D2 dark green / D6 light green
)

# ---- S42: combine Day2 + Day6 sets ----
S42_byday_list <- list(
  "HC vs Control (D2)" = S42_D2_list[["HC vs Control"]],
  "HC vs Control (D6)" = S42_D6_list[["HC vs Control"]],
  "LC vs Control (D2)" = S42_D2_list[["LC vs Control"]],
  "LC vs Control (D6)" = S42_D6_list[["LC vs Control"]],
  "LC vs HC (D2)"      = S42_D2_list[["LC vs HC"]],
  "LC vs HC (D6)"      = S42_D6_list[["LC vs HC"]]
)

make_upset_byday(
  S42_byday_list,
  byday_set_order,
  byday_colors,
  paste0(out_deg, "upset_S42_byDay_combined.png")
)
cat("✓ Saved: upset_S42_byDay_combined.png\n")

# ---- S399: combine Day2 + Day6 sets ----
S399_byday_list <- list(
  "HC vs Control (D2)" = S399_D2_list[["HC vs Control"]],
  "HC vs Control (D6)" = S399_D6_list[["HC vs Control"]],
  "LC vs Control (D2)" = S399_D2_list[["LC vs Control"]],
  "LC vs Control (D6)" = S399_D6_list[["LC vs Control"]],
  "LC vs HC (D2)"      = S399_D2_list[["LC vs HC"]],
  "LC vs HC (D6)"      = S399_D6_list[["LC vs HC"]]
)

make_upset_byday(
  S399_byday_list,
  byday_set_order,
  byday_colors,
  paste0(out_deg, "upset_S399_byDay_combined.png")
)
cat("✓ Saved: upset_S399_byDay_combined.png\n")

# PART 4C — UPSET PLOTS PER GENOTYPE: 3 SETS (comparisons), BARS SPLIT BY DAY
# One row per comparison (HC vs Control / LC vs Control / LC vs HC).
# Set-size bars AND intersection bars are stacked: Day 2 vs Day 6 (+ "Both",
# for genes significant on both days in that comparison, so bars stay exact).
# Requires: ComplexHeatmap
## orange "#ffbb00", "#ff7f00"
##purple   "#CCA6D2" "#984ea3  and 
# ------------------------------------------------------------------------

# install.packages("BiocManager"); BiocManager::install("ComplexHeatmap")  # if needed
library(ComplexHeatmap)
library(grid)

# ---- per-gene day category for one comparison ----
gene_day_category <- function(genes, d2_genes, d6_genes) {
  in_d2 <- genes %in% d2_genes
  in_d6 <- genes %in% d6_genes
  ifelse(in_d2 & in_d6, "Both",
         ifelse(in_d2, "Day 2",
                ifelse(in_d6, "Day 6", NA)))
}

make_stacked_upset <- function(D2_list, D6_list, set_order, outfile, 
                               width = 12, height = 7) {
  
  # 1. Union of sets for UpSet combination matrix
  union_list <- lapply(set_order, function(s) union(D2_list[[s]], D6_list[[s]]))
  names(union_list) <- set_order
  all_genes <- unique(unlist(union_list))
  
  m <- sapply(set_order, function(s) all_genes %in% union_list[[s]])
  rownames(m) <- all_genes
  colnames(m) <- set_order
  m01 <- 1 * m
  comb_mat <- make_comb_mat(m01)
  
  # Sort combinations by size (descending)
  combo_order <- order(-comb_size(comb_mat))
  
  # Total combined DEGs per intersection column
  column_totals <- comb_size(comb_mat)[combo_order]
  
  # 2. Build STACKED INTERSECTION matrix (Rows = Intersections, Cols = Day 2 / Day 6)
  comb_elements <- comb_name(comb_mat)[combo_order]
  
  combo_mat_ordered <- t(sapply(comb_elements, function(cc) {
    # Genes in this specific intersection combination
    genes_in_comb <- extract_comb(comb_mat, cc)
    
    # Count how many belong to Day 2 vs Day 6 lists
    d2_count <- sum(sapply(D2_list, function(x) genes_in_comb %in% x)) > 0
    d6_count <- sum(sapply(D6_list, function(x) genes_in_comb %in% x)) > 0
    
    c("Day 2" = length(intersect(genes_in_comb, unique(unlist(D2_list)))),
      "Day 6" = length(intersect(genes_in_comb, unique(unlist(D6_list)))))
  }))
  
  # 3. Build STACKED SET SIZE matrix (Rows = Sets, Cols = Day 2 / Day 6)
  set_size_mat <- t(sapply(set_order, function(s) {
    c("Day 2" = length(D2_list[[s]]), 
      "Day 6" = length(D6_list[[s]]))
  }))
  
  # 4. Top Annotation: Stacked Barplot + DEG Total Count Labels
  top_ann <- HeatmapAnnotation(
    "Total DEGs" = anno_text(
      as.character(column_totals), 
      rot = 0, 
      just = "center", 
      gp = gpar(fontsize = 9, fontface = "bold")
    ),
    "Intersection Size" = anno_barplot(
      combo_mat_ordered, 
      beside = FALSE, 
      gp = gpar(fill = c("#CCA6D2", "#984ea3")), 
      border = FALSE, 
      bar_width = 0.8, 
      height = unit(12, "cm")
    ),
    annotation_name_side = "left",
    gap = unit(2, "mm")
  )
  
  # 5. Left Annotation: Set sizes
  left_ann <- rowAnnotation(
    "Set Size" = anno_barplot(
      set_size_mat, 
      beside = FALSE, 
      gp = gpar(fill = c("#CCA6D2", "#984ea3")), 
      border = FALSE, 
      bar_width = 0.8, 
      axis_param = list(direction = "reverse"), 
      width = unit(4, "cm")
    )
  )
  
  # 6. Build Plot
  ht <- UpSet(
    comb_mat, 
    comb_order = combo_order, 
    set_order = set_order, 
    top_annotation = top_ann, 
    left_annotation = left_ann, 
    right_annotation = NULL, 
    row_names_side = "left"
  )
  
  png(outfile, width = width, height = height, units = "in", res = 300, bg = "white")
  draw(ht)
  
  # Custom Legend for Day 2 vs Day 6
  day_colors <- c("Day 2" = "#CCA6D2", "Day 6" = "#984ea3")
  lgd <- Legend(labels = names(day_colors), legend_gp = gpar(fill = day_colors), title = "Timepoint")
  draw(lgd, x = unit(1, "npc") - unit(5, "mm"), y = unit(1, "npc") - unit(20, "mm"), just = c("right", "top"))
  
  dev.off()
  cat("✓ Saved:", basename(outfile), "\n")
}
# ---- S42 ----
make_stacked_upset(
  D2_list = S42_D2_list, D6_list = S42_D6_list,
  set_order = c("LC vs HC", "LC vs Control", "HC vs Control"),
  outfile = paste0(out_deg, "upset_S42_stackedByDay.png")
)

# ---- S399 ----
make_stacked_upset(
  D2_list = S399_D2_list, D6_list = S399_D6_list,
  set_order = c("LC vs HC", "LC vs Control", "HC vs Control"),
  outfile = paste0(out_deg, "upset_S399_stackedByDay.png")
)


# PART 5 — VOLCANO PLOTS-------------------

# ---------- Legend colour scheme ----------
legend_colors <- c(
  "#b2182b",   # HC (right side / up)
  "#2166ac",   # LC (right side / up)
  "#bf812d",   # Control (right side / up)
  "gray70"     # Not significant
)
names(legend_colors) <- c(
  paste0("Upregulated in ", cond_hc),
  paste0("Upregulated in ", cond_lc),
  paste0("Upregulated in ", cond_ctrl),
  "Not significant"
)

# ---------- Shared axis limits across all 12 per-genotype-day plots ----------
all_fc <- c(
  S42_HC_D2_vs_C_D2$log2FoldChange,  S42_HC_D6_vs_C_D6$log2FoldChange,
  S42_LC_D2_vs_C_D2$log2FoldChange,  S42_LC_D6_vs_C_D6$log2FoldChange,
  S42_LC_D2_vs_HC_D2$log2FoldChange, S42_LC_D6_vs_HC_D6$log2FoldChange,
  S399_HC_D2_vs_C_D2$log2FoldChange,  S399_HC_D6_vs_C_D6$log2FoldChange,
  S399_LC_D2_vs_C_D2$log2FoldChange,  S399_LC_D6_vs_C_D6$log2FoldChange,
  S399_LC_D2_vs_HC_D2$log2FoldChange, S399_LC_D6_vs_HC_D6$log2FoldChange
)
all_p <- c(
  -log10(S42_HC_D2_vs_C_D2$padj),  -log10(S42_HC_D6_vs_C_D6$padj),
  -log10(S42_LC_D2_vs_C_D2$padj),  -log10(S42_LC_D6_vs_C_D6$padj),
  -log10(S42_LC_D2_vs_HC_D2$padj), -log10(S42_LC_D6_vs_HC_D6$padj),
  -log10(S399_HC_D2_vs_C_D2$padj),  -log10(S399_HC_D6_vs_C_D6$padj),
  -log10(S399_LC_D2_vs_C_D2$padj),  -log10(S399_LC_D6_vs_C_D6$padj),
  -log10(S399_LC_D2_vs_HC_D2$padj), -log10(S399_LC_D6_vs_HC_D6$padj)
)

x_max    <- ceiling(max(abs(all_fc), na.rm = TRUE))
#x_limits <- c(-x_max, x_max)
x_limits <- c(-10, 10)
y_limits <- c(0, max(all_p, na.rm = TRUE) * 1.1)

# ---------- Volcano function ----------
create_volcano_plot <- function(deg_data, title, condition_right, condition_left,
                                pval_threshold = 0.05,
                                log2fc_threshold = 1,
                                top_genes = 20,
                                ylim = NULL, xlim = NULL) {
  
  n_each <- ceiling(top_genes / 2)
  
  deg_data <- deg_data %>%
    mutate(
      log2FC = log2FoldChange,
      neglog10p = -log10(padj),
      significance = case_when(
        log2FC > log2fc_threshold & padj < pval_threshold ~ paste0("Upregulated in ", condition_right),
        log2FC < -log2fc_threshold & padj < pval_threshold ~ paste0("Upregulated in ", condition_left),
        TRUE ~ "Not significant"
      )
    )
  
  deg_data$significance <- factor(deg_data$significance, levels = names(legend_colors))
  
  deg_data <- deg_data %>% 
    filter(!is.na(log2FC), !is.na(neglog10p), !is.infinite(neglog10p))
  #filter two logfold 2 changes in the plotting
  deg_data <- deg_data %>%
    filter(log2FC > -10, log2FC < 10)
  
  # Select top genes
  genes_right <- deg_data %>%
    filter(significance == paste0("Upregulated in ", condition_right)) %>%
    arrange(padj, desc(abs(log2FC))) %>% head(n_each)
  
  genes_left <- deg_data %>%
    filter(significance == paste0("Upregulated in ", condition_left)) %>%
    arrange(padj, desc(abs(log2FC))) %>% head(n_each)
  
  genes_to_label <- bind_rows(genes_right, genes_left)
  
  # === ONLY USE gene_name ===
  if (!"gene_name" %in% colnames(deg_data)) {
    stop("Column 'gene_name' is required but not found in deg_data.")
  }
  
  # Make sure genes_to_label also has gene_name
  if (!"gene_name" %in% colnames(genes_to_label)) {
    genes_to_label$gene_name <- deg_data$gene_name[match(rownames(genes_to_label), rownames(deg_data))]
  }
  
  p <- ggplot(deg_data, aes(x = log2FC, y = neglog10p)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 8)+ #SIZE OF THE DOTS
    scale_color_manual(values = legend_colors, limits = names(legend_colors), drop = FALSE) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
               linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray70", linewidth = 0.3) +
    geom_text_repel(
      data = genes_to_label,
      aes(label = gene_name),      # <-- Only gene_name is used AND SIZE OF THE GENE NAME
      size = 9, box.padding = 0.3, point.padding = 0.2,
      segment.color = "gray50", segment.size = 0.4,
      min.segment.length = 0.05, force = 1.5, max.overlaps = Inf
    ) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Legend",
   #   caption = paste0("|log2FC| > ", log2fc_threshold, " & padj < ", pval_threshold)
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      plot.caption = element_text(hjust = 0.5, size = 9, color = "gray40"),
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 16, face = "bold"),
      legend.position = "left",
      legend.key.size = unit(2, "cm"),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18, face = "bold")
    )
  
  if (!is.null(ylim) || !is.null(xlim)) {
    p <- p + coord_cartesian(ylim = ylim, xlim = xlim)
  }
  
  return(p)
}
# ---------- Shared legend ----------
dummy <- ggplot(
  data.frame(x = 1:4, y = 1:4,
             label = factor(names(legend_colors), levels = names(legend_colors))),
  aes(x = x, y = y, color = label)
) +
  geom_point(size = 10) +
  scale_color_manual(values = legend_colors, drop = FALSE) +
  guides(color = guide_legend(title = "Legend", override.aes = list(size = 10))) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title    = element_text(size = 22, face = "bold"),
    legend.text     = element_text(size = 18, face = "bold"),
    legend.key.size = unit(2.2, "cm")
  )
leg <- get_legend(dummy)
stopifnot(!inherits(leg, "zeroGrob"))

# ---------- NOTE on log2FoldChange direction ----------
# For "X vs Control" files that come as "Control vs X",
# flip the sign so the named condition always sits on the RIGHT.
# Adjust the lines below if your DESeq2 output already has the correct sign.

# ---- Helper to flip sign ----
flip_fc <- function(df) { df %>% mutate(log2FoldChange = -log2FoldChange) }

# S42 — flip so HC / LC are always on the right (positive) side
S42_HC_D2 <- flip_fc(S42_HC_D2_vs_C_D2)
S42_HC_D6 <- flip_fc(S42_HC_D6_vs_C_D6)
S42_LC_D2 <- flip_fc(S42_LC_D2_vs_C_D2)
S42_LC_D6 <- flip_fc(S42_LC_D6_vs_C_D6)
# LC vs HC: flip so LC is on the right
S42_LC_vs_HC_D2 <- flip_fc(S42_LC_D2_vs_HC_D2)
S42_LC_vs_HC_D6 <- flip_fc(S42_LC_D6_vs_HC_D6)

# S399 — same logic
S399_HC_D2 <- flip_fc(S399_HC_D2_vs_C_D2)
S399_HC_D6 <- flip_fc(S399_HC_D6_vs_C_D6)
S399_LC_D2 <- flip_fc(S399_LC_D2_vs_C_D2)
S399_LC_D6 <- flip_fc(S399_LC_D6_vs_C_D6)
S399_LC_vs_HC_D2 <- flip_fc(S399_LC_D2_vs_HC_D2)
S399_LC_vs_HC_D6 <- flip_fc(S399_LC_D6_vs_HC_D6)

# ---------- Build all 12 individual volcano plots ----------
# S42 Day 2
v_S42_HC_D2 <- create_volcano_plot(S42_HC_D2, paste(cond_hc, "vs", cond_ctrl, "— S42 Day 2"),
  cond_hc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S42_LC_D2 <- create_volcano_plot(S42_LC_D2, paste(cond_lc, "vs", cond_ctrl, "— S42 Day 2"),
  cond_lc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S42_LH_D2 <- create_volcano_plot(S42_LC_vs_HC_D2, paste(cond_lc, "vs", cond_hc, "— S42 Day 2"),
  cond_lc, cond_hc, xlim = x_limits, ylim = y_limits)

# S42 Day 6
v_S42_HC_D6 <- create_volcano_plot(S42_HC_D6, paste(cond_hc, "vs", cond_ctrl, "— S42 Day 6"),
  cond_hc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S42_LC_D6 <- create_volcano_plot(S42_LC_D6, paste(cond_lc, "vs", cond_ctrl, "— S42 Day 6"),
  cond_lc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S42_LH_D6 <- create_volcano_plot(S42_LC_vs_HC_D6, paste(cond_lc, "vs", cond_hc, "— S42 Day 6"),
  cond_lc, cond_hc, xlim = x_limits, ylim = y_limits)

# S399 Day 2
v_S399_HC_D2 <- create_volcano_plot(S399_HC_D2, paste(cond_hc, "vs", cond_ctrl, "— S399 Day 2"),
  cond_hc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S399_LC_D2 <- create_volcano_plot(S399_LC_D2, paste(cond_lc, "vs", cond_ctrl, "— S399 Day 2"),
  cond_lc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S399_LH_D2 <- create_volcano_plot(S399_LC_vs_HC_D2, paste(cond_lc, "vs", cond_hc, "— S399 Day 2"),
  cond_lc, cond_hc, xlim = x_limits, ylim = y_limits)

# S399 Day 6
v_S399_HC_D6 <- create_volcano_plot(S399_HC_D6, paste(cond_hc, "vs", cond_ctrl, "— S399 Day 6"),
  cond_hc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S399_LC_D6 <- create_volcano_plot(S399_LC_D6, paste(cond_lc, "vs", cond_ctrl, "— S399 Day 6"),
  cond_lc, cond_ctrl, xlim = x_limits, ylim = y_limits)
v_S399_LH_D6 <- create_volcano_plot(S399_LC_vs_HC_D6, paste(cond_lc, "vs", cond_hc, "— S399 Day 6"),
  cond_lc, cond_hc, xlim = x_limits, ylim = y_limits)

# ---------- Strip legends ----------
nl <- function(p) p + theme(legend.position = "none")

# ---------- Combined: S42 (Day 2 top row, Day 6 bottom row) ----------
combined_S42 <- plot_grid(
  nl(v_S42_HC_D2), nl(v_S42_LC_D2), nl(v_S42_LH_D2),
  nl(v_S42_HC_D6), nl(v_S42_LC_D6), nl(v_S42_LH_D6),
  ncol = 3, nrow = 2, align = "hv",
  labels = c("Day 2", "", "", "Day 6", "", ""),
  label_size = 16, label_fontface = "bold"
)
final_S42 <- plot_grid(combined_S42, leg, ncol = 2, rel_widths = c(3, 0.8))
ggsave(paste0(out_deg, "volcano_S42_Day2_and_Day6.png"),
       final_S42, width = 26, height = 12, dpi = 600)
cat("✓ Saved: volcano_S42_Day2_and_Day6.png\n")

# ---------- Combined: S399 ----------
combined_S399 <- plot_grid(
  nl(v_S399_HC_D2), nl(v_S399_LC_D2), nl(v_S399_LH_D2),
  nl(v_S399_HC_D6), nl(v_S399_LC_D6), nl(v_S399_LH_D6),
  ncol = 3, nrow = 2, align = "hv",
  labels = c("Day 2", "", "", "Day 6", "", ""),
  label_size = 16, label_fontface = "bold"
)
final_S399 <- plot_grid(combined_S399, leg, ncol = 2, rel_widths = c(3, 0.8))
ggsave(paste0(out_deg, "volcano_S399_Day2_and_Day6.png"),
       final_S399, width = 26, height = 12, dpi = 600)
cat("✓ Saved: volcano_S399_Day2_and_Day6.png\n")

# ---------- Grand combined: S42 + S399 (4 rows × 3 cols) ----------
combined_all <- plot_grid(
  nl(v_S42_HC_D2),  nl(v_S42_LC_D2),  nl(v_S42_LH_D2),
  nl(v_S42_HC_D6),  nl(v_S42_LC_D6),  nl(v_S42_LH_D6),
  nl(v_S399_HC_D2), nl(v_S399_LC_D2), nl(v_S399_LH_D2),
  nl(v_S399_HC_D6), nl(v_S399_LC_D6), nl(v_S399_LH_D6),
  ncol = 3, nrow = 4, align = "hv",
#  labels = c("S42 Day 2", "", "", "S42 Day 6", "", "",
#             "S399 Day 2", "", "", "S399 Day 6", "", ""),
  label_size = 14, label_fontface = "bold"
)
final_all <- plot_grid(combined_all, leg, ncol = 2, rel_widths = c(3, 0.8))
ggsave(paste0(out_deg, "volcano_ALL_genotypes_days.png"),
       final_all, width = 30, height = 19, dpi = 600)
cat("✓ Saved: volcano_ALL_genotypes_days.png\n")


# ── Save LC vs HC volcano plots individually (no legend) ---------------

ggsave(paste0(out_deg, "volcano_S42_LCvsHC_Day2.png"),
       nl(v_S42_LH_D2), width = 10, height = 8, dpi = 600)
cat("✓ Saved: volcano_S42_LCvsHC_Day2.png\n")

ggsave(paste0(out_deg, "volcano_S42_LCvsHC_Day6.png"),
       nl(v_S42_LH_D6), width = 10, height = 8, dpi = 600)
cat("✓ Saved: volcano_S42_LCvsHC_Day6.png\n")

ggsave(paste0(out_deg, "volcano_S399_LCvsHC_Day2.png"),
       nl(v_S399_LH_D2), width = 10, height = 8, dpi = 600)
cat("✓ Saved: volcano_S399_LCvsHC_Day2.png\n")

ggsave(paste0(out_deg, "volcano_S399_LCvsHC_Day6.png"),
       nl(v_S399_LH_D6), width = 10, height = 8, dpi = 600)
cat("✓ Saved: volcano_S399_LCvsHC_Day6.png\n")

# PART 6 — GO ENRICHMENT DOT PLOTS-----------

# ---------- Load all GO files ----------
load_go <- function(file, comparison, direction) {
  path <- paste0(go_path, file)
  if (!file.exists(path)) {
    message("  [SKIP - not found] ", file)
    return(data.frame())
  }
  df <- tryCatch(read_csv_file(path), error = function(e) data.frame())
  if (nrow(df) == 0) {
    message("  [SKIP - empty, nothing significant] ", file)
    return(data.frame())
  }
  df$Comparison <- comparison
  df$Direction  <- direction
  return(df)
}

go_up <- bind_rows(
  load_go("HighCalorievsControl.up_GOenrich_significant.xls",       "HC vs Control",  "Up"),
  load_go("LowCalorievsControl.up_GOenrich_significant.xls",        "LC vs Control",  "Up"),
  load_go("LowCalorievsHighCalorie.up_GOenrich_significant.xls",    "LC vs HC",       "Up"),
  load_go("S42_HC_D2vsS42_C_D2.up_GOenrich_significant.xls",       "S42 HC D2 vs C", "Up"),
  load_go("S42_HC_D6vsS42_C_D6.up_GOenrich_significant.xls",       "S42 HC D6 vs C", "Up"),
  load_go("S42_LC_D2vsS42_C_D2.up_GOenrich_significant.xls",       "S42 LC D2 vs C", "Up"),
  load_go("S42_LC_D2vsS42_HC_D2.up_GOenrich_significant.xls",      "S42 LC D2 vs HC","Up"),
  load_go("S42_LC_D6vsS42_C_D6.up_GOenrich_significant.xls",       "S42 LC D6 vs C", "Up"),
  load_go("S42_LC_D6vsS42_HC_D6.up_GOenrich_significant.xls",      "S42 LC D6 vs HC","Up"),
  load_go("S399_HC_D2vsS399_C_D2.up_GOenrich_significant.xls",     "S399 HC D2 vs C","Up"),
  load_go("S399_HC_D6vsS399_C_D6.up_GOenrich_significant.xls",     "S399 HC D6 vs C","Up"),
  load_go("S399_LC_D2vsS399_C_D2.up_GOenrich_significant.xls",     "S399 LC D2 vs C","Up"),
  load_go("S399_LC_D2vsS399_HC_D2.up_GOenrich_significant.xls",    "S399 LC D2 vs HC","Up"),
  load_go("S399_LC_D6vsS399_C_D6.up_GOenrich_significant.xls",     "S399 LC D6 vs C","Up"),
  load_go("S399_LC_D6vsS399_HC_D6.up_GOenrich_significant.xls",    "S399 LC D6 vs HC","Up")
) %>% filter(padj < 0.05)

go_down <- bind_rows(
  load_go("HighCalorievsControl.down_GOenrich_significant.xls",     "HC vs Control",  "Down"),
  load_go("LowCalorievsControl.down_GOenrich_significant.xls",      "LC vs Control",  "Down"),
  load_go("LowCalorievsHighCalorie.down_GOenrich_significant.xls",  "LC vs HC",       "Down"),
  load_go("S42_HC_D2vsS42_C_D2.down_GOenrich_significant.xls",     "S42 HC D2 vs C", "Down"),
  load_go("S42_HC_D6vsS42_C_D6.down_GOenrich_significant.xls",     "S42 HC D6 vs C", "Down"),
  load_go("S42_LC_D2vsS42_C_D2.down_GOenrich_significant.xls",     "S42 LC D2 vs C", "Down"),
  load_go("S42_LC_D2vsS42_HC_D2.down_GOenrich_significant.xls",    "S42 LC D2 vs HC","Down"),
  load_go("S42_LC_D6vsS42_C_D6.down_GOenrich_significant.xls",     "S42 LC D6 vs C", "Down"),
  load_go("S42_LC_D6vsS42_HC_D6.down_GOenrich_significant.xls",    "S42 LC D6 vs HC","Down"),
  load_go("S399_HC_D2vsS399_C_D2.down_GOenrich_significant.xls",   "S399 HC D2 vs C","Down"),
  load_go("S399_HC_D6vsS399_C_D6.down_GOenrich_significant.xls",   "S399 HC D6 vs C","Down"),
  load_go("S399_LC_D2vsS399_C_D2.down_GOenrich_significant.xls",   "S399 LC D2 vs C","Down"),
  load_go("S399_LC_D2vsS399_HC_D2.down_GOenrich_significant.xls",  "S399 LC D2 vs HC","Down"),
  load_go("S399_LC_D6vsS399_C_D6.down_GOenrich_significant.xls",   "S399 LC D6 vs C","Down"),
  load_go("S399_LC_D6vsS399_HC_D6.down_GOenrich_significant.xls",  "S399 LC D6 vs HC","Down")
) %>% filter(padj < 0.05)

cat("GO up entries:  ", nrow(go_up), "\n")
cat("GO down entries:", nrow(go_down), "\n")

# ====================== HELPER FUNCTIONS ======================

prepare_dotplot_data <- function(df) {
  df %>%
    mutate(
      neg_log10_padj = -log10(padj),
      GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"), 
                                 function(x) if(length(x)==2) as.numeric(x[1])/as.numeric(x[2]) else NA_real_),
      Description_short = ifelse(nchar(Description) > 55,
                                 paste0(substr(Description, 1, 52), "..."), 
                                 Description)
    )
}

make_go_dotplot <- function(df, direction, title_label, 
                            all_comparisons = NULL) {   # <--- New argument
  
  if (is.null(all_comparisons)) {
    all_comparisons <- unique(df$Comparison)
  }
  
  if (nrow(df) == 0) {
    p <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No significant GO terms", 
               size = 6, color = "gray50") +
      theme_void() +
      labs(title = title_label)
    return(list(plot = p, height = 6))
  }
  
  df <- prepare_dotplot_data(df)
  
  # Unified GO terms from available data
  all_go_terms <- df %>%
    group_by(Description_short) %>%
    summarise(max_sig = max(neg_log10_padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_sig)) %>%
    pull(Description_short)
  
  if (length(all_go_terms) == 0) all_go_terms <- "No significant terms"
  
  # === FORCE ALL 3 COMPARISONS ===
  full_grid <- expand.grid(
    Description_short = all_go_terms,
    Comparison = all_comparisons,
    stringsAsFactors = FALSE
  )
  
  full_df <- full_grid %>%
    left_join(df, by = c("Description_short", "Comparison"))
  
  full_df <- full_df %>%
    mutate(Description_short = factor(Description_short, levels = all_go_terms))
  
  colors <- if (direction == "Up") c("#FCBBA1", "#FC4E2A", "#800026") else
    c("#BDD7E7", "#2171B5", "#08306B")
  
  p <- ggplot(full_df, aes(x = GeneRatio_numeric, y = Description_short,
                           size = Count, color = neg_log10_padj)) +
    geom_point(alpha = 0.85, na.rm = TRUE) +
    facet_wrap(~ Comparison, ncol = 3) +
    scale_color_gradientn(colors = colors, name = "-log10(padj)") +
    scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
    scale_x_continuous(name = "Gene Ratio", labels = percent_format(accuracy = 1)) +
    labs(title = title_label, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.8),
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.y = element_text(size = 10.5, face = "bold"),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position = "right"
    )
  
  n_terms <- length(all_go_terms)
  fig_h <- max(10, n_terms * 0.42 + 6)
  
  return(list(plot = p, height = fig_h))
}

# ====================== SAVE FUNCTION (Updated) ======================
save_go_plot <- function(df, direction, title_label, outfile, 
                         all_comparisons, width = 21) {
  res <- make_go_dotplot(df, direction, title_label, all_comparisons)
  ggsave(outfile, res$plot, width = width, height = res$height, dpi = 300, bg = "white")
  cat("✓ Saved:", basename(outfile), "\n")
}

# ====================== DEFINE COMPARISONS ======================

day2_s42 <- c("S42 HC D2 vs C", "S42 LC D2 vs C", "S42 LC D2 vs HC")
day6_s42 <- c("S42 HC D6 vs C", "S42 LC D6 vs C", "S42 LC D6 vs HC")

day2_s399 <- c("S399 HC D2 vs C", "S399 LC D2 vs C", "S399 LC D2 vs HC")
day6_s399 <- c("S399 HC D6 vs C", "S399 LC D6 vs C", "S399 LC D6 vs HC")

# ====================== GENERATE PLOTS ======================

# S42 - Day 2
save_go_plot(go_up_day2_s42,   "Up",   "GO Enrichment — Up-Regulated (S42 - Day 2)", 
             paste0(out_go, "dotplot_GO_UP_S42_Day2.png"), day2_s42)

save_go_plot(go_down_day2_s42, "Down", "GO Enrichment — Down-Regulated (S42 - Day 2)", 
             paste0(out_go, "dotplot_GO_DOWN_S42_Day2.png"), day2_s42)

# S42 - Day 6
save_go_plot(go_up_day6_s42,   "Up",   "GO Enrichment — Up-Regulated (S42 - Day 6)", 
             paste0(out_go, "dotplot_GO_UP_S42_Day6.png"), day6_s42)

save_go_plot(go_down_day6_s42, "Down", "GO Enrichment — Down-Regulated (S42 - Day 6)", 
             paste0(out_go, "dotplot_GO_DOWN_S42_Day6.png"), day6_s42)

# S399 - Day 2
save_go_plot(go_up_day2_s399,   "Up",   "GO Enrichment — Up-Regulated (S399 - Day 2)", 
             paste0(out_go, "dotplot_GO_UP_S399_Day2.png"), day2_s399)

save_go_plot(go_down_day2_s399, "Down", "GO Enrichment — Down-Regulated (S399 - Day 2)", 
             paste0(out_go, "dotplot_GO_DOWN_S399_Day2.png"), day2_s399)

# S399 - Day 6
save_go_plot(go_up_day6_s399,   "Up",   "GO Enrichment — Up-Regulated (S399 - Day 6)", 
             paste0(out_go, "dotplot_GO_UP_S399_Day6.png"), day6_s399)

save_go_plot(go_down_day6_s399, "Down", "GO Enrichment — Down-Regulated (S399 - Day 6)", 
             paste0(out_go, "dotplot_GO_DOWN_S399_Day6.png"), day6_s399)


# PART 7 — KEGG ENRICHMENT DOT PLOTS
# Empty files (no significant pathways) produce a blank placeholder plot.


# ---------- Load one KEGG file — direction from filename ----------
load_kegg <- function(file, comparison, direction) {
  path <- paste0(kegg_path, file)
  if (!file.exists(path)) {
    message("  [NOT FOUND] ", file)
    return(data.frame(Comparison = comparison, Direction = direction))
  }
  df <- tryCatch(read_csv_file(path), error = function(e) data.frame())
  if (nrow(df) == 0) {
    message("  [EMPTY - nothing significant] ", file)
    return(data.frame(Comparison = comparison, Direction = direction))
  }
  df$Comparison <- comparison
  df$Direction  <- direction
  return(df)
}

# ---------- All KEGG file pairs (up + down per comparison) ----------
kegg_files <- list(
  list(comp = "HC vs Control",   up = "HighCalorievsControl.up_KEGGenrich_significant.xls",    down = "HighCalorievsControl.down_KEGGenrich_significant.xls"),
  list(comp = "LC vs Control",   up = "LowCalorievsControl.up_KEGGenrich_significant.xls",     down = "LowCalorievsControl.down_KEGGenrich_significant.xls"),
  list(comp = "LC vs HC",        up = "LowCalorievsHighCalorie.up_KEGGenrich_significant.xls", down = "LowCalorievsHighCalorie.down_KEGGenrich_significant.xls"),
  list(comp = "S42 HC D2 vs C",  up = "S42_HC_D2vsS42_C_D2.up_KEGGenrich_significant.xls",    down = "S42_HC_D2vsS42_C_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S42 HC D6 vs C",  up = "S42_HC_D6vsS42_C_D6.up_KEGGenrich_significant.xls",    down = "S42_HC_D6vsS42_C_D6.down_KEGGenrich_significant.xls"),
  list(comp = "S42 LC D2 vs C",  up = "S42_LC_D2vsS42_C_D2.up_KEGGenrich_significant.xls",    down = "S42_LC_D2vsS42_C_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S42 LC D2 vs HC", up = "S42_LC_D2vsS42_HC_D2.up_KEGGenrich_significant.xls",   down = "S42_LC_D2vsS42_HC_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S42 LC D6 vs C",  up = "S42_LC_D6vsS42_C_D6.up_KEGGenrich_significant.xls",    down = "S42_LC_D6vsS42_C_D6.down_KEGGenrich_significant.xls"),
  list(comp = "S42 LC D6 vs HC", up = "S42_LC_D6vsS42_HC_D6.up_KEGGenrich_significant.xls",   down = "S42_LC_D6vsS42_HC_D6.down_KEGGenrich_significant.xls"),
  list(comp = "S399 HC D2 vs C", up = "S399_HC_D2vsS399_C_D2.up_KEGGenrich_significant.xls",  down = "S399_HC_D2vsS399_C_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S399 HC D6 vs C", up = "S399_HC_D6vsS399_C_D6.up_KEGGenrich_significant.xls",  down = "S399_HC_D6vsS399_C_D6.down_KEGGenrich_significant.xls"),
  list(comp = "S399 LC D2 vs C", up = "S399_LC_D2vsS399_C_D2.up_KEGGenrich_significant.xls",  down = "S399_LC_D2vsS399_C_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S399 LC D2 vs HC",up = "S399_LC_D2vsS399_HC_D2.up_KEGGenrich_significant.xls", down = "S399_LC_D2vsS399_HC_D2.down_KEGGenrich_significant.xls"),
  list(comp = "S399 LC D6 vs C", up = "S399_LC_D6vsS399_C_D6.up_KEGGenrich_significant.xls",  down = "S399_LC_D6vsS399_C_D6.down_KEGGenrich_significant.xls"),
  list(comp = "S399 LC D6 vs HC",up = "S399_LC_D6vsS399_HC_D6.up_KEGGenrich_significant.xls", down = "S399_LC_D6vsS399_HC_D6.down_KEGGenrich_significant.xls")
)

# Load all files and combine
kegg_combined <- bind_rows(lapply(kegg_files, function(x) {
  bind_rows(
    load_kegg(x$up,   x$comp, "Up-Regulated"),
    load_kegg(x$down, x$comp, "Down-Regulated")
  )
}))

# Ensure Direction is a factor
kegg_combined$Direction <- factor(kegg_combined$Direction,
                                  levels = c("Up-Regulated", "Down-Regulated"))

# ---- KEGG prep function ----
prepare_and_order_kegg <- function(df) {
  df <- df %>%
    filter(!is.na(padj), !is.na(GeneRatio), !is.na(Description)) %>%
    mutate(
      neg_log10_padj    = -log10(padj),
      GeneRatio_numeric = sapply(
        strsplit(as.character(GeneRatio), "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      ),
      Description_short = ifelse(
        nchar(Description) > 50,
        paste0(substr(Description, 1, 47), "..."),
        Description
      )
    )
  term_order <- df %>%
    group_by(Description_short) %>%
    summarise(max_sig = max(neg_log10_padj, na.rm = TRUE), .groups = "drop") %>%
    arrange(max_sig) %>%
    pull(Description_short)
  df %>% mutate(Description_short = factor(Description_short, levels = term_order))
}

# ---- Empty placeholder plot ----
empty_kegg_plot <- function(title_label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "No significant KEGG pathways",
             size = 6, hjust = 0.5, vjust = 0.5, color = "gray50") +
    theme_void() +
    labs(title = title_label) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13))
}

# ---- KEGG dot plot function ----
make_kegg_dotplot <- function(df, title_label, outfile, fig_width = 16) {

  # Keep only rows with actual data (padj not NA = has enrichment results)
  df_data <- df %>% filter(!is.na(padj))

  # All comparisons and directions present in this subset
  all_comps <- unique(df$Comparison)
  all_dirs  <- c("Up-Regulated", "Down-Regulated")

  fig_h <- max(8, nrow(df_data) * 0.25 + 4)

  if (nrow(df_data) == 0) {
    # Fully empty — save placeholder
    p <- empty_kegg_plot(title_label)
    ggsave(outfile, p, width = fig_width, height = 4, dpi = 300, bg = "white")
    cat("✓ Saved (no significant pathways):", basename(outfile), "\n")
    return(invisible(NULL))
  }

  df_data <- prepare_and_order_kegg(df_data) %>%
    mutate(
      Comparison = factor(Comparison, levels = all_comps),
      Direction  = factor(Direction,  levels = all_dirs)
    )

  p <- ggplot() +
    geom_point(
      data = df_data %>% filter(Direction == "Up-Regulated"),
      aes(x = GeneRatio_numeric, y = Description_short,
          size = Count, color = neg_log10_padj), alpha = 1
    ) +
    scale_color_gradientn(
      colors = c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15"),
      name   = "-log10(padj) (Up)", guide = guide_colorbar(order = 1)
    ) +
    ggnewscale::new_scale_color() +
    geom_point(
      data = df_data %>% filter(Direction == "Down-Regulated"),
      aes(x = GeneRatio_numeric, y = Description_short,
          size = Count, color = neg_log10_padj), alpha = 1
    ) +
    scale_color_gradientn(
      colors = c("#B0E2FF", "#87CEFA", "#3182BD", "#27408B"),
      name   = "-log10(padj) (Down)", guide = guide_colorbar(order = 2)
    ) +
    facet_grid(Direction ~ Comparison, scales = "free_y", space = "free_y") +
    scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
    scale_x_continuous(name = "Gene Ratio", labels = percent_format(accuracy = 1)) +
    labs(title    = title_label,
         subtitle = "Dot size = Gene Count  •  Color = −log10(padj)",
         y        = "KEGG Pathway") +
    theme_minimal(base_size = 11) +
    theme(
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.background = element_rect(color = "black", fill = "grey95", linewidth = 0.8),
      strip.text.x     = element_text(size = 9,  face = "bold"),
      strip.text.y     = element_text(size = 10, face = "bold"),
      axis.text.y      = element_text(size = 8,  face = "bold", color = "black"),
      axis.text.x      = element_text(size = 8,  face = "bold"),
      axis.title       = element_text(size = 11, face = "bold"),
      plot.title       = element_text(size = 13, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 9,  hjust = 0.5, color = "grey40"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.title     = element_text(face = "bold"),
      legend.position  = "right",
      panel.spacing.y  = unit(0.8, "cm"),
      panel.spacing.x  = unit(0.5, "cm")
    )

  # Colour strip labels: Up = red, Down = blue
  g <- ggplotGrob(p)
  strip_r_idx <- which(grepl("strip-r", g$layout$name))
  for (i in seq_along(strip_r_idx)) {
    strip_color <- if (i == 1) "#C0392B" else "#2980B9"
    strip_grob  <- g$grobs[[strip_r_idx[i]]]
    for (j in seq_along(strip_grob$grobs)) {
      inner <- strip_grob$grobs[[j]]
      for (k in seq_along(inner$children)) {
        child <- inner$children[[k]]
        if (inherits(child, "rect"))
          inner$children[[k]]$gp <- gpar(fill = strip_color, col = strip_color)
        if (inherits(child, "text") || inherits(child, "titleGrob"))
          inner$children[[k]]$gp <- modifyList(
            as.list(inner$children[[k]]$gp), list(col = "white")) %>% do.call(gpar, .)
      }
      strip_grob$grobs[[j]] <- inner
    }
    g$grobs[[strip_r_idx[i]]] <- strip_grob
  }

  png(outfile, width = fig_width, height = fig_h, units = "in", res = 300, bg = "white")
  grid.draw(g)
  dev.off()
  cat("✓ Saved:", basename(outfile), "\n")
}

# ---- Overall KEGG (HC vs Ctrl, LC vs Ctrl, LC vs HC) ----
kegg_overall <- kegg_combined %>%
  filter(Comparison %in% c("HC vs Control", "LC vs Control", "LC vs HC"))
make_kegg_dotplot(kegg_overall, "KEGG Enrichment — Overall Diet Comparisons",
                  paste0(out_kegg, "dotplot_KEGG_overall.png"), fig_width = 14)

# ---- S42 KEGG ----
kegg_S42 <- kegg_combined %>% filter(grepl("^S42", Comparison))
make_kegg_dotplot(kegg_S42, "KEGG Enrichment — S42 (All Days & Comparisons)",
                  paste0(out_kegg, "dotplot_KEGG_S42.png"), fig_width = 22)

# ---- S399 KEGG ----
kegg_S399 <- kegg_combined %>% filter(grepl("^S399", Comparison))
make_kegg_dotplot(kegg_S399, "KEGG Enrichment — S399 (All Days & Comparisons)",
                  paste0(out_kegg, "dotplot_KEGG_S399.png"), fig_width = 22)

cat("\n============================================================\n")
cat("✓ ALL ANALYSES COMPLETE\n")
cat("Outputs written to:\n")
cat("  DEG  →", out_deg,  "\n")
cat("  GO   →", out_go,   "\n")
cat("  KEGG →", out_kegg, "\n")
cat("============================================================\n")

