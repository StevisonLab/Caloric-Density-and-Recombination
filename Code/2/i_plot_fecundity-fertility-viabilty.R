setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/2")

# 1. Read and prepare main data
marker_free_viability_dataset <- read.csv("F_marker_free_viability_dataset.csv")

multiply_marked_viability_dataset<- read.csv("6_multiply_marked_viability_dataset.csv")
# mothers per vial for marker-present
mothers_per_vial <- read.csv("2_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))

# ============================================================
# LIFE-HISTORY TRAITS + ANOVA (with Data_collector)
# Marker-free & Marker-present (3 mothers per vial)
# ============================================================

library(dplyr)
library(ggplot2)

MOTHERS_PER_VIAL <- 3
ALPHA <- 0.05
traits <- c("Fecundity", "Fertility", "Viability")

# -----------------------------
# MARKER-FREE
# -----------------------------
marker_free_life_history <- marker_free_viability_dataset %>%
  mutate(
    Dataset   = "Marker-free",
    Fecundity = transformed_mean_marker_free_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count /MOTHERS_PER_VIAL ,
    Viability = pmin(pmax(derived_marker_free_total_viability * 100, 0), 100)
  ) %>%
  select(
    Data_collector = data_collector,
    Treatment = treatment,
    Genotype  = genotype,
    Brood     = brood_id,
    Dataset,
    Fecundity,
    Fertility,
    Viability
  ) %>%
  mutate(
    Genotype = as.character(Genotype),
    Treatment = as.character(Treatment),
    Brood = as.character(Brood),
    Data_collector = as.character(Data_collector))

# -----------------------------
# MARKER-PRESENT
# -----------------------------
marker_present_life_history <- multiply_marked_viability_dataset %>%
  mutate(
    Dataset   = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / female_in_vial,
    Fertility = total_adult_count / female_in_vial,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(
    Data_collector = data_collector,
    Treatment = treatment,
    Genotype  = genotype,
    Brood     = brood_id,
    Dataset,
    Fecundity,
    Fertility,
    Viability
  )  %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Brood     = factor(Brood),
    Dataset   = factor(Dataset),
    Data_collector = factor(Data_collector))
# -----------------------------
# COMBINE DATASETS
# -----------------------------

life_history_all <- bind_rows(
  marker_free_life_history,
  marker_present_life_history
) %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Brood     = factor(Brood),
    Dataset   = factor(Dataset),
    Data_collector = factor(Data_collector)
  )

# ============================================================
# PLOTTING FUNCTION
# ============================================================
plot_trait <- function(trait_name, y_label) {
  # Filter to only Marker-present
  life_history_all <- life_history_all %>%
    dplyr::filter(Dataset == "Marker-present")  # <- only these rows
  
  ggplot(life_history_all,
         aes(x = Treatment, y = .data[[trait_name]], fill = Genotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
    # geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
    facet_wrap(~ Dataset) +
    theme_classic(base_size = 24) +
    labs(
      title = trait_name,
      y = y_label,
      x = "Treatment"
    )+scale_fill_manual(values = c(
      "217" = "#8dd3c7",
      "391" = "#ffffb3",
      "491" = "#80b1d3",
      "399" = "#9e9ac8",
      "42"  = "#fb8072"
    ))
}
#reaction norm

# Define which terms are significant (TRUE = bold)
sig_B   <- TRUE
sig_G   <- TRUE
sig_E   <- TRUE
sig_GxE <- TRUE

# Helper: wrap term in bold() if significant
fmt_term <- function(term, significant) {if (significant) paste0("bold(", term, ")") else term}

# Build the label as a string expression
eq_label <- paste0( "Viability * '~' * ", fmt_term("B", sig_B), " + ", fmt_term("G", sig_G), " + ", fmt_term("E", sig_E), " + ", fmt_term("GxE", sig_GxE))

plot_trait <- function(trait_name, y_label) {
  # Filter to only Marker-present
  life_history_all_filtered <- life_history_all %>%
    dplyr::filter(Dataset == "Marker-present")
  
  ggplot(life_history_all_filtered,
         aes(
           x = Treatment,
           y = .data[[trait_name]],
           group = Genotype,
           color = Genotype
         )) +
    stat_summary(fun = mean, geom = "line", linewidth = 3) +
    stat_summary(fun = mean, geom = "point", size = 6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.55) +
    ylim(45,100)+
    # Add equation 
    annotate(
      "text",
      x = 2,          # adjust if needed
      y = 98.5,
      label = eq_label,
      parse = TRUE,
      size = 16
    ) +
    facet_wrap(~ Dataset) +
    theme_classic(base_size = 24) +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    ) +
    labs(
      #title = trait_name,
      y = y_label,
      x = "Treatment"
    ) +
    scale_color_manual(values = c(
      "217" = "#e41a1c",  # red
      "391" = "#377eb8",  # blue
      "491" = "#4daf4a",  # green
      "399" = "#984ea3",  # purple
      "42"  = "#ff7f00"   # orange
    ))
}
# -----------------------------
# GENERATE PLOTS
# -----------------------------
p_fecundity <- plot_trait("Fecundity", "Eggs per mother")
p_fertility <- plot_trait("Fertility", "Viable offspring")
p_viability <- plot_trait("Viability", "Viability (%)")
p_viability
ggsave(
  "../../Summaries and Outputs/2/18_trait_fecundity.png",
  plot = p_fecundity,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/2/19_trait_fertility.png",
  plot = p_fertility,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/2/20_trait_viability.png",
  plot = p_viability,
  width = 13,
  height = 7,
  dpi = 300
)


# ============================================================
# ANOVA — SEPARATE BY DATASET AND TRAIT (NO PLOTTING)
# ============================================================
run_anova <- function(df, trait) {
  aov(
    as.formula(paste(trait, "~ Treatment * Genotype * Brood + Data_collector")),
    data = df
  )
}

for (dataset_name in unique(as.character(life_history_all$Dataset))) {
  
  df_subset <- life_history_all %>%
    filter(Dataset == dataset_name)
  
  for (trait in traits) {
    
    # remove NA/NaN/Inf values per trait
    df_trait <- df_subset %>%
      filter(is.finite(.data[[trait]]))
    
    # only run ANOVA if there is enough data
    if (nrow(df_trait) > 1) {
      model <- run_anova(df_trait, trait)
      
      anova_table <- as.data.frame(summary(model)[[1]])
      anova_table$Term <- rownames(anova_table)
      anova_table$Dataset <- dataset_name
      anova_table$Trait <- trait
      
      rownames(anova_table) <- NULL
      
      write.csv(
        anova_table,
        paste0(
          "../../Summaries and Outputs/2/21_ANOVA_",
          dataset_name, "_",
          trait,
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# ============================================================
# SAVE CLEANED DATASET
# ============================================================
write.csv(
  life_history_all,
  "../../Summaries and Outputs/2/22_Life_History_All_Cleaned.csv",
  row.names = FALSE
)

# ============================================================
# SUMMARY MEANS BY DATASET × TREATMENT × GENOTYPE
# ============================================================
life_history_summary <- life_history_all %>%
  group_by(Dataset, Treatment, Genotype) %>%
  summarise(
    Fecundity = mean(Fecundity, na.rm = TRUE),
    Fertility = mean(Fertility, na.rm = TRUE),
    Viability = mean(Viability, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  life_history_summary,
  "../../Summaries and Outputs/2/23_life_history_summary_by_dataset.csv",
  row.names = FALSE
)

