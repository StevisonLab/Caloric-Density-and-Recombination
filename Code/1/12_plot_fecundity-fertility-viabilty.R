setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/1")

# 1. Read and prepare main data
marker_free_viability_dataset <- read.csv("F_marker_free_viability_dataset.csv")

multiply_marked_viability_dataset<- read.csv("6_multiply_marked_viability_dataset.csv")
# mothers per vial for marker-present
mothers_per_vial <- read.csv("1_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))

# ============================================================
# LIFE-HISTORY TRAITS + ANOVA (with Data_collector)
# Marker-free & Marker-present (3 mothers per vial)
# ============================================================

library(dplyr)
library(ggplot2)

MOTHERS_PER_VIAL <- 5
ALPHA <- 0.05
traits <- c("Fecundity", "Fertility", "Viability")

# -----------------------------
# MARKER-FREE
# -----------------------------
marker_free_life_history <- marker_free_viability_dataset %>%
  mutate(
    Dataset   = "Marker-free",
    Fecundity = transformed_mean_marker_free_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count /MOTHERS_PER_VIAL,
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
    Fertility = total_adult_count/female_in_vial,
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
  ggplot(
    life_history_all,
    aes(
      x = Treatment,
      y = .data[[trait_name]],
      fill = Treatment
    )
  ) +
    geom_boxplot(
      alpha = 0.7,
      outlier.shape = NA,
      width = 0.6
    ) +
    facet_wrap(~ Dataset) +
    theme_classic(base_size = 24) +
    labs(
      title = trait_name,
      y = y_label,
      x = "Treatment"
    ) +
    scale_fill_manual(values = c(
      "0.5x"  = "#dadaeb",
      "0.75x" = "#bcbddc",
      "1.0x"  = "#9e9ac8",
      "1.25x" = "#807dba",
      "1.5x"  = "#6a51a3",
      "2.0x"  = "#4a1486"
    ))
}

# -----------------------------
# GENERATE PLOTS
# -----------------------------
p_fecundity <- plot_trait("Fecundity", "Eggs per mother")
p_fertility <- plot_trait("Fertility", "Viable offspring per mother")
p_viability <- plot_trait("Viability", "Viability (%)")

ggsave(
  "../../Summaries and Outputs/1/20_trait_fecundity.png",
  plot = p_fecundity,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/1/21_trait_fertility.png",
  plot = p_fertility,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/1/22_trait_viability.png",
  plot = p_viability,
  width = 12,
  height = 6,
  dpi = 300
)


# ============================================================
# ANOVA — SEPARATE BY DATASET AND TRAIT (NO PLOTTING)
# ============================================================
run_anova <- function(df, trait) {
  aov(
    as.formula(paste(trait, "~ Treatment * Brood ")),
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
          "../../Summaries and Outputs/1/23_ANOVA_",
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
  "../../Summaries and Outputs/1/24_Life_History_All_Cleaned.csv",
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
  "../../Summaries and Outputs/1/25_life_history_summary_by_dataset.csv",
  row.names = FALSE
)

