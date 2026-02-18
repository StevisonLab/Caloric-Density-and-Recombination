setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

# 1. Read and prepare main data
marker_free_viability_dataset <- read.csv("F_marker_free_viability_dataset.csv")

multiply_marked_viability_dataset<- read.csv("6_multiply_marked_viability_dataset.csv")


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
    Fertility = total_adult_count / MOTHERS_PER_VIAL ,
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
    Fecundity = transformed_mean_multiply_marked_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count /MOTHERS_PER_VIAL ,
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
# REMOVE OUTLIERS FUNCTION
# ============================================================
#remove_outliers <- function(df, trait) {
  df %>%
    group_by(Dataset, Treatment, Genotype, Data_collector) %>%
    mutate(
      Q1  = quantile(.data[[trait]], 0.25, na.rm = TRUE),
      Q3  = quantile(.data[[trait]], 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      Lower = Q1 - 1.5 * IQR,
      Upper = Q3 + 1.5 * IQR,
      Is_Outlier = .data[[trait]] < Lower | .data[[trait]] > Upper
    ) %>%
    ungroup() %>%
    filter(!Is_Outlier) %>%
    select(-Q1, -Q3, -IQR, -Lower, -Upper, -Is_Outlier)
}

# Apply outlier removal for each trait sequentially
#for(trait in traits) {
  life_history_all <- remove_outliers(life_history_all, trait)
}

# ============================================================
# PLOTTING FUNCTION
# ============================================================
plot_trait <- function(trait_name, y_label) {
  ggplot(life_history_all,
         aes(x = Treatment, y = .data[[trait_name]], fill = Genotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
   # geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
    facet_wrap(~ Dataset) +
    theme_classic(base_size = 14) +
    labs(
      title = trait_name,
      y = y_label,
      x = "Treatment") +
    scale_fill_manual(values = c(
          "42"  = "#fb6a4a",  # orange
          "399" = "#9e9ac8"   # purple
        ))
    
}

# -----------------------------
# GENERATE PLOTS
# -----------------------------
p_fecundity <- plot_trait("Fecundity", "Eggs per mother")
p_fertility <- plot_trait("Fertility", "Viable offspring")
p_viability <- plot_trait("Viability", "Viability (%)")

ggsave(
  "../../Summaries and Outputs/3/20_trait_fecundity.png",
  plot = p_fecundity,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/3/21_trait_fertility.png",
  plot = p_fertility,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "../../Summaries and Outputs/3/22_trait_viability.png",
  plot = p_viability,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================
# ANOVA
# ============================================================
anova_results <- list()

run_anova <- function(df, trait) {
  aov(
    as.formula(paste(trait, "~ Treatment * Genotype * Brood + Data_collector")),
   # as.formula(paste(trait, "~ Treatment * Genotype * Brood ")),
    
    data = df
  )
}

has_significant_effect <- function(aov_model) {
  pvals <- summary(aov_model)[[1]][["Pr(>F)"]]
  any(pvals < ALPHA, na.rm = TRUE)
}

plot_if_significant <- function(df, trait, dataset_label) {
  
  model <- run_anova(df, trait)
  anova_results[[paste(dataset_label, trait, sep = "_")]] <<- summary(model)
  
  if (has_significant_effect(model)) {
    print(
      ggplot(df,
             aes(
               x = Data_collector,
               y = .data[[trait]],
               color = Genotype,
               fill = Treatment
             )) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
      #  geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +
        facet_wrap(~ Dataset) +
        theme_classic(base_size = 14) +
        labs(
          title = paste(dataset_label, trait, "(Significant ANOVA)"),
          x = "Data Collector",
          y = ifelse(trait == "Viability", "Viability (%)", trait)
        )
    )
  }
}

# -----------------------------
# RUN ANOVAS SEPARATELY
# -----------------------------
for (dataset_name in levels(life_history_all$Dataset)) {
  
  df_subset <- life_history_all %>% filter(Dataset == dataset_name)
  
  for (trait in traits) {
    plot_if_significant(df_subset, trait, dataset_name)
  }
}

# -----------------------------
# SAVE ALL ANOVA RESULTS
# -----------------------------
sink("../../Summaries and Outputs/3/23_life_history_ANOVA_results.csv")

for (name in names(anova_results)) {
  cat("\n=============================\n")
  cat(name, "\n")
  cat("=============================\n")
  print(anova_results[[name]])
}

sink()

# -----------------------------
# SAVE CLEANED DATASET
# -----------------------------
write.csv(life_history_all, "../../Summaries and Outputs/3/24_Life_History_All_Cleaned.csv", row.names = FALSE)

#mean per genotype x treatment
life_history_summary <- life_history_all %>%
  group_by(Dataset, Treatment, Genotype) %>%
  summarise(
    n = n(),
    mean_Fecundity = mean(Fecundity, na.rm = TRUE),
    se_Fecundity   = sd(Fecundity, na.rm = TRUE) / sqrt(n),
    
    mean_Fertility = mean(Fertility, na.rm = TRUE),
    se_Fertility   = sd(Fertility, na.rm = TRUE) / sqrt(n),
    
    mean_Viability = mean(Viability, na.rm = TRUE),
    se_Viability   = sd(Viability, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )
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
  "../../Summaries and Outputs/3/25_life_history_summary_by_dataset.csv",
  row.names = FALSE
)
 