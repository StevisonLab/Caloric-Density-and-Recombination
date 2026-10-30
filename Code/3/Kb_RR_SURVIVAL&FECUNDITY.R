# ============================================================
# LIBRARIES & SETTINGS
# ============================================================
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(lme4)
library(lmerTest)
library(broom.mixed)

# Set working directory
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")
options(contrasts = c("contr.sum", "contr.poly"))

ALPHA <- 0.05

# Define standard genotype color scheme (Orange = 42, Purple = 399)
genotype_colors <- c("42" = "#ff7f00", "399" = "#984ea3")

# ============================================================
# STEP 1: LOAD & PROCESS MARKER-PRESENT LIFE HISTORY DATA
# ============================================================
# Load raw marker-present dataset and mother count file
multiply_marked_viability_dataset <- read.csv("6_multiply_marked_viability_dataset.csv")
mothers_per_vial <- read.csv("3_mom.csv")

# Merge mothers per vial to calculate per-mother metrics
marker_present_life_history <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id")) %>%
  mutate(
    Dataset = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / female_in_vial,
    Fertility = total_adult_count / female_in_vial,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(
    VialID = cross_id,
    Data_collector = data_collector,
    Treatment = treatment,
    Genotype = genotype,
    Brood = brood_id,
    Dataset,
    Fecundity,
    Fertility,
    Viability
  ) %>%
  filter(is.finite(Fecundity), is.finite(Viability))

# Ensure factors are set
life_history_all <- marker_present_life_history %>%
  mutate(across(c(Treatment, Genotype, Brood, Dataset, Data_collector), factor))

# ============================================================
# STEP 2: LOAD & MERGE RECOMBINATION DATA
# ============================================================
df <- read.csv("23_df_RR_cleaned.csv")

recomb_vial <- df %>%
  rename(
    RR = recombination_fractions,
    Genotype = GB,
    Treatment = Diet,
    Brood = Brood,
    VialID = vialID
  ) %>%
  mutate(across(c(Treatment, Genotype, Brood), factor)) %>%
  select(VialID, Brood, Treatment, Genotype, RR)

# Convert VialIDs to character to ensure a smooth join
recomb_vial <- recomb_vial %>% mutate(VialID = as.character(VialID))
life_history_all <- life_history_all %>% mutate(VialID = as.character(VialID))

# Merge Recombination and Life History data by VialID, Treatment, Genotype, and Brood
combined_vial_data <- inner_join(
  recomb_vial,
  life_history_all,
  by = c("VialID", "Brood", "Treatment", "Genotype")
) %>%
  mutate(across(c(Treatment, Genotype, Brood, Data_collector, VialID), factor))

# Save cleaned combined dataset
write.csv(combined_vial_data, "42_MarkerPresent_Combined_RR_LifeHistory.csv", row.names = FALSE)

# ============================================================
# STEP 3: SUBGROUP GLMMs (RR ~ Viability) & (RR ~ Fecundity)
# ============================================================

# --- 3A: Viability GLMMs ---
viability_glmm_results <- combined_vial_data %>%
  group_by(Treatment, Genotype, Brood) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      lmer(RR ~ Viability + (1 | Data_collector), data = .x),
      error = function(e) NULL
    )),
    coefs = map(model, ~ if(!is.null(.x)) tidy(.x, effects = "fixed") else NULL),
    fit_stats = map(model, ~ if(!is.null(.x)) glance(.x) else NULL)
  )

viability_coefficients <- viability_glmm_results %>%
  select(Treatment, Genotype, Brood, coefs) %>%
  unnest(coefs) %>%
  filter(term == "Viability")

viability_fit_stats <- viability_glmm_results %>%
  select(Treatment, Genotype, Brood, fit_stats) %>%
  unnest(fit_stats)

# Save Viability Results
write.csv(viability_coefficients, "51_MarkerPresent_Viability_GLMM_Coefficients.csv", row.names = FALSE)
write.csv(viability_fit_stats, "52_MarkerPresent_Viability_GLMM_Fit_Statistics.csv", row.names = FALSE)

# --- 3B: Fecundity GLMMs ---
fecundity_glmm_results <- combined_vial_data %>%
  group_by(Treatment, Genotype, Brood) %>%
  nest() %>%
  mutate(
    model = map(data, ~ tryCatch(
      lmer(RR ~ Fecundity + (1 | Data_collector), data = .x),
      error = function(e) NULL
    )),
    coefs = map(model, ~ if(!is.null(.x)) tidy(.x, effects = "fixed") else NULL),
    fit_stats = map(model, ~ if(!is.null(.x)) glance(.x) else NULL)
  )

fecundity_coefficients <- fecundity_glmm_results %>%
  select(Treatment, Genotype, Brood, coefs) %>%
  unnest(coefs) %>%
  filter(term == "Fecundity")

fecundity_fit_stats <- fecundity_glmm_results %>%
  select(Treatment, Genotype, Brood, fit_stats) %>%
  unnest(fit_stats)

# Save Fecundity Results
write.csv(fecundity_coefficients, "51_MarkerPresent_Fecundity_GLMM_Coefficients.csv", row.names = FALSE)
write.csv(fecundity_fit_stats, "52_MarkerPresent_Fecundity_GLMM_Fit_Statistics.csv", row.names = FALSE)


# ============================================================
# STEP 4: PLOTTING
# ============================================================

# Helper Function for Slope Plots (Figure 1: Grouped by Genotype, Faceted by Treatment)
create_slope_plot <- function(coef_df, trait_label, filename) {
  df_plot <- coef_df %>%
    mutate(
      conf.low = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      sig_symbol = ifelse(p.value < ALPHA, "*", "")
    )
  
  p <- ggplot(df_plot, aes(x = Brood, y = estimate, color = Genotype, group = Genotype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 1) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, position = position_dodge(0.4), linewidth = 0.8) +
    geom_line(position = position_dodge(0.4), linewidth = 1.2) +
    geom_point(position = position_dodge(0.4), size = 3) +
    geom_text(aes(label = sig_symbol, y = conf.high + 0.05), position = position_dodge(0.4), size = 6, show.legend = FALSE) +
    facet_wrap(~ Treatment) +
    scale_color_manual(values = genotype_colors) +
    labs(
      title = paste("Subgroup Regression Slopes: RR vs.", trait_label, "(Marker-Present)"),
      subtitle = "Slopes (β) with 95% CIs across Broods",
      x = "Brood",
      y = paste("Effect of", trait_label, "on RR (β)"),
      color = "Genotype"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  
  ggsave(filename, plot = p, width = 12, height = 6, dpi = 300)
  return(p)
}

# Generate Slope Plots
p_viability_slopes <- create_slope_plot(viability_coefficients, "Viability (%)", "MarkerPresent_RR_vs_Viability_Slopes.png")
p_fecundity_slopes <- create_slope_plot(fecundity_coefficients, "Fecundity (Eggs/Mother)", "MarkerPresent_RR_vs_Fecundity_Slopes.png")


# Helper Function for Faceted Scatter Plots (Figure 2: Faceted by Treatment, Colored by Genotype)
create_scatter_plot <- function(data, trait_col, trait_label, filename) {
  recomb_summary <- data %>%
    group_by(Treatment, Genotype, Brood) %>%
    summarise(mean_RR = median(RR, na.rm = TRUE), .groups = "drop")
  
  lh_summary <- data %>%
    group_by(Treatment, Genotype, Brood) %>%
    summarise(mean_Trait = median(.data[[trait_col]], na.rm = TRUE), .groups = "drop")
  
  combined_group_data <- inner_join(recomb_summary, lh_summary, by = c("Brood", "Treatment", "Genotype"))
  
  p <- ggplot(combined_group_data, aes(x = mean_Trait, y = mean_RR, color = Genotype)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
    facet_wrap(~ Treatment) +
    scale_color_manual(values = genotype_colors) +
    labs(
      title = paste("Recombination Rate vs.", trait_label, "by Treatment & Genotype"),
      x = paste("Median", trait_label),
      y = "Median Recombination Rate",
      color = "Genotype"
    ) +
    theme_minimal(base_size = 18) +
    theme(
      legend.position = "right",
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  
  ggsave(filename, plot = p, width = 12, height = 6, dpi = 300)
  return(p)
}

# Generate Scatter Plots
p_viability_scatter <- create_scatter_plot(combined_vial_data, "Viability", "Viability (%)", "MarkerPresent_RR_vs_Viability_Scatter.png")
p_fecundity_scatter <- create_scatter_plot(combined_vial_data, "Fecundity", "Fecundity (Eggs/Mother)", "MarkerPresent_RR_vs_Fecundity_Scatter.png")


# ============================================================
# STEP 5: MANUSCRIPT TEXT EXPORT
# ============================================================
viability_text <- viability_coefficients %>%
  mutate(
    significance = ifelse(p.value < 0.05, "significant", "non-significant"),
    sentence = sprintf(
      "For Treatment %s, Genotype %s, Brood %s (Marker-Present): Viability had a %s effect on RR (β = %.4f ± %.4f, t = %.3f, p = %.4f).",
      Treatment, Genotype, Brood, significance, estimate, std.error, statistic, p.value
    )
  )

fecundity_text <- fecundity_coefficients %>%
  mutate(
    significance = ifelse(p.value < 0.05, "significant", "non-significant"),
    sentence = sprintf(
      "For Treatment %s, Genotype %s, Brood %s (Marker-Present): Fecundity had a %s effect on RR (β = %.4f ± %.4f, t = %.3f, p = %.4f).",
      Treatment, Genotype, Brood, significance, estimate, std.error, statistic, p.value
    )
  )

cat("\n================ MARKER-PRESENT VIABILITY RESULTS TEXT ================\n")
cat(paste(viability_text$sentence, collapse = "\n"))

cat("\n\n================ MARKER-PRESENT FECUNDITY RESULTS TEXT ================\n")
cat(paste(fecundity_text$sentence, collapse = "\n"))

# Save text outputs
writeLines(viability_text$sentence, "53_MarkerPresent_Viability_Results.txt")
writeLines(fecundity_text$sentence, "53_MarkerPresent_Fecundity_Results.txt")

