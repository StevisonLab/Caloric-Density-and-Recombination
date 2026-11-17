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

create_slope_plot <- function(coef_df, trait_label, filename) {
  
  df_plot <- coef_df %>%
    mutate(
      Brood = factor(as.character(Brood),
                     levels = LETTERS[1:7]),
      
      Treatment = factor(as.character(Treatment),
                         levels = c("0.5x", "1.0x", "2.0x")),
      
      conf.low = estimate - qt(0.975, df) * std.error,
      conf.high = estimate + qt(0.975, df) * std.error,
      
      sig_symbol = ifelse(p.value < ALPHA, "*", "")
    )
  
  pd <- position_dodge(width = 0.4)
  
  p <- ggplot(df_plot, aes(
    x = Brood,
    y = estimate,
    color = Genotype,
    group = Genotype
  )) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "gray50"
    ) +
    
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0.2,
      position = pd
    ) +
    
    geom_line(
      position = pd,
      linewidth = 1.2
    ) +
    
    geom_point(
      position = pd,
      size = 3
    ) +
    
    # Significance asterisks
    geom_text(
      aes(label = sig_symbol, y = conf.high + 0.05),
      position = pd,
      size = 6,
      show.legend = FALSE
    ) +
    
    facet_wrap(~ Treatment, nrow = 1) +
    
    # Brood labels
    scale_x_discrete(
      labels = setNames(
        paste0(LETTERS[1:7], " (", 1:7, ")"),
        LETTERS[1:7]
      )
    ) +
    
    # Genotype colors and labels
    scale_color_manual(
      values = genotype_colors,
      breaks = c("42", "399"),
      labels = c("DGRP_42", "DGRP_399")
    ) +
    scale_y_continuous(
      breaks = seq(-1.5, 1.5, by = 0.5)
    ) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    
    labs(
      title = paste("Recombination Rate vs.", trait_label),
      x = "Brood",
      y = "Regression slope (β)",
      color = "Genotype"
    ) +
    
    theme_minimal(base_size = 18) +
    
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "#f0f0f0"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
  
  ggsave(
    filename,
    plot = p,
    width = 15,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

# ============================================================
# GENERATE FIGURE 7b: SURVIVAL
# ============================================================

p_viability_slopes <- create_slope_plot(
  viability_coefficients,
  "Egg-to-Adult Survival",
  "46_RR_SURVIVAL.png"
)

p_viability_slopes

# ============================================================
# GENERATE FECUNDITY FIGURE
# ============================================================

p_fecundity_slopes <- create_slope_plot(
  fecundity_coefficients,
  "Maternal Fecundity",
  "MarkerPresent_RR_vs_Fecundity_Slopes.png"
)
