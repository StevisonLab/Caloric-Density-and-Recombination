# ============================================================
# STEP 1: LOAD DATA & CLEAN LIFE HISTORY
# ============================================================
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# Set working directory (uncomment if running locally)
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")
ALPHA <- 0.05
traits <- c("Fecundity", "Fertility", "Viability")

# Load Datasets
multiply_marked_viability_dataset <- read.csv("6_multiply_marked_viability_dataset.csv")
# mothers per vial for marker-present
mothers_per_vial <- read.csv("3_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))



df <- read.csv("23_df_RR_cleaned.csv")

# 1.2 Process Marker-Present
marker_present_life_history <- multiply_marked_viability_dataset %>%
  mutate(
    Dataset   = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / female_in_vial,
    Fertility = total_adult_count / female_in_vial,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(VialID = cross_id, Data_collector = data_collector, Treatment = treatment, 
         Genotype = genotype, Brood = brood_id, Dataset, Fecundity, Fertility, Viability)

marker_present_life_history <- marker_present_life_history %>%
  filter(is.finite(Fecundity))   # keeps only finite values

# 1.3 Combine and Outlier Removal
life_history_all <- bind_rows( marker_present_life_history) %>%
  mutate(across(c(Treatment, Genotype, Brood, Dataset, Data_collector), factor))


# ============================================================
# STEP 2: MERGE WITH RECOMBINATION DATA
# ============================================================

# Assuming 'df' is your recombination dataset already in memory
# Standardizing 'df' columns to match life_history_all
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

# 1. Convert VialID to character in the recombination data
recomb_vial <- recomb_vial %>%
  mutate(VialID = as.character(VialID))

# 2. Convert VialID to character in the life history data
life_history_all <- life_history_all %>%
  mutate(VialID = as.character(VialID))

# 3. Now the join will work perfectly
combined_vial_data <- inner_join(recomb_vial, life_history_all, 
                                 by = c("VialID", "Brood", "Treatment", "Genotype"))

# ============================================================
# STEP 3: CALCULATE MEANS AND VIAL-LEVEL SLOPES
# ============================================================
slopes_df <- combined_vial_data %>%
  filter(!is.na(RR) & !is.na(Viability)) %>%
  group_by(Dataset, Treatment, Genotype, Brood) %>%
  summarise(
    mean_RR = mean(RR, na.rm = TRUE),
    mean_Viability = mean(Viability, na.rm = TRUE),
    n_vials = n(),
    # Run model only if there's enough variation
    model = if(n() >= 4) list(lm(scale(RR) ~ scale(Viability))) else list(NULL),
    .groups = "drop"
  ) %>%
  filter(!map_lgl(model, is.null)) %>%
  mutate(
    slope = map_dbl(model, ~coef(.x)[2]),
    se = map_dbl(model, ~summary(.x)$coefficients[2, "Std. Error"]),
    p.value = map_dbl(model, ~summary(.x)$coefficients[2, "Pr(>|t|)"]),
    significant = p.value < 0.05
  )

# ============================================================
# STEP 4: GENERATE THE SLOPE FIGURE
# ============================================================
slopes_df <- slopes_df %>%
  mutate(Brood = as.numeric(factor(Brood)))

pd <- position_dodge(width = 0.5)

p5 <- ggplot(slopes_df, aes(
  x = Brood,
  y = slope,
  color = Genotype,
  group = interaction(Genotype, Dataset),
  linetype = Dataset
)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  geom_errorbar(
    aes(ymin = slope - se, ymax = slope + se),
    width = 0.55,
    linewidth = 0.8,
    position = pd
  ) +
  geom_line(linewidth = 3.2, position = pd) +
  geom_point(aes(shape = Dataset), size = 6, position = pd) +
  geom_text(
    data = filter(slopes_df, significant),
    aes(label = "*"),
    vjust = -1.2,
    size = 10,
    color = "black",
    show.legend = FALSE,
    position = pd
  ) +
  facet_wrap(~ Treatment, nrow = 1) +
  scale_x_discrete(limits = sort(unique(slopes_df$Brood))) +
  scale_color_manual(values = c(
    "217" = "#e41a1c",
    "391" = "#377eb8",
    "491" = "#4daf4a",
    "399" = "#984ea3",
    "42"  = "#ff7f00"
  )) +
  labs(
    title = "Recombination Rate vs Viability",
    subtitle = "Asterisks indicate significant",
    x = "Brood",
    y = "Standardized Slope (β)",
    color = "Genotype",
    shape = "Marker Type",
    linetype = "Marker Type"
  ) +
  theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p5 
ggsave(
  "../../Summaries and Outputs/3/46_RR_SURVIVAL.png",
  plot = p5,
  width = 20,
  height = 8,
  dpi = 300
)

# Save the final merged and cleaned data for future use
write.csv(combined_vial_data, "../../Summaries and Outputs/3/47_RR_and_Life_History_Merged.csv", row.names = FALSE)


# Summarise recombination data to group level
recomb_summary <- df %>%
  rename(
    RR = recombination_fractions,
    Genotype = GB,
    Treatment = Diet,
    Brood = Brood
  ) %>%
  mutate(across(c(Treatment, Genotype,Brood), factor)) %>%
  group_by(Treatment, Genotype,Brood) %>%
  summarise(
    mean_RR = median(RR, na.rm = TRUE), 
    .groups = "drop")

# Summarise life history data to group level
life_history_summary <- life_history_all %>%
  group_by(Brood, Treatment, Genotype) %>%
  summarise(
    mean_Viability = median(Viability, na.rm = TRUE),
    mean_Fecundity  = median(Fecundity, na.rm = TRUE),
    mean_Fertility  = median(Fertility, na.rm = TRUE),
    .groups = "drop"
  )

# Join on shared keys
combined_group_data <- inner_join(recomb_summary, life_history_summary,
                                  by = c("Brood", "Treatment", "Genotype"))

# ============================================================
# STEP 4: PLOT RR vs Viability BY TREATMENT & GENOTYPE
# ============================================================

p3_faceted <- ggplot(combined_group_data,
                     aes(x = mean_Viability, y = mean_RR, color = Genotype)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 3) +
  facet_wrap(~ Treatment) +
  labs(
    # title = "Recombination Rate vs Viability by Treatment & Genotype",
    x = "Median Viability (%)",
    y = "Median Recombination Rate",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  ) +
  scale_color_manual(values = c(
    "217" = "#e41a1c",  # red
    "391" = "#377eb8",  # blue
    "491" = "#4daf4a",  # green
    "399" = "#984ea3",  # purple
    "42"  = "#ff7f00"   # orange
  ))

# Display the plot
print(p3_faceted)
ggsave(
  filename = "../../Summaries and Outputs/3/48_RR_vs_Viability_Faceted.png",
  plot = p3_faceted,
  width = 20,
  height = 7,
  dpi = 300
)
# ============================================================
# STEP 5: RESULTS TABLES FOR MANUSCRIPT
# ============================================================

library(broom)

# ------------------------------------------------------------
# 5.1 Summary statistics by Treatment × Genotype × Brood
# ------------------------------------------------------------

results_summary <- combined_vial_data %>%
  group_by(Treatment, Genotype, Brood) %>%
  summarise(
    n_vials = n(),
    
    mean_RR = mean(RR, na.rm = TRUE),
    sd_RR   = sd(RR, na.rm = TRUE),
    se_RR   = sd_RR / sqrt(n_vials),
    
    mean_Viability = mean(Viability, na.rm = TRUE),
    sd_Viability   = sd(Viability, na.rm = TRUE),
    se_Viability   = sd_Viability / sqrt(n_vials),
    
    .groups = "drop"
  )

print(results_summary)

# Save table
write.csv(
  results_summary,
  "../../Summaries and Outputs/3/49_Results_Summary_Table.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5.2 Correlation tests for RR vs Viability
# ------------------------------------------------------------

correlation_results <- combined_vial_data %>%
  group_by(Treatment, Genotype, Brood) %>%
  summarise(
    cor_test = list(cor.test(RR, Viability, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    correlation = map_dbl(cor_test, ~ .x$estimate),
    p_value     = map_dbl(cor_test, ~ .x$p.value),
    conf_low    = map_dbl(cor_test, ~ .x$conf.int[1]),
    conf_high   = map_dbl(cor_test, ~ .x$conf.int[2])
  ) %>%
  select(-cor_test)

print(correlation_results)

write.csv(
  correlation_results,
  "../../Summaries and Outputs/3/50_Correlation_Results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5.3 Linear model statistics
# RR ~ Viability
# ------------------------------------------------------------

lm_results <- combined_vial_data %>%
  group_by(Treatment, Genotype, Brood) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(RR ~ Viability, data = .x)),
    tidy  = map(model, tidy),
    glance = map(model, glance)
  )

# Extract coefficients
lm_coefficients <- lm_results %>%
  unnest(tidy) %>%
  filter(term == "Viability") %>%
  select(
    Treatment,
    Genotype,
    Brood,
    estimate,
    std.error,
    statistic,
    p.value
  )

print(lm_coefficients)

write.csv(
  lm_coefficients,
  "../../Summaries and Outputs/3/51_LM_Coefficients.csv",
  row.names = FALSE
)

# Extract model fit statistics
lm_fit_stats <- lm_results %>%
  unnest(glance) %>%
  select(
    Treatment,
    Genotype,
    Brood,
    r.squared,
    adj.r.squared,
    sigma,
    p.value,
    AIC,
    BIC
  )

print(lm_fit_stats)

write.csv(
  lm_fit_stats,
  "../../Summaries and Outputs/3/52_LM_Fit_Statistics.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5.4 Create manuscript-ready text output
# ------------------------------------------------------------

results_text <- lm_coefficients %>%
  mutate(
    significance = ifelse(p.value < 0.05, "significant", "non-significant"),
    
    sentence = paste0(
      "In genotype ", Genotype,
      " under the ", Treatment,
      " treatment during brood ", Brood,
      ", recombination rate showed a ",
      significance,
      " association with viability (β = ",
      round(estimate, 3),
      " ± ",
      round(std.error, 3),
      ", p = ",
      signif(p.value, 3),
      ")."
    )
  )

cat("\n\n================ RESULTS TEXT ================\n")
cat(paste(results_text$sentence, collapse = "\n"))
cat("\n=============================================\n")

# Save manuscript text
writeLines(
  results_text$sentence,
  "../../Summaries and Outputs/3/53_Results_Text.txt"
)

# ------------------------------------------------------------
# 5.5 Overall treatment/genotype means
# Useful for main Results section
# ------------------------------------------------------------

overall_summary <- combined_vial_data %>%
  group_by(Treatment, Genotype) %>%
  summarise(
    n_vials = n(),
    
    RR_mean = mean(RR, na.rm = TRUE),
    RR_se   = sd(RR, na.rm = TRUE) / sqrt(n_vials),
    
    Viability_mean = mean(Viability, na.rm = TRUE),
    Viability_se   = sd(Viability, na.rm = TRUE) / sqrt(n_vials),
    
    .groups = "drop"
  )

print(overall_summary)

write.csv(
  overall_summary,
  "../../Summaries and Outputs/3/54_Overall_Group_Summary.csv",
  row.names = FALSE
)

