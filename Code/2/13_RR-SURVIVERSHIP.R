# ============================================================
# CLEANED & ORGANIZED: RR vs Viability PLOT
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Set working directory (uncomment if running locally)
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/2")

# Constants
MOTHERS_PER_VIAL <- 5
traits <- c("Fecundity", "Fertility", "Viability")

# ============================================================
# STEP 1: LOAD DATA
# ============================================================

multiply_marked_viability_dataset <- read.csv("6_multiply_marked_viability_dataset.csv")
df <- read.csv("df.csv")

# Filter extreme outliers in recombination data
df <- df %>%
  #group_by(Diet, GB) %>%
  filter(
    recombination_fractions >= quantile(recombination_fractions, 0.25, na.rm = TRUE) -
      1.5 * IQR(recombination_fractions, na.rm = TRUE),
    recombination_fractions <= quantile(recombination_fractions, 0.75, na.rm = TRUE) +
      1.5 * IQR(recombination_fractions, na.rm = TRUE)
  ) %>%
  ungroup()

# ============================================================
# STEP 2: PROCESS LIFE HISTORY DATA
# ============================================================

marker_present_life_history <- multiply_marked_viability_dataset %>%
  mutate(
    Dataset   = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count /MOTHERS_PER_VIAL ,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(VialID = cross_id, Data_collector = data_collector, Treatment = treatment, 
         Genotype = genotype, Brood = brood_id, Dataset, Fecundity, Fertility, Viability)

life_history_all <- marker_present_life_history %>%
  mutate(across(c(Treatment, Genotype, Brood, Dataset, Data_collector), factor))


# ============================================================
# STEP 3: AGGREGATE THEN MERGE (no VialID join)
# ============================================================

# Summarise recombination data to group level
recomb_summary <- df %>%
  rename(
    RR = recombination_fractions,
    Genotype = GB,
    Treatment = Diet,
    Brood = Brood
  ) %>%
  mutate(across(c(Treatment, Genotype, Brood), factor)) %>%
  group_by(Brood, Treatment, Genotype) %>%
  summarise(mean_RR = mean(RR, na.rm = TRUE), .groups = "drop")

# Summarise life history data to group level
life_history_summary <- life_history_all %>%
  group_by(Brood, Treatment, Genotype) %>%
  summarise(
    mean_Viability = mean(Viability, na.rm = TRUE),
    mean_Fecundity  = mean(Fecundity, na.rm = TRUE),
    mean_Fertility  = mean(Fertility, na.rm = TRUE),
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
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  facet_wrap(~ Treatment) +
  labs(
    title = "Recombination Rate vs Viability by Treatment & Genotype",
    x = "Mean Viability (%)",
    y = "Mean Recombination Rate",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 22)
  ) +
  scale_color_manual(values = c(
    "217" = "#8dd3c7",
    "391" = "#ffffb3",
    "491" = "#80b1d3",
    "399" = "#9e9ac8",
    "42"  = "#fb8072"
  ))
# Display the plot
print(p3_faceted)

# ============================================================
# STEP 5: SAVE PLOT
# ============================================================

ggsave(
  filename = "../../Summaries and Outputs/2/26_RR_vs_Viability_Faceted.png",
  plot = p3_faceted,
  width = 16,
  height = 8,
  dpi = 300
)

# Save merged data for reproducibility
write.csv(combined_group_data, "../../Summaries and Outputs/2/27_RR_and_Life_History_Merged.csv", row.names = FALSE)
