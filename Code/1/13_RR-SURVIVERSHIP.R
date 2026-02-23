# ============================================================
# CLEANED & ORGANIZED: RR vs Viability PLOT
# ============================================================

library(dplyr)
library(ggplot2)
library(tidyr)

# Set working directory (uncomment if running locally)
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/1")

# Constants
MOTHERS_PER_VIAL <- 5
traits <- c("Fecundity", "Fertility", "Viability")

# ============================================================
# STEP 1: LOAD DATA
# ============================================================

multiply_marked_viability_dataset <- read.csv("6_multiply_marked_viability_dataset.csv")
# mothers per vial for marker-present
mothers_per_vial <- read.csv("1_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))

df <- read.csv("df.csv")
df <- df %>%
  #group_by(Diet, GB, Brood) %>%
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
    Fecundity = transformed_mean_multiply_marked_total_egg_count / female_in_vial,
    Fertility = total_adult_count / female_in_vial,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(VialID = cross_id, Data_collector = data_collector, Treatment = treatment, 
         Genotype = genotype, Brood = brood_id, Dataset, Fecundity, Fertility, Viability)

life_history_all <- marker_present_life_history %>%
  mutate(across(c(Treatment, Genotype, Brood, Dataset, Data_collector), factor))

# ============================================================
# STEP 3: MERGE WITH RECOMBINATION DATA
# ============================================================

combined_vial_data <- df %>%
  rename(
    RR = recombination_fractions,
    Genotype = GB,
    Treatment = Diet,
    Brood = Brood,
    VialID = vialID
  ) %>%
  mutate(across(c(Treatment, Genotype, Brood), factor)) %>%
  select(VialID, Brood, Treatment, Genotype, RR) %>%
  inner_join(life_history_all, by = c("VialID", "Brood", "Treatment", "Genotype"))

# ============================================================
# STEP 4: PLOT RR vs Viability BY TREATMENT & GENOTYPE
# ============================================================

dataset_clean <- combined_vial_data %>%
  mutate(
    Genotype = as.factor(Genotype),
    Treatment = as.factor(Treatment)
  )

p3_faceted <- ggplot(dataset_clean, 
                     aes(x = Viability, y = RR, color = Genotype)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, size = 3) +
  facet_wrap(~ Treatment) +
  labs(
    title = "Recombination Rate vs Viability by Treatment & Genotype",
    x = "Viability (%)", 
    y = "Recombination Rate",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
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
  filename = "../../Summaries and Outputs/1/26_RR_vs_Viability_Faceted.png",
  plot = p3_faceted,
  width = 18,
  height = 8,
  dpi = 300
)

# Save merged data for reproducibility
write.csv(combined_vial_data, "../../Summaries and Outputs/1/27_RR_and_Life_History_Merged.csv", row.names = FALSE)
