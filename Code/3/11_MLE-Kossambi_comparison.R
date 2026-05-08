# -------------------------------
# Libraries
# -------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# -------------------------------
# Set working directory + load data
# -------------------------------
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

df_cleaned <- read_csv("df_cleaned.csv")
df_RecombinationFrequency <- read_csv("df_RecombinationFrequency.csv")

# -------------------------------
# Scale correction (MLE dataset)
# -------------------------------
df_cleaned <- df_cleaned %>%
  mutate(recombination_fractions = recombination_fractions / 100)

# -------------------------------
# Rename MLE interval columns
# -------------------------------
df_cleaned2 <- df_cleaned %>%
  rename(
    RF_y_cv   = MLE_x5,
    RF_cv_v   = MLE_x4,
    RF_v_f    = MLE_x3,
    RF_f_car  = MLE_x2,
    RF_car_yP = MLE_x1
  )

# -------------------------------
# Long format (intervals)
# -------------------------------
df_cleaned_long <- df_cleaned2 %>%
  pivot_longer(starts_with("RF_"),
               names_to = "interval",
               values_to = "MLE_value")

df_RF_long <- df_RecombinationFrequency %>%
  pivot_longer(starts_with("RF_"),
               names_to = "interval",
               values_to = "Observed_value")

# -------------------------------
# Merge datasets
# -------------------------------
df_merged <- df_cleaned_long %>%
  inner_join(df_RF_long,
             by = c("experiment", "vialID", "Brood", "Diet", "GB", "People", "interval")) %>%
  mutate(genotype_diet = paste0(GB, "_", Diet))

# -------------------------------
# Order intervals
# -------------------------------
df_merged$interval <- factor(
  df_merged$interval,
  levels = c("RF_y_cv", "RF_cv_v", "RF_v_f", "RF_f_car", "RF_car_yP")
)

# -------------------------------
# Interval differences
# -------------------------------
df_intervals <- df_merged %>%
  mutate(diff = MLE_value - Observed_value) %>%
  group_by(genotype_diet, interval) %>%
  summarise(diff = mean(diff, na.rm = TRUE), .groups = "drop")

# -------------------------------
# Total recombination difference
# -------------------------------
df_total <- df_cleaned %>%
  select(experiment, vialID, Brood, Diet, GB, People, recombination_fractions) %>%
  rename(MLE_total = recombination_fractions) %>%
  inner_join(
    df_RecombinationFrequency %>%
      select(experiment, vialID, Brood, Diet, GB, People, recombination_fractions) %>%
      rename(Observed_total = recombination_fractions),
    by = c("experiment", "vialID", "Brood", "Diet", "GB", "People")
  ) %>%
  mutate(total_diff = MLE_total - Observed_total,
         genotype_diet = paste0(GB, "_", Diet)) %>%
  group_by(genotype_diet) %>%
  summarise(total_diff = mean(total_diff, na.rm = TRUE), .groups = "drop")

# -------------------------------
# Table: interval differences
# -------------------------------
table_intervals <- df_intervals %>%
  pivot_wider(names_from = interval, values_from = diff)

# -------------------------------
# Final table
# -------------------------------
final_table <- df_total %>%
  left_join(table_intervals, by = "genotype_diet")

write.csv(final_table,
          "Recombination_comparison_MLE-Kossambi.csv",
          row.names = FALSE)

# -------------------------------
# Plot data (intervals)
# -------------------------------
df_plot <- df_merged %>%
  group_by(genotype_diet, interval) %>%
  summarise(
    MLE = mean(MLE_value, na.rm = TRUE),
    Observed = mean(Observed_value, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# Total as extra panel
# -------------------------------
df_total_plot <- df_cleaned %>%
  mutate(MLE_total = recombination_fractions) %>%
  inner_join(
    df_RecombinationFrequency %>%
      mutate(Observed_total = recombination_fractions),
    by = c("experiment", "vialID", "Brood", "Diet", "GB", "People")
  ) %>%
  mutate(genotype_diet = paste0(GB, "_", Diet),
         interval = "TOTAL") %>%
  group_by(genotype_diet, interval) %>%
  summarise(
    MLE = mean(MLE_total, na.rm = TRUE),
    Observed = mean(Observed_total, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------
# Combine + order
# -------------------------------
df_plot_all <- bind_rows(df_plot, df_total_plot)

df_plot_all$interval <- factor(
  df_plot_all$interval,
  levels = c("TOTAL",
             "RF_y_cv", "RF_cv_v", "RF_v_f", "RF_f_car", "RF_car_yP")
)

# -------------------------------
# Plot
# -------------------------------
p <- ggplot(df_plot_all, aes(x = Observed, y = MLE)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ interval) +
  labs(
    title = "MLE vs Observed Recombination Fractions",
    x = "Observed",
    y = "MLE"
  ) +
  theme_minimal()

ggsave(
  filename = "../../Summaries and Outputs/3/recombinationrate_comparison_MLE-Kossambi.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------
# 1. Global bias test
# -------------------------------
t.test(final_table$total_diff, mu = 0)

# -------------------------------
# 2. Genotype × diet effect on total bias
# -------------------------------
final_table %>%
  mutate(
    genotype = sub("_.*", "", genotype_diet),
    diet = sub(".*_", "", genotype_diet)
  ) %>%
  lm(total_diff ~ genotype * diet, data = .) %>%
  anova()

# -------------------------------
# 3. Interval-specific effects (ANOVA)
# -------------------------------
library(tidyr)

rf_long <- final_table %>%
  pivot_longer(
    cols = starts_with("RF_"),
    names_to = "interval",
    values_to = "diff"
  )

lm(diff ~ interval, data = rf_long) %>%
  anova()

# -------------------------------
# 4. Interval bias visualization
# -------------------------------
library(ggplot2)
p2=ggplot(rf_long, aes(x = interval, y = diff)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Bias across recombination intervals",
    x = "Interval",
    y = "MLE − Observed"
  )

ggsave(
  filename = "../../Summaries and Outputs/3/recombinationrate_interval_comparison_MLE-Kossambi.png",
  plot = p2,
  width = 8,
  height = 6,
  dpi = 300
)
# -------------------------------
# 5. Mixed-effects model (final model)
# -------------------------------
library(lme4)

rf_long <- final_table %>%
  mutate(
    genotype = sub("_.*", "", genotype_diet),
    diet = sub(".*_", "", genotype_diet)
  ) %>%
  pivot_longer(
    cols = starts_with("RF_"),
    names_to = "interval",
    values_to = "diff"
  )

model <- lmer(
  diff ~ interval + genotype + diet + (1 | genotype_diet),
  data = rf_long
)

summary(model)

#The MLE approach systematically overestimates recombination rates,
#but the magnitude of this bias is strongly dependent on genomic interval, with the largest deviation observed in RF_v_f.
#No significant effects of genotype or diet were detected.