#set the Directory to the Summaries and Outputs of each experiment
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)
library(multcomp)
library(lme4)
library(lmerTest)
library(broom)
library(tidyr)

# =========================================================
# 1. Read and prepare main data (interference, MLE_p)
# =========================================================
df <- read.csv("Nutrition_countingmodel_outputs_all_polled_treatments.csv")

df <- df %>%
  extract(
    original_ID,
    into = c("experiment", "vialID", "Brood"),
    regex = "^(.*)_([0-9]+)_([A-Za-z])$"
  ) %>%
  rename(interference = MLE_p)

# 2. Read and prepare treatment key
treatment_key <- read.csv("Treatmentskeys.csv")
treatment_key <- treatment_key %>%
  rename(vialID = `Vial.ID`, GB = `F1.Genotype`, Diet = `F1.developmental`)

# 3. Make sure vialID type matches
df$vialID <- as.character(df$vialID)
treatment_key$vialID <- as.character(treatment_key$vialID)

# 4. Merge
df <- df %>%
  left_join(treatment_key, by = "vialID")

# 5. Set variable types
df$vialID       <- as.factor(df$vialID)
df$Brood        <- as.factor(df$Brood)
df$Diet         <- as.factor(df$Diet)
df$GB           <- as.factor(df$GB)
df$People       <- as.factor(df$People)
df$interference <- as.numeric(df$interference)

write.csv(df, "17_df_interference.csv", row.names = FALSE)

# =========================================================
# 7. ANOVA on interference
# =========================================================
anova_result <- aov(interference ~ Diet * GB, data = df)
summary(anova_result)

anova_summary <- tidy(anova_result) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )
write.csv(anova_summary, "18_anova_model_interference.csv", row.names = FALSE)

# 8. Post hoc Tukey
TukeyHSD(anova_result)

# =========================================================
# 9. Estimated marginal means and compact letter display
# =========================================================
emm <- emmeans(anova_result, ~ Diet * GB)
cld_result <- cld(emm, Letters = letters)

plot_data <- cld_result %>%
  mutate(
    GB    = factor(GB),
    Diet  = factor(Diet, levels = c("0.5x", "1.0x", "2.0x")),
    .group = trimws(.group)
  )

p1 <- ggplot(plot_data, aes(x = Diet, y = emmean, color = GB,)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, position = position_dodge(0.4)) +
  geom_text(aes(label = .group),
            vjust = -0.8,
            position = position_dodge(0.4),
            size = 5) +
  labs(x = "Diets",
       y = "Interference (Mean ± SE)",
       color = "Genotype") +
  theme_classic(base_size = 14)
p1

ggsave(
  filename = "19_cld_interference.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)


# =========================================================
# 11. Reaction-norm style plot (group means +/- SE)
# =========================================================
se <- function(x) {
  std_dev <- sd(x, na.rm = TRUE)
  sample_size <- sum(!is.na(x))
  std_dev / sqrt(sample_size)
}

interference_summary <- df %>%
  rename(
    Genotype  = GB,
    Treatment = Diet
  ) %>%
  mutate(across(c(Treatment, Genotype), factor)) %>%
  group_by(Treatment, Genotype) %>%
  summarise(
    mean_interference = mean(interference, na.rm = TRUE),
    sd_interference   = sd(interference, na.rm = TRUE),
    se_interference   = se(interference),
    n = n(),
    .groups = "drop"
  )

print(interference_summary)
write.csv(interference_summary, "20_mean_interference_by_GB_Diet.csv", row.names = FALSE)

p2b <- ggplot(interference_summary, aes(x = as.factor(Treatment),
                                        y = mean_interference,
                                        group = as.factor(Genotype),
                                        color = as.factor(Genotype))) +
  geom_point(size = 6) +
  geom_line(linewidth = 3) +
  geom_errorbar(aes(ymin = mean_interference - se_interference,
                    ymax = mean_interference + se_interference),
                width = 0.1, linewidth = 0.5) +
  labs(
    title = "Interference",
    x = "Treatments",
    y = "Interference (MLE_p)",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(size = 35, hjust = 0.5),
    axis.line        = element_line(color = "black"),
    axis.text.x      = element_text(size = 30, face = "bold"),
    axis.text.y      = element_text(size = 30, face = "bold"),
    legend.text      = element_text(size = 25, face = "bold"),
    legend.title     = element_text(size = 25, face = "bold")
  ) +
  scale_color_manual(values = c(
    "217" = "#e41a1c",
    "391" = "#377eb8",
    "491" = "#4daf4a",
    "399" = "#984ea3",
    "42"  = "#ff7f00"
  ))
p2b

ggsave(
  filename = "21_interference_GxE.png",
  plot = p2b,
  width = 20,
  height = 6,
  dpi = 300
)
