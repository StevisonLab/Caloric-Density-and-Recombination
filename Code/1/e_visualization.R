setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/1")


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
# 1. Read and prepare main data
df <- read.csv("Nutrition_countingmodel_outputs_all.csv")

df <- df %>% extract(original_ID,
    into = c("experiment", "vialID", "Brood"),
    regex = "^(.*)_([0-9]+)_([A-Za-z])$"
  ) %>%
  rename(recombination_fractions = MLE_X_length)

# 2. Read and prepare treatment key for the first 2 
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
df$vialID <- as.factor(df$vialID)
df$Brood <- as.factor(df$Brood)
df$Diet <- as.factor(df$Diet)
df$GB <- as.factor(df$GB)
df$People <- as.factor(df$People)
df$recombination_fractions <- as.numeric(df$recombination_fractions)
write.csv(df, "14_df.csv", row.names = FALSE)

###delete the outlier (just checking)
#statisticaly
df <- df %>%
  filter(
    recombination_fractions >= quantile(recombination_fractions, 0.25, na.rm = TRUE) -
      1.5 * IQR(recombination_fractions, na.rm = TRUE),
    recombination_fractions <= quantile(recombination_fractions, 0.75, na.rm = TRUE) +
      1.5 * IQR(recombination_fractions, na.rm = TRUE)
  ) %>%
  ungroup()

###delete the outlines
#Biological
boxplot(recombination_fractions ~ Diet * GB, data = df)
write.csv(df, "14_df_cleaned.csv", row.names = FALSE)

# 3. ANOVA 
#people are 2 M,N
anova_result <- aov(recombination_fractions ~ Brood * Diet , data = df)

summary(anova_result)
# Tidy the ANOVA model and save to CSV
anova_summary <- tidy(anova_result) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    ))
write.csv(anova_summary, "../../Summaries and Outputs/1/15_anova_model_recombination.csv", row.names = FALSE)

# 4. Post hoc Tukey
TukeyHSD(anova_result)

# 5. Estimated marginal means and CLD
emm <- emmeans(anova_result, ~ Diet*Brood)
cld <- cld(emm, Letters = letters)

plot_data <- cld %>%
  mutate(
    Brood=factor(Brood),
   # GB = factor(GB),
  # People= factor(People),
    Diet = factor(Diet, levels = c("0.5x","0.75x", "1.0x","1.25x","1.5x", "2.0x")),
    .group = trimws(.group) # remove extra spaces
  )
# Plot
ggplot(plot_data, aes(x = Brood, y = emmean, color = Diet, group = Diet)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, position = position_dodge(0.4)) +
  geom_text(aes(label = .group),
            vjust = -0.8,
            position = position_dodge(0.4),
            size = 5) +
  labs(x = "Diets",
       y = "Recombination Rate (Mean ± SE)",
       color = "Genotype") +
  theme_classic(base_size = 14)

# 6. Visualization
p1=ggplot(data = df, aes(x = Diet, y = recombination_fractions, fill = Diet)) +
  geom_boxplot() +
  labs(
    title = "RR",
    x = "Treatments",
    y = "RR"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5) # size increased, centered, bold
  ) +
  scale_fill_manual(values = c(
    "0.5x"  = "#dadaeb",
    "0.75x" = "#bcbddc",
    "1.0x"  = "#9e9ac8",
    "1.25x" = "#807dba",
    "1.5x"  = "#6a51a3",
    "2.0x"  = "#4a1486"
  ))

ggsave(
  filename = "../../Summaries and Outputs/1/16_RR_GxE.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)
