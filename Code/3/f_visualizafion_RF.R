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

# 1. Read and prepare main data
df <- read.csv("16_recombination_frequency_per_vial.csv")
#df <- df %>% rename(recombination_fractions = Total_RF)
#Kosambi corrected
df <- df %>% rename(recombination_fractions = cM_Total)
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
write.csv(df, "29_df_RecombinationFrequency.csv", row.names = FALSE)

# Number of observations before filtering
n_before <- nrow(df)

# Calculate IQR thresholds across all observations
Q1 <- quantile(df$recombination_fractions, 0.25, na.rm = TRUE)
Q3 <- quantile(df$recombination_fractions, 0.75, na.rm = TRUE)
IQR_value <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_value
upper_bound <- Q3 + 1.5 * IQR_value

# Count outliers
n_outliers <- sum(
  df$recombination_fractions < lower_bound |
    df$recombination_fractions > upper_bound,
  na.rm = TRUE
)

# Count missing observations separately
n_missing <- sum(is.na(df$recombination_fractions))

# Remove outliers
df <- df %>%
  ungroup() %>%
  filter(
    recombination_fractions >= lower_bound,
    recombination_fractions <= upper_bound
  )

# Number of observations after filtering
n_after <- nrow(df)

# Summarize filtering
data.frame(
  Before = n_before,
  After = n_after,
  Outliers_removed = n_outliers,
  Missing_removed = n_missing,
  Percent_outliers_removed = round(
    100 * n_outliers / (n_before - n_missing), 2
  ),
  Lower_bound = lower_bound,
  Upper_bound = upper_bound
)

###delete the outlines
#Biological
boxplot(recombination_fractions ~ Diet * GB, data = df)
write.csv(df,"30_df_RecombinationFrequency_cleaned.csv", row.names = FALSE)

# 3. 
# Fit GLMM for Recombination Rate / Recombination Fraction
rr_glmm <- lmer(
  recombination_fractions ~ Brood * Diet * GB + (1 | vialID) + (1 | People), 
  data = df
)

# Extract Type III ANOVA
rr_anova <- anova(rr_glmm, type = 3)
print(rr_anova)

# Estimated Marginal Means across Diets grouped by Genotype (GB)
rr_emm <- emmeans(rr_glmm, ~ Diet | GB)

# Pairwise Post-Hoc Comparisons
rr_pairs <- pairs(rr_emm, adjust = "tukey")

# Detailed Percent Change calculation relative to 1.0x Diet
rr_pct <- as.data.frame(rr_emm) %>%
  group_by(GB) %>%
  mutate(
    Control_Mean = emmean[Diet == "1.0x"],
    Percent_Change = ((emmean - Control_Mean) / Control_Mean) * 100
  )

# Save Outputs
write.csv(as.data.frame(rr_anova), "31_recombinationfrequency_glmm_anova.csv")
write.csv(as.data.frame(rr_pairs), "31_recombinationfrequency_glmm_pairs.csv")
write.csv(rr_pct, "31_recombinationfrequency_glmm_emmeans_pct.csv")

#anova_result <- aov(recombination_fractions ~ Brood * Diet * GB , data = df)
#summary(anova_result)
#anova_summary <- tidy(anova_result) %>%
#  mutate(
 #   significance = case_when(
  #    p.value < 0.001 ~ "***",
   #   p.value < 0.01  ~ "**",
    #  p.value < 0.05  ~ "*",
     # p.value < 0.1   ~ ".",
      #TRUE ~ ""
    #))
#write.csv(anova_summary, "../../Summaries and Outputs/3/31_anova_model_recombinationfrequency.csv", row.names = FALSE)

# 6. Visualization
p2=ggplot(df, aes(x = Diet, y = recombination_fractions, fill = as.factor(GB))) +
  geom_boxplot() +
  facet_wrap(~GB, scales = "fixed", nrow = 1) +  # one panel per genotype
  labs(
    title = "RR G X E",
    x = "Treatments",
    y = "RR"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5)
  ) +
  scale_fill_manual(values = c(
    "217" = "#8dd3c7",
    "391" = "#ffffb3",
    "491" = "#80b1d3",
    "399" = "#9e9ac8",
    "42"  = "#fb8072"
  ))
p2
#reaction norm 
# Define which terms are significant (TRUE = bold)
sig_B   <- TRUE
sig_G   <- TRUE
sig_E   <- FALSE
sig_GxE <- TRUE

# Helper: wrap term in bold() if significant
fmt_term <- function(term, significant) {if (significant) paste0("bold(", term, ")") else term}

# Build the label as a string expression
eq_label <- paste0( "RR * '~' * ", fmt_term("B", sig_B), " + ", fmt_term("G", sig_G), " + ", fmt_term("E", sig_E), " + ", fmt_term("GxE", sig_GxE))

# Define the standard error function
se <- function(x) {
  # Calculate standard deviation (sd) and sample size (length)
  std_dev <- sd(x)
  sample_size <- length(x)
  # Return standard error: sd / sqrt(n)
  standard_error <- std_dev / sqrt(sample_size)
  return(standard_error)
}

# Summarise recombination data to group level
recomb_summary <- df %>%
  rename(
    RR = recombination_fractions,
    Genotype = GB,
    Treatment = Diet,
    #   Brood = Brood
  ) %>%
  mutate(across(c(Treatment, Genotype), factor)) %>%
  group_by(Treatment, Genotype) %>%
  summarise(
    mean_RR = median(RR, na.rm = TRUE), 
    sd_RR = sd(RR,na.rm=TRUE),
    # se_RR = stderr(RR,na.rm=TRUE),
    se_RR = se(RR),
    
    .groups = "drop")

p2b <- ggplot(recomb_summary, aes(x = as.factor(Treatment),
                                  y = mean_RR,
                                  group = as.factor(Genotype),
                                  color = as.factor(Genotype))) +
  geom_point(size = 6) +
  geom_line(linewidth = 3) +
  geom_errorbar(aes(ymin = mean_RR - se_RR, ymax = mean_RR + se_RR), 
                width = 0.1, linewidth = 0.5) +
  
  # annotate(
  #   "text",
  #   x = 2,
  #   y = 99.5,
  #   label = eq_label,
  #   parse = TRUE,
  #   size = 16
  # ) +
  
  labs(
    title = "Recombination Fraction",
    x = "Treatments",
    y = "RF",
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
  )   +
  scale_color_manual(values = c(
    "217" = "#e41a1c",
    "391" = "#377eb8",
    "491" = "#4daf4a",
    "399" = "#984ea3",
    "42"  = "#ff7f00"
  ))

p2b


ggsave(
  filename = "../../Summaries and Outputs/3/31_RF_GxE.png",
  plot = p2b,
  width = 20,
  height = 6,
  dpi = 300
)


# =========================================================
# Mean recombination fraction for each Genotype × Treatment
# (gives 6 means total if you have 2 genotypes × 3 treatments)
# =========================================================

mean_recombination <- df %>%
  group_by(GB, Diet) %>%
  summarise(
    mean_RF = mean(recombination_fractions, na.rm = TRUE),
    sd_RF   = sd(recombination_fractions, na.rm = TRUE),
    se_RF   = sd_RF / sqrt(n()),
    n       = n(),
    .groups = "drop"
  )

# View results
print(mean_recombination)

write.csv(
mean_recombination,
"../../Summaries and Outputs/3/32_mean_recombination_by_GB_Diet_RECOMBINATIONFREQUENCY.csv",
row.names = FALSE
)
