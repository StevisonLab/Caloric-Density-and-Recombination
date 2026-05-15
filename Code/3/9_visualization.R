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
df <- read.csv("Nutrition_countingmodel_outputs_all.csv")
#df<- read.csv("Nutrition_countingmodel_outputs_crosspooled.csv")

df <- df %>% extract(original_ID,
    into = c("experiment", "vialID", "Brood"),
    regex = "^(.*)_([0-9]+)_([A-Za-z])$"
  ) %>%
  rename(recombination_fractions = MLE_X_length)

#_______
#df <- read.csv("14_recombination_frequency_per_vial.csv")

df <- df %>% 
  rename(recombination_fractions = Total_RF)

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
write.csv(df, "df.csv", row.names = FALSE)
write.csv(df, "df_RecombinationFrequency.csv", row.names = FALSE)

###delete the outlier (just checking)
#statisticaly
df <- df %>%
  #group_by(Diet, GB) %>%
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
#df <- df %>%filter(recombination_fractions < 1000)  
write.csv(df, "df_cleaned.csv", row.names = FALSE)

# 3. ANOVA
anova_result <- aov(recombination_fractions ~ Brood + GB + Diet + People, data = df)
anova_result <- aov(recombination_fractions ~ Brood * Diet * GB , data = df)
summary(anova_result)
anova_summary <- tidy(anova_result) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    ))
write.csv(anova_summary, "../../Summaries and Outputs/3/15_anova_model_recombination.csv", row.names = FALSE)
write.csv(anova_summary, "../../Summaries and Outputs/3/15_anova_model_recombinationfrequency.csv", row.names = FALSE)

# 4. Post hoc Tukey
TukeyHSD(anova_result)

# 5. Estimated marginal means and CLD
emm <- emmeans(anova_result, ~ Diet *GB*Brood)
cld <- cld(emm, Letters = letters)

plot_data <- cld %>%
  mutate(
   Brood=factor(Brood),
    GB = factor(GB),
   # People= factor(People),
    Diet = factor(Diet, levels = c("0.5x", "1.0x","2.0x")),
    .group = trimws(.group) # remove extra spaces
  )
# Plot
p1=ggplot(plot_data, aes(x = Brood, y = emmean, color = GB, group = Diet)) +
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
p1
ggsave(
  filename = "../../Summaries and Outputs/1/16_cld_recombinationrate.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)
#excluding the data
df <- df %>%filter(People %in% c("M","C","S"))

# 3. ANOVA
anova_result <- aov(recombination_fractions ~ Brood + Diet +GB+ People, data = df)
anova_result <- aov(recombination_fractions ~ Brood * Diet *GB, data = df)
summary(anova_result)
anova_summary <- tidy(anova_result) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    ))
#write.csv(anova_summary, "../../Summaries and Outputs/1/15_anova_model_recombination.csv", row.names = FALSE)

# 4. Post hoc Tukey
TukeyHSD(anova_result)

# 5. Estimated marginal means and CLD
emm <- emmeans(anova_result, ~ Diet *Brood)
cld <- cld(emm, Letters = letters)

plot_data <- cld %>%
  mutate(
    Brood=factor(Brood),
    #GB = factor(GB),
    # People= factor(People),
    Diet = factor(Diet, levels = c("0.5x","0.75x", "1.0x","1.25x" ,"1.5x","2.0x")),
    .group = trimws(.group) # remove extra spaces
  )
# Plot
p1=ggplot(plot_data, aes(x = Brood, y = emmean, color = GB, group = Diet)) +
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

#ggsave(
  filename = "../../Summaries and Outputs/1/16_cld_recombinationrate.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

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
sig_B   <- FALSE
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

p2b <- ggplot(recomb_summary,aes( x =as.factor(Treatment),
                                  y = mean_RR,
                                  group = as.factor(Genotype),
                                  color = as.factor(Genotype)))+geom_point(size = 6)+ geom_line(linewidth = 3)+
  geom_errorbar(aes(ymin = mean_RR - se_RR, ymax = mean_RR + se_RR), width = 0.55) +
#  ylim(65,110)+
# Add equation 
annotate(
  "text",
  x = 2,          # adjust if needed
  y = 99.5,
  label = eq_label,
  parse = TRUE,
  size = 16
) +
  
  labs(
    # title = "RR G X E",
    x = "Treatments",
    y = "RR",
    color = "Genotype"
  ) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5)
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
  filename = "../../Summaries and Outputs/3/18_RR_GxE.png",
  plot = p2b,
  width = 13,
  height = 6,
  dpi = 300
)
#per strain brood 
df$Diet_GB <- paste0(df$GB, "_", df$Diet)

p3=ggplot(df, aes(x = Brood, y = recombination_fractions, group = Diet_GB, color = Diet_GB)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, size = 0.8) +  # error bars
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(title = "RR by G × E × Age",
       x = "Age",
       y = "RR",
       color = "Diet"
  ) +
  facet_wrap(~GB) +
  theme_minimal(base_size = 25) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 14),
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  ) +
  scale_color_manual(
    values = c(
      "42_0.5x"  = "#fcbba1",
      "42_1.0x"  = "#fb6a4a",
      "42_2.0x"  = "#99000d",
      "399_0.5x" = "#dadaeb",
      "399_1.0x" = "#9e9ac8",
      "399_2.0x" = "#4a1486"
    ),
    labels = c(
      "42_0.5x"  = "0.5x",
      "42_1.0x"  = "1.0x",
      "42_2.0x"  = "2.0x",
      "399_0.5x" = "0.5x",
      "399_1.0x" = "1.0x",
      "399_2.0x" = "2.0x"
    )
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1, shape = 16)))
p3
ggsave(
  filename = "../../Summaries and Outputs/3/19_RR_GxExBrood.png",
  plot = p3,
  width = 12,
  height = 6,
  dpi = 300
)
