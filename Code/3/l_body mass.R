#setting working directory
setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3")


# Load necessary libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(car)
library(emmeans)
library(ggplot2)
library(multcompView)
library(multcomp)
library(readr)
library(stringr)
library(lme4)
library(lmerTest) # For p-values in mixed models
library(broom)


# Import the merged_data
data <- read.csv("BodyMass.csv")
Key <- read.csv("BodyMasscross.csv")

# Merge key  + data 
merged_data <- data %>%
  left_join(Key, by = c("ID"))
# Remove rows with missing data in important columns
merged_data <- merged_data %>% drop_na(Scale)

# Ensure categorical variables are factors
merged_data$ID <- as.character(merged_data$ID)
merged_data$Scale <- as.numeric(merged_data$Scale)
merged_data$treatment <- as.character(merged_data$treatment)
merged_data$Genotype <- as.character(merged_data$Genotype)

#=================================================
#                ANOVA + Post-hoc
#=================================================

anova_result <- aov(Scale ~ Genotype * treatment * Brood, data = merged_data)
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
write.csv(anova_summary,"../../Summaries and Outputs/3/55_bodymass_Anova_table.csv",row.names = FALSE, quote = FALSE)

#-------------------------
# Post-hoc with emmeans
#-------------------------

# --- Genotype post-hoc ---
genotype_emm <- emmeans(anova_result, ~ Genotype | Brood )  # consider Brood because of interaction
genotype_posthoc <- pairs(genotype_emm, adjust = "tukey")
summary(genotype_posthoc)
write.csv(as.data.frame(genotype_posthoc),
          "../../Summaries and Outputs/3/56_Tukey_Genotype.csv", row.names = FALSE)

# Get letters for plotting
genotype_cld <- cld(genotype_emm, Letters = letters, adjust = "tukey")

# --- Treatment post-hoc ---
treatment_emm <- emmeans(anova_result, ~ treatment | Brood)
treatment_posthoc <- pairs(treatment_emm, adjust = "tukey")
summary(treatment_posthoc)
write.csv(as.data.frame(treatment_posthoc),
          "../../Summaries and Outputs/3/57_Tukey_Treatment.csv", row.names = FALSE)

# Get letters for plotting
treatment_cld <- cld(treatment_emm, Letters = letters, adjust = "tukey")

#-------------------------
# Boxplots with letters
#-------------------------
# Treatment boxplot faceted by Brood
abc <- ggplot(merged_data, aes(x = treatment, y = Scale, fill = Genotype)) +
  geom_boxplot(na.rm = TRUE) + 
  labs(
    title = "Body Mass",
    x = "Treatment", 
    y = "Body Mass (mg)"
  ) +  
  facet_wrap(~ Genotype) +
  theme_minimal(base_size = 28) +  # Increased base size for overall elements
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 38, face = "bold", hjust = 0.5),
    
    # --- BLACK & BIGGER AXES STYLING ---
    axis.line = element_line(color = "black", linewidth = 1),         # Solid black x & y axis lines
    axis.ticks = element_line(color = "black", linewidth = 1),        # Black axis tick marks
    axis.ticks.length = unit(0.3, "cm"),                              # Longer/bigger tick marks
    axis.text = element_text(size = 22, color = "black"),             # Black & larger tick labels
    axis.title = element_text(size = 26, color = "black", face = "bold") # Black & larger x/y titles
  ) +
  scale_fill_manual(values = c(
    "399" = "#984ea3",  # purple
    "42"  = "#ff7f00"   # orange
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA))

# Display the plot
abc

ggsave(
  filename = "../../Summaries and Outputs/3/58_bodymass.png",
  plot = abc,
  width = 20,
  height = 7,
  dpi = 300
)

