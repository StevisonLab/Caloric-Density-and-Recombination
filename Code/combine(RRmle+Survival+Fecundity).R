df <- read.csv("~/GitHub/Caloric-Density-and-Recombination/Raw Data/Pilots/Pilot_MLE_analyzed_dataset_diets.csv")
merged_data <- read.csv("~/GitHub/Caloric-Density-and-Recombination/Raw Data/Pilots/Recombination Rate/merged_data_output.csv")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(multcompView)
library(tidyverse)
library(emmeans)    # For post-hoc tests
library(ggpubr)     # For publication-ready plots
library(car)        # For ANOVA diagnostics
library(MASS)       # For robust regression (if needed)
library(multcomp)
library(broom)  # for tidy model summaries



# First, standardize the column names between the two datasets

df <- df %>%
  rename(
    VialID = vialID,
    F1.developmental = diet,
    F1.Genotype = GB,
    recombination_fraction = recombination_fractions
  )

# Create survival data from merged_data
survival_summary <- merged_data %>%
  group_by(VialID, Brood, F1.developmental, F1.Genotype) %>%
  summarise(
    # Sum Count where Sex is not NA (alive)
    FemaleMale = sum(Count[!is.na(Sex)], na.rm = TRUE),
    # Sum Count where Sex is NA (dead)
    dead = sum(Count[is.na(Sex)], na.rm = TRUE),
    # Egg count per vial
    Egg_Count = first(as.numeric(Egg_sumcount)),
    .groups = 'drop'
  ) %>%
  mutate(
    Total_Count = FemaleMale + dead,
    Survival = Total_Count / Egg_Count,
    Survival = if_else(Survival > 1, 1, Survival)  # Cap Survival at 1
  )

# Create fecundity data from merged_data
fecundity_data <- merged_data %>%
  group_by(VialID, Brood) %>%
  summarise(
    F1.developmental = first(F1.developmental),
    F1.Genotype = first(F1.Genotype),
    fecundity = first(Egg_sumcount) / first(Female),
    .groups = 'drop'
  )

# Now combine all three datasets
combined_data <- survival_summary %>%
  full_join(fecundity_data, by = c("VialID", "Brood", "F1.developmental", "F1.Genotype")) %>%
  full_join(df, by = c("VialID", "Brood", "F1.developmental", "F1.Genotype"))

# Convert genotype to factor for analysis
combined_data <- combined_data %>%
  mutate(F1.Genotype = as.factor(F1.Genotype),
         F1.developmental = as.factor(F1.developmental))

# Create analysis subsets
selected_genotypes <- c("42", "217", "391", "399", "491")
selected_treatments <- c("0.5x", "1.0x", "2.0x")

# 1. Full dataset
#full_dataset <- combined_data

# 2. Genotype 399 only
#dataset_399 <- combined_data %>% 
# filter(F1.Genotype == "399")

# 3. Three diets and five genotypes
dataset_3diet_5geno <- combined_data %>%
  filter(F1.Genotype %in% selected_genotypes,
         F1.developmental %in% selected_treatments)

# 4. 1.0x treatment only
#dataset_1x <- combined_data %>% 
#  filter(F1.developmental == "1.0x")

#______________________________________________________________________________
# Convert factors if needed
dataset_3diet_5geno$F1.developmental <- as.factor(dataset_3diet_5geno$F1.developmental)
dataset_3diet_5geno$F1.Genotype <- as.factor(dataset_3diet_5geno$F1.Genotype)

#Model 1: Main effects only
anova_model1 <- aov(recombination_fraction ~ F1.developmental + F1.Genotype + Survival + fecundity,
                    data = dataset_3diet_5geno)
summary(anova_model1)
#Model 2: With interactions
anova_model2 <- aov(recombination_fraction ~ F1.developmental * F1.Genotype * Survival * fecundity,
                    data = dataset_3diet_5geno)
summary(anova_model2)

# For developmental diets
diet_means <- emmeans(anova_model2, ~ F1.developmental)
diet_cld <- cld(diet_means, Letters = letters, adjust = "tukey")
print(diet_cld)

#For genotypes:
geno_means <- emmeans(anova_model2, ~ F1.Genotype)
geno_cld <- cld(geno_means, Letters = letters, adjust = "tukey")
print(geno_cld)
#For interaction between diet and genotype:
diet_geno_means <- emmeans(anova_model2, ~ F1.developmental * F1.Genotype)
diet_geno_cld <- cld(diet_geno_means, Letters = letters, adjust = "tukey")
print(diet_geno_cld)

# --------------------------------------------------------------------------
# STEP 1: Data Preparation & Exploration
# --------------------------------------------------------------------------

# Convert factors
dataset_clean <- dataset_3diet_5geno %>%
  mutate(
    F1.developmental = factor(F1.developmental),
    F1.Genotype = factor(F1.Genotype)
  )

# Check distribution
hist(dataset_clean$recombination_fraction, 
     main = "Distribution of Recombination Fraction",
     xlab = "Recombination Fraction")

# Check for zeros (log-transform if needed)
if(any(dataset_clean$recombination_fraction <= 0)){
  dataset_clean$log_rec <- log(dataset_clean$recombination_fraction + 0.001)
}

# --------------------------------------------------------------------------
# STEP 2: Main Analysis - Linear Model
# --------------------------------------------------------------------------

# Full model with interactions
model <- lm(recombination_fraction ~ F1.developmental * F1.Genotype + 
              Survival + fecundity,
            data = dataset_clean)

# Model summary
summary(model)

# ANOVA table (Type III SS)
Anova(model, type = 3) 

# --------------------------------------------------------------------------
# STEP 3: Model Diagnostics
# --------------------------------------------------------------------------

# Check assumptions
par(mfrow = c(2,2))
plot(model)

# Formal tests
shapiro.test(resid(model))                       # Normality
ncvTest(model)                                   # Homoscedasticity
durbinWatsonTest(model)                          # Independence

# If assumptions violated:
if(ncvTest(model)$p < 0.05 | shapiro.test(resid(model))$p < 0.05){
  robust_model <- rlm(recombination_fraction ~ F1.developmental * F1.Genotype + 
                        Survival + fecundity,
                      data = dataset_clean)
  print(summary(robust_model))
}

# --------------------------------------------------------------------------
# STEP 4: Post-Hoc Analysis
# --------------------------------------------------------------------------

# Get marginal means
emm_diet <- emmeans(model, ~ F1.developmental)
emm_geno <- emmeans(model, ~ F1.Genotype)
emm_interaction <- emmeans(model, ~ F1.developmental * F1.Genotype)

# Tukey-adjusted comparisons
diet_cld <- cld(emm_diet, Letters = letters, adjust = "tukey")
geno_cld <- cld(emm_geno, Letters = letters, adjust = "tukey")
interaction_cld <- cld(emm_interaction, Letters = letters, adjust = "tukey")

print(diet_cld)
print(geno_cld)
print(interaction_cld)

# --------------------------------------------------------------------------
# STEP 5: Visualization
# --------------------------------------------------------------------------

# Main effects plot

p1 <- ggplot(dataset_clean, aes(F1.developmental, recombination_fraction)) +
  geom_boxplot(aes(fill = F1.developmental)) +
  geom_text(data = diet_cld, aes(y = max(dataset_clean$recombination_fraction)*1.1, 
                                 label = .group)) +
  labs(x = "Diet", y = "Recombination Fraction")

# Interaction plot
p2 <- ggplot(dataset_clean, aes(F1.developmental, recombination_fraction, 
                                color = F1.Genotype, group = F1.Genotype)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  labs(x = "Diet", y = "Mean Recombination Fraction")

# Survival relationship
p3 <- ggplot(dataset_clean, aes(Survival, recombination_fraction, color = F1.Genotype, group = F1.Genotype)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "Survival", y = "Recombination Fraction")

#p3 <- ggplot(dataset_clean, aes(Survival, recombination_fraction)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(x = "Survival", y = "Recombination Fraction")

# Arrange plots
ggarrange(p1, p2, p3, ncol = 3, widths = c(1, 1.2, 1))


#___________________
#Color by Diet 
p3_diet <- ggplot(dataset_clean, 
                  aes(x = Survival, 
                      y = recombination_fraction, 
                      color = F1.developmental,    # Color by diet
                      group = F1.developmental)) + # Group by diet
  geom_point(alpha = 0.6) +                       # Semi-transparent points
  geom_smooth(method = "lm", se = TRUE) +          # Add regression lines with CI
  labs(x = "Survival", 
       y = "Recombination Fraction",
       color = "Diet") +                          # Better legend title
  theme_minimal() +
  scale_color_viridis_d()                          # Colorblind-friendly palette

print(p3_diet)

#Facet by Diet 
p3_faceted <- ggplot(dataset_clean, 
                     aes(x = Survival, 
                         y = recombination_fraction, 
                         color = F1.Genotype)) +   # Color by genotype
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +         # Simpler lines without CI
  facet_wrap(~ F1.developmental) +                 # Separate panels by diet
  labs(x = "Survival", 
       y = "Recombination Fraction",
       color = "Genotype") +
  theme_bw()

print(p3_faceted)

#Diet-Specific Interactions
p3_interaction <- ggplot(dataset_clean, 
                         aes(x = Survival, 
                             y = recombination_fraction,
                             color = F1.Genotype,
                             linetype = F1.developmental)) +  # Diet as linetype
 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Survival", 
       y = "Recombination Fraction",
       color = "Genotype",
       linetype = "Diet") +
  scale_color_brewer(palette = "Set1") +           # Distinct colors
  theme_minimal()

print(p3_interaction)


#__________________________
#slope of RR vs survival plotted against mean survival and fecundity per strain
#__________________
# Step 1: Fit separate linear models for each genotype and extract slope

slopes_by_genotype <- dataset_clean %>%
  group_by(F1.Genotype) %>%
  do(tidy(lm(recombination_fraction ~ Survival, data = .))) %>%
  ungroup() %>%
  filter(term == "Survival") %>%
  rename(Slope = estimate) %>%
  dplyr::select(F1.Genotype, Slope)  # explicitly call dplyr::select

# Step 2: Get mean survival and fecundity per genotype
means_by_genotype <- dataset_clean %>%
  group_by(F1.Genotype) %>%
  summarise(
    Mean_Survival = mean(Survival, na.rm = TRUE),
    Mean_Fecundity = mean(fecundity, na.rm = TRUE)
  )

# Step 3: Merge slopes and means into one table
summary_table <- left_join(slopes_by_genotype, means_by_genotype, by = "F1.Genotype")

# View the summary table
print(summary_table)

# A. Slope vs. Mean Survival
ggplot(summary_table, aes(x = Mean_Survival, y = Slope)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = F1.Genotype), vjust = -0.6, size = 4) +
  labs(
    title = "Slope of RR ~ Survival vs Mean Survival",
    x = "Mean Survival",
    y = "Slope (Recombination ~ Survival)"
  ) +
  theme_minimal()


#🔹 B. Slope vs. Mean Fecundity
ggplot(summary_table, aes(x = Mean_Fecundity, y = Slope)) +
  geom_point(size = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  geom_text(aes(label = F1.Genotype), vjust = -0.6, size = 4) +
  labs(
    title = "Slope of RR ~ Survival vs Mean Fecundity",
    x = "Mean Fecundity",
    y = "Slope (Recombination ~ Survival)"
  ) +
  theme_minimal()


#🔹 C. Optional: Color by one and plot the other
lm_summary_surv <- lm(Slope ~ Mean_Survival, data = summary_table)
summary(lm_summary_surv)

lm_summary_fec <- lm(Slope ~ Mean_Fecundity, data = summary_table)
summary(lm_summary_fec)


