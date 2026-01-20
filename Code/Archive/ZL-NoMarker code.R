# Set working directory
setwd("C:/Users/Melika Shiran/Box/Plasticity and Nutrition NEW (Melika)/03. Raw Data and Analysis/1.ZL-NoMarker Piolot")
# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(readxl)
library(emmeans)
library(multcomp)
library(multcompView)
library(lme4)
library(lmerTest)
library(car)

# Load dataset
data <- read_excel("ZL-NoMarker PILOT DATA.xlsx")
colnames(data)

# Convert variables to factors
data <- data %>%
  mutate(Genotype = as.factor(Genotype),
         F1develop_env = as.factor(F1develop_env),
         F2develop_env = as.factor(F2develop_env),
         Brood = as.factor(Brood),
         Day = as.character(Day))

# Remove rows where any of the relevant columns have NA values(clearing the crosses that i have not done)
data_cleaned <- data %>%
  filter(!is.na(FemaleCount) & !is.na(MaleCount))

# Convert EggCount NA to 0
data_cleaned <- data_cleaned %>%
  mutate(EggCount = replace(EggCount, is.na(EggCount), 0))

data_cleaned <- data_cleaned %>%
  mutate(EggCount_piont = replace(Piont_click_eggcount, is.na(EggCount), 0))

# Ensure numeric conversion for all relevant columns
data_cleaned <- data_cleaned %>%
  mutate(
    EggCount = as.numeric(EggCount),
    EggCount_piont = as.numeric(Piont_click_eggcount),
    FemaleCount = as.numeric(FemaleCount),
    MaleCount = as.numeric(MaleCount),
    TotalOffspring = FemaleCount + MaleCount  # Compute TotalOffspring
  )
#adding new column for f1 and f2env
data_cleaned <- data_cleaned %>%
  mutate(Development_Env = paste(F1develop_env, F2develop_env, sep = "|"))

#Aggregate Data by VialID & Brood

data_summary <- data_cleaned %>%
  group_by(VialID, Brood, Genotype, F1develop_env, F2develop_env, Development_Env) %>%
  summarise(
    EggCount = sum(EggCount, na.rm = TRUE), 
    EggCount_piont = sum(EggCount_piont, na.rm = TRUE),       # Sum eggs per vial per brood
    FemaleCount = sum(FemaleCount, na.rm = TRUE), # Sum females per vial per brood
    MaleCount = sum(MaleCount, na.rm = TRUE),     # Sum males per vial per brood
    TotalOffspring = sum(FemaleCount, na.rm = TRUE) + sum(MaleCount, na.rm = TRUE)  # Compute total offspring
  ) %>%
  ungroup()

# Calculate Zygotic Lethality (ZL)
data_summary <- data_summary %>%
  mutate(ZygoticLethality = 100 - (TotalOffspring /EggCount) * 100)

data_summary <- data_summary %>%
  mutate(ZygoticLethality_piont = 100 - (TotalOffspring / EggCount_piont) * 100)

# Check and Save the Cleaned Data
summary(data_summary)  # Summary statistics
write.csv(data_summary, "Data_summary.csv", row.names = FALSE)

# Fit the model again with the updated summary data
###model <- lm(ZygoticLethality ~ Brood + Genotype + F1develop_env + Development_Env , data = data_summary)
##model <- lmer(ZygoticLethality ~ F2develop_env:F1develop_env + (1|Brood) + Genotype, data = data_summary)
#model <- lmer(ZygoticLethality ~ Development_Env + (1|Brood) + Genotype, data = data_summary)
#model_nested <- lm(ZygoticLethality ~ Brood + Genotype + F1develop_env / F2develop_env, data = data_summary)
                  # View the model summary
#summary(model)
#Anova(model)
#summary(model_nested)
#Anova(model_nested)

##piont and click: Fit the model again with the updated summary data
#model <- lmer(ZygoticLethality_piont ~ F2develop_env:F1develop_env + (1|Brood) + Genotype, data = data_summary)
#model <- lmer(ZygoticLethality_piont~ Development_Env + (1|Brood) + Genotype, data = data_summary)
#model_nested <- lm(ZygoticLethality_piont ~ Brood + Genotype + F1develop_env / F2develop_env, data = data_summary)
#F1develop_env:F2develop_env has a variance of 0, you should remove it from the model:
model_fixed <- lmer(ZygoticLethality_piont ~ F1develop_env + (1 | Brood) + Genotype, data = data_summary)

#summary(model_nested)
#Anova(model_nested)
summary(model_fixed)
Anova(model_fixed)

# Perform ANOVA on the model
anova_results <- Anova(model_fixed)

# View ANOVA results
print(anova_results)
write.csv(anova_results,file="Anova_table.csv")

################em_means
#em_means <- emmeans(model_nested, ~ Development_Env)
#cont_emm=contrast(em_means,method="pairwise")

em_means <- emmeans(model_fixed, ~ F1develop_env)
cont_emm=contrast(em_means,method="pairwise")

# View the EMMs results
output=summary(cont_emm)
write.csv(output,file="Posthoc_table.csv")

################contrast for Genotype
Genotype_contrast <- emmeans(model_fixed, "Genotype")
summary(Genotype_contrast)
pairs(Genotype_contrast)


#Boxplot of Zygotic Lethality by Genotype
ggplot(data_summary, aes(x = Genotype, y = ZygoticLethality_piont, fill = Genotype)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Zygotic Lethality by Genotype", 
       x = "Genotype", y = "Zygotic Lethality (%)")
#Boxplot of Zygotic Lethality by Developmental Environment
ggplot(data_summary, aes(x = F1develop_env, y = ZygoticLethality_piont, fill = F1develop_env)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Zygotic Lethality by F1 Developmental Environment", 
       x = "F1 Developmental Environment", y = "Zygotic Lethality (%)")


#Post Hoc Letters (Compact Letter Display)
# Post hoc letters for Genotype
out_put= Genotype_letters <- cld(Genotype_contrast, Letters = letters)
print(Genotype_letters)
write.csv(out_put,file="Genotype_posthoc.csv")


############Adult count/ engg count_pint and click.
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure necessary calculations
data_summary <- data_summary %>%
  mutate(
    AdultRatio = TotalOffspring / EggCount_piont  # Compute Adult Count / Egg Count
  )

# Create the boxplot
ggplot(data_summary, aes(x = F1develop_env, y = AdultRatio, fill = F1develop_env)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +  # Boxplots without outlier points
  geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.7) +  # Jittered raw data points
  stat_summary(fun = "mean", geom = "point", shape = 4, size = 3, color = "black") +  # Mean markers
  scale_fill_manual(values = c("purple", "orange", "brown")) +  # Adjust colors as needed
  theme_minimal() +
  theme(
    legend.position = "none",  # Hide legend since x-axis is F1develop_env
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Adult Count / Egg Count by F1 Development Environment",
    x = "F1 Development Environment",
    y = "Adult Count / Egg Count"
  )
