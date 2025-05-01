#_______________________________________________
# 0. Reading the data
#_______________________________________________
df <- read.csv("~/GitHub/Caloric-Density-and-Recombination/Raw Data/Pilots/Pilot_MLE_analyzed_dataset_diets.csv")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)
library(multcomp)
library(lme4)
library(lmerTest)
 # 2. Set as factors
 df$GB<- as.factor(df$GB)
 
 # 3. ANOVA
 anova_model<- aov(recombination_fractions ~ Brood * diet* GB, data = df)
 summary(anova_model)
 
 # 4. Post hoc Tukey
 TukeyHSD(anova_model)
 
 # 5. Estimated marginal means and CLD
 emm <- emmeans(anova_model, ~ Brood * diet* GB)
 cld <- cld(emm, Letters = letters)
 
 # 6. Visualization
 ggplot(data=df,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()
 ggplot(data=df,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()+ facet_wrap(~Brood)


#_______________________________________________
#                        399 with  6 treatmemt
#_______________________________________________

# 1. Filter genotype 399
df_399 <- df %>% filter(GB == 399)

# 2. Set Brood and diet as factors
df_399$Brood <- as.factor(df_399$Brood)
df_399$diet <- as.factor(df_399$diet)

# 3. ANOVA
anova_model_399 <- aov(recombination_fractions ~ Brood * diet, data = df_399)
summary(anova_model_399)

# 4. Post hoc Tukey
TukeyHSD(anova_model_399)

# 5. Estimated marginal means and CLD
emm_399 <- emmeans(anova_model_399, ~ Brood * diet)
cld_399 <- cld(emm_399, Letters = letters)

# 6. Visualization
ggplot(data=df_399,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()
ggplot(data=df_399,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()+ facet_wrap(~Brood)


#____________________________________________
#                       all that have 0.5x 1x 2x diets only 42,217,391,399,491
#_____________________________________________

# 1. Filter genotypes and diets of interest
df_1X <- df %>%
  filter(
    GB %in% c(42, 217, 391, 399, 491),
    diet %in% c("0.5x", "1.0x", "2.0x")
    )

# 2. Set Brood and diet as factors
df_1X$Brood <- as.factor(df_1X$Brood)
df_1X$diet <- as.factor(df_1X$diet)

# 3. Run ANOVA
anova_model <- aov(recombination_fractions ~ Brood * diet, data = df_1X)
summary(anova_model)

# 4. Post hoc Tukey test
TukeyHSD(anova_model)

# 5. Estimated marginal means and compact letter display (CLD)
emm <- emmeans(anova_model, ~ Brood * diet)
cld_results <- cld(emm, Letters = letters)

# 6. Visualization

# Boxplot by diet
ggplot(data = df_1X, aes(x = diet, y = recombination_fractions, color = diet)) +
  geom_boxplot() +
  theme_minimal()

# Boxplot faceted by Brood
ggplot(data = df_1X, aes(x = diet, y = recombination_fractions, color = diet)) +
  geom_boxplot() +
  facet_wrap(~Brood) +
  theme_minimal()

#box plot of genotype
#ggplot(data = df_1X, aes(x = diet, y = recombination_fractions, color = diet)) +
  geom_boxplot() +
  facet_wrap(~GB) +
  theme_minimal() +
  labs(title = "Recombination Fractions by Diet and Genotype",
       x = "Diet", y = "Recombination Fraction")

ggplot(data = df_1X, aes(x = as.factor(GB), y = recombination_fractions, fill = diet)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "Recombination Fractions by Genotype and Diet",
       x = "Genotype (GB)", y = "Recombination Fraction")



#____________________________________________________________________
#    only 1x
#___________________________________________________________________
# 1. Filter 1X  
df_1X <- df %>%
  filter(
    diet %in% c("1.0x")
  )

# 2. Set Brood and diet as factors
df_1X$Brood <- as.factor(df_1X$Brood)
df_1X$diet <- as.factor(df_1X$diet)
df_1X$GB <- as.factor(df_1X$GB)


# 3. Run ANOVA
anova_model <- aov(recombination_fractions ~ Brood * GB, data = df_1X)
summary(anova_model)

# 4. Post hoc Tukey test
TukeyHSD(anova_model)

# 5. Estimated marginal means and compact letter display (CLD)
emm <- emmeans(anova_model, ~ Brood * GB)
cld_results <- cld(emm, Letters = letters)

# 6. Visualization

# Boxplot by gb
ggplot(data = df_1X, aes(x = GB, y = recombination_fractions, color = GB)) +
  geom_boxplot() +
  theme_minimal()

# Boxplot faceted by Brood
ggplot(data = df_1X, aes(x = GB, y = recombination_fractions, color = GB)) +
  geom_boxplot() +
  facet_wrap(~Brood) +
  theme_minimal()

