#_______________________________________________
# 0. Reading the data
#_______________________________________________
df <- read.csv("~/GitHub/Caloric-Density-and-Recombination/Raw Data/Pilots/Pilot_MLE_analyzed_dataset_diets.csv")
observer <- read.csv("~/GitHub/Caloric-Density-and-Recombination/Raw Data/Pilots/Recombination Rate/observer.csv")

df <- df %>%
  left_join(observer, by = "VialID")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)
library(multcomp)
library(lme4)
library(lmerTest)
library(broom)

 # 2. Set as factors
 df$GB<- as.factor(df$GB)
 
 # 3. ANOVA
 anova_model<- aov(recombination_fractions ~ Brood * diet* GB, data = df)
 summary(anova_model)
 
# Tidy the ANOVA model and save to CSV
anova_summary <- tidy(anova_model)
write.csv(anova_summary, "~/GitHub/Caloric-Density-and-Recombination/Statistical Analysis/anova_model.csv", row.names = FALSE)

 # 4. Post hoc Tukey
 TukeyHSD(anova_model)
 
 # 5. Estimated marginal means and CLD
 emm <- emmeans(anova_model, ~ Brood * diet* GB)
 cld <- cld(emm, Letters = letters)
 
 # 6. Visualization
 ggplot(data=df,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()
 ggplot(data=df,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()+ facet_wrap(~GB)
 ggplot(df, aes(x = Brood, y = recombination_fractions, group = diet, color = diet)) +
   stat_summary(fun = mean, geom = "line", size = 1) +             # Line plot with mean values
   stat_summary(fun = mean, geom = "point", size = 2) +            # Points on the line
   labs(
     title = "Recombination Fractions by Brood and Diet",
     x = "Brood",
     y = "Recombination Fraction",
     color = "Diet"
   ) +
   theme_minimal()
 


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
anova_summary <- tidy(anova_model_399)
write.csv(anova_summary, "~/GitHub/Caloric-Density-and-Recombination/Statistical Analysis/RR-399anova_model.csv", row.names = FALSE)

# 4. Post hoc Tukey
TukeyHSD(anova_model_399)

# 5. Estimated marginal means and CLD
emm_399 <- emmeans(anova_model_399, ~ Brood * diet)
cld_399 <- cld(emm_399, Letters = letters)

# 6. Visualization
ggplot(data=df_399,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()
ggplot(data=df_399,aes(x=diet,y=recombination_fractions, color=diet))+geom_boxplot()+ facet_wrap(~Brood)

ggplot(df_399, aes(x = Brood, y = recombination_fractions, group = diet, color = diet)) +
  stat_summary(fun = mean, geom = "line", size = 1) +             # Line plot with mean values
  stat_summary(fun = mean, geom = "point", size = 2) +            # Points on the line
  labs(
    title = "Recombination Fractions by Brood and Diet",
    x = "Brood",
    y = "Recombination Fraction",
    color = "Diet"
  ) +
  theme_minimal()

#____________________________________________
#                       all that have 0.5x 1x 2x diets only 42,217,391,399,491
#_____________________________________________

# 1. Filter genotypes and diets of interest

df_3diet <- df %>%
  filter(
    GB %in% c(42, 217, 391, 399, 491),
    diet %in% c("0.5x", "1.0x", "2.0x")
    )

# 2. Set Brood and diet as factors
df_3diet$Brood <- as.factor(df_3diet$Brood)
df_3diet$diet <- as.factor(df_3diet$diet)
df_3diet$GB = as.factor(df_3diet$GB)

# 3. Run ANOVA
anova_model <- aov(recombination_fractions ~ Brood * diet* GB, data = df_3diet)
summary(anova_model)
write.csv(anova_summary, "~/GitHub/Caloric-Density-and-Recombination/Statistical Analysis/anova_model5genotype.csv",row.names = FALSE)


# 4. Post hoc Tukey test
TukeyHSD(anova_model)

# 5. Estimated marginal means and compact letter display (CLD)
emm <- emmeans(anova_model, ~ GB* diet)
cld_results <- cld(emm, Letters = letters)

# 6. Visualization

# Boxplot by diet
ggplot(data = df_3diet, aes(x = diet, y = recombination_fractions, color = diet)) +
  geom_boxplot() +
  theme_minimal()

# Boxplot faceted by Brood
ggplot(data = df_3diet, aes(x = diet, y = recombination_fractions, color = diet)) +
  geom_boxplot() +
  facet_wrap(~GB) +
  theme_minimal()

ggplot(data = df_3diet, aes(x = as.factor(GB), y = recombination_fractions, fill = diet)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "Recombination Fractions by Genotype and Diet",
       x = "Genotype (GB)", y = "Recombination Fraction")


ggplot(df_3diet, aes(x = Brood, y = recombination_fractions, group = diet, color = diet)) +
  stat_summary(fun = mean, geom = "line", size = 1) +             # Line plot with mean values
  stat_summary(fun = mean, geom = "point", size = 2) +            # Points on the line
  labs(
    title = "Recombination Fractions by Brood and Diet",
    x = "Brood",
    y = "Recombination Fraction",
    color = "Diet"
  ) +
  facet_wrap(~GB) +
  theme_minimal()

#____________________________________________________________________
#    only 1x
#___________________________________________________________________
# 1. Filter 1X  
df_1x <- df %>%
  filter(
    diet %in% c("1.0x")
  )

# 2. Set Brood and diet as factors
df_1x$Brood <- as.factor(df_1x$Brood)
df_1x$diet <- as.factor(df_1x$diet)
df_1x$GB <- as.factor(df_1x$GB)


# 3. Run ANOVA
anova_model1X <- aov(recombination_fractions ~ Brood * GB, data = df_1x)
summary(anova_model1X)
anova_summary <- tidy(anova_model1X)
write.csv(anova_summary, "~/GitHub/Caloric-Density-and-Recombination/Statistical Analysis/RR-1Xanova_model.csv", row.names = FALSE)

# 4. Post hoc Tukey test
TukeyHSD(anova_model)

# 5. Estimated marginal means and compact letter display (CLD)
emm <- emmeans(anova_model, ~ Brood * GB)
cld_results <- cld(emm, Letters = letters)

# 6. Visualization

# Boxplot by gb
ggplot(data = df_1x, aes(x = GB, y = recombination_fractions, color = GB)) +
  geom_boxplot() +
  theme_minimal()

# Boxplot faceted by Brood
ggplot(data = df_1x, aes(x = GB, y = recombination_fractions, color = GB)) +
  geom_boxplot() +
  facet_wrap(~Brood) +
  theme_minimal()


ggplot(df_1x, aes(x = Brood, y = recombination_fractions, group = diet, color = diet)) +
  stat_summary(fun = mean, geom = "line", size = 1) +             # Line plot with mean values
  stat_summary(fun = mean, geom = "point", size = 2) +            # Points on the line
  labs(
    title = "Recombination Fractions by Brood and Diet",
    x = "Brood",
    y = "Recombination Fraction",
    color = "Diet"
  ) +
  facet_wrap(~GB) +
  theme_minimal()


