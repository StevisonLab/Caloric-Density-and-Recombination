setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs")

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

df<- read.csv("Nutrition_countingmodel_outputs_crosspooled.csv")

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
df$Brood <- as.factor(df$Brood)
df$Diet <- as.factor(df$Diet)
df$GB <- as.factor(df$GB)
df$People <- as.factor(df$People)
df$recombination_fractions <- as.numeric(df$recombination_fractions)

###delete the outlier (just checking)
#statisticaly
#df <- df %>%
  group_by(Diet, GB) %>%
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

df <- df %>%filter(recombination_fractions < 1000)  # example threshold
df <- df %>%filter(People %in% c("M","L","S")) # example threshold

df <- df %>%filter(!(Diet == "0.5X" & People == "C"))


# 3. ANOVA
anova_result <- aov(recombination_fractions ~  Brood * Diet * GB * People, data = df)
anova_result <- aov(recombination_fractions ~ People, data = df)

anova_result <- aov(recombination_fractions ~ Brood * Diet * GB, data = df)
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
#write.csv(anova_summary, "../Statistical Analysis/anova_model_recombination.csv", row.names = FALSE)

# 4. Post hoc Tukey
TukeyHSD(anova_result)

# 5. Estimated marginal means and CLD
emm <- emmeans(anova_result, ~ People* GB)
cld <- cld(emm, Letters = letters)

plot_data <- cld %>%
  mutate(
    GB = factor(GB),
    People= factor(People),
    #Diet = factor(Diet, levels = c("0.5x", "1.0x", "2.0x")),
    .group = trimws(.group) # remove extra spaces
  )
# Plot
ggplot(plot_data, aes(x = People, y = emmean, color = GB, group = GB)) +
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
ggplot(data=df,aes(x=Diet,y=recombination_fractions, color=Diet))+geom_boxplot()
ggplot(data=df,aes(x=Diet,y=recombination_fractions, color=Diet))+geom_boxplot()+ facet_wrap(~GB)
ggplot(df, aes(x = Brood, y = recombination_fractions, group = Diet, color = Diet)) +
  stat_summary(fun = mean, geom = "line", size = 1) +             # Line plot with mean values
  stat_summary(fun = mean, geom = "point", size = 2) +            # Points on the line
  labs(
    title = "Recombination Fractions by Brood and Diet",
    x = "Brood",
    y = "Recombination Fraction",
    color = "Diet"
  ) +
  theme_minimal()

#per strain brood 
df$Diet_GB <- paste0(df$GB, "_", df$Diet)

ggplot(df, aes(x = Brood, y = recombination_fractions, group = Diet_GB, color = Diet_GB)) +
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

#____________________________________________
#                       only 0.5x 2x
#_____________________________________________

# 1. Filter Diets of interest

df_2Diet <- df %>%
  filter(
    Diet %in% c("0.5x", "2.0x")
  )

# 2. Set Brood and Diet as factors
df_2Diet$Brood <- as.factor(df_2Diet$Brood)
df_2Diet$Diet <- as.factor(df_2Diet$Diet)
df_2Diet$GB = as.factor(df_2Diet$GB)

# 3. Run ANOVA
anova_model <- aov(recombination_fractions ~ Brood * Diet* GB, data = df_2Diet)
summary(anova_model)
#write.csv(anova_summary, "~/GitHub/Caloric-Density-and-Recombination/Statistical Analysis/anova_model5genotype.csv",row.names = FALSE)


# 4. Post hoc Tukey test
TukeyHSD(anova_model)

# 5. Estimated marginal means and compact letter display (CLD)
emm <- emmeans(anova_model, ~ GB* Diet)
cld_results1 <- cld(emm, Letters = letters)

# 6. Visualization

# Boxplot by Diet
ggplot(data = df_2Diet, aes(x = Diet, y = recombination_fractions, color = Diet)) +
  geom_boxplot() +
  theme_minimal()

# Boxplot faceted by Brood

ggplot(data = df_2Diet, aes(x = as.factor(GB), y = recombination_fractions, fill = Diet)) +
  geom_boxplot(position = position_dodge()) +
  theme_minimal() +
  labs(title = "Recombination Fractions by Genotype and Diet",
       x = "Genotype (GB)", y = "Recombination Fraction")


ggplot(df_2Diet, aes(x = Brood, y = recombination_fractions, group = Diet, color = Diet)) +
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







