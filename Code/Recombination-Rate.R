#set the working Directory
#setwd("C:/Users/Melika Shiran/Box/Plasticity and Nutrition NEW (Melika)/03. Raw Data and Analysis/2_RR and ZL Marker present/Data_Analysed(csv)")
setwd("~/GitHub/Caloric-Density-and-Recombination/Code")

# Read and combine all files
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(lme4)
library(lmerTest) # For p-values in mixed models
library(emmeans)  # For post-hoc tests

# SetWD
#Laurie version
#setwd("/Users/lss0021/Documents/GitHub/Caloric-Density-and-Recombination/Code")
#Melika version


# Combining all the phenotypes
phenotype_files <- list.files(path = "../Raw Data/Pilots/Recombination Rate/Phenotypes", pattern = "\\.csv$", full.names = TRUE)

all_phenotypes <- phenotype_files %>%
  lapply(function(file) {
    df <- read_csv(file, col_types = cols(.default = "c"))
    df$Source <- basename(file)
    return(df)
  }) %>%
  bind_rows()

# Check the combined data
head(all_phenotypes)
colnames(all_phenotypes)

# Clean phenotype data
all_phenotypes <- all_phenotypes %>%
  rename(VialID = `Vial ID`)
#Deleting the rows that have Na in initials (deleting extras)
all_phenotypes <- all_phenotypes %>%
  filter(!is.na(Initials) & Initials != "")

#deleting rows that have na or 0 in theirs counts all_phenotypes 
all_phenotypes <- all_phenotypes %>%
  filter(
    !is.na(Count),                # Remove NAs
    !str_detect(Count, "o")       # Remove rows with "o" in Count
  )

#making sure their numerics and characters
all_phenotypes <- all_phenotypes %>%
  mutate(
    Count = as.numeric(Count),
    VialID = as.character(VialID)
  )

# Read your treatment key
treatment_key <- read_csv("../Raw Data/Pilots/Recombination Rate/Treatmentskeys.csv")
treatment_key <- treatment_key %>%
  rename(VialID = `Vial ID`)

#make sure that it is chracter
treatment_key <- treatment_key %>%
  mutate(VialID = as.character(VialID), `F1 Genotype`= as.character(`F1 Genotype`) )

#reading egg count
egg_count <- read_csv("../Raw Data/Pilots/Recombination Rate/All-eggcount.csv")
colnames(egg_count)
#str(egg_count)
egg_count <- egg_count %>%
  rename(VialID = `Vial ID`)  # Fix if there's a space or typo

egg_count <- egg_count %>%
  mutate(
    Egg_sumcount = as.numeric(Egg_sumcount),
    VialID = as.character(VialID)
  )

# Merge egg count + data 
merged_data <- all_phenotypes %>%
  left_join(egg_count, by = c("VialID", "Brood"))
 
#treatment key + merged_data
merged_data <- merged_data %>%
  left_join(treatment_key, by = "VialID")

#sum the total of same class (same Genotype and sex)of each vialId and brood in one 
#(REMOVING DATES,NOTES, source)

merged_data <- merged_data %>%
  select(-Date,-Source, -Notes, ) %>%  # Remove the date column
  group_by(across(-Count)) %>%  # Group by all columns except Count
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

# Save merged_data as CSV
write_csv(merged_data, "../Raw Data/Pilots/Recombination Rate/merged_data_output.csv")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

######tOotal scored flies
total_scored <- sum(merged_data$Count, na.rm = TRUE)
print(total_scored)

write.csv(total_scored, "../Summary Tables/Totalfliesscored.csv")

#### how many each person scored?
# Group by Initials and sum Count
initials_count <- merged_data %>%
  group_by(Initials) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop")

# Print the result
print(initials_count)
####how many crossover each individual has counted??
total_crossovers_per_initial <- merged_data %>%
  mutate(Count = as.numeric(Count)) %>%
  mutate(Crossover_Events_Per_Fly = case_when(
    `crossover Classification` == "Non-Crossover" ~ 0,
    `crossover Classification` == "Single Crossover" ~ 1,
    `crossover Classification` == "Double Crossover" ~ 2,
    `crossover Classification` == "Triple Crossover" ~ 3,
    `crossover Classification` == "Quadruple Crossover"~ 4,
    TRUE ~ NA
  )) %>%
  mutate(Total_Crossover_Events_Per_Row = Crossover_Events_Per_Fly * Count) %>%
  group_by(Initials) %>%
  summarise(Total_Crossover_Events = sum(Total_Crossover_Events_Per_Row, na.rm = TRUE), .groups = "drop")

print(total_crossovers_per_initial)

####how may crossover each individual counted per group 

crossover_counts_per_initial <- merged_data %>%
  group_by(Initials, `crossover Classification`) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop")

print(crossover_counts_per_initial)
#################################################################
###how many egg in each treatment:
# Prepare and join treatment info
egg_with_treatment<- egg_count %>%
  left_join(treatment_key, by = "VialID") %>%
  mutate(
    Egg_sumcount = as.numeric(Egg_sumcount),
    Brood = as.factor(Brood),
    `F1 developmental` = as.factor(`F1 developmental`)
  )

# Create faceted boxplot across Broods
ggplot(egg_with_treatment, aes(x = `F1 developmental`, y = Egg_sumcount, fill = `F1 developmental`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.5) +
  facet_wrap(~ Brood, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Egg Counts by Food Treatment Across Broods",
    x = "Food Treatment",
    y = "Egg Count"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )
################################################################
##### day of ecolsion in different broods 
# Summarize data
summary_counts <- merged_data %>%
  group_by(Date, Brood, `F1 developmental`) %>%
  summarise(mean_count = mean(Count, na.rm = TRUE), .groups = "drop")

# Line plot of mean counts over time
ggplot(summary_counts, aes(x = Date, y = mean_count, color = `F1 developmental`, group = `F1 developmental`)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ Brood, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean Progeny Count Over Time by Brood and Food Treatment",
    x = "Scoring Date",
    y = "Mean Progeny Count"
  ) +
  theme(strip.text = element_text(face = "bold"))

geom_smooth(method = "loess", se = FALSE)

# Convert to Date if needed
all$ScoringDate <- as.Date(all$ScoDate)

# Create a day number column
all <- all %>%
  group_by(Brood) %>%  # Optional: reset Day 1 per brood
  mutate(DayNum = as.integer(ScoringDate - min(ScoringDate)) + 1) %>%
  ungroup()
ggplot(all, aes(x = DayNum, y = ProgenyCount, color = `F1 developmental`)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ Brood) +
  labs(x = "Days After First Scoring", y = "Progeny Count")


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
########++++++++++++++++++Statistic RR+++++++++##################

# Data processing
rr_summary <- merged_data %>%
  filter(`crossover Classification` %in% c("Single Crossover", "Double Crossover", "Triple Crossover", "Non-Crossover")) %>% 
  group_by(Brood, VialID, `F1 developmental`, `F1 Genotype`, `crossover Classification`) %>%
  summarise(TotalCount = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = `crossover Classification`, values_from = TotalCount, values_fill = 0) %>%
  mutate(
    Total = `Single Crossover` + `Double Crossover` + `Non-Crossover` + `Triple Crossover`, 
    Recombinant = `Single Crossover` + `Double Crossover` + `Triple Crossover`,  
    RR = Recombinant / Total
  )

# Convert variables to factors
rr_summary <- rr_summary %>%
  mutate(
    `F1.developmental` = as.factor(`F1.developmental`),
    `F1.Genotype` = as.factor(`F1.Genotype`),
    Brood = as.factor(Brood)
  )

# Make column names syntactically valid
names(rr_summary) <- make.names(names(rr_summary))

# Mixed-effects model with Brood as random effect
mixed_model <- lmer(RR ~ F1.developmental * F1.Genotype + (1 | Brood), data = rr_summary)

# Summary of the model
summary(mixed_model)

# Post-hoc tests for fixed effects (if significant)
# For F1.Genotype
emmeans(mixed_model, pairwise ~ F1.Genotype, adjust = "tukey")

# For F1.developmental (uncomment if needed)
# emmeans(mixed_model, pairwise ~ F1.developmental, adjust = "tukey")

# Visualization
ggplot(rr_summary, aes(x = F1.developmental, y = RR, fill = Brood)) +
  geom_boxplot(position = position_dodge()) +
  facet_wrap(~ F1.Genotype) +
  theme_minimal(base_size = 14) +
  labs(title = "Recombination Rate by Treatment, Genotype, and Brood")



############RR###################
#grouping based on the class for RR

all_phenotypes_summary <- all_phenotypes %>%
  group_by(Brood, VialID, Class) %>%
  summarise(TotalCount = sum(as.numeric(Count), na.rm = TRUE), .groups = "drop")

############ZL###################
#Calculate survival
survival_summary <- merged_data %>%
  group_by(Brood, VialID, `F1 developmental`, `F1 Genotype`) %>%
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
    Survival = Total_Count / Egg_Count
  )
# Plot by treatment!

ggplot(survival_summary, aes(x = `F1 developmental`, y = Survival, fill = `F1 developmental`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Survival Rate by Treatment",
    x = "Treatment",
    y = "Survival Rate"
  )
###plot genotype and treatment
ggplot(survival_summary, aes(x = `F1 developmental`, y = Survival, fill = `F1 developmental`)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # cleaner boxplot
  geom_jitter(width = 0.2, alpha = 0.5, color = "black") +  # adds points
  facet_wrap(~ `F1 Genotype`, scales = "free_y") +  # one plot per genotype
  theme_minimal(base_size = 14) +
  labs(
    title = "Survival Rate by Treatment and Genotype",
    x = "Treatment",
    y = "Survival Rate"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )
#########zl by sex f+half dead ############## M+HALF OF DEAD

survival_sex_summary <- merged_data %>%
  group_by(Brood, VialID, `F1 developmental`, `F1 Genotype`) %>%
  summarise(
    live_F = sum(Count[Sex == "F"], na.rm = TRUE),
    live_M = sum(Count[Sex == "M"], na.rm = TRUE),
    dead = sum(Count[is.na(Sex)], na.rm = TRUE),
    Egg_Count = first(as.numeric(Egg_sumcount)),
    .groups = 'drop'
  ) %>%
  mutate(
    # Split the dead equally
    dead_per_sex = dead / 2,
    
    adjusted_F = live_F + dead_per_sex,
    adjusted_M = live_M + dead_per_sex,
    
    expected_sex = Egg_Count / 2,
    
    Survival_F = adjusted_F / expected_sex,
    Survival_M = adjusted_M / expected_sex
  )

######PLOT

# Convert to long format for ggplot
survival_long <- survival_sex_summary %>%
  select(`F1 developmental`, `F1 Genotype`, Survival_F, Survival_M) %>%
  pivot_longer(cols = c(Survival_F, Survival_M), names_to = "Sex", values_to = "Survival") %>%
  mutate(Sex = ifelse(Sex == "Survival_F", "Female", "Male"))

# Plot
ggplot(survival_long, aes(x = `F1 developmental`, y = Survival, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ `F1 Genotype`, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Survival by Treatment, Genotype, and Sex",
    x = "Treatment",
    y = "Survival Rate"
  ) +
  scale_fill_manual(values = c("Female" = "#F0899C", "Male" = "lightblue")) +
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )

#######################survival ship of only in 399 f+half of the dead //////////// M+half of the dead
# Filter for genotype 399
survival_filtered <- merged_data %>%
  filter(`F1 Genotype` == "399")

# Calculate survival with half dead assigned to each sex
survival_399 <- survival_filtered %>%
  group_by(Brood, VialID, `F1 developmental`, `F1 Genotype`) %>%
  summarise(
    live_F = sum(Count[Sex == "F"], na.rm = TRUE),
    live_M = sum(Count[Sex == "M"], na.rm = TRUE),
    dead = sum(Count[is.na(Sex)], na.rm = TRUE),
    Egg_Count = first(as.numeric(Egg_sumcount)),
    .groups = 'drop'
  ) %>%
  mutate(
    dead_half = dead / 2,
    adjusted_F = live_F + dead_half,
    adjusted_M = live_M + dead_half,
    expected_per_sex = Egg_Count / 2,
    female_survival = adjusted_F / expected_per_sex,
    male_survival = adjusted_M / expected_per_sex
  )
# For general use with ggplot and facetting by sex
survival_long_399 <- survival_399 %>%
  select(`F1 developmental`, female_survival, male_survival) %>%
  pivot_longer(cols = c(female_survival, male_survival), names_to = "Sex", values_to = "Survival") %>%
  mutate(Sex = ifelse(Sex == "female_survival", "Female", "Male"))
ggplot(survival_long_399, aes(x = `F1 developmental`, y = Survival, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(
    title = "Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Survival Rate"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Female" = "#F0899C", "Male" = "lightblue")) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )
###plot seprate F 399
ggplot(survival_399, aes(x = `F1 developmental`, y = female_survival, fill = `F1 developmental`)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Female Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Female Survival Rate"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "RdPu") +
  theme(legend.position = "none")

### plot seprate M
ggplot(survival_399, aes(x = `F1 developmental`, y = male_survival, fill = `F1 developmental`)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Male Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Male Survival Rate"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Blues") +
  theme(legend.position = "none")


# Female survival ANOVA
anova_f <- aov(female_survival ~ `F1 developmental`, data = survival_399)
summary(anova_f)
TukeyHSD(anova_f)

# Male survival ANOVA
anova_m <- aov(male_survival ~ `F1 developmental`, data = survival_399)
summary(anova_m)
TukeyHSD(anova_m)

