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
library(emmeans)
library(multcompView)  # For cld()


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
#reading the female in cross (for egg per female)
femalescross <- read_csv("../Raw Data/Pilots/Recombination Rate/Femaleinvial.csv")
colnames(femalescross)

femalescross <- femalescross %>%
  rename(VialID = `Vial ID`)

femalescross <- femalescross %>%
  mutate( VialID = as.character(VialID) )

# Merge egg count + data 
merged_data <- all_phenotypes %>%
  left_join(egg_count, by = c("VialID", "Brood"))
 
#treatment key + merged_data
merged_data <- merged_data %>%
  left_join(treatment_key, by = "VialID")

#femalescross+ merged_data
merged_data <- merged_data %>%
  left_join(femalescross, c("VialID", "Brood"))

#sum the total of same class (same Genotype and sex)of each vialId and brood in one 
#(REMOVING DATES,NOTES, source)

merged_data <- merged_data %>%
  select(-Date,-Source, -Notes, ) %>%  # Remove the date column
  group_by(across(-Count)) %>%  # Group by all columns except Count
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")

# creating the egg per female column
merged_data$Eggperfemale <- merged_data$Egg_sumcount / merged_data$Female

# Save merged_data as CSV
write_csv(merged_data, "../Raw Data/Pilots/Recombination Rate/merged_data_output.csv")


#_____________________________________________________________________________
#                                             General
#______________________________________________________________________________
######total scored flies
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
                
#_______________________________________________________________________________________
#                                               Fecundity
#_______________________________________________________________________________________

###how many egg per female in each treatment

fecundity_data <- merged_data %>%
  group_by(VialID, Brood) %>%
  mutate(fecundity = Egg_sumcount / Female) %>%
  slice(1) %>%
  dplyr::select(VialID, Brood,`F1 Genotype`, `F1 developmental`, fecundity) %>%
  ungroup()

# 3. ANOVA
anova_model_Fecundity<- aov(fecundity ~ `F1 Genotype` * `F1 developmental`, data = fecundity_data)
summary(anova_model_Fecundity)

#Estimated marginal means and CLD
emm_Fecundity <- emmeans(anova_model_Fecundity, ~ `F1 Genotype` * `F1 developmental`)
cld_Fecundity <- cld(emm_Fecundity, Letters = letters)
print(cld_Fecundity)

#Box plot $ diet
ggplot(fecundity_data, aes(x = `F1 developmental`, y = fecundity, fill = `F1 Genotype`)) +
  geom_boxplot() +
  labs(title = "Fecundity by Treatment and Genotype",
       x = "Treatment",
       y = "Fecundity (Eggs per Female)") +
  theme_minimal()

#Box plot $ GB
ggplot(fecundity_data, aes(x = `F1 Genotype`, y = fecundity, fill = `F1 developmental`)) +
  geom_boxplot() +
  labs(title = "Fecundity by Treatment and Genotype",
       x = "Treatment",
       y = "Fecundity (Eggs per Female)") +
  theme_minimal()
#_______________________
#         399
#________________________
# Genotype 399 
fecundity_data_399 <- fecundity_data %>%
  filter(`F1 Genotype` == 399 )

# Box plot: fecundity for Genotype 399 with 0.5x, 1x, and 2x treatments
ggplot(fecundity_data_399, aes(x = `F1 developmental`, y = fecundity, fill = `F1 developmental`)) +
  geom_boxplot() +
  labs(title = "Fecundity for Genotype 399 by Treatment",
       x = "Treatment (F1 Developmental)",
       y = "Fecundity (Eggs per Female)") +
  theme_minimal()
#_______________________
#   0.5x, 1.0x, and 2.0x
#________________________
# Subset data for Treatments 0.5x, 1.0x, and 2.0x
fecundity_data_3diet <- fecundity_data %>%
  filter(`F1 developmental` %in% c("0.5x", "1.0x", "2.0x") & `F1 Genotype` %in% c(42, 217, 391, 399, 491),)

# Box plot: with 0.5x, 1x, and 2x treatments
ggplot(fecundity_data_3diet, aes(x = `F1 developmental`, y = fecundity, fill = `F1 Genotype`)) +
  geom_boxplot() +
  labs(title = "Fecundity for 3 diet by Treatment",
       x = "Treatment (F1 Developmental)",
       y = "Fecundity (Eggs per Female)") +
  theme_minimal()
#_______________________
#            1x
#________________________
# Subset data for Treatment 1x
fecundity_data_1x <- fecundity_data %>%
  filter(`F1 developmental` == "1.0x")

# Box plot: only for 1.0x treatment across selected genotypes
ggplot(fecundity_data_1x, aes(x = as.factor(`F1 Genotype`), y = fecundity, fill = as.factor(`F1 Genotype`))) +
  geom_boxplot() +
  labs(title = "Fecundity for Genotypes under 1.0x Treatment",
       x = "F1 Genotype",
       y = "Fecundity (Eggs per Female)") +
  theme_minimal()
#________________________________________________________________________________
#                             survival>1    create the total adult count
#_______________________________________________________________________________

merged_data<- merged_data %>%
  group_by(VialID, Brood) %>%
  mutate(total_adult = sum(Count, na.rm = TRUE)) %>%
  ungroup()

###Filter rows where egg-sumcount is less than total_adult 

merged_problematic <- merged_data %>%
  filter(Egg_sumcount < total_adult) %>%
  dplyr::select(VialID, Brood, Egg_sumcount, total_adult) %>%
  distinct()

#write_csv(merged_data, "../Raw Data/Pilots/Recombination Rate/merged_problematic.csv")


#_____________________________________________________________________________________
#                                          Survival
#______________________________________________________________________________________

#Calculate survival alive +dead

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
    Survival = Total_Count / Egg_Count,
    Survival = if_else(Survival > 1, 1, Survival)  # Cap Survival at 1
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
#___________________________________________________
#########zl by sex f+half dead ############## M+HALF OF DEAD
#__________________________________________________
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
    Survival_F = pmin(adjusted_F / expected_sex, 1), #not bigger than 1
    Survival_M = pmin(adjusted_M / expected_sex, 1)
  )
# Convert to long format for plotting
survival_long <- survival_sex_summary %>%
  dplyr::select(`F1 developmental`, `F1 Genotype`, Survival_F, Survival_M) %>%
  pivot_longer(cols = c(Survival_F, Survival_M), names_to = "Sex", values_to = "Survival") %>%
  mutate(Sex = ifelse(Sex == "Survival_F", "Female", "Male"))


# Plot
ggplot(survival_long, aes(x = `F1 developmental`, y = Survival, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_wrap(~ `F1 Genotype`) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Survival by Treatment, Genotype, and Sex",
    x = "Treatment",
    y = "Survival Rate"
  ) +
  scale_fill_manual(values = c("Female" = "#F0899C", "Male" = "lightblue")) +
  ylim(0, 1) +  # this line sets the y-axis from 0 to 1
  theme(
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )

#________________________
# survival 399 => half F and M
#________________________

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
    female_survival = pmin(adjusted_F / expected_per_sex, 1),
    male_survival = pmin(adjusted_M / expected_per_sex, 1)#not bigger than 1
  )

# Convert to long format for combined plot
survival_long_399 <- survival_399 %>%
  select(`F1 developmental`, female_survival, male_survival) %>%
  pivot_longer(cols = c(female_survival, male_survival),
               names_to = "Sex", values_to = "Survival") %>%
  mutate(Sex = ifelse(Sex == "female_survival", "Female", "Male"))

# Combined plot (female + male)
ggplot(survival_long_399, aes(x = `F1 developmental`, y = Survival, fill = Sex)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(
    title = "Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Survival Rate"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Female" = "#F0899C", "Male" = "lightblue")) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "#f0f0f0"),
    strip.text = element_text(face = "bold")
  )

# Female-only plot
ggplot(survival_399, aes(x = `F1 developmental`, y = female_survival, fill = `F1 developmental`)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Female Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Female Survival Rate"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "RdPu") +
  theme(legend.position = "none")

# Male-only plot
ggplot(survival_399, aes(x = `F1 developmental`, y = male_survival, fill = `F1 developmental`)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Male Survival Rate by Treatment (Genotype 399)",
    x = "Food Treatment",
    y = "Male Survival Rate"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Blues") +
  theme(legend.position = "none")

# ANOVA - Female
anova_f <- aov(female_survival ~ `F1 developmental`, data = survival_399)
summary(anova_f)
TukeyHSD(anova_f)

# ANOVA - Male
anova_m <- aov(male_survival ~ `F1 developmental`, data = survival_399)
summary(anova_m)
TukeyHSD(anova_m)


