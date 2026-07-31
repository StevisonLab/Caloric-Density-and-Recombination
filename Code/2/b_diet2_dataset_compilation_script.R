
#setwd("Raw Data/") # Please use Raw Data as working directory
setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/2")

library(dplyr)
library(dbplyr)
library(emmeans)
library(multcompView)  # for cld() support if needed
library(multcomp)
#markerfree----
# Read in csv files of raw count data from the marker free cross
CH_marker_free_egg_count <- read.csv("A_CH_raw_marker_free_egg_counts.csv", header=TRUE)
LM_marker_free_egg_count <- read.csv("B_LM_raw_marker_free_egg_counts.csv", header=TRUE)
MS_marker_free_adult_count <- read.csv("C_MS_raw_marker_free_adult_counts.csv", header=TRUE)

# Compile the raw data from different observers and life history stages
CH_raw_marker_free_total_egg_count <- CH_marker_free_egg_count$total_egg_count
LM_raw_marker_free_total_egg_count <- LM_marker_free_egg_count$total_egg_count
compiled_raw_marker_free_dataset <- cbind(MS_marker_free_adult_count, CH_raw_marker_free_total_egg_count, LM_raw_marker_free_total_egg_count)

# Split treatment_id into treatment and genotype for marker-free dataset
treatment_genotype <- strsplit(as.character(compiled_raw_marker_free_dataset$treatment_id), "_")
compiled_raw_marker_free_dataset$treatment <- sapply(treatment_genotype, function(x) x[1])
compiled_raw_marker_free_dataset$genotype <- sapply(treatment_genotype, function(x) x[2])

write.csv(compiled_raw_marker_free_dataset, "../../Summaries and Outputs/2/D_compiled_raw_marker_free_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting the average and round these two counts to the nearest whole number
CH_transformed_marker_free_total_egg_count <- (CH_raw_marker_free_total_egg_count - 6.479979) / 0.9224861
LM_transformed_marker_free_total_egg_count <- (LM_raw_marker_free_total_egg_count - 5.346608 ) / 0.9101858        
transformed_mean_marker_free_total_egg_count <- round((CH_transformed_marker_free_total_egg_count + LM_transformed_marker_free_total_egg_count) / 2)

derived_marker_free_dataset <- cbind(compiled_raw_marker_free_dataset, transformed_mean_marker_free_total_egg_count)
write.csv(derived_marker_free_dataset, "../../Summaries and Outputs/2/E_derived_marker_free_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_marker_free_female_viability <- derived_marker_free_dataset$female_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_male_viability <- derived_marker_free_dataset$male_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_total_viability <- derived_marker_free_dataset$total_adult_count / derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count
marker_free_viability_dataset <- cbind(derived_marker_free_dataset, derived_marker_free_female_viability, derived_marker_free_male_viability, derived_marker_free_total_viability)
write.csv(marker_free_viability_dataset, "../../Summaries and Outputs/2/F_marker_free_viability_dataset.csv", row.names = TRUE)

#$#
marker_free_viability_long <- data.frame(
  Cross = rep(marker_free_viability_dataset$cross_id, times = 2),
  Brood = rep(marker_free_viability_dataset$brood_id, times = 2),
  Treatment = rep(marker_free_viability_dataset$treatment, times = 2),
  Data_collector = rep(marker_free_viability_dataset$data_collector, times = 2),
  Genotype = rep(marker_free_viability_dataset$genotype, times = 2),
  Sex = factor(
    rep(c("Female", "Male"), each = nrow(marker_free_viability_dataset))
  ),
  Viability = rep(marker_free_viability_dataset$derived_marker_free_total_viability, times = 2))

marker_free_viability_long$Cross <- factor(marker_free_viability_long$Cross)
marker_free_viability_long$Brood <- factor(marker_free_viability_long$Brood)
marker_free_viability_long$Treatment <- factor(marker_free_viability_long$Treatment)
marker_free_viability_long$Genotype <- factor(marker_free_viability_long$Genotype)
marker_free_viability_long$Sex <- factor(marker_free_viability_long$Sex)
marker_free_viability_long$Data_collector <- factor(marker_free_viability_long$Data_collector)
write.csv(marker_free_viability_long,"../../Summaries and Outputs/2/G_marker_free_viability_long.csv",row.names = FALSE)
#anova
marked_free_viability_anova <- aov(Viability ~ Treatment  * Sex * Genotype *Brood * Data_collector,data = marker_free_viability_long)
summary(marked_free_viability_anova)

# emmeans for Treatment × Genotype (ignoring Sex if not significant)
marker_free_emm <- emmeans(marked_free_viability_anova, ~ Treatment * Genotype*Brood )

# Get pairwise comparisons if needed
marker_free_pairs <- pairs(marker_free_emm)
marker_free_pairs

# Compute compact letter display (CLD)
marker_free_cld <- cld(marker_free_emm, Letters = letters, adjust = "tukey")
marker_free_cld_df <- as.data.frame(marker_free_cld)
write.csv(marker_free_cld_df,
          "../../Summaries and Outputs/2/I_marker_free_viability_emmeans_cld.csv",
          row.names = FALSE)

#2
viability_means_tgs <- marker_free_viability_long %>%
group_by(Treatment,Brood) %>%
  summarise(
    Mean_Viability = mean(Viability, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(viability_means_tgs,"../../Summaries and Outputs/2/H_marker_free_viability_sumarry.csv",row.names = FALSE)

#Marker present------
# Read in csv files of raw count data from the multiply marked cross
CH_multiply_marked_egg_count <- read.csv("1_CH_raw_multiply_marked_egg_counts.csv", header=TRUE)
LH_multiply_marked_egg_count <- read.csv("2_LH_raw_multiply_marked_egg_counts.csv", header=TRUE)
MS_multiply_marked_adult_phenotypic_class_counts <- read.csv("3_MS_raw_multiply_marked_adult_phenotypic_class_2experiment.csv", header=TRUE)


# Compile the raw data from different observers and life history stages
CH_raw_multiply_marked_total_egg_count <- CH_multiply_marked_egg_count$total_egg_count
LH_raw_multiply_marked_total_egg_count <- LH_multiply_marked_egg_count$total_egg_count
compiled_raw_multiply_marked_dataset <- cbind(MS_multiply_marked_adult_phenotypic_class_counts, CH_raw_multiply_marked_total_egg_count, LH_raw_multiply_marked_total_egg_count)

# Split treatment_id into treatment and genotype for multiply-marked dataset
treatment_genotype_mm <- strsplit(as.character(compiled_raw_multiply_marked_dataset$treatment_id), "_")
compiled_raw_multiply_marked_dataset$treatment <- sapply(treatment_genotype_mm, function(x) x[1])
compiled_raw_multiply_marked_dataset$genotype <- sapply(treatment_genotype_mm, function(x) x[2])

write.csv(compiled_raw_multiply_marked_dataset, "../../Summaries and Outputs/2/4_compiled_raw_multiply_marked_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting
CH_transformed_multiply_marked_total_egg_count <- (CH_raw_multiply_marked_total_egg_count - 6.479979) / 0.9224861    
LH_transformed_multiply_marked_total_egg_count <- (LH_raw_multiply_marked_total_egg_count - 2.5569307 ) / 0.8528595 
transformed_mean_multiply_marked_total_egg_count <- round((CH_transformed_multiply_marked_total_egg_count + LH_transformed_multiply_marked_total_egg_count) / 2)
derived_multiply_marked_dataset <- cbind(compiled_raw_multiply_marked_dataset, transformed_mean_multiply_marked_total_egg_count)
derived_multiply_marked_dataset <-derived_multiply_marked_dataset %>% filter(!is.na(total_adult_count))
write.csv(derived_multiply_marked_dataset, "../../Summaries and Outputs/2/5_derived_multiply_marked_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_multiply_marked_female_viability <- derived_multiply_marked_dataset$female_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_male_viability <- derived_multiply_marked_dataset$male_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_total_viability <- derived_multiply_marked_dataset$total_adult_count / derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count
multiply_marked_viability_dataset <- cbind(derived_multiply_marked_dataset, derived_multiply_marked_female_viability, derived_multiply_marked_male_viability, derived_multiply_marked_total_viability)
write.csv(multiply_marked_viability_dataset, "../../Summaries and Outputs/2/6_multiply_marked_viability_dataset.csv", row.names = TRUE)

# Compute the overall mean viabilities for reporting in main text 
overall_multiply_marked_viability <- mean(derived_multiply_marked_total_viability, na.rm = TRUE)
female_multiply_marked_viability <- mean(derived_multiply_marked_female_viability, na.rm = TRUE)
male_multiply_marked_viability <- mean(derived_multiply_marked_male_viability, na.rm = TRUE)

# Print overall viability statistics
cat("Overall multiply marked viability:", overall_multiply_marked_viability, "\n")
cat("Female multiply marked viability:", female_multiply_marked_viability, "\n")
cat("Male multiply marked viability:", male_multiply_marked_viability, "\n")

# Generate a long-format data frame for sex-specific viability effects in the multiply marked crosses
multiply_marked_viability_long <- data.frame(
  Cross = rep(multiply_marked_viability_dataset$cross_id, times = 2),
  Brood = rep(multiply_marked_viability_dataset$brood_id, times = 2),
  Treatment = rep(multiply_marked_viability_dataset$treatment, times = 2),
  Genotype = rep(multiply_marked_viability_dataset$genotype, times = 2),
  Data_Collector = rep(multiply_marked_viability_dataset$data_collector, times = 2),
  Sex = factor(rep(c("Female", "Male"), each = nrow(multiply_marked_viability_dataset))),
  Viability = rep(multiply_marked_viability_dataset$derived_multiply_marked_total_viability, times = 2))

# Convert variables to appropriate types
multiply_marked_viability_long$Cross <- factor(multiply_marked_viability_long$Cross)
multiply_marked_viability_long$Brood <- factor(multiply_marked_viability_long$Brood)
multiply_marked_viability_long$Treatment <- factor(multiply_marked_viability_long$Treatment)
multiply_marked_viability_long$Genotype <- factor(multiply_marked_viability_long$Genotype)
multiply_marked_viability_long$Data_Collector <- factor(multiply_marked_viability_long$Data_Collector)
multiply_marked_viability_long$Sex <- factor(multiply_marked_viability_long$Sex)

# Save the long-format dataset
write.csv(multiply_marked_viability_long, 
          "../../Summaries and Outputs/2/7_multiply_marked_viability_long.csv", 
          row.names = FALSE)

# Run ANOVA on viability
multiply_marked_viability_anova <- aov(Viability ~ Treatment * Genotype * Sex * Brood * Data_Collector, 
                                       data = multiply_marked_viability_long)
summary(multiply_marked_viability_anova)

# Compute emmeans for Treatment × Genotype × Data_Collector
multiply_marked_emm <- emmeans(multiply_marked_viability_anova, ~ Treatment * Genotype)

# Pairwise comparisons
multiply_marked_pairs <- pairs(multiply_marked_emm)
multiply_marked_pairs

# Compute CLD
multiply_marked_cld <- cld(multiply_marked_emm, Letters = letters, adjust = "tukey")
multiply_marked_cld_df <- as.data.frame(multiply_marked_cld)
write.csv(multiply_marked_cld_df,
          "../../Summaries and Outputs/2/J_multiply_marked_viability_emmeans_cld.csv",
          row.names = FALSE)
          
# Compute Treatment × Genotype means (matching the marker-free approach)
 multiply_marked_viability_means_tg <- multiply_marked_viability_long %>%
   group_by(Treatment, Genotype) %>%
   summarise(
     Mean_Viability = mean(Viability, na.rm = TRUE),
     SD_Viability = sd(Viability, na.rm = TRUE),
     N = n(),
     .groups = "drop"
   )
 
# Save the summary statistics
 write.csv(multiply_marked_viability_means_tg, 
           "../../Summaries and Outputs/2/8_multiply_marked_viability_summary.csv", 
           row.names = FALSE)
 

# Input for Counting Model ------------------------------------------------
#individual vials-------------
# Create input data file for fitting Cx(Co)M model to flies pooled by individual vials
# Calculate the number of lethal zygotes by subtracting total adults from mean transformed egg count
# Add marker-free viabilities as classes 130, 131, 132, rename columns to conform to standard format
# 1) Compute lethal zygotes PER VIAL
derived_multiply_marked_dataset <- derived_multiply_marked_dataset %>%
  mutate(
    count_lethal_zygotes =
      transformed_mean_multiply_marked_total_egg_count - total_adult_count
  )

# 2) Split treatment_id into Treatment and Genotype
treatment_genotype <- strsplit(
  as.character(derived_multiply_marked_dataset$treatment_id), "_"
)

derived_multiply_marked_dataset$Treatment <- sapply(treatment_genotype, `[`, 1)
derived_multiply_marked_dataset$Genotype  <- sapply(treatment_genotype, `[`, 2)
derived_multiply_marked_dataset$Brood <- derived_multiply_marked_dataset$brood_id

#2
derived_multiply_marked_dataset <- derived_multiply_marked_dataset %>%
  left_join(viability_means_tgs, by = c("Treatment","Brood"))
# 4) Prune columns
pruned_derived_multiply_marked_dataset <- derived_multiply_marked_dataset[  , -c(1, 2, 3, 5, 6, 135, 136, 137, 138, 139, 140, 141, 142, 144, 145,146) ]

# 5) Rename for multilocus format
female_viability_vector <- pruned_derived_multiply_marked_dataset$Mean_Viability
male_viability_vector <- pruned_derived_multiply_marked_dataset$Mean_Viability
multi_locus_individual_vials_dataset <- cbind(pruned_derived_multiply_marked_dataset,female_viability_vector,male_viability_vector)
names(multi_locus_individual_vials_dataset) <- c("ID", paste0("class_", 1:132))

# 6) Remove rows with missing final class
multi_locus_individual_vials_dataset <-multi_locus_individual_vials_dataset %>% filter(!is.na(class_132))

# 7) Write output
write.csv(multi_locus_individual_vials_dataset,"../../Summaries and Outputs/2/13_multi_locus_individual_vials_dataset.csv", row.names = FALSE)

