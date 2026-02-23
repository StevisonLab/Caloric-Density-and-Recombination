#setwd("Raw Data/") # Please use Raw Data as working directory
library(dplyr)
library(dbplyr)
library(emmeans)
library(multcomp) 
library(multcompView)  # for cld() support if needed

#Marker free----
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

write.csv(compiled_raw_marker_free_dataset, "../../Summaries and Outputs/3/D_compiled_raw_marker_free_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting the average and round these two counts to the nearest whole number
CH_transformed_marker_free_total_egg_count <- (CH_raw_marker_free_total_egg_count - 6.479979) / 0.9224861
LM_transformed_marker_free_total_egg_count <- (LM_raw_marker_free_total_egg_count - 4.866269) / 0.8351445      
transformed_mean_marker_free_total_egg_count <- round((CH_transformed_marker_free_total_egg_count + LM_transformed_marker_free_total_egg_count) / 2)

derived_marker_free_dataset <- cbind(compiled_raw_marker_free_dataset, transformed_mean_marker_free_total_egg_count)
write.csv(derived_marker_free_dataset, "../../Summaries and Outputs/3/E_derived_marker_free_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_marker_free_female_viability <- derived_marker_free_dataset$female_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_male_viability <- derived_marker_free_dataset$male_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_total_viability <- derived_marker_free_dataset$total_adult_count / derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count
marker_free_viability_dataset <- cbind(derived_marker_free_dataset, derived_marker_free_female_viability, derived_marker_free_male_viability, derived_marker_free_total_viability)
write.csv(marker_free_viability_dataset, "../../Summaries and Outputs/3/F_marker_free_viability_dataset.csv", row.names = TRUE)

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
write.csv(marker_free_viability_long,"../../Summaries and Outputs/3/G_marker_free_viability_long.csv",row.names = FALSE)
#anova
marked_free_viability_anova <- aov(Viability ~ Treatment  * Sex * Genotype *Brood * Data_collector,data = marker_free_viability_long)
anova_table <- summary(marked_free_viability_anova)[[1]]
write.csv(anova_table,"../../Summaries and Outputs/3/H_marker_free_viability_anova.csv",row.names = FALSE)

# emmeans for Treatment × Genotype (ignoring Sex if not significant)
marker_free_emm <- emmeans(marked_free_viability_anova, ~ Treatment * Genotype)

# Get pairwise comparisons if needed
marker_free_pairs <- pairs(marker_free_emm)
marker_free_pairs

# Compute compact letter display (CLD)
marker_free_cld <- cld(marker_free_emm, Letters = letters, adjust = "tukey")
marker_free_cld_df <- as.data.frame(marker_free_cld)
write.csv(marker_free_cld_df,
          "../../Summaries and Outputs/3/I_marker_free_viability_emmeans_cld.csv",
          row.names = FALSE)

#Compute Treatment × Genotype × means => because sex was not significant we are not using #Sex
#3
viability_means_tgs <- marker_free_viability_long %>%
  group_by(Treatment, Genotype,Brood) %>%
  summarise(
    Mean_Viability = mean(Viability, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(viability_means_tgs,"../../Summaries and Outputs/3/H_marker_free_viability_sumarry.csv",row.names = FALSE)

#Marker present------
# Read in csv files of raw count data from the multiply marked cross
CH_multiply_marked_egg_count <- read.csv("1_CH_raw_multiply_marked_egg_counts.csv", header=TRUE)
LM_multiply_marked_egg_count <- read.csv("2_LM_raw_multiply_marked_egg_counts.csv", header=TRUE)
MS_multiply_marked_adult_phenotypic_class_counts <- read.csv("3_MS_raw_multiply_marked_adult_phenotypic_class_counts.csv", header=TRUE)

# Compile the raw data from different observers and life history stages
CH_raw_multiply_marked_total_egg_count <- CH_multiply_marked_egg_count$total_egg_count
LM_raw_multiply_marked_total_egg_count <- LM_multiply_marked_egg_count$total_egg_count
compiled_raw_multiply_marked_dataset <- cbind(MS_multiply_marked_adult_phenotypic_class_counts, CH_raw_multiply_marked_total_egg_count, LM_raw_multiply_marked_total_egg_count)

# Split treatment_id into treatment and genotype for multiply-marked dataset
treatment_genotype_mm <- strsplit(as.character(compiled_raw_multiply_marked_dataset$treatment_id), "_")
compiled_raw_multiply_marked_dataset$treatment <- sapply(treatment_genotype_mm, function(x) x[1])
compiled_raw_multiply_marked_dataset$genotype <- sapply(treatment_genotype_mm, function(x) x[2])

write.csv(compiled_raw_multiply_marked_dataset, "../../Summaries and Outputs/3/4_compiled_raw_multiply_marked_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting
CH_transformed_multiply_marked_total_egg_count <- (CH_raw_multiply_marked_total_egg_count - 6.479979) / 0.9224861    
LM_transformed_multiply_marked_total_egg_count <- (LM_raw_multiply_marked_total_egg_count - 4.866269) / 0.8351445
transformed_mean_multiply_marked_total_egg_count <- round((CH_transformed_multiply_marked_total_egg_count + LM_transformed_multiply_marked_total_egg_count) / 2)
derived_multiply_marked_dataset <- cbind(compiled_raw_multiply_marked_dataset, transformed_mean_multiply_marked_total_egg_count)
derived_multiply_marked_dataset <-derived_multiply_marked_dataset %>% filter(!is.na(total_adult_count))
write.csv(derived_multiply_marked_dataset, "../../Summaries and Outputs/3/5_derived_multiply_marked_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_multiply_marked_female_viability <- derived_multiply_marked_dataset$female_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_male_viability <- derived_multiply_marked_dataset$male_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_total_viability <- derived_multiply_marked_dataset$total_adult_count / derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count
multiply_marked_viability_dataset <- cbind(derived_multiply_marked_dataset, derived_multiply_marked_female_viability, derived_multiply_marked_male_viability, derived_multiply_marked_total_viability)
write.csv(multiply_marked_viability_dataset, "../../Summaries and Outputs/3/6_multiply_marked_viability_dataset.csv", row.names = TRUE)

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
          "../../Summaries and Outputs/3/7_multiply_marked_viability_long.csv", 
          row.names = FALSE)

# Run ANOVA on viability
multiply_marked_viability_anova <- aov(Viability ~ Treatment * Genotype * Sex * Brood * Data_Collector, 
                                       data = multiply_marked_viability_long)
anova_table2 <- summary(multiply_marked_viability_anova)[[1]]
write.csv(anova_table2,"../../Summaries and Outputs/3/8_multiply_marked_viability_anova.csv",row.names = TRUE)

# Compute emmeans for Treatment × Genotype × Data_Collector
multiply_marked_emm <- emmeans(multiply_marked_viability_anova, ~ Treatment * Genotype)

# Pairwise comparisons
multiply_marked_pairs <- pairs(multiply_marked_emm)
multiply_marked_pairs

# Compute CLD
multiply_marked_cld <- cld(multiply_marked_emm, Letters = letters, adjust = "tukey")
multiply_marked_cld_df <- as.data.frame(multiply_marked_cld)
write.csv(multiply_marked_cld_df,
          "../../Summaries and Outputs/3/9_multiply_marked_viability_emmeans_cld.csv",
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
           "../../Summaries and Outputs/3/10_multiply_marked_viability_summary.csv", 
           row.names = FALSE)
 
 #single locus -----------------
 #Create 6 single locus datasets as vectors
 phenotypic_class_dataset <- multiply_marked_viability_dataset[,c(7:134)]
 phenotypic_class_sums_vector <- colSums(phenotypic_class_dataset)
 scute_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                     phenotypic_class_sums_vector[3],
                                     phenotypic_class_sums_vector[5],
                                     phenotypic_class_sums_vector[7],
                                     phenotypic_class_sums_vector[9],
                                     phenotypic_class_sums_vector[11],
                                     phenotypic_class_sums_vector[13],
                                     phenotypic_class_sums_vector[15],
                                     phenotypic_class_sums_vector[17],
                                     phenotypic_class_sums_vector[19],
                                     phenotypic_class_sums_vector[21],
                                     phenotypic_class_sums_vector[23],
                                     phenotypic_class_sums_vector[25],
                                     phenotypic_class_sums_vector[27],
                                     phenotypic_class_sums_vector[29],
                                     phenotypic_class_sums_vector[31],
                                     phenotypic_class_sums_vector[33],
                                     phenotypic_class_sums_vector[35],
                                     phenotypic_class_sums_vector[37],
                                     phenotypic_class_sums_vector[39],
                                     phenotypic_class_sums_vector[41],
                                     phenotypic_class_sums_vector[43],
                                     phenotypic_class_sums_vector[45],
                                     phenotypic_class_sums_vector[47],
                                     phenotypic_class_sums_vector[49],
                                     phenotypic_class_sums_vector[51],
                                     phenotypic_class_sums_vector[53],
                                     phenotypic_class_sums_vector[55],
                                     phenotypic_class_sums_vector[57],
                                     phenotypic_class_sums_vector[59],
                                     phenotypic_class_sums_vector[61],
                                     phenotypic_class_sums_vector[63]),
                                 sum(phenotypic_class_sums_vector[65],
                                     phenotypic_class_sums_vector[67],
                                     phenotypic_class_sums_vector[69],
                                     phenotypic_class_sums_vector[71],
                                     phenotypic_class_sums_vector[73],
                                     phenotypic_class_sums_vector[75],
                                     phenotypic_class_sums_vector[77],
                                     phenotypic_class_sums_vector[79],
                                     phenotypic_class_sums_vector[81],
                                     phenotypic_class_sums_vector[83],
                                     phenotypic_class_sums_vector[85],
                                     phenotypic_class_sums_vector[87],
                                     phenotypic_class_sums_vector[89],
                                     phenotypic_class_sums_vector[91],
                                     phenotypic_class_sums_vector[93],
                                     phenotypic_class_sums_vector[95],
                                     phenotypic_class_sums_vector[97],
                                     phenotypic_class_sums_vector[99],
                                     phenotypic_class_sums_vector[101],
                                     phenotypic_class_sums_vector[103],
                                     phenotypic_class_sums_vector[105],
                                     phenotypic_class_sums_vector[107],
                                     phenotypic_class_sums_vector[109],
                                     phenotypic_class_sums_vector[111],
                                     phenotypic_class_sums_vector[113],
                                     phenotypic_class_sums_vector[115],
                                     phenotypic_class_sums_vector[117],
                                     phenotypic_class_sums_vector[119],
                                     phenotypic_class_sums_vector[121],
                                     phenotypic_class_sums_vector[123],
                                     phenotypic_class_sums_vector[125],
                                     phenotypic_class_sums_vector[127]),
                                 sum(phenotypic_class_sums_vector[2],
                                     phenotypic_class_sums_vector[4],
                                     phenotypic_class_sums_vector[6],
                                     phenotypic_class_sums_vector[8],
                                     phenotypic_class_sums_vector[10],
                                     phenotypic_class_sums_vector[12],
                                     phenotypic_class_sums_vector[14],
                                     phenotypic_class_sums_vector[16],
                                     phenotypic_class_sums_vector[18],
                                     phenotypic_class_sums_vector[20],
                                     phenotypic_class_sums_vector[22],
                                     phenotypic_class_sums_vector[24],
                                     phenotypic_class_sums_vector[26],
                                     phenotypic_class_sums_vector[28],
                                     phenotypic_class_sums_vector[30],
                                     phenotypic_class_sums_vector[32],
                                     phenotypic_class_sums_vector[34],
                                     phenotypic_class_sums_vector[36],
                                     phenotypic_class_sums_vector[38],
                                     phenotypic_class_sums_vector[40],
                                     phenotypic_class_sums_vector[42],
                                     phenotypic_class_sums_vector[44],
                                     phenotypic_class_sums_vector[46],
                                     phenotypic_class_sums_vector[48],
                                     phenotypic_class_sums_vector[50],
                                     phenotypic_class_sums_vector[52],
                                     phenotypic_class_sums_vector[54],
                                     phenotypic_class_sums_vector[56],
                                     phenotypic_class_sums_vector[58],
                                     phenotypic_class_sums_vector[60],
                                     phenotypic_class_sums_vector[62],
                                     phenotypic_class_sums_vector[64]),
                                 sum(phenotypic_class_sums_vector[66],
                                     phenotypic_class_sums_vector[68],
                                     phenotypic_class_sums_vector[70],
                                     phenotypic_class_sums_vector[72],
                                     phenotypic_class_sums_vector[74],
                                     phenotypic_class_sums_vector[76],
                                     phenotypic_class_sums_vector[78],
                                     phenotypic_class_sums_vector[80],
                                     phenotypic_class_sums_vector[82],
                                     phenotypic_class_sums_vector[84],
                                     phenotypic_class_sums_vector[86],
                                     phenotypic_class_sums_vector[88],
                                     phenotypic_class_sums_vector[90],
                                     phenotypic_class_sums_vector[92],
                                     phenotypic_class_sums_vector[94],
                                     phenotypic_class_sums_vector[96],
                                     phenotypic_class_sums_vector[98],
                                     phenotypic_class_sums_vector[100],
                                     phenotypic_class_sums_vector[102],
                                     phenotypic_class_sums_vector[104],
                                     phenotypic_class_sums_vector[106],
                                     phenotypic_class_sums_vector[108],
                                     phenotypic_class_sums_vector[110],
                                     phenotypic_class_sums_vector[112],
                                     phenotypic_class_sums_vector[114],
                                     phenotypic_class_sums_vector[116],
                                     phenotypic_class_sums_vector[118],
                                     phenotypic_class_sums_vector[120],
                                     phenotypic_class_sums_vector[122],
                                     phenotypic_class_sums_vector[124],
                                     phenotypic_class_sums_vector[126],
                                     phenotypic_class_sums_vector[128]),
                                 4795)
 
 # Write csv for the single locus dataset
 names(scute_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(scute_single_locus_dataset,"../../Summaries and Outputs/3/7_scute_single_locus_dataset.csv",row.names = TRUE)
 
 crossveinless_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                             phenotypic_class_sums_vector[3],
                                             phenotypic_class_sums_vector[5],
                                             phenotypic_class_sums_vector[7],
                                             phenotypic_class_sums_vector[9],
                                             phenotypic_class_sums_vector[12],
                                             phenotypic_class_sums_vector[13],
                                             phenotypic_class_sums_vector[15],
                                             phenotypic_class_sums_vector[17],
                                             phenotypic_class_sums_vector[20],
                                             phenotypic_class_sums_vector[21],
                                             phenotypic_class_sums_vector[23],
                                             phenotypic_class_sums_vector[26],
                                             phenotypic_class_sums_vector[27],
                                             phenotypic_class_sums_vector[30],
                                             phenotypic_class_sums_vector[32],
                                             phenotypic_class_sums_vector[33],
                                             phenotypic_class_sums_vector[35],
                                             phenotypic_class_sums_vector[38],
                                             phenotypic_class_sums_vector[39],
                                             phenotypic_class_sums_vector[42],
                                             phenotypic_class_sums_vector[44],
                                             phenotypic_class_sums_vector[45],
                                             phenotypic_class_sums_vector[48],
                                             phenotypic_class_sums_vector[50],
                                             phenotypic_class_sums_vector[52],
                                             phenotypic_class_sums_vector[53],
                                             phenotypic_class_sums_vector[56],
                                             phenotypic_class_sums_vector[58],
                                             phenotypic_class_sums_vector[60],
                                             phenotypic_class_sums_vector[62],
                                             phenotypic_class_sums_vector[64]),
                                         sum(phenotypic_class_sums_vector[65],
                                             phenotypic_class_sums_vector[67],
                                             phenotypic_class_sums_vector[69],
                                             phenotypic_class_sums_vector[71],
                                             phenotypic_class_sums_vector[73],
                                             phenotypic_class_sums_vector[76],
                                             phenotypic_class_sums_vector[77],
                                             phenotypic_class_sums_vector[79],
                                             phenotypic_class_sums_vector[81],
                                             phenotypic_class_sums_vector[84],
                                             phenotypic_class_sums_vector[85],
                                             phenotypic_class_sums_vector[87],
                                             phenotypic_class_sums_vector[90],
                                             phenotypic_class_sums_vector[91],
                                             phenotypic_class_sums_vector[94],
                                             phenotypic_class_sums_vector[96],
                                             phenotypic_class_sums_vector[97],
                                             phenotypic_class_sums_vector[99],
                                             phenotypic_class_sums_vector[102],
                                             phenotypic_class_sums_vector[103],
                                             phenotypic_class_sums_vector[106],
                                             phenotypic_class_sums_vector[108],
                                             phenotypic_class_sums_vector[109],
                                             phenotypic_class_sums_vector[112],
                                             phenotypic_class_sums_vector[114],
                                             phenotypic_class_sums_vector[116],
                                             phenotypic_class_sums_vector[117],
                                             phenotypic_class_sums_vector[120],
                                             phenotypic_class_sums_vector[122],
                                             phenotypic_class_sums_vector[124],
                                             phenotypic_class_sums_vector[126],
                                             phenotypic_class_sums_vector[128]),
                                         sum(phenotypic_class_sums_vector[2],
                                             phenotypic_class_sums_vector[4],
                                             phenotypic_class_sums_vector[6],
                                             phenotypic_class_sums_vector[8],
                                             phenotypic_class_sums_vector[10],
                                             phenotypic_class_sums_vector[11],
                                             phenotypic_class_sums_vector[14],
                                             phenotypic_class_sums_vector[16],
                                             phenotypic_class_sums_vector[18],
                                             phenotypic_class_sums_vector[19],
                                             phenotypic_class_sums_vector[22],
                                             phenotypic_class_sums_vector[24],
                                             phenotypic_class_sums_vector[25],
                                             phenotypic_class_sums_vector[28],
                                             phenotypic_class_sums_vector[29],
                                             phenotypic_class_sums_vector[31],
                                             phenotypic_class_sums_vector[34],
                                             phenotypic_class_sums_vector[36],
                                             phenotypic_class_sums_vector[37],
                                             phenotypic_class_sums_vector[40],
                                             phenotypic_class_sums_vector[41],
                                             phenotypic_class_sums_vector[43],
                                             phenotypic_class_sums_vector[46],
                                             phenotypic_class_sums_vector[47],
                                             phenotypic_class_sums_vector[49],
                                             phenotypic_class_sums_vector[51],
                                             phenotypic_class_sums_vector[54],
                                             phenotypic_class_sums_vector[55],
                                             phenotypic_class_sums_vector[57],
                                             phenotypic_class_sums_vector[59],
                                             phenotypic_class_sums_vector[61],
                                             phenotypic_class_sums_vector[63]),
                                         sum(phenotypic_class_sums_vector[66],
                                             phenotypic_class_sums_vector[68],
                                             phenotypic_class_sums_vector[70],
                                             phenotypic_class_sums_vector[72],
                                             phenotypic_class_sums_vector[74],
                                             phenotypic_class_sums_vector[75],
                                             phenotypic_class_sums_vector[78],
                                             phenotypic_class_sums_vector[80],
                                             phenotypic_class_sums_vector[82],
                                             phenotypic_class_sums_vector[83],
                                             phenotypic_class_sums_vector[86],
                                             phenotypic_class_sums_vector[88],
                                             phenotypic_class_sums_vector[89],
                                             phenotypic_class_sums_vector[92],
                                             phenotypic_class_sums_vector[93],
                                             phenotypic_class_sums_vector[95],
                                             phenotypic_class_sums_vector[98],
                                             phenotypic_class_sums_vector[100],
                                             phenotypic_class_sums_vector[101],
                                             phenotypic_class_sums_vector[104],
                                             phenotypic_class_sums_vector[105],
                                             phenotypic_class_sums_vector[107],
                                             phenotypic_class_sums_vector[110],
                                             phenotypic_class_sums_vector[111],
                                             phenotypic_class_sums_vector[113],
                                             phenotypic_class_sums_vector[115],
                                             phenotypic_class_sums_vector[118],
                                             phenotypic_class_sums_vector[119],
                                             phenotypic_class_sums_vector[121],
                                             phenotypic_class_sums_vector[123],
                                             phenotypic_class_sums_vector[125],
                                             phenotypic_class_sums_vector[127]),
                                         4795)
 names(crossveinless_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(crossveinless_single_locus_dataset,"../../Summaries and Outputs/3/8_crossveinless_single_locus_dataset.csv",row.names = TRUE)
 
 vermilion_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                         phenotypic_class_sums_vector[3],
                                         phenotypic_class_sums_vector[5],
                                         phenotypic_class_sums_vector[7],
                                         phenotypic_class_sums_vector[10],
                                         phenotypic_class_sums_vector[12],
                                         phenotypic_class_sums_vector[13],
                                         phenotypic_class_sums_vector[15],
                                         phenotypic_class_sums_vector[18],
                                         phenotypic_class_sums_vector[20],
                                         phenotypic_class_sums_vector[21],
                                         phenotypic_class_sums_vector[24],
                                         phenotypic_class_sums_vector[26],
                                         phenotypic_class_sums_vector[28],
                                         phenotypic_class_sums_vector[30],
                                         phenotypic_class_sums_vector[31],
                                         phenotypic_class_sums_vector[33],
                                         phenotypic_class_sums_vector[36],
                                         phenotypic_class_sums_vector[38],
                                         phenotypic_class_sums_vector[40],
                                         phenotypic_class_sums_vector[42],
                                         phenotypic_class_sums_vector[43],
                                         phenotypic_class_sums_vector[46],
                                         phenotypic_class_sums_vector[48],
                                         phenotypic_class_sums_vector[49],
                                         phenotypic_class_sums_vector[51],
                                         phenotypic_class_sums_vector[54],
                                         phenotypic_class_sums_vector[56],
                                         phenotypic_class_sums_vector[57],
                                         phenotypic_class_sums_vector[59],
                                         phenotypic_class_sums_vector[61],
                                         phenotypic_class_sums_vector[63]),
                                     sum(phenotypic_class_sums_vector[65],
                                         phenotypic_class_sums_vector[67],
                                         phenotypic_class_sums_vector[69],
                                         phenotypic_class_sums_vector[71],
                                         phenotypic_class_sums_vector[74],
                                         phenotypic_class_sums_vector[76],
                                         phenotypic_class_sums_vector[77],
                                         phenotypic_class_sums_vector[79],
                                         phenotypic_class_sums_vector[82],
                                         phenotypic_class_sums_vector[84],
                                         phenotypic_class_sums_vector[85],
                                         phenotypic_class_sums_vector[88],
                                         phenotypic_class_sums_vector[90],
                                         phenotypic_class_sums_vector[92],
                                         phenotypic_class_sums_vector[94],
                                         phenotypic_class_sums_vector[95],
                                         phenotypic_class_sums_vector[97],
                                         phenotypic_class_sums_vector[100],
                                         phenotypic_class_sums_vector[102],
                                         phenotypic_class_sums_vector[104],
                                         phenotypic_class_sums_vector[106],
                                         phenotypic_class_sums_vector[107],
                                         phenotypic_class_sums_vector[110],
                                         phenotypic_class_sums_vector[112],
                                         phenotypic_class_sums_vector[113],
                                         phenotypic_class_sums_vector[115],
                                         phenotypic_class_sums_vector[118],
                                         phenotypic_class_sums_vector[120],
                                         phenotypic_class_sums_vector[121],
                                         phenotypic_class_sums_vector[123],
                                         phenotypic_class_sums_vector[125],
                                         phenotypic_class_sums_vector[127]), 
                                     sum(phenotypic_class_sums_vector[2],
                                         phenotypic_class_sums_vector[4],
                                         phenotypic_class_sums_vector[6],
                                         phenotypic_class_sums_vector[8],
                                         phenotypic_class_sums_vector[9],
                                         phenotypic_class_sums_vector[11],
                                         phenotypic_class_sums_vector[14],
                                         phenotypic_class_sums_vector[16],
                                         phenotypic_class_sums_vector[17],
                                         phenotypic_class_sums_vector[19],
                                         phenotypic_class_sums_vector[22],
                                         phenotypic_class_sums_vector[23],
                                         phenotypic_class_sums_vector[25],
                                         phenotypic_class_sums_vector[27],
                                         phenotypic_class_sums_vector[29],
                                         phenotypic_class_sums_vector[32],
                                         phenotypic_class_sums_vector[34],
                                         phenotypic_class_sums_vector[35],
                                         phenotypic_class_sums_vector[37],
                                         phenotypic_class_sums_vector[39],
                                         phenotypic_class_sums_vector[41],
                                         phenotypic_class_sums_vector[44],
                                         phenotypic_class_sums_vector[45],
                                         phenotypic_class_sums_vector[47],
                                         phenotypic_class_sums_vector[50],
                                         phenotypic_class_sums_vector[52],
                                         phenotypic_class_sums_vector[53],
                                         phenotypic_class_sums_vector[55],
                                         phenotypic_class_sums_vector[58],
                                         phenotypic_class_sums_vector[60],
                                         phenotypic_class_sums_vector[62],
                                         phenotypic_class_sums_vector[64]),
                                     sum(phenotypic_class_sums_vector[66],
                                         phenotypic_class_sums_vector[68],
                                         phenotypic_class_sums_vector[70],
                                         phenotypic_class_sums_vector[72],
                                         phenotypic_class_sums_vector[73],
                                         phenotypic_class_sums_vector[75],
                                         phenotypic_class_sums_vector[78],
                                         phenotypic_class_sums_vector[80],
                                         phenotypic_class_sums_vector[81],
                                         phenotypic_class_sums_vector[83],
                                         phenotypic_class_sums_vector[86],
                                         phenotypic_class_sums_vector[87],
                                         phenotypic_class_sums_vector[89],
                                         phenotypic_class_sums_vector[91],
                                         phenotypic_class_sums_vector[93],
                                         phenotypic_class_sums_vector[96],
                                         phenotypic_class_sums_vector[98],
                                         phenotypic_class_sums_vector[99],
                                         phenotypic_class_sums_vector[101],
                                         phenotypic_class_sums_vector[103],
                                         phenotypic_class_sums_vector[105],
                                         phenotypic_class_sums_vector[108],
                                         phenotypic_class_sums_vector[109],
                                         phenotypic_class_sums_vector[111],
                                         phenotypic_class_sums_vector[114],
                                         phenotypic_class_sums_vector[116],
                                         phenotypic_class_sums_vector[117],
                                         phenotypic_class_sums_vector[119],
                                         phenotypic_class_sums_vector[122],
                                         phenotypic_class_sums_vector[124],
                                         phenotypic_class_sums_vector[126],
                                         phenotypic_class_sums_vector[128]),
                                     4795)
 names(vermilion_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(vermilion_single_locus_dataset,"../../Summaries and Outputs/3/9_vermilion_single_locus_dataset.csv",row.names = TRUE)
 
 forked_single_locus_dataset <-c(sum(phenotypic_class_sums_vector[1],
                                     phenotypic_class_sums_vector[3],
                                     phenotypic_class_sums_vector[5],
                                     phenotypic_class_sums_vector[8],
                                     phenotypic_class_sums_vector[10],
                                     phenotypic_class_sums_vector[12],
                                     phenotypic_class_sums_vector[13],
                                     phenotypic_class_sums_vector[16],
                                     phenotypic_class_sums_vector[18],
                                     phenotypic_class_sums_vector[20],
                                     phenotypic_class_sums_vector[22],
                                     phenotypic_class_sums_vector[24],
                                     phenotypic_class_sums_vector[26],
                                     phenotypic_class_sums_vector[27],
                                     phenotypic_class_sums_vector[29],
                                     phenotypic_class_sums_vector[31],
                                     phenotypic_class_sums_vector[34],
                                     phenotypic_class_sums_vector[36],
                                     phenotypic_class_sums_vector[38],
                                     phenotypic_class_sums_vector[39],
                                     phenotypic_class_sums_vector[41],
                                     phenotypic_class_sums_vector[43],
                                     phenotypic_class_sums_vector[45],
                                     phenotypic_class_sums_vector[47],
                                     phenotypic_class_sums_vector[49],
                                     phenotypic_class_sums_vector[52],
                                     phenotypic_class_sums_vector[53],
                                     phenotypic_class_sums_vector[55],
                                     phenotypic_class_sums_vector[57],
                                     phenotypic_class_sums_vector[60],
                                     phenotypic_class_sums_vector[62],
                                     phenotypic_class_sums_vector[64]),
                                 sum(phenotypic_class_sums_vector[65],
                                     phenotypic_class_sums_vector[67],
                                     phenotypic_class_sums_vector[69],
                                     phenotypic_class_sums_vector[72],
                                     phenotypic_class_sums_vector[74],
                                     phenotypic_class_sums_vector[76],
                                     phenotypic_class_sums_vector[77],
                                     phenotypic_class_sums_vector[80],
                                     phenotypic_class_sums_vector[82],
                                     phenotypic_class_sums_vector[84],
                                     phenotypic_class_sums_vector[86],
                                     phenotypic_class_sums_vector[88],
                                     phenotypic_class_sums_vector[90],
                                     phenotypic_class_sums_vector[91],
                                     phenotypic_class_sums_vector[93],
                                     phenotypic_class_sums_vector[95],
                                     phenotypic_class_sums_vector[98],
                                     phenotypic_class_sums_vector[100],
                                     phenotypic_class_sums_vector[102],
                                     phenotypic_class_sums_vector[103],
                                     phenotypic_class_sums_vector[105],
                                     phenotypic_class_sums_vector[107],
                                     phenotypic_class_sums_vector[109],
                                     phenotypic_class_sums_vector[111],
                                     phenotypic_class_sums_vector[113],
                                     phenotypic_class_sums_vector[116],
                                     phenotypic_class_sums_vector[117],
                                     phenotypic_class_sums_vector[119],
                                     phenotypic_class_sums_vector[121],
                                     phenotypic_class_sums_vector[124],
                                     phenotypic_class_sums_vector[126],
                                     phenotypic_class_sums_vector[128]),
                                 sum(phenotypic_class_sums_vector[2],
                                     phenotypic_class_sums_vector[4],
                                     phenotypic_class_sums_vector[6],
                                     phenotypic_class_sums_vector[7],
                                     phenotypic_class_sums_vector[9],
                                     phenotypic_class_sums_vector[11],
                                     phenotypic_class_sums_vector[14],
                                     phenotypic_class_sums_vector[15],
                                     phenotypic_class_sums_vector[17],
                                     phenotypic_class_sums_vector[19],
                                     phenotypic_class_sums_vector[21],
                                     phenotypic_class_sums_vector[23],
                                     phenotypic_class_sums_vector[25],
                                     phenotypic_class_sums_vector[28],
                                     phenotypic_class_sums_vector[30],
                                     phenotypic_class_sums_vector[32],
                                     phenotypic_class_sums_vector[33],
                                     phenotypic_class_sums_vector[35],
                                     phenotypic_class_sums_vector[37],
                                     phenotypic_class_sums_vector[40],
                                     phenotypic_class_sums_vector[42],
                                     phenotypic_class_sums_vector[44],
                                     phenotypic_class_sums_vector[46],
                                     phenotypic_class_sums_vector[48],
                                     phenotypic_class_sums_vector[50],
                                     phenotypic_class_sums_vector[51],
                                     phenotypic_class_sums_vector[54],
                                     phenotypic_class_sums_vector[56],
                                     phenotypic_class_sums_vector[58],
                                     phenotypic_class_sums_vector[59],
                                     phenotypic_class_sums_vector[61],
                                     phenotypic_class_sums_vector[63]),
                                 sum(phenotypic_class_sums_vector[66],
                                     phenotypic_class_sums_vector[68],
                                     phenotypic_class_sums_vector[70],
                                     phenotypic_class_sums_vector[71],
                                     phenotypic_class_sums_vector[73],
                                     phenotypic_class_sums_vector[75],
                                     phenotypic_class_sums_vector[78],
                                     phenotypic_class_sums_vector[79],
                                     phenotypic_class_sums_vector[81],
                                     phenotypic_class_sums_vector[83],
                                     phenotypic_class_sums_vector[85],
                                     phenotypic_class_sums_vector[87],
                                     phenotypic_class_sums_vector[89],
                                     phenotypic_class_sums_vector[92],
                                     phenotypic_class_sums_vector[94],
                                     phenotypic_class_sums_vector[96],
                                     phenotypic_class_sums_vector[97],
                                     phenotypic_class_sums_vector[99],
                                     phenotypic_class_sums_vector[101],
                                     phenotypic_class_sums_vector[104],
                                     phenotypic_class_sums_vector[106],
                                     phenotypic_class_sums_vector[108],
                                     phenotypic_class_sums_vector[110],
                                     phenotypic_class_sums_vector[112],
                                     phenotypic_class_sums_vector[114],
                                     phenotypic_class_sums_vector[115],
                                     phenotypic_class_sums_vector[118],
                                     phenotypic_class_sums_vector[120],
                                     phenotypic_class_sums_vector[122],
                                     phenotypic_class_sums_vector[123],
                                     phenotypic_class_sums_vector[125],
                                     phenotypic_class_sums_vector[127]),
                                 4795)
 names(forked_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(forked_single_locus_dataset,"../../Summaries and Outputs/3/10_forked_single_locus_dataset.csv",row.names = TRUE)
 
 carnation_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                         phenotypic_class_sums_vector[3],
                                         phenotypic_class_sums_vector[6],
                                         phenotypic_class_sums_vector[8],
                                         phenotypic_class_sums_vector[10],
                                         phenotypic_class_sums_vector[12],
                                         phenotypic_class_sums_vector[14],
                                         phenotypic_class_sums_vector[16],
                                         phenotypic_class_sums_vector[18],
                                         phenotypic_class_sums_vector[20],
                                         phenotypic_class_sums_vector[21],
                                         phenotypic_class_sums_vector[23],
                                         phenotypic_class_sums_vector[25],
                                         phenotypic_class_sums_vector[27],
                                         phenotypic_class_sums_vector[29],
                                         phenotypic_class_sums_vector[31],
                                         phenotypic_class_sums_vector[33],
                                         phenotypic_class_sums_vector[35],
                                         phenotypic_class_sums_vector[37],
                                         phenotypic_class_sums_vector[39],
                                         phenotypic_class_sums_vector[41],
                                         phenotypic_class_sums_vector[43],
                                         phenotypic_class_sums_vector[46],
                                         phenotypic_class_sums_vector[48],
                                         phenotypic_class_sums_vector[50],
                                         phenotypic_class_sums_vector[52],
                                         phenotypic_class_sums_vector[54],
                                         phenotypic_class_sums_vector[56],
                                         phenotypic_class_sums_vector[58],
                                         phenotypic_class_sums_vector[60],
                                         phenotypic_class_sums_vector[61],
                                         phenotypic_class_sums_vector[63]),
                                     sum(phenotypic_class_sums_vector[65],
                                         phenotypic_class_sums_vector[67],
                                         phenotypic_class_sums_vector[70],
                                         phenotypic_class_sums_vector[72],
                                         phenotypic_class_sums_vector[74],
                                         phenotypic_class_sums_vector[76],
                                         phenotypic_class_sums_vector[78],
                                         phenotypic_class_sums_vector[80],
                                         phenotypic_class_sums_vector[82],
                                         phenotypic_class_sums_vector[84],
                                         phenotypic_class_sums_vector[85],
                                         phenotypic_class_sums_vector[87],
                                         phenotypic_class_sums_vector[89],
                                         phenotypic_class_sums_vector[91],
                                         phenotypic_class_sums_vector[93],
                                         phenotypic_class_sums_vector[95],
                                         phenotypic_class_sums_vector[97],
                                         phenotypic_class_sums_vector[99],
                                         phenotypic_class_sums_vector[101],
                                         phenotypic_class_sums_vector[103],
                                         phenotypic_class_sums_vector[105],
                                         phenotypic_class_sums_vector[107],
                                         phenotypic_class_sums_vector[110],
                                         phenotypic_class_sums_vector[112],
                                         phenotypic_class_sums_vector[114],
                                         phenotypic_class_sums_vector[116],
                                         phenotypic_class_sums_vector[118],
                                         phenotypic_class_sums_vector[120],
                                         phenotypic_class_sums_vector[122],
                                         phenotypic_class_sums_vector[124],
                                         phenotypic_class_sums_vector[125],
                                         phenotypic_class_sums_vector[127]),
                                     sum(phenotypic_class_sums_vector[2],
                                         phenotypic_class_sums_vector[4],
                                         phenotypic_class_sums_vector[5],
                                         phenotypic_class_sums_vector[7],
                                         phenotypic_class_sums_vector[9],
                                         phenotypic_class_sums_vector[11],
                                         phenotypic_class_sums_vector[13],
                                         phenotypic_class_sums_vector[15],
                                         phenotypic_class_sums_vector[17],
                                         phenotypic_class_sums_vector[19],
                                         phenotypic_class_sums_vector[22],
                                         phenotypic_class_sums_vector[24],
                                         phenotypic_class_sums_vector[26],
                                         phenotypic_class_sums_vector[28],
                                         phenotypic_class_sums_vector[30],
                                         phenotypic_class_sums_vector[32],
                                         phenotypic_class_sums_vector[34],
                                         phenotypic_class_sums_vector[36],
                                         phenotypic_class_sums_vector[38],
                                         phenotypic_class_sums_vector[40],
                                         phenotypic_class_sums_vector[42],
                                         phenotypic_class_sums_vector[44],
                                         phenotypic_class_sums_vector[45],
                                         phenotypic_class_sums_vector[47],
                                         phenotypic_class_sums_vector[49],
                                         phenotypic_class_sums_vector[51],
                                         phenotypic_class_sums_vector[53],
                                         phenotypic_class_sums_vector[55],
                                         phenotypic_class_sums_vector[57],
                                         phenotypic_class_sums_vector[59],
                                         phenotypic_class_sums_vector[62],
                                         phenotypic_class_sums_vector[64]),
                                     sum(phenotypic_class_sums_vector[66],
                                         phenotypic_class_sums_vector[68],
                                         phenotypic_class_sums_vector[69],
                                         phenotypic_class_sums_vector[71],
                                         phenotypic_class_sums_vector[73],
                                         phenotypic_class_sums_vector[75],
                                         phenotypic_class_sums_vector[77],
                                         phenotypic_class_sums_vector[79],
                                         phenotypic_class_sums_vector[81],
                                         phenotypic_class_sums_vector[83],
                                         phenotypic_class_sums_vector[86],
                                         phenotypic_class_sums_vector[88],
                                         phenotypic_class_sums_vector[90],
                                         phenotypic_class_sums_vector[92],
                                         phenotypic_class_sums_vector[94],
                                         phenotypic_class_sums_vector[96],
                                         phenotypic_class_sums_vector[98],
                                         phenotypic_class_sums_vector[100],
                                         phenotypic_class_sums_vector[102],
                                         phenotypic_class_sums_vector[104],
                                         phenotypic_class_sums_vector[106],
                                         phenotypic_class_sums_vector[108],
                                         phenotypic_class_sums_vector[109],
                                         phenotypic_class_sums_vector[111],
                                         phenotypic_class_sums_vector[113],
                                         phenotypic_class_sums_vector[115],
                                         phenotypic_class_sums_vector[117],
                                         phenotypic_class_sums_vector[119],
                                         phenotypic_class_sums_vector[121],
                                         phenotypic_class_sums_vector[123],
                                         phenotypic_class_sums_vector[126],
                                         phenotypic_class_sums_vector[128]),
                                     4795)
 names(carnation_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(carnation_single_locus_dataset,"../../Summaries and Outputs/3/11_carnation_single_locus_dataset.csv",row.names = TRUE)
 
 yellow_plus_single_locus_dataset <- c(sum(phenotypic_class_sums_vector[1],
                                           phenotypic_class_sums_vector[4],
                                           phenotypic_class_sums_vector[6],
                                           phenotypic_class_sums_vector[8],
                                           phenotypic_class_sums_vector[10],
                                           phenotypic_class_sums_vector[12],
                                           phenotypic_class_sums_vector[13],
                                           phenotypic_class_sums_vector[15],
                                           phenotypic_class_sums_vector[17],
                                           phenotypic_class_sums_vector[19],
                                           phenotypic_class_sums_vector[21],
                                           phenotypic_class_sums_vector[23],
                                           phenotypic_class_sums_vector[25],
                                           phenotypic_class_sums_vector[27],
                                           phenotypic_class_sums_vector[29],
                                           phenotypic_class_sums_vector[31],
                                           phenotypic_class_sums_vector[34],
                                           phenotypic_class_sums_vector[36],
                                           phenotypic_class_sums_vector[38],
                                           phenotypic_class_sums_vector[40],
                                           phenotypic_class_sums_vector[42],
                                           phenotypic_class_sums_vector[44],
                                           phenotypic_class_sums_vector[46],
                                           phenotypic_class_sums_vector[48],
                                           phenotypic_class_sums_vector[50],
                                           phenotypic_class_sums_vector[52],
                                           phenotypic_class_sums_vector[53],
                                           phenotypic_class_sums_vector[55],
                                           phenotypic_class_sums_vector[57],
                                           phenotypic_class_sums_vector[59],
                                           phenotypic_class_sums_vector[61],
                                           phenotypic_class_sums_vector[64]),
                                       sum(phenotypic_class_sums_vector[65],
                                           phenotypic_class_sums_vector[67],
                                           phenotypic_class_sums_vector[70],
                                           phenotypic_class_sums_vector[72],
                                           phenotypic_class_sums_vector[74],
                                           phenotypic_class_sums_vector[76],
                                           phenotypic_class_sums_vector[78],
                                           phenotypic_class_sums_vector[80],
                                           phenotypic_class_sums_vector[82],
                                           phenotypic_class_sums_vector[84],
                                           phenotypic_class_sums_vector[85],
                                           phenotypic_class_sums_vector[87],
                                           phenotypic_class_sums_vector[89],
                                           phenotypic_class_sums_vector[91],
                                           phenotypic_class_sums_vector[93],
                                           phenotypic_class_sums_vector[95],
                                           phenotypic_class_sums_vector[97],
                                           phenotypic_class_sums_vector[99],
                                           phenotypic_class_sums_vector[101],
                                           phenotypic_class_sums_vector[103],
                                           phenotypic_class_sums_vector[105],
                                           phenotypic_class_sums_vector[107],
                                           phenotypic_class_sums_vector[110],
                                           phenotypic_class_sums_vector[112],
                                           phenotypic_class_sums_vector[114],
                                           phenotypic_class_sums_vector[116],
                                           phenotypic_class_sums_vector[118],
                                           phenotypic_class_sums_vector[120],
                                           phenotypic_class_sums_vector[122],
                                           phenotypic_class_sums_vector[124],
                                           phenotypic_class_sums_vector[125],
                                           phenotypic_class_sums_vector[127]),
                                       sum(phenotypic_class_sums_vector[2],
                                           phenotypic_class_sums_vector[3],
                                           phenotypic_class_sums_vector[5],
                                           phenotypic_class_sums_vector[7],
                                           phenotypic_class_sums_vector[9],
                                           phenotypic_class_sums_vector[11],
                                           phenotypic_class_sums_vector[14],
                                           phenotypic_class_sums_vector[16],
                                           phenotypic_class_sums_vector[18],
                                           phenotypic_class_sums_vector[20],
                                           phenotypic_class_sums_vector[22],
                                           phenotypic_class_sums_vector[24],
                                           phenotypic_class_sums_vector[26],
                                           phenotypic_class_sums_vector[28],
                                           phenotypic_class_sums_vector[30],
                                           phenotypic_class_sums_vector[32],
                                           phenotypic_class_sums_vector[33],
                                           phenotypic_class_sums_vector[35],
                                           phenotypic_class_sums_vector[37],
                                           phenotypic_class_sums_vector[39],
                                           phenotypic_class_sums_vector[41],
                                           phenotypic_class_sums_vector[43],
                                           phenotypic_class_sums_vector[45],
                                           phenotypic_class_sums_vector[47],
                                           phenotypic_class_sums_vector[49],
                                           phenotypic_class_sums_vector[51],
                                           phenotypic_class_sums_vector[54],
                                           phenotypic_class_sums_vector[56],
                                           phenotypic_class_sums_vector[58],
                                           phenotypic_class_sums_vector[60],
                                           phenotypic_class_sums_vector[62],
                                           phenotypic_class_sums_vector[63]),
                                       sum(phenotypic_class_sums_vector[66],
                                           phenotypic_class_sums_vector[67],
                                           phenotypic_class_sums_vector[69],
                                           phenotypic_class_sums_vector[71],
                                           phenotypic_class_sums_vector[73],
                                           phenotypic_class_sums_vector[75],
                                           phenotypic_class_sums_vector[78],
                                           phenotypic_class_sums_vector[80],
                                           phenotypic_class_sums_vector[82],
                                           phenotypic_class_sums_vector[84],
                                           phenotypic_class_sums_vector[86],
                                           phenotypic_class_sums_vector[88],
                                           phenotypic_class_sums_vector[90],
                                           phenotypic_class_sums_vector[92],
                                           phenotypic_class_sums_vector[94],
                                           phenotypic_class_sums_vector[96],
                                           phenotypic_class_sums_vector[97],
                                           phenotypic_class_sums_vector[99],
                                           phenotypic_class_sums_vector[101],
                                           phenotypic_class_sums_vector[103],
                                           phenotypic_class_sums_vector[105],
                                           phenotypic_class_sums_vector[107],
                                           phenotypic_class_sums_vector[109],
                                           phenotypic_class_sums_vector[111],
                                           phenotypic_class_sums_vector[113],
                                           phenotypic_class_sums_vector[115],
                                           phenotypic_class_sums_vector[118],
                                           phenotypic_class_sums_vector[120],
                                           phenotypic_class_sums_vector[122],
                                           phenotypic_class_sums_vector[124],
                                           phenotypic_class_sums_vector[126],
                                           phenotypic_class_sums_vector[127]),
                                       4795)
 names(yellow_plus_single_locus_dataset) <- c("female_wildtype", "male_wildtype", "female_mutant", "male_mutant", "lethal_zygote")
 write.csv(yellow_plus_single_locus_dataset,"../../Summaries and Outputs/3/12_yellow_plus_single_locus_dataset.csv",row.names = TRUE)
 

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

# 3) Join Mean_Viability ONCE (markerfree viability)
derived_multiply_marked_dataset <- derived_multiply_marked_dataset %>%
  left_join(viability_means_tgs, by = c("Treatment", "Genotype", "Brood"))

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
write.csv(multi_locus_individual_vials_dataset,"../../Summaries and Outputs/3/13_multi_locus_individual_vials_dataset.csv", row.names = FALSE)

#pooled by replicate cross------
# Create input data file for fitting Cx(Co)M model to flies pooled by replicate cross using individual_vials_dataset
pooled_dataset <- derived_multiply_marked_dataset %>%
  arrange(data_collector, treatment_id,) %>%
  group_by(data_collector, treatment_id,brood_id) %>%
  summarise(
    # keep identifiers from first vial
    experiment_id  = first(experiment_id),
    cross_id       = first(cross_id),
   # brood_id       = first(brood_id),
    vial_unique_id = first(vial_unique_id),
    
    # sum phenotypes
    across(starts_with("phenotypic_class_"), sum, na.rm = TRUE),
    
    # sum adults
    female_adult_count = sum(female_adult_count, na.rm = TRUE),
    male_adult_count   = sum(male_adult_count, na.rm = TRUE),
    total_adult_count  = sum(total_adult_count, na.rm = TRUE),
    
    # sum eggs and lethals
    CH_raw_multiply_marked_total_egg_count =
      sum(CH_raw_multiply_marked_total_egg_count, na.rm = TRUE),
    LM_raw_multiply_marked_total_egg_count =
      sum(LM_raw_multiply_marked_total_egg_count, na.rm = TRUE),
    transformed_mean_multiply_marked_total_egg_count =
      sum(transformed_mean_multiply_marked_total_egg_count, na.rm = TRUE),
    count_lethal_zygotes =
      sum(count_lethal_zygotes, na.rm = TRUE),
    
    # carry forward marker-free viability
    Mean_Viability = first(Mean_Viability),
    
    .groups = "drop"
  )

# 5) Prune columns for multilocus format #keep the vial_unique_id &count_lethal_zygotes& Mean_Viability

pruned_pooled_dataset <- pooled_dataset[ , -c(1:5, 135:140)]

# 6) Rename + duplicate viability for sex-specific slots

female_viability_vector <- pruned_pooled_dataset$Mean_Viability
male_viability_vector   <- pruned_pooled_dataset$Mean_Viability

multi_locus_cross_pooled_dataset <-cbind(pruned_pooled_dataset, female_viability_vector,male_viability_vector)

names(multi_locus_cross_pooled_dataset) <- c("ID", paste0("class_", 1:132))

# 7) Remove rows with missing final class

multi_locus_cross_pooled_dataset <- multi_locus_cross_pooled_dataset %>%filter(!is.na(class_132))

# 8) Write output
write.csv( multi_locus_cross_pooled_dataset, "../../Summaries and Outputs/3/14_multi_locus_pooled_by_collector_cross.csv", row.names = FALSE)


