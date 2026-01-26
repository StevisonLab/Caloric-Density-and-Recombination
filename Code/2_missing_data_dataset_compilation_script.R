#setwd("Raw Data/") # Please use Raw Data as working directory
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

write.csv(compiled_raw_marker_free_dataset, "../../Summaries and Outputs/D_compiled_raw_marker_free_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting the average and round these two counts to the nearest whole number
#Lauren=> Y=0.83X-4.8
#Connor=> Y=0.92X-6.4
CH_transformed_marker_free_total_egg_count <- (CH_raw_marker_free_total_egg_count - 6.479979) / 0.9224861
LM_transformed_marker_free_total_egg_count <- (LM_raw_marker_free_total_egg_count - 4.866269) / 0.8351445
transformed_mean_marker_free_total_egg_count <- round((CH_transformed_marker_free_total_egg_count + LM_transformed_marker_free_total_egg_count) / 2)
derived_marker_free_dataset <- cbind(compiled_raw_marker_free_dataset, transformed_mean_marker_free_total_egg_count)
write.csv(derived_marker_free_dataset, "../../Summaries and Outputs/E_derived_marker_free_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_marker_free_female_viability <- derived_marker_free_dataset$female_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_male_viability <- derived_marker_free_dataset$male_adult_count / (derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count / 2)
derived_marker_free_total_viability <- derived_marker_free_dataset$total_adult_count / derived_marker_free_dataset$transformed_mean_marker_free_total_egg_count
marker_free_viability_dataset <- cbind(derived_marker_free_dataset, derived_marker_free_female_viability, derived_marker_free_male_viability, derived_marker_free_total_viability)
write.csv(marker_free_viability_dataset, "../../Summaries and Outputs/F_marker_free_viability_dataset.csv", row.names = TRUE)

# Compute the overall viabilities for reporting in main text
overall_marker_free_viability <- mean(derived_marker_free_total_viability)
female_marker_free_viability <- mean(derived_marker_free_female_viability)
male_marker_free_viability <- mean(derived_marker_free_male_viability)

# Generate a data frame for sex-specific viability effects in the marker free crosses
female_vector <- rep("Female", nrow(marker_free_viability_dataset))
male_vector <- rep("Male", nrow(marker_free_viability_dataset))
marker_free_female_viabilities <- data.frame(
  Sex = female_vector,
  Cross = marker_free_viability_dataset$cross_id,
  Brood = marker_free_viability_dataset$brood_id,
  Treatment = marker_free_viability_dataset$treatment,
  Genotype = marker_free_viability_dataset$genotype,
  Viability = marker_free_viability_dataset$derived_marker_free_female_viability
)
marker_free_male_viabilities <- data.frame(
  Sex = male_vector,
  Cross = marker_free_viability_dataset$cross_id,
  Brood = marker_free_viability_dataset$brood_id,
  Treatment = marker_free_viability_dataset$treatment,
  Genotype = marker_free_viability_dataset$genotype,
  Viability = marker_free_viability_dataset$derived_marker_free_male_viability
)
marker_free_viability_anova_input <- rbind(marker_free_female_viabilities, marker_free_male_viabilities)

# Convert variable "Cross" to a categorical factor and "Brood" to a numeric value
marker_free_viability_anova_input$Cross <- factor(marker_free_viability_anova_input$Cross)
marker_free_viability_anova_input$Brood <- factor(marker_free_viability_anova_input$Brood)
marker_free_viability_anova_input$Brood <- as.numeric(marker_free_viability_anova_input$Brood)
marker_free_viability_anova_input$Treatment <- factor(marker_free_viability_anova_input$Treatment)
marker_free_viability_anova_input$Genotype <- factor(marker_free_viability_anova_input$Genotype)

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

write.csv(compiled_raw_multiply_marked_dataset, "../../Summaries and Outputs/4_compiled_raw_multiply_marked_dataset.csv", row.names = TRUE)

# Transform individual raw egg counts to address undercounting
CH_transformed_multiply_marked_total_egg_count <- (CH_raw_multiply_marked_total_egg_count - 6.479979) / 0.9224861
LM_transformed_multiply_marked_total_egg_count <- (LM_raw_multiply_marked_total_egg_count - 4.866269) / 0.8351445
transformed_mean_multiply_marked_total_egg_count <- round((CH_transformed_multiply_marked_total_egg_count + LM_transformed_multiply_marked_total_egg_count) / 2)
derived_multiply_marked_dataset <- cbind(compiled_raw_multiply_marked_dataset, transformed_mean_multiply_marked_total_egg_count)
write.csv(derived_multiply_marked_dataset, "../../Summaries and Outputs/5_derived_multiply_marked_dataset.csv", row.names = TRUE)

# Calculate the per vial viability of females, males, and total under the 1:1 primary sex ratio assumption
derived_multiply_marked_female_viability <- derived_multiply_marked_dataset$female_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_male_viability <- derived_multiply_marked_dataset$male_adult_count / (derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count / 2)
derived_multiply_marked_total_viability <- derived_multiply_marked_dataset$total_adult_count / derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count
multiply_marked_viability_dataset <- cbind(derived_multiply_marked_dataset, derived_multiply_marked_female_viability, derived_multiply_marked_male_viability, derived_multiply_marked_total_viability)
write.csv(multiply_marked_viability_dataset, "../../Summaries and Outputs/6_multiply_marked_viability_dataset.csv", row.names = TRUE)

# Compute the overall mean viabilities for reporting in main text 
overall_multiply_marked_viability <- mean(derived_multiply_marked_total_viability)
female_multiply_marked_viability <- mean(derived_multiply_marked_female_viability)
male_multiply_marked_viability <- mean(derived_multiply_marked_male_viability)

# Generate a data frame for sex-specific viability effects in the multiply marked crosses
female_vector <- rep("Female", nrow(multiply_marked_viability_dataset))
male_vector <- rep("Male", nrow(multiply_marked_viability_dataset))
multiply_marked_female_viabilities <- data.frame(
  Sex = female_vector,
  Cross = multiply_marked_viability_dataset$cross_id,
  Brood = multiply_marked_viability_dataset$brood_id,
  Treatment = multiply_marked_viability_dataset$treatment,
  Genotype = multiply_marked_viability_dataset$genotype,
  Viability = multiply_marked_viability_dataset$derived_multiply_marked_female_viability
)
multiply_marked_male_viabilities <- data.frame(
  Sex = male_vector,
  Cross = multiply_marked_viability_dataset$cross_id,
  Brood = multiply_marked_viability_dataset$brood_id,
  Treatment = multiply_marked_viability_dataset$treatment,
  Genotype = multiply_marked_viability_dataset$genotype,
  Viability = multiply_marked_viability_dataset$derived_multiply_marked_male_viability
)
multiply_marked_viability_anova_input <- rbind(multiply_marked_female_viabilities, multiply_marked_male_viabilities)

# Convert variable "Cross" to a categorical factor and "Brood" to a numeric value
multiply_marked_viability_anova_input$Cross <- factor(multiply_marked_viability_anova_input$Cross)
multiply_marked_viability_anova_input$Brood <- factor(multiply_marked_viability_anova_input$Brood)
multiply_marked_viability_anova_input$Brood <- as.numeric(multiply_marked_viability_anova_input$Brood)
multiply_marked_viability_anova_input$Treatment <- factor(multiply_marked_viability_anova_input$Treatment)
multiply_marked_viability_anova_input$Genotype <- factor(multiply_marked_viability_anova_input$Genotype)

# Rows 61:90 and 161:190 were G, H, I broods
# cleaned_multiply_marked_viability_anova_input <- multiply_marked_viability_anova_input[-c(61:90,161:190), ]
# For now, keeping all rows:
cleaned_multiply_marked_viability_anova_input <- multiply_marked_viability_anova_input


# Input for Counting Model ------------------------------------------------
#individual vials-------------
# Create input data file for fitting Cx(Co)M model to flies pooled by individual vials
# Calculate the number of lethal zygotes by subtracting total adults from mean transformed egg count
# Add marker-free viabilities as classes 130, 131, 132, rename columns to conform to standard format
# Create vector of IDs to remove
ids_to_remove <- c(
  "MS_5_88_A",
  paste0("MS_5_86_", LETTERS[1:7]),
  paste0("MS_5_87_", LETTERS[1:7]),
  paste0("MS_5_88_", LETTERS[1:7])
)
# Remove rows with IDs in ids_to_remove
multi_locus_individual_vials_dataset <- multi_locus_individual_vials_dataset[!multi_locus_individual_vials_dataset$ID %in% ids_to_remove, ]
# Create cleaned dataset
cleaned_derived_multiply_marked_dataset <- multi_locus_individual_vials_dataset
# Calculate lethal zygotes
count_lethal_zygotes <- cleaned_derived_multiply_marked_dataset$transformed_mean_multiply_marked_total_egg_count - cleaned_derived_multiply_marked_dataset$total_adult_count
# Add lethal zygotes column
cleaned_derived_multiply_marked_dataset <- cbind(cleaned_derived_multiply_marked_dataset, count_lethal_zygotes)
# Remove specified columns
pruned_derived_multiply_marked_dataset <- cleaned_derived_multiply_marked_dataset[, -c(1,2,3,5,6,135,136,137,138,139,140)]

#change
total_viability_vector <- rep(0.921, 672)
female_viability_vector <- rep(0.921, 672)
male_viability_vector <- rep(0.921, 672)

multi_locus_individual_vials_dataset <- cbind(pruned_derived_multiply_marked_dataset,total_viability_vector,female_viability_vector,male_viability_vector)
names(multi_locus_individual_vials_dataset) <- c("ID", paste0("class_", 1:132))
ids_to_remove <- c(
  paste0("MS_5_86_", LETTERS[1:7]),
  paste0("MS_5_88_", LETTERS[1:7])
)

# Remove those rows
multi_locus_individual_vials_dataset <-
  multi_locus_individual_vials_dataset[
    !multi_locus_individual_vials_dataset$ID %in% ids_to_remove, ]

write.csv(multi_locus_individual_vials_dataset,"../../Summaries and Outputs/13_multi_locus_individual_vials_dataset.csv",row.names = FALSE)


######???
#pooled by replicate cross------
# Create input data file for fitting Cx(Co)M model to flies pooled by replicate cross using individual_vials_dataset
cross_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
MS_14_1_rows <- c(1, 11, 21, 31, 41, 51, 61)
MS_14_1 <- colSums(cross_pooling_matrix[MS_14_1_rows, ])
MS_14_3_rows <- c(2, 12, 22, 32, 42, 52, 62)
MS_14_3 <- colSums(cross_pooling_matrix[MS_14_3_rows, ])
MS_14_5_rows <- c(3, 13, 23, 33, 43, 53, 63)
MS_14_5 <- colSums(cross_pooling_matrix[MS_14_5_rows, ])
MS_14_6_rows <- c(4, 14, 24, 34, 44, 54, 64)
MS_14_6 <- colSums(cross_pooling_matrix[MS_14_6_rows, ])
MS_14_8_rows <- c(5, 15, 25, 35, 45, 55, 65)
MS_14_8 <- colSums(cross_pooling_matrix[MS_14_8_rows, ])
MS_14_12_rows <- c(6, 16, 26, 36, 46, 56, 66)
MS_14_12 <- colSums(cross_pooling_matrix[MS_14_12_rows, ])
MS_14_14_rows <- c(7, 17, 27, 37, 47, 57, 67)
MS_14_14 <- colSums(cross_pooling_matrix[MS_14_14_rows, ])
MS_14_17_rows <- c(8, 18, 28, 38, 48, 58, 68)
MS_14_17 <- colSums(cross_pooling_matrix[MS_14_17_rows, ])
MS_14_18_rows <- c(9, 19, 29, 39, 49, 59, 69)
MS_14_18 <- colSums(cross_pooling_matrix[MS_14_18_rows, ])
MS_14_19_rows <- c(10, 20, 30, 40, 50, 60, 70)
MS_14_19 <- colSums(cross_pooling_matrix[MS_14_19_rows, ])
multi_locus_cross_pooled_dataset <- rbind(MS_14_1, MS_14_3, MS_14_5, MS_14_6, MS_14_8, MS_14_12, MS_14_14, MS_14_17, MS_14_18, MS_14_19)
unique_cross_id <- c("MS_14_1", "MS_14_3", "MS_14_5", "MS_14_6", "MS_14_8", "MS_14_12", "MS_14_14", "MS_14_17", "MS_14_18", "MS_14_19")
cross_viability_vector <- c(0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921)
multi_locus_cross_pooled_dataset <- as.data.frame(cbind(unique_cross_id, multi_locus_cross_pooled_dataset, cross_viability_vector, cross_viability_vector, cross_viability_vector))
names(multi_locus_cross_pooled_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_cross_pooled_dataset,"../../Summaries and Outputs/14_multi_locus_cross_pooled_dataset.csv",row.names = TRUE)

#pooled by brooding------
# Create input data file for fitting Cx(Co)M model to flies pooled by brooding period using individual_vials_dataset
brood_pooling_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
MS_14_A_rows <- c(1:10)
MS_14_A <- colSums(brood_pooling_matrix[MS_14_A_rows, ])
MS_14_B_rows <- c(11:20)
MS_14_B <- colSums(brood_pooling_matrix[MS_14_B_rows, ])
MS_14_C_rows <- c(21:30)
MS_14_C <- colSums(brood_pooling_matrix[MS_14_C_rows, ])
MS_14_D_rows <- c(31:40)
MS_14_D <- colSums(brood_pooling_matrix[MS_14_D_rows, ])
MS_14_E_rows <- c(41:50)
MS_14_E <- colSums(brood_pooling_matrix[MS_14_E_rows, ])
MS_14_F_rows <- c(51:60)
MS_14_F <- colSums(brood_pooling_matrix[MS_14_F_rows, ])
MS_14_J_rows <- c(61:70)
MS_14_J <- colSums(brood_pooling_matrix[MS_14_J_rows, ])
multi_locus_brood_pooled_dataset <- rbind(MS_14_A, MS_14_B, MS_14_C, MS_14_D, MS_14_E, MS_14_F, MS_14_J)
unique_brood_id <- c("MS_14_A", "MS_14_B", "MS_14_C", "MS_14_D", "MS_14_E", "MS_14_F", "MS_14_J")
brood_viability_vector <- c(0.921, 0.921, 0.921, 0.921, 0.921, 0.921, 0.921)
multi_locus_brood_pooled_dataset <- as.data.frame(cbind(unique_brood_id, multi_locus_brood_pooled_dataset, brood_viability_vector, brood_viability_vector, brood_viability_vector))
names(multi_locus_brood_pooled_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_brood_pooled_dataset,"../../Summaries and Outputs/15_multi_locus_brood_pooled_dataset.csv",row.names = TRUE)

#individual_vials_dataset--------
# Create input data file for fitting Cx(Co)M model to flies pooled for the full experiment using individual_vials_dataset
full_experiment_matrix <- multi_locus_individual_vials_dataset[,-c(1,131,132,133)]
MS_14_rows <- c(1:70)
MS_14 <- colSums(full_experiment_matrix[MS_14_rows, ])
multi_locus_full_experiment_dataset <- rbind(MS_14, MS_14)
full_experiment_id <- c("full_experiment", "full_experiment")
full_experiment_vector <- c(0.921, 0.921)
multi_locus_full_experiment_dataset <- as.data.frame(cbind(full_experiment_id, multi_locus_full_experiment_dataset, full_experiment_vector, full_experiment_vector, full_experiment_vector))
multi_locus_full_experiment_dataset <- multi_locus_full_experiment_dataset[-2,]
names(multi_locus_full_experiment_dataset) <- c("ID", paste0("class_", 1:132))
write.csv(multi_locus_full_experiment_dataset,"../../Summaries and Outputs/16_multi_locus_full_experiment_dataset.csv",row.names = TRUE)

# Create a vector of observed values for fitting CxCoM model
single_experimental_unit <- as.numeric(multi_locus_full_experiment_dataset)
observed_count <- c(single_experimental_unit[130],single_experimental_unit[2:129])


