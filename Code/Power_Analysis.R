# Power Analysis

#First, we need to read in the recombination data and summarize it by vial
#then we run it through the counting model script
#then we use total genetic map length to get the effect size for each possible comparison we might want to do in the main experiment
#Goal 1: Pick the diets
#All six
#Goal 2: Pick the strains
#399,42,491
#Goal 3: Run power analysis based on choices to get sample size
#Mean squared error: 354.01 #from Spencer's experiments
#Effect sizes:
#17 for 42 
#14.32 for 399
#5.31 for 491

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(emmeans)
library(lme4)
library(lmerTest)

# install.packages("Superpower")
library(Superpower)

df <- read.csv("Raw Data/Pilots/Pilot_MLE_analyzed_dataset_diets.csv")

#start with a subset of strains and diets
df_3diet <- df %>%
  filter(
    GB %in% c(42, 217, 391, 399, 491),
    diet %in% c("0.5x", "1.0x", "2.0x")
  )


# Step 0: Ensure factors are set correctly
df_3diet$diet <- as.factor(df_3diet$diet)
df_3diet$GB <- as.factor(df_3diet$GB)
df_3diet$GB=droplevels(df_3diet$GB,exclude=c("57","508"))
df_3diet$Brood <- as.factor(df_3diet$Brood)
df_3diet$vialID <- as.factor(df_3diet$vialID)

# Step 1: Expand to ensure balanced repeated measures
all_combos <- expand.grid(vialID = unique(df_3diet$vialID),
                          Brood = unique(df_3diet$Brood))

df_fixed <- merge(all_combos, df_3diet, all = TRUE)

df_fixed <- df_fixed %>%
  group_by(vialID) %>%
  tidyr::fill(diet, GB, .direction = "downup") %>%
  ungroup()


# run model to get effect sizes:
model <- glm(formula =  recombination_fractions ~ diet*GB, data = df_fixed)
anova(model)
emmeans_table <- emmeans(model, specs = pairwise ~ diet, by=c("GB"), mode="kenward-roger")
emmeans_table$emmeans
ES=as.data.frame(summary(emmeans_table))[c('emmean')]
ES=subset(ES,!is.na(ES$emmean))

ES_42=ES[1:3,"emmean"]
ES_217=ES[4:6,"emmean"]
ES_391=ES[7:9,"emmean"]
ES_399=ES[10:12,"emmean"]
ES_491=ES[13:15,"emmean"]

## First run the power analysis using only 3 strains and 3 diets: ----

mu_vector = c(ES_399, #estimates for 399 from model
       ES_491, #estimates for 491 from model
       ES_42) #estimates for 42 from model
sd_value = sqrt(354.01)    # from Spencer's experiments

sample_sizes <- seq(2, 15, by = 1)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "3b*3b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2", "g3", "treatment", "t1", "t2", "t3")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  interaction_powers[i] <- result$main_results$power[3]
  cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
}

# Clean up: remove NAs just in case
valid_idx <- !is.na(interaction_powers)
plot(sample_sizes[valid_idx], interaction_powers[valid_idx],
     type = "b", pch = 19, col = "blue",
     xlab = "Sample Size per Cell",
     ylab = "Power (Interaction: group × treatment)",
     main = "Power Curve for Interaction Effect")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

loess_fit <- loess(interaction_powers ~ sample_sizes, span = 0.4)
lines(sample_sizes, predict(loess_fit), col = "darkgreen", lwd = 2)

# GxE with only 3 diets and three strains is achievable with at least 8 vials per combo which is 82% power (72 vials). 

##  So, what if we were to increase # strains or the number of diets. Let's test out number of diets first: -----

mu_vector = c(86.7,101.0, 98.3,102.3,92.1,96.0, #estimates for 399 from model
              102.5,104.8,104.8,104.8,104.8,95.7, #estimates for 491 from model with control values for intermediates
              109.5,100.95,92.4,85.27,78.13,71.8) #estimates for 42 from model with control values for intermediates
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(2, 15, by = 1)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "3b*6b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2", "g3", "treatment", "t1", "t2", "t3","t4","t5","t6")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 100, verbose = FALSE)
  
    interaction_powers[i] <- result$main_results$power[3]
    cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")

}

# Clean up: remove NAs just in case
valid_idx <- !is.na(interaction_powers)
plot(sample_sizes[valid_idx], interaction_powers[valid_idx],
     type = "b", pch = 19, col = "blue",
     xlab = "Sample Size per Cell",
     ylab = "Power (Interaction: group × treatment)",
     main = "Power Curve for Interaction Effect")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

loess_fit <- loess(interaction_powers ~ sample_sizes, span = 0.4)
lines(sample_sizes, predict(loess_fit), col = "darkgreen", lwd = 2)


# GxE with 6 diets and three strains is achievable with at least 9 vials per combo (162 vials). 


## Ok, what if we just used pilot results: ----

mu_vector = c(86.7,98.3,96.0, #estimates for 399 from model
              102.5,104.8,95.7, #estimates for 491 from model 
              109.5,92.4,71.8, #estimates for 42 from model 
              96.7,100.1,86.9, #estimates for 391 from model
              119.9,105.8,100.2) #estimates for 217 from model
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(2, 15, by = 1)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "5b*3b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2", "g3", "g4","g5",  "treatment", "t1", "t2", "t3")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  interaction_powers[i] <- result$main_results$power[3]
  cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  

}


#Ok, no luck there! What if we went big and did MORE gb and more diets?

mu_vector = c(86.7,101.0, 98.3,102.3,92.1,96.0, #estimates for 399 from model
              102.5,104.8,104.8,104.8,104.8,95.7, #estimates for 491 from model with control values for intermediates 
              109.5,100.95,92.4,85.27,78.13,71.8, #estimates for 42 from model assuming a linear step function
              96.7,98.4,100.1,95.7,91.3,86.9, #estimates for 391 from model assuming a linear step function
              119.9,112.85,105.8,103.93,102.1,100.2) #estimates for 217 from model assuming a linear step function
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(5, 50, by = 5)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "5b*6b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2", "g3", "g4","g5",  "treatment", "t1", "t2", "t3","t4","t5","t6")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  interaction_powers[i] <- result$main_results$power[3]
  cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  
}


## Ok, what if we did only 2 genotypes and 6 diets? ----

mu_vector = c(86.7,101.0, 98.3,102.3,92.1,96.0, #estimates for 399 from model
              109.5,100.95,92.4,85.27,78.13,71.8) #estimates for 42 from model assuming a linear step function
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(2, 15, by = 1)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "2b*6b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2",  "treatment", "t1", "t2", "t3","t4","t5","t6")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  interaction_powers[i] <- result$main_results$power[3]
  cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  
}

# Clean up: remove NAs just in case
valid_idx <- !is.na(interaction_powers)
plot(sample_sizes[valid_idx], interaction_powers[valid_idx],
     type = "b", pch = 19, col = "blue",
     xlab = "Sample Size per Cell",
     ylab = "Power (Interaction: group × treatment)",
     main = "Power Curve for Interaction Effect")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

loess_fit <- loess(interaction_powers ~ sample_sizes, span = 0.4)
lines(sample_sizes, predict(loess_fit), col = "darkgreen", lwd = 2)

# GxE is for 2 strains and 6 diets achievable with 8 replicates which gets to 87.2% power (96 vials).

## Now what if we reduced the number of diets?: -----

# remove 2.0x:
mu_vector = c(86.7,101.0, 98.3,102.3,92.1, #estimates for 399 from model
              102.5,104.8,104.8,102.5431,100.2681, #estimates for 491 from model with control values for intermediates 
              109.5,100.95,92.4,85.27,78.13) #estimates for 42 from model assuming a linear step function
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(2, 15, by = 1)
interaction_powers <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  current_n <- sample_sizes[i]
  
  design <- ANOVA_design(
    design = "3b*5b",
    n = current_n,
    mu = mu_vector,
    sd = sd_value,
    labelnames = c("group", "g1", "g2", "g3", "treatment", "t1", "t2", "t3","t4","t5")
  )
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  interaction_powers[i] <- result$main_results$power[3]
  cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  
}

# Clean up: remove NAs just in case
valid_idx <- !is.na(interaction_powers)
plot(sample_sizes[valid_idx], interaction_powers[valid_idx],
     type = "b", pch = 19, col = "blue",
     xlab = "Sample Size per Cell",
     ylab = "Power (Interaction: group × treatment)",
     main = "Power Curve for Interaction Effect")
abline(h = 0.8, col = "red", lty = 2)  # reference line for 80% power

loess_fit <- loess(interaction_powers ~ sample_sizes, span = 0.4)
lines(sample_sizes, predict(loess_fit), col = "darkgreen", lwd = 2)

# GxE is with 4 diets (no 2.0x or 0.5x) and 2 strains is NOT achievable with 15 replicates.
# GxE is with 5 diets (no 2.0x) and 2 strains is achievable with 9 replicates getting 81.8% power (90 vials).
# GxE is with 5 diets (no 2.0x) and 3 strains is achievable with 10 replicates getting 83% power (150 vials).


# Sim R Power Analysis (does NOT work!!) ----------------------------------------------------

# Load required libraries
library(simr)
library(lme4)
library(dplyr)
library(tidyr)

#start with a subset of strains and diets
df_3diet <- df %>%
  filter(
    GB %in% c(42, 217, 391, 399, 491),
    diet %in% c("0.5x", "1.0x", "2.0x")
  )


# Step 0: Ensure factors are set correctly
df_3diet$diet <- as.factor(df_3diet$diet)
df_3diet$GB <- as.factor(df_3diet$GB)
df_3diet$GB=droplevels(df_3diet$GB,exclude=c("57","508"))
df_3diet$Brood <- as.factor(df_3diet$Brood)
df_3diet$vialID <- as.factor(df_3diet$vialID)

# Step 1: Expand to ensure balanced repeated measures
all_combos <- expand.grid(vialID = unique(df_3diet$vialID),
                          Brood = unique(df_3diet$Brood))

df_fixed <- merge(all_combos, df_3diet, all = TRUE)

df_fixed <- df_fixed %>%
  group_by(vialID) %>%
  tidyr::fill(diet, GB, .direction = "downup") %>%
  ungroup()

# Optional: Check for balance
# table(df_fixed$vialID, df_fixed$Brood) #should be NO zeros!!!
# Check levels
# table(df_fixed$diet)
# table(df_fixed$GB)
# table(df_fixed$vialID)

# Step 2: Fit model on cleaned data
model <- lmer(recombination_fractions ~ diet * GB + (1 | vialID), data = df_fixed, method = "REML")
anova(model)

# Step 3: Prepare for simr
#model_sim <- makeLmer(model)
model_sim <- model

# Step 4: Power analysis for varying numbers of vials
vial_counts <- seq(420, 3100, by = 50)
vial_counts <- seq(10, 30, by = 5)
interaction_powers <- numeric(length(vial_counts))

for (i in seq_along(vial_counts)) {
  current_n <- vial_counts[i]
  mod_ext <- try(extend(model, along = "vialID", n = current_n), silent = TRUE)
  
  if (inherits(mod_ext, "try-error")) {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Model extension failed.\n")
    next
  }
  
  sim <- try(powerSim(mod_ext, test = fixed("diet:GB", method = "f"), nsim = 50), silent = TRUE)
  
  if (inherits(sim, "try-error")) {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Simulation failed.\n")
  } else {
    interaction_powers[i] <- sim$power
    cat("n =", current_n, "→ Power =", sim$power, "\n")
  }
}

# Step 5: Plot the power curve
plot(vial_counts, interaction_powers, type = "b", pch = 19, col = "blue",
     xlab = "Number of Vials", ylab = "Power (Interaction: diet × GB)",
     main = "Power Curve for diet × GB Interaction")
abline(h = 0.8, col = "red", lty = 2)



#Goal 4: Design the specifics of the experiment and downstream data analysis
