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

# install.packages("Superpower")
library(Superpower)

mu_vector = c(86.7,98.3,96.0, #estimates for 399 from model
       102.5,104.8,95.7, #estimates for 491 from model
       109.5,92.4,71.8) #estimates for 42 from model
sd_value = sqrt(354.01)    # from Spencer's experiments

#base_design <- ANOVA_design(
#  design = "3b*3b",  # 3 groups x 3 treatments, both between-subjects
#  n = 2,            # initial sample size per condition
#  mu = c(86.7,98.3,96.0, #estimates for 399 from model
#         102.5,104.8,95.7, #estimates for 491 from model
#         109.5,92.4,71.8), #estimates for 42 from model       
  # specify means for each of the 9 cells
#  sd = sqrt(354.01),    # from Spencer's experiments
#  labelnames = c("group", "S399", "S491", "S42", "treatment", "0.5x", "1.0x", "2.0x")
#)

#result <- ANOVA_power(ANOVA_design, alpha_level = 0.05, nsims = 1000)
#print(result)

sample_sizes <- seq(5, 50, by = 5)
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
  
  # Check if the interaction term is present
  interaction_index <- which(result$power$Effect == "group:treatment")
  
  if (length(interaction_index) > 0) {
    interaction_powers[i] <- result$power$power[interaction_index]
    cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  } else {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Interaction Power = NA (effect not found)\n")
  }
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

# Ok, not able to get GxE with only 3 diets and three strains. So, either we increase # strains or the number of diets. Let's test out number of diets first:

mu_vector = c(86.7,101.0, 98.3,102.3,92.1,96.0, #estimates for 399 from model
              102.5,104.8,104.8,104.8,104.8,95.7, #estimates for 491 from model with control values for intermediates
              109.5,100.95,92.4,85.27,78.13,71.8) #estimates for 42 from model with control values for intermediates
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(5, 50, by = 5)
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
  
  result <- ANOVA_power(design, alpha_level = 0.05, nsims = 500, verbose = FALSE)
  
  # Check if the interaction term is present
  interaction_index <- which(result$power$Effect == "group:treatment")
  
  if (length(interaction_index) > 0) {
    interaction_powers[i] <- result$power$power[interaction_index]
    cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  } else {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Interaction Power = NA (effect not found)\n")
  }
}

# Even a sample size as large as 50 would not yield a significant result


# Ok, not able to get GxE with 6 diets and three strains. So, now we increase # strains based on pilot results

mu_vector = c(86.7,98.3,96.0, #estimates for 399 from model
              102.5,104.8,95.7, #estimates for 491 from model 
              109.5,92.4,71.8, #estimates for 42 from model 
              96.7,100.1,86.9, #estimates for 391 from model
              119.9,105.8,100.2) #estimates for 217 from model
sd_value = sqrt(354.01)    # from Spencer's experiments


sample_sizes <- seq(5, 50, by = 5)
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
  
  # Check if the interaction term is present
  interaction_index <- which(result$power$Effect == "group:treatment")
  
  if (length(interaction_index) > 0) {
    interaction_powers[i] <- result$power$power[interaction_index]
    cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  } else {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Interaction Power = NA (effect not found)\n")
  }
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
  
  # Check if the interaction term is present
  interaction_index <- which(result$power$Effect == "group:treatment")
  
  if (length(interaction_index) > 0) {
    interaction_powers[i] <- result$power$power[interaction_index]
    cat("n =", current_n, "→ Interaction Power =", interaction_powers[i], "\n")
  } else {
    interaction_powers[i] <- NA
    cat("n =", current_n, "→ Interaction Power = NA (effect not found)\n")
  }
}



#Goal 4: Design the specifics of the experiment and downstream data analysis
