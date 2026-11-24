# --- Set working directory ---
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

library(dplyr)
library(tidyr)
library(emmeans)
library(broom)

# --- Load your data ---
df <- read.csv("23_df_RR_cleaned.csv")

# Ensure MLE columns are numeric
df <- df %>% mutate(across(starts_with("MLE_"), as.numeric))

# --- Define markers and intervals ---
#markers <- c("sc","cv","v","f","car","y+")
markers <- c("y+","car","f","v","cv","sc")


intervals <- list()
for(i in 1:(length(markers)-1)){
  intervals[[i]] <- list(
    name = paste(markers[i], markers[i+1], sep="-"),
    xcol = paste0("MLE_x", i)
  )
}

# --- Table 6: Recombination fractions (Mean ± SE) ---

table6_summary <- lapply(intervals, function(intv){
  df %>%
    group_by(GB, Diet) %>%
    summarise(
      n = sum(!is.na(.data[[intv$xcol]])),
      mean_val = mean(.data[[intv$xcol]], na.rm = TRUE),
      sd_val = sd(.data[[intv$xcol]], na.rm = TRUE),
      se_val = sd_val / sqrt(n),
      min_val = min(.data[[intv$xcol]], na.rm = TRUE),
      max_val = max(.data[[intv$xcol]], na.rm = TRUE),
      above_05 = sum(.data[[intv$xcol]] > 0.5, na.rm = TRUE),
      below_0 = sum(.data[[intv$xcol]] < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(interval = intv$name)
}) %>% bind_rows()

# Check the original values, SD, SE and sample sizes
print(table6_summary, n = Inf, width = Inf)

# Create the manuscript table
table6 <- table6_summary %>%
  mutate(
    Mean_SE = sprintf("%.3f ± %.3f", mean_val, se_val)
  ) %>%
  select(interval, GB, Diet, Mean_SE) %>%
  pivot_wider(
    names_from = c(GB, Diet),
    values_from = Mean_SE,
    names_sep = "_"
  )

# View corrected table
print(table6, width = Inf)

View(table6)

# Save corrected table
write.csv(
  table6,
  "33_Table6_RecombinationRates_intervals_SE.csv",
  row.names = FALSE
)

# Save detailed diagnostic results
write.csv(
  table6_summary,
  "33_Table6_Diagnostic.csv",
  row.names = FALSE
)

# --- Table 3: Crossover interference (1 - coefficient of coincidence) ---
double_intervals <- list(
  list(name="c12", xcols1="MLE_x1", xcols2="MLE_x2"),
  list(name="c23", xcols1="MLE_x2", xcols2="MLE_x3"),
  list(name="c34", xcols1="MLE_x3", xcols2="MLE_x4"),
  list(name="c45", xcols1="MLE_x4", xcols2="MLE_x5")
)

table3 <- lapply(double_intervals, function(intv){
  df %>%
    mutate(
      obs_double = rowSums(select(., all_of(c(intv$xcols1,intv$xcols2)))),
      exp_double = .data[[intv$xcols1]] * .data[[intv$xcols2]],
      # Avoid division by zero or extremely small numbers
      interference = ifelse(exp_double < 1e-6, NA, 1 - (obs_double / exp_double))
    ) %>%
    group_by(GB, Diet) %>%
    summarise(
      mean_val = mean(interference, na.rm=TRUE),
      sd_val   = sd(interference, na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(interval = intv$name)
}) %>% bind_rows() %>%
  pivot_wider(names_from = c(GB, Diet), values_from = c(mean_val, sd_val), names_sep="_")

# Combine mean ± sd into one string
mean_cols <- grep("^mean_val", names(table3), value = TRUE)
for(col in mean_cols){
  sd_col <- sub("mean_val", "sd_val", col)
  table3[[col]] <- ifelse(
    is.na(table3[[col]]), 
    NA, 
    paste0(round(table3[[col]],2), " ± ", round(table3[[sd_col]],2))
  )
}

# Drop SD columns
table3 <- table3 %>% select(-starts_with("sd_val"))

write.csv(table3, "34_Table3_CrossoverInterference.csv", row.names=FALSE)

# --- Table 4: ANOVA p-values per interval ---
table4 <- lapply(intervals, function(intv){
  df$interval_val <- df[[intv$xcol]]
  model <- aov(interval_val ~ Brood * Diet * GB * People, data=df)
  tidy(model) %>%
    select(term, df, p.value) %>%
    mutate(interval=intv$name)
}) %>% bind_rows() %>%
  pivot_wider(names_from = interval, values_from = p.value)

write.csv(table4, "35_Table4_ANOVA_pvalues.csv", row.names=FALSE)
#people effect:

library(lme4)
library(emmeans)
library(multcompView)
library(dplyr)

# Make sure your interval columns are numeric
df <- df %>% mutate(across(starts_with("MLE_x"), as.numeric))

# Define intervals (sum the markers for each interval)

intervals <- list(
"car-y+"= c("MLE_x1","MLE_x2"),
"f-car" = c("MLE_x2","MLE_x3"),
"v-f" = c("MLE_x3","MLE_x4"),
"cv-v" = c("MLE_x4","MLE_x5"),
"sc-cv"= c("MLE_x5") 
)


cld_people_list <- lapply(names(intervals), function(intv_name){
  
  # Sum columns for this interval
  df$interval_sum <- rowSums(select(df, all_of(intervals[[intv_name]])), na.rm = TRUE)
  
  # Fit model with People as main effect (you can include other factors if desired)
  model <- aov(interval_sum ~ People, data = df)
  
  # Estimated marginal means
  emm <- emmeans(model, ~ People)
  
  # Compact Letter Display (significance grouping)
  cld_df <- cld(emm, Letters = letters)
  
  cld_df %>%
    select(People, .group) %>%
    mutate(interval = intv_name)
  
}) %>% bind_rows()

# View the CLD for People
cld_people_list

#plot
library(ggplot2)
library(dplyr)

# Prepare mean recombination per People for each interval
interval_means <- lapply(names(intervals), function(intv_name){
  df$interval_sum <- rowSums(select(df, all_of(intervals[[intv_name]])), na.rm = TRUE)
  
  df %>%
    group_by(People) %>%
    summarise(mean = mean(interval_sum), sd = sd(interval_sum), .groups = "drop") %>%
    mutate(interval = intv_name)
}) %>% bind_rows()

# Join with CLD letters
plot_df <- interval_means %>%
  left_join(cld_people_list, by = c("People", "interval"))

# Plot with fixed y-axis
ggplot(plot_df, aes(x = People, y = mean, fill = People)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, position = position_dodge(0.6)) +
  geom_text(aes(label = .group, y = mean + sd + 0.02), size = 5) +
  facet_wrap(~interval, scales = "fixed") +  # fixed y-axis for all intervals
  labs(title = "Recombination Rate by People with CLD",
       y = "Recombination rate (sum of markers)",
       x = "People") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
