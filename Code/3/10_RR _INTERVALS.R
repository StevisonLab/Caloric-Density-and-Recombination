# --- Set working directory ---
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

library(dplyr)
library(tidyr)
library(emmeans)
library(broom)

# --- Load your data ---
df <- read.csv("df.csv")

# Ensure MLE columns are numeric
df <- df %>% mutate(across(starts_with("MLE_"), as.numeric))

# --- Define markers and intervals ---
markers <- c("y","sc","cv","v","f","car")
intervals <- list()
for(i in 1:(length(markers)-1)){
  intervals[[i]] <- list(
    name = paste(markers[i], markers[i+1], sep="-"),
    xcol = paste0("MLE_x", i)
  )
}

# --- Table 2: Recombination rates ---
table2 <- lapply(intervals, function(intv){
  df %>%
    group_by(GB, Diet) %>%
    summarise(
      mean_val = mean(.data[[intv$xcol]], na.rm=TRUE),
      sd_val   = sd(.data[[intv$xcol]], na.rm=TRUE),
      .groups = "drop"
    ) %>%
    mutate(interval = intv$name)
}) %>% bind_rows() %>%
  pivot_wider(names_from = c(GB, Diet), values_from = c(mean_val, sd_val), names_sep="_")

# Combine mean ± sd into one string column
mean_cols <- grep("^mean_val", names(table2), value = TRUE)
for(col in mean_cols){
  sd_col <- sub("mean_val", "sd_val", col)
  table2[[col]] <- paste0(round(table2[[col]],2), " ± ", round(table2[[sd_col]],2))
}

# Drop SD columns
table2 <- table2 %>% select(-starts_with("sd_val"))
write.csv(table2, "Table2_RecombinationRates.csv", row.names=FALSE)

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

write.csv(table3, "Table3_CrossoverInterference.csv", row.names=FALSE)

# --- Table 4: ANOVA p-values per interval ---
table4 <- lapply(intervals, function(intv){
  df$interval_val <- df[[intv$xcol]]
  model <- aov(interval_val ~ Brood * Diet * GB * People, data=df)
  tidy(model) %>%
    select(term, df, p.value) %>%
    mutate(interval=intv$name)
}) %>% bind_rows() %>%
  pivot_wider(names_from = interval, values_from = p.value)

write.csv(table4, "Table4_ANOVA_pvalues.csv", row.names=FALSE)
