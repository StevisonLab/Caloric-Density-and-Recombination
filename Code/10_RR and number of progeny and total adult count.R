setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

# Recombination data (the df you use in ANOVAs)
df <- read.csv("Nutrition_countingmodel_outputs_all.csv", stringsAsFactors = FALSE)

df <- df %>% extract(original_ID,
                     into = c("experiment", "vialID", "Brood"),
                     regex = "^(.*)_([0-9]+)_([A-Za-z])$"
) %>%
  rename(recombination_fractions = MLE_X_length)

df <- df %>%filter(recombination_fractions < 1000)  

df$vialID  <- as.character(df$vialID)
df$Brood  <- as.factor(df$Brood)

# 2. Read and prepare treatment key for the first 2 
treatment_key <- read.csv("Treatmentskeys.csv")
treatment_key <- treatment_key %>%
  rename(vialID = `Vial.ID`, GB = `F1.Genotype`, Diet = `F1.developmental`)

# Egg / progeny data
progeny_df <- read.csv("5_derived_multiply_marked_dataset.csv", 
                       stringsAsFactors = FALSE)

progeny_df$vialID  <- progeny_df$cross_id 
progeny_df$Brood <- progeny_df$brood_id

progeny_df$vialID <- as.character(progeny_df$vialID)
progeny_df$Brood <- as.factor(progeny_df$Brood)

#subseting the data to ones that we want 
progeny_subset <- progeny_df[, c(
  "vialID",
  "Brood",
  "transformed_mean_multiply_marked_total_egg_count",
  "total_adult_count"
)]

#merge df +key
df <- merge(
  df,
  treatment_key,
  by = c("vialID"),
  all.x = TRUE
)

merged_df <- merge(
  df,
  progeny_subset,
  by = c("vialID", "Brood"),
  all.x = TRUE
)

#making numeric 
merged_df$transformed_mean_multiply_marked_total_egg_count <-
  as.numeric(merged_df$transformed_mean_multiply_marked_total_egg_count)

merged_df$total_adult_count <-
  as.numeric(merged_df$total_adult_count)

merged_df$People <-
  as.factor(merged_df$People)

#model
lm1a <- lm(
  recombination_fractions ~ 
    transformed_mean_multiply_marked_total_egg_count +
    GB + Diet + People,
  data = merged_df
)

lm1b <- lm(
  recombination_fractions ~ 
    transformed_mean_multiply_marked_total_egg_count +
    Brood + Diet + People,
  data = merged_df
)

summary(lm1a)
summary(lm1b)

lm_no_prog <- lm(
  recombination_fractions ~ GB + Diet + Brood + People,
  data = merged_df
)

anova(lm_no_prog, lm1a)
anova(lm_no_prog, lm1b)

plot(
  merged_df$transformed_mean_multiply_marked_total_egg_count,
  merged_df$recombination_fractions,
  pch = 16,
  col = as.factor(merged_df$People),
  xlab = "Progeny (transformed egg count)",
  ylab = "Recombination fraction"
)

abline(
  lm(recombination_fractions ~ transformed_mean_multiply_marked_total_egg_count,
     data = merged_df),
  lwd = 2
)
# adult counts 
lm_adults <- lm(
  recombination_fractions ~ 
    total_adult_count + GB + Diet + Brood + People,
  data = merged_df
)

summary(lm_adults)

#plot
merged_df$People <- factor(
  merged_df$People,
  levels = c("C", "L", "M", "S")
)

cols <- c(
  "C" = "#000000",  # black (clean, neutral)
  "L" = "#D55E00",  # orange
  "M" = "#0072B2",  # blue
  "S" = "#009E73"   # green
)

plot(
  merged_df$transformed_mean_multiply_marked_total_egg_count,
  merged_df$recombination_fractions,
  pch = 16,
  col = cols[merged_df$People],
  xlab = "Progeny number (transformed egg count)",
  ylab = "Recombination fraction",
  main = "Recombination increases with progeny number"
)

abline(
  lm(recombination_fractions ~ transformed_mean_multiply_marked_total_egg_count,
     data = merged_df),
  lwd = 2
)

legend(
  "topleft",
  legend = names(cols),
  col = cols,
  pch = 16,
  title = "Collector"
)


