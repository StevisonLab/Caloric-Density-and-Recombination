#setting working directory
setwd("C:/Users/Melika Shiran/Box/Plasticity and Nutrition NEW (Melika)/03. Raw Data and Analysis/2.Body Mass (ZL-NoMarker Piolot)")

#install.packages("emmeans")

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyverse)

# Import the data
data <- read_csv("Body_Mass.csv")

# Remove rows with missing data in important columns
data <- data %>% drop_na(BodyMass)

# Ensure categorical variables are factors
data$Food <- as.factor(data$Food)
data$Genetics <- as.factor(data$Genetics)
data$MatingStatus <- as.factor(data$MatingStatus)
data$grams=data$BodyMass*10^3

# Calculate the mean body mass for each group
summary_stats <- data %>%
  group_by(Genetics, Food, MatingStatus) %>%
  summarise(
    Mean_Body_Mass = mean(grams, na.rm = TRUE),
    SD_Body_Mass = sd(grams, na.rm = TRUE),
    N = n()
  )

# View the results
print(summary_stats)
write.csv(summary_stats,file="Summary_stats.csv",quote=F,row.names = F)


# Perform ANOVA on the raw data
##anova_result <- aov(BodyMass ~ Genetics * Food + MatingStatus, data = data)
anova_result <- aov(grams ~ Genetics * Food, data = data)

# View the summary of the ANOVA
summary(anova_result)
library(car)
output=Anova(anova_result)

#write anova table to CSV
write.csv(output,file="Anova_table.csv",row.names = T,quote=F)

# Load necessary libraries
library(emmeans)
library(ggplot2)

#Calculate estimated marginal means (EMMs)
em_means <- emmeans(anova_result, specs = pairwise ~ Food, by=c("Genetics"))
cont_emm=contrast(em_means,method="pairwise")

#only do posthoc on significant model parameters
em_means <- emmeans(anova_result, specs = pairwise ~ Food)
cont_emm=contrast(em_means,method="pairwise")

# View the EMMs results
summary(cont_emm)

#Create a basic plot using emmeans' plot function
plot(em_means, comparisons = TRUE, horizontal = FALSE)

ggplot(data,aes(x=Food,y=grams)) + geom_boxplot()

##posthoc Letters
# Load required libraries
library(multcompView)
library(ggplot2)

# Perform Tukey post hoc test
posthoc <- TukeyHSD(anova_result, "Food")

# Extract significance letters
letters_df <- multcompLetters(posthoc$Food[, "p adj"])$Letters
letters_df <- data.frame(Food = names(letters_df), Letters = letters_df)

# Find max values for proper label positioning
boxplot_stats <- data %>%
  group_by(Food) %>%
  summarise(y_pos = max(grams, na.rm = TRUE) + 0.02)  # Adjust as needed

# Merge y positions with letters
letters_df <- merge(letters_df, boxplot_stats, by = "Food")

# Create improved boxplot with correctly positioned letters
ggplot(data, aes(x = Food, y = grams)) + 
  geom_boxplot() +
  geom_text(data = letters_df, aes(x = Food, y = y_pos+1, label = Letters),
            size = 6, fontface = "bold") +  # Increase font size
  labs(title = "Post Hoc Analysis of Body Mass by Food Treatment",
       x = "Food Treatment",
       y = "Body Mass (mg)") +
  theme_minimal()

ggsave("Posthoc_tukeys.png")
