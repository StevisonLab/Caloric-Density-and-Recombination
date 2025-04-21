#adapting Melika data to counting model code

library(doBy)
library(tidyr)
library(dplyr)


merged_data$Treatment=merged_data$`F1 developmental`
merged_data$Stock=merged_data$`F1 Genotype`

summary_data <- summaryBy(
  Count ~ VialID + Treatment + Stock +
    Scute + Crossveinless + Vermilion + Fork + Carnation + Sex,
  data = merged_data,
  FUN = sum,
  keep.names = TRUE
)

summary_data$Scute=ifelse(summary_data$Scute==TRUE,1,0)
summary_data$Crossveinless=ifelse(summary_data$Crossveinless==TRUE,1,0)
summary_data$Vermilion=ifelse(summary_data$Vermilion==TRUE,1,0)
summary_data$Fork=ifelse(summary_data$Fork==TRUE,1,0)
summary_data$Carnation=ifelse(summary_data$Carnation==TRUE,1,0)

#New Marker Class Column
summary_data$Class <- paste0(
  summary_data$Carnation,
  summary_data$Fork,
  summary_data$Vermilion,
  summary_data$Crossveinless,
  summary_data$Scute
)

#Need to also add in the number of egg count per vial that did NOT survive to adulthood as the dead class here!


# Step 1: Separate rows with NA in Sex
sex_known <- summary_data %>%
  filter(!is.na(Sex))

sex_unknown <- summary_data %>%
  filter(is.na(Sex)) %>%
  mutate(Class = "dead")

# Step 2: Pivot the known-sex data
pivoted_known <- sex_known %>%
  pivot_wider(
    names_from = Sex,
    values_from = Count,
    values_fill = 0
  )

# Step 3: Combine back the two sets
pivoted_data <- bind_rows(pivoted_known, sex_unknown)

set.seed(123)  # for reproducibility

# Loop through rows with Markers == "dead"
for (i in which(pivoted_data$Class == "dead")) {
 # total <- pivoted_data$M[i] + pivoted_data$F[i] + pivoted_data$Count[i]
  #old_male=pivoted_data$M[i]
  #old_female=pivoted_data$F[i]

  # Skip if there's no total to split
 # if (total == 0) next
  
  half <- floor(pivoted_data$Count[i] / 2)
  remainder <- pivoted_data$Count[i] %% 2
  
  # Randomly assign the extra 1 (if odd total)
  if (remainder == 1) {
    if (sample(c(TRUE, FALSE), 1)) {
      male_new <- half + 1
      female_new <- half
    } else {
      male_new <- half
      female_new <- half + 1
    }
  } else {
    male_new <- female_new <- half
  }
  
  pivoted_data$M[i] <- male_new
  pivoted_data$F[i] <- female_new
}

pivoted_data$Females=pivoted_data$F
pivoted_data$Males=pivoted_data$M

vial_IDs=unique(sort(pivoted_data$VialID))
my_df <- data.frame(
  VialID = vial_IDs,
  genetic_map_length = NA_real_  # numeric NA, assuming this will hold numbers
)

for (vial in 1:length(vial_IDs)) {
  v=vial_IDs[vial]
  temp=pivoted_data[pivoted_data$VialID==v,]
  temp$Count=temp$F+temp$M
  observed_dataset=temp[,c("Class","Carnation","Fork","Vermilion","Crossveinless","Scute","Females","Males","Count")]
  source("counting_model.R")
  tml=df2[1,9]+df2[1,10]+df2[1,11]+df2[1,12]
  my_df$genetic_map_length[output$VialID==v]=tml
}

