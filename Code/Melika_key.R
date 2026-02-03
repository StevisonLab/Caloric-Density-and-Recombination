# create_key.R - Run this ONCE before submitting array jobs
setwd("/home/mzg0146/Stevison_shared/Diet_Plasticity/New/")

full_experiment_dataset <- read.csv("13_multi_locus_individual_vials_dataset.csv", 
                                    stringsAsFactors = FALSE)

Key <- data.frame(
  task_id = 1:nrow(full_experiment_dataset),
  original_ID = full_experiment_dataset[, 1]
)

write.csv(Key, "Key.csv", row.names = FALSE)
cat("Key created with", nrow(Key), "entries\n")
