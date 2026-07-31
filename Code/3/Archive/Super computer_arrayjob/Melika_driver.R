setwd("/home/mzg0146/Stevison_shared/Diet_Plasticity/New/")

# Load packages
list.of.packages <- c("doParallel", "foreach", "iterators", "dfoptim", "parallel")
lapply(list.of.packages, library, character.only = TRUE)

# Load dataset
full_experiment_dataset <- read.csv("13_multi_locus_individual_vials_dataset.csv", stringsAsFactors = FALSE)
#full_experiment_dataset <- read.csv("14_multi_locus_pooled_by_collector_cross.csv", stringsAsFactors = FALSE)

# Get SLURM array task ID
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (is.na(task_id)) {
  stop("SLURM_ARRAY_TASK_ID is not set")
}

if (task_id > nrow(full_experiment_dataset)) {
  stop(paste("SLURM_ARRAY_TASK_ID", task_id, 
             "exceeds number of rows in dataset", nrow(full_experiment_dataset)))
}

# Making the single experiment
single_experimental_unit <- full_experiment_dataset[task_id, ]

# Call the source script
source("Nutrition_counting_model_v3.R")

# The output already has task_id in the ID column from source script
# Just save it
output_file <- paste0("Nutrition_countingmodel_outputs_", task_id, ".csv")
write.csv(MLE_Uniform_CxCoM, 
          paste0("/home/mzg0146/Stevison_shared/Diet_Plasticity/New/Outputs/", output_file), 
          row.names = FALSE)

cat("Task", task_id, "completed\n")
