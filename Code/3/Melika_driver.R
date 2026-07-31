# Fetch settings from environment variables passed by Bash
work_dir   <- Sys.getenv("WORK_DIR")
input_file <- Sys.getenv("INPUT_FILE")

# Set working directory
setwd(work_dir)

# Load packages
list.of.packages <- c("doParallel", "foreach", "iterators", "dfoptim", "parallel")
lapply(list.of.packages, library, character.only = TRUE)

output_dir <- file.path(work_dir, "Outputs")

# Ensure output folder exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load dataset dynamically
full_experiment_dataset <- read.csv(input_file, stringsAsFactors = FALSE)

# Get SLURM array task ID
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) { 
  stop("SLURM_ARRAY_TASK_ID is not set") 
} 
if (task_id > nrow(full_experiment_dataset)) { 
  stop(paste("SLURM_ARRAY_TASK_ID", task_id, "exceeds number of rows in dataset", nrow(full_experiment_dataset))) 
}

# Process single experiment
single_experimental_unit <- full_experiment_dataset[task_id, ]

# Run source model script
source("Nutrition_counting_model_v3.R")

# Save task output
output_file <- paste0("Nutrition_countingmodel_outputs_", task_id, ".csv")
write.csv(MLE_Uniform_CxCoM, file.path(output_dir, output_file), row.names = FALSE)

cat("Task", task_id, "completed successfully for row", task_id, "in", input_file, "\n")