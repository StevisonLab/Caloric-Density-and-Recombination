# Read environment variables
work_dir   <- Sys.getenv("WORK_DIR")
data_dir   <- Sys.getenv("DATA_DIR")
input_file <- Sys.getenv("INPUT_FILE")

setwd(work_dir)

# Load required packages
list.of.packages <- c("doParallel", "foreach", "iterators", "dfoptim", "parallel")
lapply(list.of.packages, library, character.only = TRUE)

output_dir <- file.path(work_dir, "Outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read input dataset from DATA_DIR
input_path <- file.path(data_dir, input_file)
full_experiment_dataset <- read.csv(input_path, stringsAsFactors = FALSE)

# Get SLURM task ID
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) { 
  stop("SLURM_ARRAY_TASK_ID is not set") 
} 
if (task_id > nrow(full_experiment_dataset)) { 
  stop(paste("SLURM_ARRAY_TASK_ID", task_id, "exceeds rows in dataset", nrow(full_experiment_dataset))) 
}

# Process current row
single_experimental_unit <- full_experiment_dataset[task_id, ]

# Run analysis script
source("Nutrition_counting_model_v3.R")

# Save individual output in WORK_DIR/Outputs
output_file <- paste0("Nutrition_countingmodel_outputs_", task_id, ".csv")
write.csv(MLE_Uniform_CxCoM, file.path(output_dir, output_file), row.names = FALSE)

cat("Task", task_id, "completed successfully.\n")
