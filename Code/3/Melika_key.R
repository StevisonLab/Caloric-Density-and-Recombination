# Read environment variables
work_dir   <- Sys.getenv("WORK_DIR")
data_dir   <- Sys.getenv("DATA_DIR")
input_file <- Sys.getenv("INPUT_FILE")

setwd(work_dir)

# Read input CSV from DATA_DIR
input_path <- file.path(data_dir, input_file)
full_experiment_dataset <- read.csv(input_path, stringsAsFactors = FALSE)

# Generate key mapping
Key <- data.frame(
  task_id     = 1:nrow(full_experiment_dataset),
  original_ID = full_experiment_dataset[, 1]
)

# Save Key.csv in WORK_DIR
write.csv(Key, "Key.csv", row.names = FALSE)
cat("Key created successfully for file:", input_path, "\n")
