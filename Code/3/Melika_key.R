# Fetch settings from environment variables passed by Bash
work_dir   <- Sys.getenv("WORK_DIR")
input_file <- Sys.getenv("INPUT_FILE")

# Set working directory
setwd(work_dir)

# Read input data dynamically
full_experiment_dataset <- read.csv(input_file, stringsAsFactors = FALSE)

Key <- data.frame(
  task_id = 1:nrow(full_experiment_dataset),
  original_ID = full_experiment_dataset[, 1]
)

write.csv(Key, "Key.csv", row.names = FALSE)
cat("Key created with", nrow(Key), "entries for file:", input_file, "\n")