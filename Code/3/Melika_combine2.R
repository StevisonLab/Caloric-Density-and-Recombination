# Fetch settings from environment variables passed by Bash
work_dir          <- Sys.getenv("WORK_DIR")
final_output_name <- Sys.getenv("FINAL_OUTPUT_NAME")

# Set working directory
setwd(work_dir)

output_folder <- file.path(work_dir, "Outputs")
key_path      <- file.path(work_dir, "Key.csv")

# 1. Load the Key
key <- read.csv(key_path, stringsAsFactors = FALSE)

# 2. Get all individual output files
csv_files <- list.files(
  path = output_folder, 
  pattern = "^Nutrition_countingmodel_outputs_.*\\.csv$", 
  full.names = TRUE
)

# Filter out any pre-existing combined file
csv_files <- csv_files[!grep(final_output_name, csv_files)]

# 3. Combine output files into one data frame
df_list <- lapply(csv_files, read.csv)
all_outputs <- do.call(rbind, df_list)

# 4. Merge with the Key to restore original_ID
colnames(all_outputs)[1] <- "task_id"
final_combined <- merge(key, all_outputs, by = "task_id")

# 5. Remove task_id column (keeping original_ID)
final_combined$task_id <- NULL

# 6. Save the final result
write.csv(final_combined, file.path(work_dir, final_output_name), row.names = FALSE)

cat("Successfully combined", length(csv_files), "files into:", final_output_name, "\n")