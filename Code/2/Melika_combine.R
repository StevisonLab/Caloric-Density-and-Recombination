# Read environment variables
work_dir          <- Sys.getenv("WORK_DIR")
data_dir          <- Sys.getenv("DATA_DIR")
final_output_name <- Sys.getenv("FINAL_OUTPUT_NAME")

setwd(work_dir)

output_folder <- file.path(work_dir, "Outputs")
key_path      <- file.path(work_dir, "Key.csv")

# 1. Load Key
key <- read.csv(key_path, stringsAsFactors = FALSE)

# 2. Find all individual task output CSVs
csv_files <- list.files(
  path       = output_folder, 
  pattern    = "^Nutrition_countingmodel_outputs_.*\\.csv$", 
  full.names = TRUE
)

# 3. Read and combine all outputs
df_list     <- lapply(csv_files, read.csv)
all_outputs <- do.call(rbind, df_list)

# 4. Merge outputs with Key to preserve original IDs
colnames(all_outputs)[1] <- "task_id"
final_combined           <- merge(key, all_outputs, by = "task_id")
final_combined$task_id   <- NULL

# 5. Save final result into DATA_DIR ("Summaries and Outputs/3")
final_output_path <- file.path(data_dir, final_output_name)
write.csv(final_combined, final_output_path, row.names = FALSE)

cat("Successfully combined", length(csv_files), "files.\n")
cat("Final merged file written to:", final_output_path, "\n")
