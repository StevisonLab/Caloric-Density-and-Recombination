# combine_all_outputs.R
output_folder <- "/home/mzg0146/Stevison_shared/Diet_Plasticity/New/Outputs/"
key_path <- "/home/mzg0146/Stevison_shared/Diet_Plasticity/New/Key.csv"

# 1. Load the Key
key <- read.csv(key_path, stringsAsFactors = FALSE)

# 2. Get all individual output files
csv_files <- list.files(path = output_folder,
                        pattern = "^Nutrition_countingmodel_outputs_.*\\.csv$",
                        full.names = TRUE)

# 3. Combine them into one data frame
df_list <- lapply(csv_files, read.csv)
all_outputs <- do.call(rbind, df_list)

# 4. Merge with the Key to get the original_ID
# We assume the first column of your output is the task_id
colnames(all_outputs)[1] <- "task_id"
final_combined <- merge(key, all_outputs, by = "task_id")

# 5. (Optional) Remove the task_id column and keep only original_ID
# This makes original_ID the new first column
final_combined$task_id <- NULL

# 6. Save the final result
write.csv(final_combined, 
          file.path("/home/mzg0146/Stevison_shared/Diet_Plasticity/New/", "Nutrition_countingmodel_outputs_all.csv"), 
          row.names = FALSE)

cat("Combined", length(csv_files), "files and mapped original IDs successfully.\n")
