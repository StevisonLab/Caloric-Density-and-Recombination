# ============================================================
# Driver 1: Data preparation (Run locally)
# ============================================================

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3")

# Step 0. Egg count observer correction
source("../../Code/3/a_observer_correction_egg_counts.R")

# Step 1. Compile all datasets
source("../../Code/3/b_diet_dataset_compilation_script.R")

cat("\nData preparation complete.\n")
cat("Next step: submit 2_Driverscript.sh to the supercomputer.\n")
