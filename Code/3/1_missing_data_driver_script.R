setwd("Raw Data/") # Please use Raw Data as working directory
setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3")

# Compute the egg correction
source(file = "../../Code/0_observer_correction_egg_counts.R")

# Compile, format, transform all raw data into appropriate input files
source(file = "../../Code/2_missing_data_dataset_compilation_script.R")


# Perform organismal viability analyses
source(file = "../../Code/9_.R")

# Perform multi-locus viability analyses
# This script will take a LONG time!
source(file = "../../Code/7_missing_data_multi_locus_analysis.R")

# Perform alternate pooling viability analyses
# This script will also take a LONG time!
#wiat for now!
source(file = "../../Code/8_missing_data_pooling_analysis.R")

