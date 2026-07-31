# ============================================================
# Driver 2: Downstream analyses (Run after cluster finishes)
# ============================================================

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3")

# Recombination frequencies
source("../../Code/c_RecombinationFrequency.R")

# Interference
source("../../Code/d_Interference.R")

# Figures
source("../../Code/e_visualization.R")
source("../../Code/f_visualizafion_RF.R")

# Confidence intervals
source("../../Code/g_RR_Intervals.R")

# MLE vs Kosambi comparison
source("../../Code/h_MLE-Kossambi_comparison.R")

# Fecundity, fertility and viability
source("../../Code/i_plot_fecundity-fertility-viability.R")

# Marker viability
source("../../Code/j_marker present viability.R")

# RR survivorship
source("../../Code/k_RR-SURVIVERSHIP.R")

# Body mass
source("../../Code/l_body mass.R")

cat("\nPipeline complete.\n")