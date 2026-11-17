# ============================================================
# Driver 2: Downstream analyses (Run after cluster finishes)
# ============================================================

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3")

# Recombination frequencies
source("../../Code/3/c_RecombinationFrequency.R")

# Interference
source("../../Code/3/d_Interference.R")

# Figures
source("../../Code/3/e_visualization.R")
source("../../Code/3/f_visualizafion_RF.R")

# Confidence intervals
source("../../Code/3/g_RR_Intervals.R")

# MLE vs Kosambi comparison
source("../../Code/3/h_MLE-Kossambi_comparison.R")

# Fecundity, fertility and viability
source("../../Code/3/i_plot_fecundity-fertility-viabilty.R")

# Marker viability
source("../../Code/3/j_marker present viabilty.R")

# RR survivorship
source("../../Code/3/Kb_RR_SURVIVAL&FECUNDITY.R")

# Body mass
source("../../Code/3/l_body mass.R")

cat("\nPipeline complete.\n")
