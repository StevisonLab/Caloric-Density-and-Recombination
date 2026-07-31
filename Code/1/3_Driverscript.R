# ============================================================
# Driver 2: Downstream analyses (Run after cluster finishes)
# ============================================================

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/2")

# Figures
source("../../Code/3/e_visualization.R")

# Fecundity, fertility and viability
source("../../Code/3/i_plot_fecundity-fertility-viabilty.R")

# RR survivorship
source("../../Code/3/k_RR-SURVIVERSHIP.R")

cat("\nPipeline complete.\n")
