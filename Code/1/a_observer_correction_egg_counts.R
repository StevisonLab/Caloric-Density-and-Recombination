#Need to adjust lines 30 and 31 and 121 and 122 of script #b for observers in THIS dataset

##pick connor1 and lauren2 because they are more closest to the time of experiment

library(lmodel2)
library(ggplot2)

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3") #use whatever working directory you like

pac_comparison <- read.csv("pac_averages 1.csv")
#pac_comparison <- subset(pac_comparison, strain != "BDSC 6326")


#lauren1 <- lmodel2(PAC_lm ~ microscope_count, data=pac_comparison)
#lauren1

lauren2 <- lmodel2(PAC_lm2 ~ microscope_count, data=pac_comparison)
lauren2

Connor1 <- lmodel2(PAC_ch ~ microscope_count, data=pac_comparison)
Connor1

#Connor2 <- lmodel2(PAC_ch2 ~ microscope_count, data=pac_comparison)
#Connor2

laney<-lmodel2(PAC_lh ~ microscope_count, data=pac_comparison)
laney