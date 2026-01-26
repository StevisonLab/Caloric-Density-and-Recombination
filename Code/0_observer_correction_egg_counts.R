#Need to adjust lines 18 and 19 of script #2 for observers in THIS dataset

#Read in Connor and Lauren egg counts from Protocol study
#compare to microscope counts
#Calculate adjustment
#edit numbers in script as objects!
#add to driver script a source line for this script at the top!

setwd("~/GitHub/Caloric-Density-and-Recombination/Raw Data/3") #use whatever working directory you like
pac_comparison <- read.csv("pac_averages 1.csv")
library(lmodel2) #make sure to install.packages("lmodel2")
lauren <- lmodel2(PAC_lm ~ microscope_count, data=pac_comparison)
lauren

Connor <- lmodel2(PAC_ch ~ microscope_count, data=pac_comparison)
Connor
#Always use the 2 MA regression results:
#lauren's y=mx+b is 
#y= Number of point and clicks
#x= real number of the vial (macroscope)


#Lauren=> Y=0.83X-4.8
#Connor=> Y=0.92X-6.4
