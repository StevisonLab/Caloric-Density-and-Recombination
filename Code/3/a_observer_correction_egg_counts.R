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
#Always use the 2 MA regression results:
#lauren's y=mx+b is 
#y= Number of point and clicks
#x= real number of the vial (macroscope)

#Lauren=> Y=0.83X-4.8
#Connor=> Y=0.92X-6.4

# --- Extract MA slope and intercept for each observer ---
get_ma <- function(model) model$regression.results[model$regression.results$Method == "MA", c("Slope", "Intercept")]

#L1 <- get_ma(lauren1);
L2 <- get_ma(lauren2)
C1 <- get_ma(Connor1);
#C2 <- get_ma(Connor2)
La <- get_ma(laney)

# --- Build predicted lines over the range of real counts ---
x_range <- range(pac_comparison$microscope_count, na.rm = TRUE)
xs <- seq(x_range[1], x_range[2], length.out = 100)

make_line <- function(coefs, observer) {
  data.frame(
    x         = xs,
    y         = coefs$Slope * xs + coefs$Intercept,
    observer  = observer
  )
}

lines_df <- rbind(
#  make_line(L1, "Lauren1"),
  make_line(L2, "Lauren"),
  make_line(C1, "Connor"),
 # make_line(C2, "Connor2"),
  make_line(La, "Laney")
)

# --- Plot ---
a = ggplot() +
  # True 1:1 line
  geom_abline(slope = 1, intercept = 0,
              color = "black", linewidth = 1, linetype = "solid") +
  annotate("text", x = x_range[2] * 0.95, y = x_range[2],
           label = "True count (1:1)", hjust = 1, size = 4) +
  
  # Observer regression lines
  geom_line(data = lines_df,
            aes(x = x, y = y, color = observer),
            linewidth = 1) +
  
  # Raw observed points (optional — shows scatter around each line)
  geom_point(data = pac_comparison,
             aes(x = microscope_count, y = PAC_lm),  color = "#3266ad", size = 1.5) +
  geom_point(data = pac_comparison,
             aes(x = microscope_count, y = PAC_lm2), color = "#5ba3d9", size = 1.5) +
  geom_point(data = pac_comparison,
             aes(x = microscope_count, y = PAC_ch),  color = "#c0392b", size = 1.5) +
  geom_point(data = pac_comparison,
             aes(x = microscope_count, y = PAC_ch2), color = "#e67e22", size = 1.5) +
  geom_point(data = pac_comparison,
             aes(x = microscope_count, y = PAC_lh),  color = "#27ae60", size = 1.5) +
  
  scale_color_manual(values = c(
   # "Lauren" = "#3266ad",
    "Lauren" = "#5ba3d9",
    "Connor" = "#c0392b", 
    #"Connor after" = "#e67e22",
    "Laney"    = "#27ae60"
  )) +

  labs(
    title    = "Observer egg counts vs. true microscope count",
    x        = "Microscope count (true)",
    y        = "Observer count (predicted by MA regression)",
    color    = "Observer", linetype = "Observer"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

ggsave(
  filename = "../../Summaries and Outputs/3/1_Observerslope.png",
  plot = a,
  width = 6,
  height = 6,
  dpi = 300
)



