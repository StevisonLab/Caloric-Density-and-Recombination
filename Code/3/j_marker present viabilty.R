setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")


#Prepare marker-present life history---------

multiply_marked_viability_dataset<- read.csv("6_multiply_marked_viability_dataset.csv")
# mothers per vial for marker-present
mothers_per_vial <- read.csv("3_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))


marker_present_life_history <- multiply_marked_viability_dataset %>%
  mutate(
    Dataset   = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / female_in_vial,
    Fertility = total_adult_count / female_in_vial,
    Viability = pmin(pmax(derived_multiply_marked_total_viability * 100, 0), 100)
  ) %>%
  select(
    Data_collector = data_collector,
    Treatment = treatment,
    Genotype  = genotype,
    Brood     = brood_id,
    Dataset,
    Fecundity,
    Fertility,
    Viability
  )  %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Brood     = factor(Brood),
    Dataset   = factor(Dataset),
    Data_collector = factor(Data_collector))



marker_present_life_history <- marker_present_life_history %>%
  filter(is.finite(Fecundity))   # keeps only finite values

# 2. Run ANOVA for Viability

viability_anova <- aov(
  Viability ~ Brood * Treatment * Data_collector * Genotype,
  data = marker_present_life_history
)
summary(viability_anova)

Fecundity_anova <- aov(
  Fecundity ~ Brood * Treatment *Genotype * Data_collector,
  data = marker_present_life_history
)
summary(Fecundity_anova)

Fertility_anova <- aov(
  Fertility ~ Brood * Treatment *Genotype * Data_collector,
  data = marker_present_life_history
)
summary(Fecundity_anova)
# 4. Estimated Marginal Means + CLD

emm_viability <- emmeans(viability_anova, ~  Data_collector)
pairs(emm_viability)

emm_Fecundity <- emmeans(Fecundity_anova , ~  Data_collector)
pairs(emm_Fecundity)

emm_Fertility <- emmeans(Fertility_anova, ~  Data_collector)
pairs(emm_Fertility)

fec_df  <- as.data.frame(emm_Fecundity)
fert_df <- as.data.frame(emm_Fertility)

png("45_marker_present_fertility_fecundity_emmeans.png", width = 10, height = 6, units = "in", res = 300)
plot(
  fec_df$emmean,
  fert_df$emmean,
  xlab = "Fecundity (estimated marginal mean)",
  ylab = "Fertility (estimated marginal mean)",
  pch = 19
)
dev.off()

