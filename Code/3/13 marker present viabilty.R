setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")


#Prepare marker-present life history---------

multiply_marked_viability_dataset<- read.csv("6_multiply_marked_viability_dataset.csv")

MOTHERS_PER_VIAL <- 3

marker_present_life_history <- multiply_marked_viability_dataset %>%
  mutate(
    Dataset   = "Marker-present",
    Fecundity = transformed_mean_multiply_marked_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count / MOTHERS_PER_VIAL,
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

png("marker_present_fertility_fecundity_emmeans.png", width = 10, height = 6, units = "in", res = 300)
plot(
  fec_df$emmean,
  fert_df$emmean,
  xlab = "Fecundity (estimated marginal mean)",
  ylab = "Fertility (estimated marginal mean)",
  pch = 19
)
dev.off()


#Prepare marker-free life history---------

marker_free_viability_dataset <- read.csv("F_marker_free_viability_dataset.csv")

MOTHERS_PER_VIAL <- 3

marker_free_life_history <- marker_free_viability_dataset %>%
  mutate(
    Dataset   = "Marker-free",
    Fecundity = transformed_mean_marker_free_total_egg_count / MOTHERS_PER_VIAL,
    Fertility = total_adult_count / MOTHERS_PER_VIAL,
    Viability = pmin(pmax(derived_marker_free_total_viability * 100, 0), 100)
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
  ) %>%
  mutate(
    Genotype = as.character(Genotype),
    Treatment = as.character(Treatment),
    Brood = as.character(Brood),
    Data_collector = as.character(Data_collector)
  )


# 2. Run ANOVA for Viability

viability_anova <- aov(
  Viability ~ Brood * Treatment * Genotype * Data_collector,
  data = marker_free_life_history
)

# Tidy ANOVA summary with significance stars
viability_anova_summary <- tidy(viability_anova) %>%
  mutate(
    significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      p.value < 0.1   ~ ".",
      TRUE ~ ""
    )
  )

# Save ANOVA summary
write.csv(
  viability_anova_summary,
  "/marker_free_viability_anova_summary.csv",
  row.names = FALSE
)


# 3. Post hoc Tukey

tukey_results <- TukeyHSD(viability_anova)
print(tukey_results)


# 4. Estimated Marginal Means + CLD

emm_viability <- emmeans(viability_anova, ~ Treatment * Genotype * Data_collector)

# CLD letters
cld_viability <- multcomp::cld(emm_viability, Letters = letters)

plot_data <- cld_viability %>%
  mutate(
    Genotype = factor(Genotype),
    Data_collector = factor(Data_collector),
    Treatment = factor(Treatment),
    .group = trimws(.group)  # remove extra spaces
  )


# 5. Plot Viability with CLD

p_viability <- ggplot(plot_data, aes(x = Data_collector, y = emmean, color = Genotype, group = Treatment)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.2, position = position_dodge(0.4)) +
  geom_text(aes(label = .group),
            vjust = -0.8,
            position = position_dodge(0.4),
            size = 5) +
  labs(
    x = "Data Collector",
    y = "Viability (%)",
    color = "Genotype"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave(
  "../../Summaries and Outputs/3/marker_free_viability_cld.png",
  plot = p_viability,
  width = 10,
  height = 6,
  dpi = 300
)
