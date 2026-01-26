#setwd("Raw Data/") # Please use Raw Data as working directory

# Perform Regression and write ANOVA table for supplemental material
marker_free_viability_anova.lm <- lm(Viability ~ Brood + Sex + Treatment + Genotype, data=marker_free_viability_anova_input)
marker_free_viability_anova.table <- anova(marker_free_viability_anova.lm)
marker_free_viability_regression_coefficients <- summary(marker_free_viability_anova.lm)
write.csv(marker_free_viability_anova.table, "../../Summaries and Outputs/table_S1_marker_free_viability_regression.csv", row.names = TRUE)

# Perform Regression and write ANOVA table for supplemental material
multiply_marked_viability_anova.lm <- lm(Viability ~ Brood + Sex + Treatment + Genotype, data=cleaned_multiply_marked_viability_anova_input)
multiply_marked_viability_anova.table <- anova(multiply_marked_viability_anova.lm)
multiply_marked_viability_regression_coefficients <- summary(multiply_marked_viability_anova.lm)
write.csv(multiply_marked_viability_anova.table, "../../Summaries and Outputs/table_S2_multiply_marked_viability_regression.csv", row.names = TRUE)

# With interactions between Treatment and other variables
marker_free_viability_anova.lm <- lm(Viability ~ Brood + Sex + Treatment + Genotype + Treatment:Sex + Treatment:Genotype, data=marker_free_viability_anova_input)

# Or with all two-way interactions
marker_free_viability_anova.lm <- lm(Viability ~ (Brood + Sex + Cross + Treatment + Genotype)^2, data=marker_free_viability_anova_input)
