# Load required libraries
library(emmeans)
library(multcomp)
library(ggplot2)
library(dplyr)

# Run ANOVA (assuming multiply_marked_viability_long is already created)
multiply_marked_viability_anova <- aov(Viability ~ Treatment * Genotype * Sex * Brood * Data_Collector, 
                                       data = multiply_marked_viability_long)

# Display ANOVA summary
cat("ANOVA Summary:\n")
print(summary(multiply_marked_viability_anova))

# ============================================================================
# COMPACT LETTER DISPLAY (CLD) FOR SIGNIFICANT EFFECTS
# ===========================================================================

# 3. Data_Collector effect (highly significant: p < 2e-16)
cat("\n\n========================================\n")
cat("3. Data_Collector Main Effect\n")
cat("========================================\n")
emm_collector <- emmeans(multiply_marked_viability_anova, ~ Data_Collector)
print(emm_collector)
cld_collector <- cld(emm_collector, Letters = letters, adjust = "tukey")
print(cld_collector)

# 4. Treatment × Data_Collector interaction (significant: p = 0.016251)
cat("\n\n========================================\n")
cat("4. Treatment × Data_Collector Interaction\n")
cat("========================================\n")
emm_trt_collector <- emmeans(multiply_marked_viability_anova, ~ Treatment * Data_Collector * Genotype)
print(emm_trt_collector)
cld_trt_collector <- cld(emm_trt_collector, Letters = letters, adjust = "tukey")
print(cld_trt_collector)

# 5. Sex × Data_Collector interaction (highly significant: p = 0.000709)
cat("\n\n========================================\n")
cat("5. Sex × Data_Collector Interaction\n")
cat("========================================\n")
emm_sex_collector <- emmeans(multiply_marked_viability_anova, ~ Sex * Data_Collector)
print(emm_sex_collector)
cld_sex_collector <- cld(emm_sex_collector, Letters = letters, adjust = "tukey")
print(cld_sex_collector)


# ============================================================================
# VISUALIZATIONS
# ============================================================================


# Plot 3: Data_Collector Effect
p3 <- ggplot(as.data.frame(cld_collector), aes(x = Data_Collector, y = emmean)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = .group, y = upper.CL), vjust = -0.5, size = 5) +
  labs(title = "Data Collector Effect on Viability",
       subtitle = "Letters indicate significant differences (Tukey HSD, α=0.05)",
       x = "Data Collector",
       y = "Estimated Marginal Mean Viability") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p3)
ggsave("../../Summaries and Outputs/Plot3_Data_Collector_Effect.png", p3, width = 8, height = 6, dpi = 300)

# Plot 4: Treatment × Data_Collector Interaction
p4 <- ggplot(as.data.frame(cld_trt_collector), aes(x = Treatment, y = emmean, color = Data_Collector, group = Genotype )) +
  geom_point(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label = .group, y = upper.CL), 
            position = position_dodge(0.9), vjust = -0.5, size = 3) +
  labs(title = "Treatment × Data Collector Interaction",
       x = "Treatment",
       y = "Estimated Marginal Mean Viability",
       fill = "Data Collector") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

print(p4)
ggsave("../../Summaries and Outputs/Plot4_Treatment_DataCollector_Interaction.png", p4, width = 12, height = 6, dpi = 300)

# Plot 5: Sex × Data_Collector Interaction
p5 <- ggplot(as.data.frame(cld_sex_collector), aes(x = Sex, y = emmean, fill = Data_Collector)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), alpha = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                position = position_dodge(0.9), width = 0.25) +
  geom_text(aes(label = .group, y = upper.CL), 
            position = position_dodge(0.9), vjust = -0.5, size = 3.5) +
  labs(title = "Sex × Data Collector Interaction",
       subtitle = "Letters indicate significant differences (Tukey HSD, α=0.05)",
       x = "Sex",
       y = "Estimated Marginal Mean Viability",
       fill = "Data Collector") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "right")

print(p5)
ggsave("../../Summaries and Outputs/Plot5_Sex_DataCollector_Interaction.png", p5, width = 10, height = 6, dpi = 300)


# Save all CLD results
write.csv(as.data.frame(cld_treatment), 
          "../../Summaries and Outputs/CLD_Treatment.csv", row.names = FALSE)
