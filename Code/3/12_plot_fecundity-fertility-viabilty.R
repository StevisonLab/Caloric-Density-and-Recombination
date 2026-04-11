# ============================================================
# SET WORKING DIRECTORY
# ============================================================
setwd("~/GitHub/Caloric-Density-and-Recombination/Summaries and Outputs/3")

# ============================================================
# READ DATA
# ============================================================
marker_free_viability_dataset <- read.csv("F_marker_free_viability_dataset.csv")
multiply_marked_viability_dataset <- read.csv("6_multiply_marked_viability_dataset.csv")

# mothers per vial for marker-present
mothers_per_vial <- read.csv("3_mom.csv")

multiply_marked_viability_dataset <- multiply_marked_viability_dataset %>%
  left_join(mothers_per_vial, by = c("cross_id", "brood_id"))

# ============================================================
# LIBRARIES & CONSTANTS
# ============================================================
library(dplyr)
library(ggplot2)
library(broom)

MOTHERS_PER_VIAL <- 3
traits <- c("Fecundity", "Fertility", "Viability")

# ============================================================
# MARKER-FREE DATA
# ============================================================
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
    vialID      = cross_id,   # ✅ ADD THIS LINE
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

# ============================================================
# MARKER-PRESENT DATA
# ============================================================
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
    vialID      = cross_id,   
    Dataset,
    Fecundity,
    Fertility,
    Viability
  ) %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Brood     = factor(Brood),
    Dataset   = factor(Dataset),
    Data_collector = factor(Data_collector)
  )

# ============================================================
# COMBINE DATASETS
# ============================================================
life_history_all <- bind_rows(
  marker_free_life_history,
  marker_present_life_history
) %>%
  mutate(
    Treatment = factor(Treatment),
    Genotype  = factor(Genotype),
    Brood     = factor(Brood),
    Dataset   = factor(Dataset),
    Data_collector = factor(Data_collector)
  )

# ============================================================
# TRAIT PLOTS (UNCHANGED)
# ============================================================
plot_trait <- function(trait_name, y_label) {
  ggplot(
    life_history_all,
    aes(x = Treatment, y = .data[[trait_name]], fill = Genotype)
  ) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA,
                 position = position_dodge(width = 0.8)) +
    facet_wrap(~ Dataset) +
    theme_classic(base_size = 14) +
    labs(
      title = trait_name,
      y = y_label,
      x = "Treatment"
    ) +
    scale_fill_manual(values = c(
   #   "42"  = "#ff7f00",
  #    "399" = "#984ea3"
      "399" = "#4a1486",  # purple
      "42"  = "#99000d"   # orange
    ))
}
#reaction norm-----------
get_significance <- function(anova_table) {
  
  # Trim whitespace from term names
  anova_table$Term <- trimws(anova_table$Term)
  
  get_p <- function(term) {
    row <- anova_table[anova_table$Term == term, ]
    if (nrow(row) == 0) return(NA)
    return(row$`Pr(>F)`)
  }
  
  list(
    sig_E   = isTRUE(get_p("Treatment") < 0.05),              # Environment
    sig_G   = isTRUE(get_p("Genotype") < 0.05),
    sig_B   = isTRUE(get_p("Brood") < 0.05),
    sig_GxE = isTRUE(get_p("Treatment:Genotype") < 0.05)
  )
}
# Format terms (bold if significant) — no quotes inside bold()
fmt_term <- function(term, significant) {
  if (isTRUE(significant)) {
    paste0("bold(", term, ")")
  } else {
    term
  }
}

# Build equation label: plain trait name + tilde + bold significant terms
build_eq_label <- function(sig_list, trait_name) {
  paste0(
    "paste('", trait_name, "', ' ~ ', ",
    if (sig_list$sig_G) "bold(G)" else "G", " + ",
    if (sig_list$sig_E) "bold(E)" else "E", " + ",
    if (sig_list$sig_B) "bold(B)" else "B", " + ",
    if (sig_list$sig_GxE) "bold(G %*% E)" else "G %*% E",
    ")"
  )
}
# ============================================================
# RUN ANOVA + STORE LABELS
# ============================================================

anova_labels <- list()

for (dataset_name in unique(as.character(life_history_all$Dataset))) {
  
  df_subset <- life_history_all %>%
    filter(Dataset == dataset_name)
  
  for (trait in traits) {
    
    df_trait <- df_subset %>%
      filter(is.finite(.data[[trait]]))
    
    if (nrow(df_trait) > 1) {
      
      model <- aov(
        as.formula(paste(trait, "~ Treatment * Genotype * Brood + Data_collector")),
        data = df_trait
      )
      
      anova_table <- as.data.frame(summary(model)[[1]])
      anova_table$Term <- rownames(anova_table)
      rownames(anova_table) <- NULL
      
      sig_list <- get_significance(anova_table)
      
      eq_label <- build_eq_label(sig_list, trait)
      
      anova_labels[[paste(dataset_name, trait, sep = "_")]] <- eq_label
    }
  }
}

# ============================================================
# UPDATED PLOT FUNCTION
# ============================================================

plot_trait <- function(trait_name, y_label) {
  
  label_df <- data.frame(
    Dataset = unique(life_history_all$Dataset),
    x = 2,
    y = 95,
    label = c(
      anova_labels[[paste("Marker-free", trait_name, sep = "_")]],
      anova_labels[[paste("Marker-present", trait_name, sep = "_")]]
    )
  )
  
  ggplot(life_history_all,
         aes(
           x = Treatment,
           y = .data[[trait_name]],
           group = Genotype,
           color = Genotype
         )) +
    
    stat_summary(fun = mean, geom = "line", linewidth = 3) +
    stat_summary(fun = mean, geom = "point", size = 6) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.55) +
    
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      parse = TRUE,
      inherit.aes = FALSE,
      size = 8
    ) +
    
    ylim(0, 100) +
    facet_wrap(~ Dataset) +
    
    theme_classic(base_size = 24) +
    theme(
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    ) +
    
    labs(
      title = trait_name,
      y = y_label,
      x = "Treatment"
    ) +
    
    scale_color_manual(values = c(
      "217" = "#e41a1c",
      "391" = "#377eb8",
      "491" = "#4daf4a",
      "399" = "#984ea3",
      "42"  = "#ff7f00"
    ))
}

# ============================================================
# GENERATE PLOTS
# ============================================================

p_fecundity <- plot_trait("Fecundity", "Eggs per mother")
p_fertility <- plot_trait("Fertility", "Viable offspring")
p_viability <- plot_trait("Viability", "Viability (%)")

ggsave("../../Summaries and Outputs/3/20_trait_fecundity.png",
       plot = p_fecundity, width = 12, height = 8, dpi = 300)

ggsave("../../Summaries and Outputs/3/21_trait_fertility.png",
       plot = p_fertility, width = 12, height = 8, dpi = 300)

ggsave("../../Summaries and Outputs/3/22_trait_viability.png",
       plot = p_viability, width = 12, height = 8, dpi = 300)

# ============================================================
# ANOVA — SEPARATE BY DATASET AND TRAIT (NO PLOTTING)
# ============================================================
run_anova <- function(df, trait) {
  aov(
    as.formula(paste(trait, "~ Treatment * Genotype * Brood + Data_collector")),
    data = df
  )
}

for (dataset_name in unique(as.character(life_history_all$Dataset))) {
  
  df_subset <- life_history_all %>%
    filter(Dataset == dataset_name)
  
  for (trait in traits) {
    
    # remove NA/NaN/Inf values per trait
    df_trait <- df_subset %>%
      filter(is.finite(.data[[trait]]))
    
    # only run ANOVA if there is enough data
    if (nrow(df_trait) > 1) {
      model <- run_anova(df_trait, trait)
      
      anova_table <- as.data.frame(summary(model)[[1]])
      anova_table$Term <- rownames(anova_table)
      anova_table$Dataset <- dataset_name
      anova_table$Trait <- trait
      
      rownames(anova_table) <- NULL
      
      write.csv(
        anova_table,
        paste0(
          "../../Summaries and Outputs/3/23_ANOVA_",
          dataset_name, "_",
          trait,
          ".csv"
        ),
        row.names = FALSE
      )
    }
  }
}

# ============================================================
# SAVE CLEANED DATASET
# ============================================================
write.csv(
  life_history_all,
  "../../Summaries and Outputs/3/24_Life_History_All_Cleaned.csv",
  row.names = FALSE
)

# ============================================================
# SUMMARY MEANS BY DATASET × TREATMENT × GENOTYPE
# ============================================================
life_history_summary <- life_history_all %>%
  group_by(Dataset, Treatment, Genotype) %>%
  summarise(
    Fecundity = mean(Fecundity, na.rm = TRUE),
    Fertility = mean(Fertility, na.rm = TRUE),
    Viability = mean(Viability, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  life_history_summary,
  "../../Summaries and Outputs/3/25_life_history_summary_by_dataset.csv",
  row.names = FALSE
)
