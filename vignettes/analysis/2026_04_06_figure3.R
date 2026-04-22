pacman::p_load(tidyverse)


plot_data <- readRDS("~/2026_fedGee/analysis/2026_04_03_wnar_output.RDS")

# 2. Variable Mapping (Including Reference Groups in Parentheses)
clean_data <- plot_data %>%
  filter(Term != "(Intercept)") %>%
  mutate(
    OR      = exp(Estimate),
    CI_low  = exp(Estimate - 1.96 * SE),
    CI_high = exp(Estimate + 1.96 * SE),
    
    Clean_Term = case_when(
      # Demographics
      str_detect(Term, "age_cat22") ~ "Age (ref: 18-39): 40-64",
      str_detect(Term, "age_cat23") ~ "Age (ref: 18-39): 65+",
      str_detect(Term, "sex1")      ~ "Sex (ref: Male): Female",
      str_detect(Term, "race_cat1") ~ "Race/Ethnicity (ref: White): Am. Indian/Alaska Native",
      str_detect(Term, "race_cat3") ~ "Race/Ethnicity (ref: White): Black or African American",
      str_detect(Term, "race_cat4") ~ "Race/Ethnicity (ref: White): Hispanic",
      str_detect(Term, "rural_cat_ORH2") ~ "Rurality (ref: Urban): Rural",
      str_detect(Term, "adi_cat2")  ~ "ADI (ref: 1-25): 26-50",
      str_detect(Term, "adi_cat3")  ~ "ADI (ref: 1-25): 51-75",
      str_detect(Term, "adi_cat4")  ~ "ADI (ref: 1-25): 76-100",
      
      # Substance Use History
      #str_detect(Term, "prior_AUD_dx1") ~ "Prior Ambulatory AUD Diagnosis (ref: No)",
      str_detect(Term, "remote_MAUD1")  ~ "Any MAUD Use in Prior Year (ref: No)",
      str_detect(Term, "oud_dx1")       ~ "OUD Diagnosis (ref: No)",
      str_detect(Term, "auditc_baseline_cat_d1") ~ "AUDIT-C (ref: Unhealthy): No Use",
      str_detect(Term, "auditc_baseline_cat_d2") ~ "AUDIT-C (ref: Unhealthy): Low Use",
      str_detect(Term, "auditc_baseline_cat_d4") ~ "AUDIT-C (ref: Unhealthy): High Risk Use",
      str_detect(Term, "auditc_baseline_cat_d5") ~ "AUDIT-C (ref: Unhealthy): Very High Risk",
      
      # Medical Complexity
      str_detect(Term, "CAN_quint2") ~ "VA CAN Score (ref: Q1): Quintile 2",
      str_detect(Term, "CAN_quint3") ~ "VA CAN Score (ref: Q1): Quintile 3",
      str_detect(Term, "CAN_quint4") ~ "VA CAN Score (ref: Q1): Quintile 4",
      str_detect(Term, "CAN_quint5") ~ "VA CAN Score (ref: Q1): Quintile 5",
      str_detect(Term, "frail_cat22") ~ "Frailty (ref: Fit): Prefrail",
      str_detect(Term, "frail_cat23") ~ "Frailty (ref: Fit): Mildly Frail",
      str_detect(Term, "frail_cat24") ~ "Frailty (ref: Fit): Moderately Frail",
      str_detect(Term, "frail_cat25") ~ "Frailty (ref: Fit): Severely Frail",
      str_detect(Term, "icu1")       ~ "ICU Care (ref: No)",
      
      # Hospitalization Factors
      str_detect(Term, "speciality_cat1") ~ "Speciality (ref: Med/Surg): Psychiatry", 
      str_detect(Term, "speciality_cat3") ~ "Speciality (ref: Med/Surg): Surgery",
      str_detect(Term, "discharge_cat2")  ~ "Discharge (ref: Routine): Against Medical Advice",
      str_detect(Term, "discharge_cat3")  ~ "Discharge (ref: Routine): Domiciliary",
      str_detect(Term, "hospital_addiction_consult1") ~ "Addiction Consult (ref: None): Clinician",
      str_detect(Term, "hospital_addiction_consult2") ~ "Addiction Consult (ref: None): Non-Clinician",
      TRUE ~ "Exclude"
    ),
    
    Category = case_when(
      str_detect(Clean_Term, "Age|Sex|Race|Ethnicity|Rurality|ADI") ~ "Demographics",
      str_detect(Clean_Term, "AUD|MAUD|OUD|AUDIT-C") ~ "Substance Use History",
      str_detect(Clean_Term, "CAN|Frailty|ICU") ~ "Medical Complexity",
      str_detect(Clean_Term, "Consult|Speciality|Discharge") ~ "Hospitalization\nFactors",
      TRUE ~ "Exclude"
    )
  ) %>%
  filter(Category != "Exclude")

# 3. Ordering logic (Table 1 Blueprint)
ordered_levels <- rev(c(
  "Age (ref: 18-39): 40-64", "Age (ref: 18-39): 65+",
  "Sex (ref: Male): Female",
  "Race/Ethnicity (ref: White): Am. Indian/Alaska Native", "Race/Ethnicity (ref: White): Black or African American", "Race/Ethnicity (ref: White): Hispanic",
  "Rurality (ref: Urban): Rural",
  "ADI (ref: 1-25): 26-50", "ADI (ref: 1-25): 51-75", "ADI (ref: 1-25): 76-100",
  "Prior Ambulatory AUD Diagnosis (ref: No)", "Any MAUD Use in Prior Year (ref: No)", "OUD Diagnosis (ref: No)",
  "AUDIT-C (ref: Unhealthy): No Use", "AUDIT-C (ref: Unhealthy): Low Use", "AUDIT-C (ref: Unhealthy): High Risk Use", "AUDIT-C (ref: Unhealthy): Very High Risk",
  "VA CAN Score (ref: Q1): Quintile 2", "VA CAN Score (ref: Q1): Quintile 3", "VA CAN Score (ref: Q1): Quintile 4", "VA CAN Score (ref: Q1): Quintile 5",
  "Frailty (ref: Fit): Prefrail", "Frailty (ref: Fit): Mildly Frail", "Frailty (ref: Fit): Moderately Frail", "Frailty (ref: Fit): Severely Frail",
  "ICU Care (ref: No)",
  "Speciality (ref: Med/Surg): Psychiatry", "Speciality (ref: Med/Surg): Surgery",
  "Discharge (ref: Routine): Against Medical Advice", "Discharge (ref: Routine): Domiciliary",
  "Addiction Consult (ref: None): Clinician", "Addiction Consult (ref: None): Non-Clinician"
))

clean_data$Clean_Term <- factor(clean_data$Clean_Term, levels = ordered_levels)

# 5. Generate the Faceted Forest Plot
p_forest <- ggplot(clean_data, aes(x = OR, y = Clean_Term, color = Model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha = 0.6) +
  
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high), 
                  position = position_dodge(width = 0.6), 
                  size = 0.5, fatten = 2.5) +
  
  scale_color_manual(values = c("Pooled GEE" = "#0072B2", "FedGEE" = "#D55E00")) +
  
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), 
                labels = c("0.25", "0.50", "1.0", "2.0", "4.0")) +
  
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  
  labs(
    x = "Adjusted Odds Ratio (95% CI) of MAUD initiation",
    y = NULL
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"),
    panel.spacing = unit(1, "lines"),
    
    # Clean up the facet headers to act as Category titles
    strip.text.y = element_text(angle = 270, face = "bold", hjust = 0.5, size = 12, color = "white"),
    strip.background = element_rect(fill = "black", color = NA),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.y = element_text(color = "black")
  )

# Display Plot
print(p_forest)
ggsave("~/2026_fedGee/analysis/2026_04_06_figure3.pdf", 
       p_forest, 
       height=11, 
       width=8.5, 
       units="in",
       device=cairo_pdf)

