############################################################
############################################################
# TABLE GENERATION SCRIPT
# Generates:
#   - Table 1: Baseline Characteristics (with SMD)
#   - Table 2: Event Counts & Incidence Rates
# COMPLIANCE: STROBE
############################################################
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

setDT(dt)

log_msg("📊 GENERATING BASELINE & EVENT TABLES")

############################################################
# STANDARDIZED MEAN DIFFERENCE (SMD) CALCULATION
############################################################
# SMD measures covariate balance between exposure groups
# Threshold: SMD > 0.1 indicates meaningful imbalance

calc_smd <- function(x, group, ref_level = "A1") {
  if (is.numeric(x)) {
    # Continuous variables: mean difference / pooled SD
    m_ref <- mean(x[group == ref_level], na.rm = TRUE)
    sd_ref <- sd(x[group == ref_level], na.rm = TRUE)
    
    smd_list <- list()
    for (lvl in unique(group)[unique(group) != ref_level]) {
      m_lvl <- mean(x[group == lvl], na.rm = TRUE)
      sd_lvl <- sd(x[group == lvl], na.rm = TRUE)
      sd_pooled <- sqrt((sd_ref^2 + sd_lvl^2) / 2)
      smd <- abs(m_ref - m_lvl) / sd_pooled
      smd_list[[paste0("A1_vs_", lvl)]] <- smd
    }
    return(max(unlist(smd_list), na.rm = TRUE))
    
  } else {
    # Binary/categorical: proportion difference / pooled SD
    p_ref <- mean(x[group == ref_level], na.rm = TRUE)
    
    smd_list <- list()
    for (lvl in unique(group)[unique(group) != ref_level]) {
      p_lvl <- mean(x[group == lvl], na.rm = TRUE)
      p_pooled <- (p_ref + p_lvl) / 2
      sd_pooled <- sqrt(p_pooled * (1 - p_pooled))
      
      if (sd_pooled > 0) {
        smd <- abs(p_ref - p_lvl) / sd_pooled
      } else {
        smd <- 0
      }
      smd_list[[paste0("A1_vs_", lvl)]] <- smd
    }
    return(max(unlist(smd_list), na.rm = TRUE))
  }
}

############################################################
# TABLE 1: BASELINE CHARACTERISTICS BY ACR LEVEL
############################################################

generate_table1 <- function(data, acr_var = "acr.level", 
                            output_file = "Table1_Baseline.csv") {
  
  log_msg("  Generating Table 1: Baseline Characteristics...")
  
  # Formatting helpers
  fmt_cont <- function(x) sprintf("%.2f (%.2f–%.2f)", 
                                   median(x, na.rm = TRUE),
                                   quantile(x, 0.25, na.rm = TRUE),
                                   quantile(x, 0.75, na.rm = TRUE))
  
  fmt_cat <- function(x) sprintf("%.1f", mean(x, na.rm = TRUE) * 100)
  
  acr_levels <- sort(unique(data[[acr_var]]))
  n_total <- nrow(data)
  table1 <- data.table()
  
  # Row 1: N (%)
  n_by_acr <- data[, .N, by = get(acr_var)]
  table1 <- rbind(table1, data.table(
    Variable = "N (%)",
    Overall = sprintf("%d", n_total),
    A1 = sprintf("%d (%.1f%%)", n_by_acr[get(acr_var) == "A1", N], 
                 n_by_acr[get(acr_var) == "A1", N]/n_total*100),
    A2 = sprintf("%d (%.1f%%)", n_by_acr[get(acr_var) == "A2", N], 
                 n_by_acr[get(acr_var) == "A2", N]/n_total*100),
    A3 = sprintf("%d (%.1f%%)", n_by_acr[get(acr_var) == "A3", N], 
                 n_by_acr[get(acr_var) == "A3", N]/n_total*100),
    SMD = ""
  ))
  
  # Demographics section
  table1 <- rbind(table1, data.table(
    Variable = "**Demographics**",
    Overall = "", A1 = "", A2 = "", A3 = "", SMD = ""
  ))
  
  demographics <- list(
    list(var = "age", label = "Age, years", type = "cont"),
    list(var = "sex", label = "Female, %", type = "cat", cat_val = 1),
    list(var = "follow_up", label = "Follow-up duration, years", type = "cont"),
    list(var = "iyear", label = "Index year category", type = "cat_multi")
  )
  
  for (dem in demographics) {
    if (dem$type == "cont") {
      overall_val <- fmt_cont(data[[dem$var]])
      a1_val <- fmt_cont(data[get(acr_var) == "A1", get(dem$var)])
      a2_val <- fmt_cont(data[get(acr_var) == "A2", get(dem$var)])
      a3_val <- fmt_cont(data[get(acr_var) == "A3", get(dem$var)])
      smd <- calc_smd(data[[dem$var]], data[[acr_var]])
      
    } else if (dem$type == "cat") {
      overall_val <- fmt_cat(data[[dem$var]] == dem$cat_val)
      a1_val <- fmt_cat(data[get(acr_var) == "A1", get(dem$var)] == dem$cat_val)
      a2_val <- fmt_cat(data[get(acr_var) == "A2", get(dem$var)] == dem$cat_val)
      a3_val <- fmt_cat(data[get(acr_var) == "A3", get(dem$var)] == dem$cat_val)
      smd <- calc_smd(data[[dem$var]] == dem$cat_val, data[[acr_var]])
      
    } else if (dem$type == "cat_multi") {
      overall_val <- "See distribution"
      a1_val <- "-"
      a2_val <- "-"
      a3_val <- "-"
      smd <- calc_smd(as.numeric(data[[dem$var]]), data[[acr_var]])
    }
    
    table1 <- rbind(table1, data.table(
      Variable = dem$label,
      Overall = overall_val,
      A1 = a1_val, A2 = a2_val, A3 = a3_val,
      SMD = sprintf("%.3f%s", smd, ifelse(smd > 0.1, "*", ""))
    ))
  }
  
  # iyear distribution (Era 1/2/3)
  if ("iyear" %in% names(data)) {
    for (era in c("Era 1", "Era 2", "Era 3")) {
      era_label <- switch(era,
                          "Era 1" = "  Era 1: 2005–2009, %",
                          "Era 2" = "  Era 2: 2010–2014, %",
                          "Era 3" = "  Era 3: 2015–2020, %")
      
      table1 <- rbind(table1, data.table(
        Variable = era_label,
        Overall = fmt_cat(data$iyear == era),
        A1 = fmt_cat(data[get(acr_var) == "A1", iyear] == era),
        A2 = fmt_cat(data[get(acr_var) == "A2", iyear] == era),
        A3 = fmt_cat(data[get(acr_var) == "A3", iyear] == era),
        SMD = sprintf("%.3f", calc_smd(data$iyear == era, data[[acr_var]]))
      ))
    }
  }
  
  # Physical Measurements
  table1 <- rbind(table1, data.table(
    Variable = "**Physical Measurements, median (IQR)**",
    Overall = "", A1 = "", A2 = "", A3 = "", SMD = ""
  ))
  
  phys_vars <- list(
    list(var = "b.bmi", label = "BMI, kg/m²"),
    list(var = "b.DBP", label = "Diastolic BP, mmHg"),
    list(var = "b.SBP", label = "Systolic BP, mmHg")
  )
  
  for (pv in phys_vars) {
    overall_val <- fmt_cont(data[[pv$var]])
    a1_val <- fmt_cont(data[get(acr_var) == "A1", get(pv$var)])
    a2_val <- fmt_cont(data[get(acr_var) == "A2", get(pv$var)])
    a3_val <- fmt_cont(data[get(acr_var) == "A3", get(pv$var)])
    smd <- calc_smd(data[[pv$var]], data[[acr_var]])
    
    table1 <- rbind(table1, data.table(
      Variable = pv$label,
      Overall = overall_val,
      A1 = a1_val, A2 = a2_val, A3 = a3_val,
      SMD = sprintf("%.3f%s", smd, ifelse(smd > 0.1, "*", ""))
    ))
  }
  
  # Laboratory Measurements
  table1 <- rbind(table1, data.table(
    Variable = "**Laboratory Measurements, median (IQR)**",
    Overall = "", A1 = "", A2 = "", A3 = "", SMD = ""
  ))
  
  lab_vars <- list(
    list(var = "b.A1c", label = "HbA1c (NGSP %)"),
    list(var = "acr", label = "uACR, mg/mmol"),
    list(var = "b.eGFR", label = "eGFR, mL/min/1.73m²"),
    list(var = "b.LDL", label = "LDL-Cholesterol"),
    list(var = "b.HDL", label = "HDL-Cholesterol")
  )
  
  for (lv in lab_vars) {
    overall_val <- fmt_cont(data[[lv$var]])
    a1_val <- fmt_cont(data[get(acr_var) == "A1", get(lv$var)])
    a2_val <- fmt_cont(data[get(acr_var) == "A2", get(lv$var)])
    a3_val <- fmt_cont(data[get(acr_var) == "A3", get(lv$var)])
    smd <- calc_smd(data[[lv$var]], data[[acr_var]])
    
    table1 <- rbind(table1, data.table(
      Variable = lv$label,
      Overall = overall_val,
      A1 = a1_val, A2 = a2_val, A3 = a3_val,
      SMD = sprintf("%.3f%s", smd, ifelse(smd > 0.1, "*", ""))
    ))
  }
  
  # Comorbidities
  table1 <- rbind(table1, data.table(
    Variable = "**Baseline Comorbidities, %**",
    Overall = "", A1 = "", A2 = "", A3 = "", SMD = ""
  ))
  
  comorb_vars <- c(
    "ckd" = "Chronic Kidney Disease (CKD)",
    "ht" = "Hypertension",
    "lipid" = "Dyslipidaemia",
    "ane" = "Anaemia",
    "cancer" = "Cancer",
    "copd" = "COPD",
    "ast" = "Asthma",
    "dem" = "Dementia",
    "gout" = "Gout",
    "pad" = "Peripheral Artery Disease (PAD)",
    "lupus" = "Systemic lupus erythematosus (SLE)"
  )
  
  for (cv in names(comorb_vars)) {
    overall_val <- fmt_cat(data[[cv]])
    a1_val <- fmt_cat(data[get(acr_var) == "A1", get(cv)])
    a2_val <- fmt_cat(data[get(acr_var) == "A2", get(cv)])
    a3_val <- fmt_cat(data[get(acr_var) == "A3", get(cv)])
    smd <- calc_smd(data[[cv]], data[[acr_var]])
    
    table1 <- rbind(table1, data.table(
      Variable = comorb_vars[cv],
      Overall = paste0(overall_val, "%"),
      A1 = paste0(a1_val, "%"),
      A2 = paste0(a2_val, "%"),
      A3 = paste0(a3_val, "%"),
      SMD = sprintf("%.3f%s", smd, ifelse(smd > 0.1, "*", ""))
    ))
  }
  
  # Medications
  table1 <- rbind(table1, data.table(
    Variable = "**Baseline Medications, %**",
    Overall = "", A1 = "", A2 = "", A3 = "", SMD = ""
  ))
  
  med_vars <- c(
    "m.Metformin" = "Metformin",
    "m.Sulfony" = "Sulphonylureas",
    "m.Insulins" = "Insulin",
    "m.DPP4" = "SGLT2i, GLP-1 RAs, DPP-4i, and others",
    "m.ACEi" = "ACE Inhibitors",
    "m.ARBs" = "Angiotensin II receptor blockers (ARBs)",
    "m.BetaB" = "Beta-blockers",
    "m.CCBlocker" = "Calcium channel blockers (CCBs)",
    "m.Diuretics" = "Loop diuretics",
    "m.Thiazide" = "Thiazide diuretics",
    "m.MRAs" = "Mineralocorticoid receptor antagonists (MRAs)",
    "m.Antiplatelet" = "Antiplatelet agents",
    "m.OAC" = "Oral anticoagulants (OACs)",
    "m.Statins" = "Lipid-lowering therapies"
  )
  
  for (mv in names(med_vars)) {
    overall_val <- fmt_cat(data[[mv]])
    a1_val <- fmt_cat(data[get(acr_var) == "A1", get(mv)])
    a2_val <- fmt_cat(data[get(acr_var) == "A2", get(mv)])
    a3_val <- fmt_cat(data[get(acr_var) == "A3", get(mv)])
    smd <- calc_smd(data[[mv]], data[[acr_var]])
    
    table1 <- rbind(table1, data.table(
      Variable = med_vars[mv],
      Overall = paste0(overall_val, "%"),
      A1 = paste0(a1_val, "%"),
      A2 = paste0(a2_val, "%"),
      A3 = paste0(a3_val, "%"),
      SMD = sprintf("%.3f%s", smd, ifelse(smd > 0.1, "*", ""))
    ))
  }
  
  # Export
  fwrite(table1, output_file, sep = ",")
  log_msg("  ✓ Table 1 saved: ", output_file)
  log_msg("  📌 SMD > 0.1 marked with * (meaningful imbalance)")
  
  return(table1)
}

############################################################
# TABLE 2: EVENT COUNTS AND INCIDENCE RATES
############################################################

generate_table2 <- function(data, output_file = "Table2_Events.csv") {
  
  log_msg("  Generating Table 2: Event Counts and Incidence Rates...")
  
  outcomes <- list(
    cvd = list(time = "time_cvd", status = "status_cvd", label = "CVD composite"),
    acs = list(time = "time_acs", status = "status_acs", label = "ACS"),
    hf = list(time = "time_hf", status = "status_hf", label = "Heart Failure"),
    af = list(time = "time_af", status = "status_af", label = "Atrial Fibrillation"),
    stroke = list(time = "time_stroke", status = "status_stroke", label = "Stroke"),
    death = list(time = "time_death", status = "status_death", label = "Mortality")
  )
  
  acr_levels <- c("A1", "A2", "A3")
  table2 <- data.table()
  
  for (out_name in names(outcomes)) {
    out <- outcomes[[out_name]]
    
    for (acr_lvl in acr_levels) {
      d_sub <- data[acr.level == acr_lvl]
      
      t_vec <- d_sub[[out$time]]
      s_vec <- d_sub[[out$status]]
      
      n_events <- sum(s_vec == 1, na.rm = TRUE)
      
      follow_up_times <- t_vec[!is.na(t_vec)]
      median_fu <- median(follow_up_times, na.rm = TRUE)
      q1_fu <- quantile(follow_up_times, 0.25, na.rm = TRUE)
      q3_fu <- quantile(follow_up_times, 0.75, na.rm = TRUE)
      fu_str <- sprintf("%.2f (%.2f–%.2f)", median_fu, q1_fu, q3_fu)
      
      # Auto-detect time units (days vs years)
      if (median_fu > 20) {
        person_years <- sum(follow_up_times, na.rm = TRUE) / 365.25
      } else {
        person_years <- sum(follow_up_times, na.rm = TRUE)
      }
      
      ir <- (n_events / person_years) * 100
      
      if (n_events > 0) {
        ci_lower <- qpois(0.025, n_events) / person_years * 100
        ci_upper <- qpois(0.975, n_events + 1) / person_years * 100
      } else {
        ci_lower <- 0
        ci_upper <- 0
      }
      ir_str <- sprintf("%.2f (%.2f–%.2f)", ir, ci_lower, ci_upper)
      
      table2 <- rbind(table2, data.table(
        Outcome = out$label,
        uACR_Level = acr_lvl,
        N = nrow(d_sub),
        Events = n_events,
        `Follow-up time Median [IQR]` = fu_str,
        `Person-years` = sprintf("%.2f", person_years),
        `Incident Rate per 100 PY (95% CI)` = ir_str
      ))
    }
  }
  
  fwrite(table2, output_file, sep = ",")
  log_msg("  ✓ Table 2 saved: ", output_file)
  
  return(table2)
}

############################################################
# RUN TABLE GENERATION
############################################################

table1 <- generate_table1(dt, acr_var = "acr.level", 
                          output_file = "acr_prognosis_v3/tables/Table1_Baseline.csv")
table2 <- generate_table2(dt, 
                          output_file = "acr_prognosis_v3/tables/Table2_Events.csv")

log_msg("\n✅ TABLE GENERATION COMPLETE")
log_msg("📋 Table 1 variables: ", nrow(table1))
log_msg("📋 Table 2 outcomes: ", nrow(table2))