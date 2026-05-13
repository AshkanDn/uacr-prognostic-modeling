############################################################
# FINAL RESULTS PIPELINE 
# Input: .RData bundles from the final analysis pipeline
#
# PURPOSE
#   - Validate analysis bundle integrity before exporting results
#   - Extract HR tables from stored fitted models
#   - Produce absolute risk, risk difference, and risk ratio tables
#   - Produce unadjusted CIF curves for competing-risk endpoints
#   - Produce adjusted curves (whole cohort and sex-stratified)
#   - Produce grouped forest plots for A2 and A3 on log scale
#   - Produce spline-band figures (whole cohort and sex-stratified)
#   - Include imputation diagnostics outputs in the report index
#   - Compile a manuscript-ready Word report
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(flextable)
  library(officer)
  library(riskRegression)
  library(survival)
  library(prodlim)
})

############################################################
# 0) LOAD ANALYSIS ARTIFACTS + STRUCTURE CHECKS
############################################################
base_dir <- "acr_prognosis_final"
dir_table  <- file.path(base_dir, "tables")
dir_plot   <- file.path(base_dir, "figures")
dir_report <- file.path(base_dir, "report")
dir_score  <- file.path(base_dir, "scores")
dir.create(dir_table,  recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plot,   recursive = TRUE, showWarnings = FALSE)
dir.create(dir_report, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_score,   recursive = TRUE, showWarnings = FALSE)

load_bundle <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists("bundle", envir = e)) stop("No 'bundle' object found in: ", path)
  get("bundle", envir = e)
}

manifest_path <- file.path(base_dir, "manifest.RData")
if (!file.exists(manifest_path)) stop("manifest.RData not found in: ", base_dir)
load(manifest_path)

if (!exists("manifest") || !all(c("ep", "file") %in% names(manifest))) {
  stop("manifest.RData exists but does not contain a valid manifest with ep/file columns.")
}

bundles <- lapply(manifest$file, load_bundle)
names(bundles) <- manifest$ep

msg <- function(...) message(paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = "")))
msg("Result pipeline initialised | endpoints = ", paste(names(bundles), collapse = ", "))

############################################################
# 1) VISUAL THEME / LABELS
############################################################
COL <- list(
  A1 = "#1B5E20", A2 = "#0D47A1", A3 = "#C62828",
  txt = "#37474F", grid = "#ECEFF1",
  male = "#0D47A1", female = "#C62828"
)

thm <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_line(color = COL$grid, linewidth = 0.2),
    axis.line = element_line(color = COL$txt, linewidth = 0.4),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    legend.position = "bottom",
    legend.key = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

endpoint_labels <- c(
  cvd = "CVD",
  acs = "ACS",
  hf = "HF",
  af = "AF",
  stroke = "Stroke",
  pad = "PAD",
  CVDp = "CVD+P",
  death = "Mortality"
)

outcome_order <- c("CVD", "ACS", "HF", "AF", "Stroke", "PAD", "CVD+P", "Mortality")

pretty_endpoint <- function(x) {
  ifelse(x %in% names(endpoint_labels), unname(endpoint_labels[x]), x)
}

safe_endpoint_file <- function(x) gsub("[^A-Za-z0-9]+", "", x)

############################################################
# 2) SMALL UTILITY FUNCTIONS
############################################################
`%||%` <- function(x, y) if (is.null(x)) y else x

get_bundle_models <- function(bundle) bundle$models %||% NULL
get_sex_bundle <- function(bundle) bundle$sex_analysis %||% bundle$sex %||% NULL

assert_bundle_structure <- function(bundle, ep = "unknown") {
  required_top <- c("meta", "data", "ate", "models", "spline", "sex_analysis", "rhs", "pooled")
  miss <- setdiff(required_top, names(bundle))
  if (length(miss) > 0) {
    stop("Bundle '", ep, "' is missing required top-level fields: ", paste(miss, collapse = ", "))
  }
  if (is.null(bundle$meta$ep) || is.null(bundle$meta$type) || is.null(bundle$meta$label)) {
    stop("Bundle '", ep, "' has incomplete meta information.")
  }
  if (is.null(bundle$pooled$mean) || is.null(bundle$pooled$diff) || is.null(bundle$pooled$ratio)) {
    stop("Bundle '", ep, "' lacks pooled ATE outputs needed for results.")
  }
  TRUE
}

# Robust coefficient table extraction from Cox / CSC / related fits.
coef_table_from_fit <- function(mod) {
  s <- summary(mod)
  tab <- NULL
  if (!is.null(s$coefficients)) tab <- as.data.frame(s$coefficients)
  if (is.null(tab) && !is.null(s$coef)) tab <- as.data.frame(s$coef)
  if (is.null(tab) && !is.null(s$coefmat)) tab <- as.data.frame(s$coefmat)
  if (is.null(tab) && !is.null(s$mat)) tab <- as.data.frame(s$mat)
  if (is.null(tab)) stop("Could not extract coefficient table from model summary.")

  beta <- if ("coef" %in% names(tab)) tab[["coef"]] else
    if ("Estimate" %in% names(tab)) tab[["Estimate"]] else
      if ("exp(coef)" %in% names(tab)) log(tab[["exp(coef)"]]) else tab[[1]]
  se <- if ("se(coef)" %in% names(tab)) tab[["se(coef)"]] else
    if ("Std. Error" %in% names(tab)) tab[["Std. Error"]] else
      if ("se" %in% names(tab)) tab[["se"]] else rep(NA_real_, length(beta))
  p <- if ("Pr(>|z|)" %in% names(tab)) tab[["Pr(>|z|)"]] else
    if ("p" %in% names(tab)) tab[["p"]] else rep(NA_real_, length(beta))

  data.table(
    Variable = rownames(tab),
    HR = exp(beta),
    LCL = exp(beta - 1.96 * se),
    UCL = exp(beta + 1.96 * se),
    p = p,
    beta = beta,
    se = se
  )
}

fmt_hr <- Vectorize(function(hr, l, u, p) {
  x <- sprintf("%.2f(%.2f-%.2f)", hr, l, u)
  if (!is.na(p) && p < 0.001) paste0(x, "
*p<0.001") else x
})

fmt_risk <- Vectorize(function(est, l, u) {
  sprintf("%.2f(%.2f-%.2f)", est * 100, l * 100, u * 100)
})

standardize_ate_tbl <- function(x) {
  dt <- as.data.table(x)
  if (nrow(dt) == 0) return(dt)

  old <- names(dt)
  clean <- tolower(gsub("[^a-z0-9]+", "", old))
  setnames(dt, old, clean)

  rename_first <- function(candidates, new) {
    hit <- intersect(names(dt), candidates)
    if (length(hit) > 0) setnames(dt, hit[1], new)
  }

  rename_first(c("time", "times"), "time")
  rename_first(c("estimate", "est", "risk", "value"), "estimate")
  rename_first(c("se", "stderr", "stderror", "standarderror"), "se")
  rename_first(c("lcl", "lower", "cilower", "lowerci", "ci.lower"), "lcl")
  rename_first(c("ucl", "upper", "ciupper", "upperci", "ci.upper"), "ucl")
  rename_first(c("pvalue", "p", "pvalueadj"), "p")
  rename_first(c("treatment", "trt", "level", "acrlevel", "acr"), "treatment")
  rename_first(c("reference", "ref", "base"), "reference")
  rename_first(c("comparison", "contrast", "comp", "trt2"), "comparison")

  if (!all(c("lcl", "ucl") %in% names(dt)) && all(c("se", "estimate") %in% names(dt))) {
    dt[, `:=`(lcl = estimate - 1.96 * se, ucl = estimate + 1.96 * se)]
  }
  dt[]
}

calc_evalue_hr <- function(hr, lcl) {
  hr_use <- ifelse(hr < 1, 1 / hr, hr)
  lcl_use <- ifelse(lcl < 1, 1 / lcl, lcl)
  point <- hr_use + sqrt(hr_use * (hr_use - 1))
  ci_e <- ifelse(lcl_use > 1, lcl_use + sqrt(lcl_use * (lcl_use - 1)), 1)
  data.table(evalue = point, evalue_ci = ci_e)
}

############################################################
# 3) PERFORMANCE / QC HELPERS
############################################################
run_calibration <- function(model, data, time_var, status_var, cause, times) {
  if (is.na(cause)) {
    Score(
      list(full = model),
      formula  = as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ 1")),
      data     = data,
      times    = times,
      plots    = "Calibration",
      metrics  = c("AUC", "Brier"),
      conf.int = TRUE,
      se.fit   = TRUE
    )
  } else {
    Score(
      list(full = model),
      formula  = as.formula(paste0("Hist(", time_var, ", ", status_var, ") ~ 1")),
      data     = data,
      cause    = cause,
      times    = times,
      plots    = "Calibration",
      metrics  = c("AUC", "Brier"),
      conf.int = TRUE,
      se.fit   = TRUE
    )
  }
}

extract_performance_summary <- function(score_obj, model_name = "full") {
  if (is.null(score_obj)) return(NULL)
  out <- list()
  for (metric in c("AUC", "Brier")) {
    obj <- score_obj[[metric]]
    if (is.null(obj)) next
    dtm <- if (!is.null(obj$score)) as.data.table(obj$score) else as.data.table(obj)
    if (nrow(dtm) == 0) next

    old <- names(dtm)
    clean <- tolower(gsub("[^a-z0-9]+", "", old))
    setnames(dtm, old, clean)

    model_col <- intersect(names(dtm), c("model", "models"))
    time_col  <- intersect(names(dtm), c("time", "times"))
    est_col   <- intersect(names(dtm), c("auc", "brier", "score", "estimate"))
    se_col    <- intersect(names(dtm), c("se", "stderr", "stderror"))

    if (length(model_col) > 0) {
      keep <- dtm[[model_col[1]]] %in% c(model_name, "full", 1)
      if (any(keep, na.rm = TRUE)) dtm <- dtm[keep]
    }
    if (length(est_col) == 0) next

    out[[metric]] <- data.table(
      metric = metric,
      time = if (length(time_col) > 0) dtm[[time_col[1]]] else NA_real_,
      estimate = dtm[[est_col[1]]],
      se = if (length(se_col) > 0) dtm[[se_col[1]]] else NA_real_
    )
  }
  rbindlist(out, fill = TRUE)
}

pool_scalar_rubin <- function(est, se) {
  m <- length(est)
  qbar <- mean(est, na.rm = TRUE)
  ubar <- mean(se^2, na.rm = TRUE)
  b <- stats::var(est, na.rm = TRUE)
  tvar <- ubar + (1 + 1/m) * b
  list(est = qbar, se = sqrt(tvar), lcl = qbar - 1.96 * sqrt(tvar), ucl = qbar + 1.96 * sqrt(tvar))
}

pool_performance_metrics <- function(perf_list) {
  dat <- rbindlist(lapply(seq_along(perf_list), function(i) {
    x <- perf_list[[i]]
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x[, imp := i]
    x
  }), fill = TRUE)
  if (nrow(dat) == 0) return(NULL)
  dat[, {
    pooled <- pool_scalar_rubin(estimate, ifelse(is.na(se), 0, se))
    .(
      mean_est = mean(estimate, na.rm = TRUE),
      min_est = min(estimate, na.rm = TRUE),
      max_est = max(estimate, na.rm = TRUE),
      pooled_est = pooled$est,
      pooled_lcl = pooled$lcl,
      pooled_ucl = pooled$ucl,
      n_imp = .N
    )
  }, by = .(metric, time)]
}

summarise_missingness_diagnostics <- function(base_dir) {
  list(
    trace_pdf = file.path(base_dir, "figures", "Imputation_TracePlots.pdf"),
    density_pdf = file.path(base_dir, "figures", "Imputation_DensityPlots.pdf"),
    strip_pdf = file.path(base_dir, "figures", "Imputation_StripPlots.pdf"),
    missingness_table = file.path(base_dir, "tables", "Table_Missingness.csv"),
    predictor_matrix = file.path(base_dir, "tables", "Imputation_PredictorMatrix.csv"),
    visit_sequence = file.path(base_dir, "tables", "Imputation_VisitSequence.txt"),
    imputation_model = file.path(base_dir, "imputation_model.RData")
  )
}

############################################################
# 4) PLOT HELPERS
############################################################
plot_grouped_forest <- function(df, file_png, title_text, group_levels, group_cols,
                                outcome_order = outcome_order,
                                x_limits = c(0.8, 4.5), label_x = 3.7) {
  x_breaks <- c(1, 1.5, 2, 2.5, 3, 4)
  x_labels <- c("1.0", "1.5", "2.0", "2.5", "3.0", "4.0")

  df[, outcome := factor(outcome, levels = rev(outcome_order))]
  df[, group := factor(group, levels = group_levels)]

  p <- ggplot(df, aes(x = est, y = outcome, color = group)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(
      aes(xmin = low, xmax = high),
      height = 0.18,
      linewidth = 0.7,
      position = position_dodge(width = 0.55)
    ) +
    geom_point(size = 3.0, position = position_dodge(width = 0.55)) +
    geom_text(
      aes(x = label_x, label = label_text),
      hjust = 0,
      size = 3.15,
      position = position_dodge(width = 0.55)
    ) +
    scale_color_manual(values = group_cols) +
    scale_x_log10(limits = x_limits, breaks = x_breaks, labels = x_labels) +
    coord_cartesian(clip = "off") +
    labs(x = "Hazard ratio (95% CI)", y = NULL, color = NULL, title = title_text) +
    thm +
    theme(panel.grid.major.y = element_blank())

  ggsave(file_png, p, width = 8.2, height = 5.8, dpi = 300)
  invisible(p)
}

############################################################
# 5) TABLE / FIGURE WRITERS
############################################################
write_grouped_forest_outputs <- function(bundles, out_table, out_plot,
                                         terms = c("acr.levelA2", "acr.levelA3")) {
  out <- list()
  for (term in terms) {
    safe_term <- gsub("[^A-Za-z0-9]+", "", term)

    main_df <- rbindlist(lapply(names(bundles), function(nm) {
      b <- bundles[[nm]]
      mods <- get_bundle_models(b)
      mod_un <- NULL; mod_fu <- NULL
      if (!is.null(mods$csc)) {
        mod_un <- mods$csc$un %||% NULL
        mod_fu <- mods$csc$fu %||% NULL
      }
      if (is.null(mod_un) && !is.null(mods$cox)) mod_un <- mods$cox$un %||% NULL
      if (is.null(mod_fu) && !is.null(mods$cox)) mod_fu <- mods$cox$fu %||% NULL

      extract_row <- function(mod, model_label) {
        if (is.null(mod)) return(NULL)
        ct <- coef_table_from_fit(mod)
        row <- ct[Variable == term]
        if (nrow(row) == 0) return(NULL)
        data.table(
          outcome = nm,
          group = model_label,
          est = row$HR,
          low = row$LCL,
          high = row$UCL,
          p = row$p,
          label_text = sprintf("%.2f (%.2f-%.2f)", row$HR, row$LCL, row$UCL)
        )
      }

      rbind(extract_row(mod_un, "Unadjusted"), extract_row(mod_fu, "Fully adjusted"), fill = TRUE)
    }), fill = TRUE)

    fwrite(main_df, file.path(out_table, paste0("Forest_MainPopulation_", safe_term, ".csv")))
    plot_grouped_forest(
      df = main_df,
      file_png = file.path(out_plot, paste0("Forest_MainPopulation_", safe_term, ".png")),
      title_text = paste0("Main population: ", term),
      group_levels = c("Unadjusted", "Fully adjusted"),
      group_cols = c("Unadjusted" = "#455A64", "Fully adjusted" = "#0D47A1")
    )

    sex_df <- rbindlist(lapply(names(bundles), function(nm) {
      sx <- get_sex_bundle(bundles[[nm]])
      if (is.null(sx)) return(NULL)
      strat <- sx$strat_models %||% sx$strat %||% NULL
      if (is.null(strat)) return(NULL)

      extract_row <- function(mod, sex_label) {
        if (is.null(mod)) return(NULL)
        ct <- coef_table_from_fit(mod)
        row <- ct[Variable == term]
        if (nrow(row) == 0) return(NULL)
        data.table(
          outcome = nm,
          group = sex_label,
          est = row$HR,
          low = row$LCL,
          high = row$UCL,
          p = row$p,
          label_text = sprintf("%.2f (%.2f-%.2f)", row$HR, row$LCL, row$UCL)
        )
      }

      rbind(extract_row(strat$male, "Male"), extract_row(strat$female, "Female"), fill = TRUE)
    }), fill = TRUE)

    fwrite(sex_df, file.path(out_table, paste0("Forest_SexStratified_", safe_term, ".csv")))
    plot_grouped_forest(
      df = sex_df,
      file_png = file.path(out_plot, paste0("Forest_SexStratified_", safe_term, ".png")),
      title_text = paste0("Sex-stratified: ", term),
      group_levels = c("Male", "Female"),
      group_cols = c("Male" = COL$male, "Female" = COL$female)
    )

    out[[term]] <- list(main = main_df, sex = sex_df)
  }
  invisible(out)
}

write_adjusted_curve_outputs <- function(bundle, out_table, out_plot, c_times) {
  ep <- bundle$meta$ep
  cause <- if (bundle$meta$type == "competing") bundle$meta$cause else NA_integer_

  if (is.null(bundle$ate$fu)) return(invisible(FALSE))
  df_overall <- standardize_ate_tbl(bundle$ate$fu$curve)
  trt_col <- intersect(names(df_overall), c("treatment", "acr", "level", "acrlevel"))
  df_overall[, ACR := if (length(trt_col) > 0) as.character(get(trt_col[1])) else NA_character_]
  df_overall[, `:=`(
    Risk = estimate * 100,
    LCL  = if ("lcl" %in% names(df_overall)) lcl * 100 else NA_real_,
    UCL  = if ("ucl" %in% names(df_overall)) ucl * 100 else NA_real_
  )]
  fwrite(df_overall, file.path(out_table, paste0("AdjustedCurve_Overall_", ep, ".csv")))

  p_overall <- ggplot(df_overall, aes(x = time, y = Risk, color = ACR, fill = ACR)) +
    geom_line(linewidth = 1.1) +
    geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.12, color = NA) +
    scale_color_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
    scale_fill_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
    labs(
      title = paste0("Fully adjusted risk curves: ", ep),
      x = "Years from diagnosis",
      y = "Adjusted absolute risk (%)",
      color = "uACR",
      fill = "uACR"
    ) + thm
  ggsave(file.path(out_plot, paste0("AdjustedCurve_Overall_", ep, ".png")), p_overall, width = 8, height = 6, dpi = 300)

  sx <- get_sex_bundle(bundle)
  strat <- sx$strat_models %||% sx$strat %||% NULL
  if (!is.null(strat)) {
    d_m <- bundle$data[sex == "Male"]
    d_f <- bundle$data[sex == "Female"]
    male_fit <- strat$male
    female_fit <- strat$female

    if (!is.null(male_fit) && !is.null(female_fit)) {
      male_curve <- run_ate_local(male_fit, d_m, "acr.level", c_times, cause)$curve
      female_curve <- run_ate_local(female_fit, d_f, "acr.level", c_times, cause)$curve
      male_curve[, sex := "Male"]
      female_curve[, sex := "Female"]
      sex_curve <- rbind(male_curve, female_curve, fill = TRUE)

      trt_col2 <- intersect(names(sex_curve), c("treatment", "acr", "level", "acrlevel"))
      sex_curve[, ACR := if (length(trt_col2) > 0) as.character(get(trt_col2[1])) else NA_character_]
      sex_curve[, `:=`(
        Risk = estimate * 100,
        LCL  = if ("lcl" %in% names(sex_curve)) lcl * 100 else NA_real_,
        UCL  = if ("ucl" %in% names(sex_curve)) ucl * 100 else NA_real_
      )]

      fwrite(sex_curve, file.path(out_table, paste0("AdjustedCurve_SexStratified_", ep, ".csv")))

      p_sex <- ggplot(sex_curve, aes(x = time, y = Risk, color = ACR, fill = ACR)) +
        geom_line(linewidth = 1.0) +
        geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.10, color = NA) +
        facet_wrap(~ sex) +
        scale_color_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
        scale_fill_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
        labs(
          title = paste0("Sex-stratified adjusted risk curves: ", ep),
          x = "Years from diagnosis",
          y = "Adjusted absolute risk (%)",
          color = "uACR",
          fill = "uACR"
        ) + thm
      ggsave(file.path(out_plot, paste0("AdjustedCurve_SexStratified_", ep, ".png")),
             p_sex, width = 9, height = 6, dpi = 300)
    }
  }
  invisible(TRUE)
}

write_spline_outputs <- function(bundle, out_table, out_plot, B = 100, seed = 1117) {
  ep <- bundle$meta$ep
  cause <- if (bundle$meta$type == "competing") bundle$meta$cause else NA_integer_
  if (is.null(bundle$spline$grid)) return(invisible(FALSE))
  grid <- bundle$spline$grid

  predict_spline_point <- function(fit, data, grid, sex_value = NULL, cause = NA_integer_) {
    rbindlist(lapply(grid, function(v) {
      nd <- copy(data)
      nd[, log_acr := log(v + 0.01)]
      if (!is.null(sex_value) && "sex" %in% names(nd)) {
        nd[, sex := factor(sex_value, levels = levels(data$sex))]
      }
      pr <- if (is.na(cause)) predictRisk(fit, newdata = nd, times = 5) else predictRisk(fit, newdata = nd, times = 5, cause = cause)
      data.table(acr = v, Risk = mean(as.matrix(pr)[, 1], na.rm = TRUE) * 100)
    }))
  }

  fit_overall <- bundle$spline$model %||% bundle$models$spl %||% bundle$models$fu %||% bundle$models$cox$spl %||% bundle$models$csc$spl
  if (!is.null(fit_overall)) {
    bootstrap_overall <- function(data, grid, B, seed) {
      set.seed(seed)
      n <- nrow(data)
      boot <- matrix(NA_real_, nrow = B, ncol = length(grid))
      for (b in seq_len(B)) {
        idx <- sample.int(n, size = n, replace = TRUE)
        db <- data[idx]
        fit_b <- if (!is.null(bundle$rhs$spl)) {
          if (bundle$meta$type == "competing") CSC(as.formula(paste0("Hist(t,s)~", bundle$rhs$spl)), data = db, cause = cause) else coxph(as.formula(paste0("Surv(t,s)~", bundle$rhs$spl)), data = db, x = TRUE, y = TRUE)
        } else {
          fit_overall
        }
        boot[b, ] <- sapply(grid, function(v) {
          nd <- copy(db)
          nd[, log_acr := log(v + 0.01)]
          pr <- if (is.na(cause)) predictRisk(fit_b, newdata = nd, times = 5) else predictRisk(fit_b, newdata = nd, times = 5, cause = cause)
          mean(as.matrix(pr)[, 1], na.rm = TRUE) * 100
        })
        if (b %% 10 == 0) msg("Spline bootstrap overall: ", ep, " ", b, "/", B)
      }
      data.table(acr = grid, risk = colMeans(boot, na.rm = TRUE), lcl = apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE), ucl = apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE))
    }

    point_overall <- predict_spline_point(fit_overall, bundle$data, grid, cause = cause)
    band_overall <- bootstrap_overall(bundle$data, grid, B, seed)
    spline_overall <- merge(point_overall, band_overall, by = "acr", all.x = TRUE)
    fwrite(spline_overall, file.path(out_table, paste0("Spline_Overall_", ep, ".csv")))

    p_overall <- ggplot(spline_overall, aes(x = acr, y = Risk)) +
      geom_line(color = COL$A2, linewidth = 1.2) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), fill = COL$A2, alpha = 0.18) +
      labs(
        title = paste0("Whole cohort 5-year spline risk: ", ep),
        x = "Continuous uACR",
        y = "Adjusted 5-year risk (%)"
      ) + thm
    ggsave(file.path(out_plot, paste0("Spline_Overall_", ep, ".png")), p_overall, width = 8, height = 6, dpi = 300)
  }

  sx <- get_sex_bundle(bundle)
  sex_fit <- sx$spline_model %||% NULL
  if (!is.null(sex_fit) && !is.null(bundle$rhs$spl_sex)) {
    point_sex <- rbindlist(lapply(c("Male", "Female"), function(sex_value) {
      p <- predict_spline_point(sex_fit, bundle$data, grid, sex_value = sex_value, cause = cause)
      p[, sex := sex_value]
      p
    }), fill = TRUE)

    band_sex <- rbindlist(lapply(c("Male", "Female"), function(sex_value) {
      set.seed(seed)
      n <- nrow(bundle$data)
      boot <- matrix(NA_real_, nrow = B, ncol = length(grid))
      for (b in seq_len(B)) {
        idx <- sample.int(n, size = n, replace = TRUE)
        db <- bundle$data[idx]
        fit_b <- if (bundle$meta$type == "competing") CSC(as.formula(paste0("Hist(t,s)~", bundle$rhs$spl_sex)), data = db, cause = cause) else coxph(as.formula(paste0("Surv(t,s)~", bundle$rhs$spl_sex)), data = db, x = TRUE, y = TRUE)
        boot[b, ] <- sapply(grid, function(v) {
          nd <- copy(db)
          nd[, log_acr := log(v + 0.01)]
          nd[, sex := factor(sex_value, levels = levels(db$sex))]
          pr <- if (bundle$meta$type == "competing") predictRisk(fit_b, newdata = nd, times = 5, cause = cause) else predictRisk(fit_b, newdata = nd, times = 5)
          mean(as.matrix(pr)[, 1], na.rm = TRUE) * 100
        })
      }
      data.table(acr = grid, sex = sex_value, risk = colMeans(boot, na.rm = TRUE), lcl = apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE), ucl = apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE))
    }), fill = TRUE)

    spline_sex <- merge(point_sex, band_sex, by = c("acr", "sex"), all.x = TRUE)
    fwrite(spline_sex, file.path(out_table, paste0("Spline_SexStratified_", ep, ".csv")))

    p_sex <- ggplot(spline_sex, aes(x = acr, y = Risk, color = sex, fill = sex)) +
      geom_line(linewidth = 1.1) +
      geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.15, color = NA) +
      facet_wrap(~ sex) +
      scale_color_manual(values = c(Male = COL$male, Female = COL$female)) +
      scale_fill_manual(values = c(Male = COL$male, Female = COL$female)) +
      labs(
        title = paste0("Sex-stratified 5-year spline risk: ", ep),
        x = "Continuous uACR",
        y = "Adjusted 5-year risk (%)"
      ) + thm
    ggsave(file.path(out_plot, paste0("Spline_SexStratified_", ep, ".png")), p_sex, width = 9, height = 6, dpi = 300)
  }
  invisible(TRUE)
}

write_sex_ate_tables <- function(bundles, out_table) {
  risk_rows <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    sx <- get_sex_bundle(b)
    if (is.null(sx) || is.null(sx$ate)) return(NULL)

    make_block <- function(x, sex_label, metric_name) {
      if (is.null(x) || is.null(x[[metric_name]])) return(NULL)
      m <- standardize_ate_tbl(x[[metric_name]])
      if (nrow(m) == 0) return(NULL)
      data.table(
        Outcome = pretty_endpoint(nm),
        Sex = sex_label,
        Metric = metric_name,
        Time = m$time,
        Estimate = m$estimate,
        LCL = m$lcl,
        UCL = m$ucl
      )
    }

    rbind(
      make_block(sx$ate$male,   "Male",   "mean"),
      make_block(sx$ate$female, "Female", "mean"),
      make_block(sx$ate$male,   "Male",   "diff"),
      make_block(sx$ate$female, "Female", "diff"),
      make_block(sx$ate$male,   "Male",   "ratio"),
      make_block(sx$ate$female, "Female", "ratio"),
      fill = TRUE
    )
  }), fill = TRUE)
  if (nrow(risk_rows) > 0) fwrite(risk_rows, file.path(out_table, "Table_SexStratified_ATE.csv"))
}

write_performance_tables <- function(bundles, out_table) {
  perf_rows <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    x <- b$pooled$performance %||% b$mi$performance_summary %||% NULL
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x[, Outcome := pretty_endpoint(nm)]
    x
  }), fill = TRUE)
  if (nrow(perf_rows) > 0) fwrite(perf_rows, file.path(out_table, "Table_Performance_Summary.csv"))
}

write_evalues_table <- function(bundles, out_table) {
  ev <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    x <- b$pooled$evalues %||% b$mi$pooled_evalues %||% NULL
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x[, Outcome := pretty_endpoint(nm)]
    x
  }), fill = TRUE)
  if (nrow(ev) > 0) fwrite(ev, file.path(out_table, "Table_EValues.csv"))
}

write_fgr_sensitivity_tables <- function(bundles, out_table) {
  fgr_rows <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    sens <- b$pooled$sensitivity_fgr %||% b$mi$sensitivity_fgr %||% NULL
    if (is.null(sens)) return(NULL)

    collect <- function(x, label) {
      if (is.null(x) || nrow(x) == 0) return(NULL)
      x[, Group := label]
      x[, Outcome := pretty_endpoint(nm)]
      x
    }

    rbind(
      collect(sens$fu, "Fully adjusted"),
      collect(sens$male, "Male"),
      collect(sens$female, "Female"),
      fill = TRUE
    )
  }), fill = TRUE)
  if (nrow(fgr_rows) > 0) fwrite(fgr_rows, file.path(out_table, "Table_FGR_Sensitivity.csv"))
}

write_interaction_table <- function(bundles, out_table) {
  int_rows <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    sx <- get_sex_bundle(b)
    p_int <- NA_real_
    if (!is.null(sx)) p_int <- sx$p_int %||% sx$p_interaction %||% NA_real_

    mods <- get_bundle_models(b)
    mod_int <- NULL
    if (!is.null(mods$csc) && !is.null(mods$csc$int)) mod_int <- mods$csc$int
    if (is.null(mod_int) && !is.null(mods$cox) && !is.null(mods$cox$int)) mod_int <- mods$cox$int
    if (is.null(mod_int) && !is.null(mods$int)) mod_int <- mods$int

    if (is.null(mod_int)) {
      data.table(Endpoint = pretty_endpoint(nm), Variable = "sex interaction", HR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = p_int)
    } else {
      ct <- coef_table_from_fit(mod_int)
      ct[, Endpoint := pretty_endpoint(nm)]
      ct[, .(Endpoint, Variable, HR, LCL, UCL, p)]
    }
  }), fill = TRUE)
  fwrite(int_rows, file.path(out_table, "Table_SexInteraction.csv"))
}

write_imputation_diagnostics_index <- function(out_table) {
  diag_files <- summarise_missingness_diagnostics(base_dir)
  idx <- data.table(
    Item = c("Imputation model", "Missingness table", "Predictor matrix", "Visit sequence", "Trace plots", "Density plots", "Strip plots"),
    File = c(diag_files$imputation_model, diag_files$missingness_table, diag_files$predictor_matrix, diag_files$visit_sequence, diag_files$trace_pdf, diag_files$density_pdf, diag_files$strip_pdf)
  )
  fwrite(idx, file.path(out_table, "Table_Imputation_Diagnostics_Index.csv"))
  idx
}

############################################################
# 6) WRITING TABLES
############################################################
msg("Generating tables...")

# 6A. Hazard ratio tables for the main population
for (adj in c("un", "de", "fu")) {
  tab_list <- lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    mods <- get_bundle_models(b)
    mod <- NULL
    if (!is.null(mods$csc) && !is.null(mods$csc[[adj]])) mod <- mods$csc[[adj]]
    if (is.null(mod) && !is.null(mods$cox) && !is.null(mods$cox[[adj]])) mod <- mods$cox[[adj]]
    if (is.null(mod) && !is.null(mods[[adj]])) mod <- mods[[adj]]
    if (is.null(mod)) return(NULL)

    ct <- coef_table_from_fit(mod)
    ct[, Est := fmt_hr(HR, LCL, UCL, p)]
    ct[, .(Variable, Est)][, setNames(.SD, nm)]
  })
  tab_list <- Filter(Negate(is.null), tab_list)
  if (length(tab_list) > 0) {
    wide <- Reduce(function(x, y) merge(x, y, by = "Variable", all = TRUE), tab_list)
    fwrite(wide, file.path(dir_table, paste0("Table_HR_", adj, ".csv")))
  }
}

# 6B. Fine-Gray sensitivity tables
write_fgr_sensitivity_tables(bundles, dir_table)

# 6C. Absolute risk tables (pooled ATEs from analysis pipeline)
for (adj in c("un", "de", "fu")) {
  out_all <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    if (is.null(b$ate[[adj]])) return(NULL)
    m <- standardize_ate_tbl(b$ate[[adj]]$mean)
    if (!all(c("time", "estimate") %in% names(m))) return(NULL)

    trt_col <- intersect(names(m), c("treatment", "acr", "level", "acrlevel"))
    trt <- if (length(trt_col) > 0) as.character(m[[trt_col[1]]]) else NA_character_

    data.table(
      Outcome = pretty_endpoint(nm),
      ACR_Level = trt,
      Time = m$time,
      Estimate = m$estimate,
      LCL = if ("lcl" %in% names(m)) m$lcl else NA_real_,
      UCL = if ("ucl" %in% names(m)) m$ucl else NA_real_
    )
  }), fill = TRUE)

  if (nrow(out_all) > 0) {
    out_all[, RiskPct := Estimate * 100]
    out_all[, LCLPct := LCL * 100]
    out_all[, UCLPct := UCL * 100]
    fwrite(out_all, file.path(dir_table, paste0("Table_AbsRisk_", adj, ".csv")))
  }
}

# 6D. Risk difference and risk ratio tables (fully adjusted only)
for (metric in c("diff", "ratio")) {
  out_all <- rbindlist(lapply(names(bundles), function(nm) {
    b <- bundles[[nm]]
    if (is.null(b$ate$fu)) return(NULL)
    m <- standardize_ate_tbl(b$ate$fu[[metric]])
    if (!all(c("time", "estimate") %in% names(m))) return(NULL)

    data.table(
      Outcome = pretty_endpoint(nm),
      Reference = if ("reference" %in% names(m)) as.character(m$reference) else NA_character_,
      Comparison = if ("comparison" %in% names(m)) as.character(m$comparison) else NA_character_,
      Time = m$time,
      Estimate = m$estimate,
      LCL = if ("lcl" %in% names(m)) m$lcl else NA_real_,
      UCL = if ("ucl" %in% names(m)) m$ucl else NA_real_
    )
  }), fill = TRUE)

  if (nrow(out_all) > 0) fwrite(out_all, file.path(dir_table, paste0("Table_", metric, "_Full.csv")))
}

# 6E. Sex-specific ATE table
write_sex_ate_tables(bundles, dir_table)

# 6F. Sex interaction table
write_interaction_table(bundles, dir_table)

# 6G. Performance, E-values, imputation QC index
write_performance_tables(bundles, dir_table)
write_evalues_table(bundles, dir_table)
imp_diag_idx <- write_imputation_diagnostics_index(dir_table)

############################################################
# 7) FIGURES
############################################################
msg("Generating figures...")
plot_data <- list()

# 7A. Grouped forest plots for A2 and A3 (log scale)
forest_out <- write_grouped_forest_outputs(
  bundles = bundles,
  out_table = dir_table,
  out_plot = dir_plot,
  terms = c("acr.levelA2", "acr.levelA3")
)

# 7B. Unadjusted CIF, adjusted curves, calibration, splines
for (nm in names(bundles)) {
  b <- bundles[[nm]]
  ep_lbl <- pretty_endpoint(nm)

  # Unadjusted CIF curves (competing-risk endpoints only)
  if (!is.null(b$cuminc)) {
    df_ci <- cuminc_to_df(b$cuminc$obj, times = seq(0.05, 5, 0.05))
    df_ci[, Risk := est * 100]
    df_ci[, LCL := lcl * 100]
    df_ci[, UCL := ucl * 100]

    p1 <- ggplot(df_ci, aes(x = time, y = Risk, color = ACR, fill = ACR)) +
      geom_step(linewidth = 1.0) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.12, color = NA) +
      scale_color_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
      scale_fill_manual(values = c(A1 = COL$A1, A2 = COL$A2, A3 = COL$A3)) +
      labs(
        title = paste0("Unadjusted cumulative incidence: ", ep_lbl),
        x = "Years from diagnosis",
        y = "Cumulative incidence (%)",
        color = "uACR",
        fill = "uACR"
      ) + thm

    ggsave(file.path(dir_plot, paste0("Unadjusted_CIF_", nm, ".png")), p1, width = 8, height = 6, dpi = 300)
    plot_data[[paste0(nm, "_cif")]] <- df_ci
  }

  # Adjusted risk curves: whole cohort and sex-stratified
  write_adjusted_curve_outputs(bundle = b, out_table = dir_table, out_plot = dir_plot, c_times = seq(0.05, 5, 0.05))

  # Calibration: fully adjusted model only (recomputed from stored fitted model)
  mods <- get_bundle_models(b)
  fit_full <- NULL
  if (!is.null(mods$csc) && !is.null(mods$csc$fu)) fit_full <- mods$csc$fu
  if (is.null(fit_full) && !is.null(mods$cox) && !is.null(mods$cox$fu)) fit_full <- mods$cox$fu
  if (is.null(fit_full) && !is.null(mods$fu)) fit_full <- mods$fu

  if (!is.null(fit_full)) {
    score_obj <- run_calibration(
      model = fit_full,
      data = b$data,
      time_var = "t",
      status_var = "s",
      cause = if (b$meta$type == "competing") b$meta$cause else NA_integer_,
      times = c(1, 3, 5)
    )
    save(score_obj, file = file.path(dir_score, paste0("Score_", nm, ".RData")))

    for (tt in c(1, 3, 5)) {
      p_cal <- plotCalibration(
        score_obj,
        models = "full",
        times = tt,
        method = "nne",
        q = 10,
        bars = FALSE,
        percent = TRUE
      ) +
        ggtitle(paste0("Calibration: ", ep_lbl, " at ", tt, "-year horizon")) + thm

      ggsave(file.path(dir_plot, paste0("Calibration_", nm, "_", tt, "y.png")),
             p_cal, width = 7.5, height = 5.5, dpi = 300)
    }
    plot_data[[paste0(nm, "_calibration")]] <- score_obj
  }

  # Spline-based 5-year absolute risk
  write_spline_outputs(bundle = b, out_table = dir_table, out_plot = dir_plot, B = 100)
}

save(plot_data, file = file.path(base_dir, "plot_data.RData"))

############################################################
# 8) WORD REPORT
############################################################
msg("Compiling Word report...")

doc <- read_docx()
doc <- body_add_par(doc, "ACR Prognosis Analysis Report", style = "heading 1")
doc <- body_add_par(doc, paste("Generated:", Sys.time()), style = "Normal")
doc <- body_add_par(doc, "", style = "Normal")

# Imputation diagnostics section
if (file.exists(file.path(dir_table, "Table_Imputation_Diagnostics_Index.csv"))) {
  doc <- body_add_par(doc, "Imputation diagnostics", style = "heading 2")
  doc <- body_add_par(doc, "Trace, density, and strip plots were generated after MICE. The imputation model, missingness table, predictor matrix, and visit sequence are indexed below.", style = "Normal")
  doc <- body_add_flextable(doc, flextable(read.csv(file.path(dir_table, "Table_Imputation_Diagnostics_Index.csv"))) %>% autofit())
}

# Performance and E-values
perf_file <- file.path(dir_table, "Table_Performance_Summary.csv")
if (file.exists(perf_file)) {
  doc <- body_add_par(doc, "Model performance summary", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(perf_file)) %>% autofit())
}

ev_file <- file.path(dir_table, "Table_EValues.csv")
if (file.exists(ev_file)) {
  doc <- body_add_par(doc, "E-values", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(ev_file)) %>% autofit())
}

# Main tables
for (adj in c("un", "de", "fu")) {
  hr_file <- file.path(dir_table, paste0("Table_HR_", adj, ".csv"))
  if (file.exists(hr_file)) {
    doc <- body_add_par(doc, paste0("Hazard ratios (", adj, ")"), style = "heading 2")
    doc <- body_add_flextable(doc, flextable(read.csv(hr_file)) %>% autofit())
  }

  ar_file <- file.path(dir_table, paste0("Table_AbsRisk_", adj, ".csv"))
  if (file.exists(ar_file)) {
    doc <- body_add_par(doc, paste0("Absolute risks (", adj, ")"), style = "heading 2")
    doc <- body_add_flextable(doc, flextable(read.csv(ar_file)) %>% autofit())
  }
}

# Fine-Gray sensitivity
fgr_file <- file.path(dir_table, "Table_FGR_Sensitivity.csv")
if (file.exists(fgr_file)) {
  doc <- body_add_par(doc, "Fine-Gray sensitivity", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(fgr_file)) %>% autofit())
}

# Risk differences and ratios
rd_file <- file.path(dir_table, "Table_diff_Full.csv")
rr_file <- file.path(dir_table, "Table_ratio_Full.csv")
if (file.exists(rd_file)) {
  doc <- body_add_par(doc, "Risk differences (fully adjusted)", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(rd_file)) %>% autofit())
}
if (file.exists(rr_file)) {
  doc <- body_add_par(doc, "Risk ratios (fully adjusted)", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(rr_file)) %>% autofit())
}

sex_ate_file <- file.path(dir_table, "Table_SexStratified_ATE.csv")
if (file.exists(sex_ate_file)) {
  doc <- body_add_par(doc, "Sex-stratified ATE outputs", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(sex_ate_file)) %>% autofit())
}

int_file <- file.path(dir_table, "Table_SexInteraction.csv")
if (file.exists(int_file)) {
  doc <- body_add_par(doc, "Sex interaction coefficients", style = "heading 2")
  doc <- body_add_flextable(doc, flextable(read.csv(int_file)) %>% autofit())
}

# Forest plots (A2 and A3) and key figure panels
for (term in c("acrlevelA2", "acrlevelA3")) {
  main_f <- file.path(dir_plot, paste0("Forest_MainPopulation_", term, ".png"))
  sex_f  <- file.path(dir_plot, paste0("Forest_SexStratified_", term, ".png"))
  if (file.exists(main_f)) {
    doc <- body_add_par(doc, paste0("Main population forest plot: ", term), style = "heading 2")
    doc <- body_add_img(doc, src = main_f, width = 6.9, height = 4.8)
  }
  if (file.exists(sex_f)) {
    doc <- body_add_par(doc, paste0("Sex-stratified forest plot: ", term), style = "heading 2")
    doc <- body_add_img(doc, src = sex_f, width = 6.9, height = 4.8)
  }
}

for (nm in names(bundles)) {
  ep_lbl <- pretty_endpoint(nm)
  adj_fig <- file.path(dir_plot, paste0("AdjustedCurve_Overall_", nm, ".png"))
  sex_adj_fig <- file.path(dir_plot, paste0("AdjustedCurve_SexStratified_", nm, ".png"))
  cal_fig <- file.path(dir_plot, paste0("Calibration_", nm, "_5y.png"))
  spl_fig <- file.path(dir_plot, paste0("Spline_Overall_", nm, ".png"))
  spl_sex_fig <- file.path(dir_plot, paste0("Spline_SexStratified_", nm, ".png"))
  cif_fig <- file.path(dir_plot, paste0("Unadjusted_CIF_", nm, ".png"))

  if (file.exists(cif_fig)) {
    doc <- body_add_par(doc, paste("Unadjusted CIF:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = cif_fig, width = 6.5, height = 4.6)
  }
  if (file.exists(adj_fig)) {
    doc <- body_add_par(doc, paste("Fully adjusted risk curve:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = adj_fig, width = 6.5, height = 4.6)
  }
  if (file.exists(sex_adj_fig)) {
    doc <- body_add_par(doc, paste("Sex-stratified adjusted risk curve:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = sex_adj_fig, width = 6.5, height = 4.6)
  }
  if (file.exists(cal_fig)) {
    doc <- body_add_par(doc, paste("Calibration:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = cal_fig, width = 6.5, height = 4.6)
  }
  if (file.exists(spl_fig)) {
    doc <- body_add_par(doc, paste("Spline risk curve:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = spl_fig, width = 6.5, height = 4.6)
  }
  if (file.exists(spl_sex_fig)) {
    doc <- body_add_par(doc, paste("Sex-stratified spline risk curve:", ep_lbl), style = "heading 2")
    doc <- body_add_img(doc, src = spl_sex_fig, width = 6.5, height = 4.6)
  }
}

# Final save
print(doc, target = file.path(dir_report, "ACR_Prognosis_Report.docx"))

msg("Result pipeline complete")
msg("Tables written to: ", dir_table)
msg("Figures written to: ", dir_plot)
msg("Report written to: ", file.path(dir_report, "ACR_Prognosis_Report.docx"))
