############################################################
# FINAL ANALYSIS PIPELINE
# Study: Baseline Albuminuria and Systemic Vascular Risk
# Design: Incident type 2 diabetes cohort (CPRD-GOLD linked to HES & ONS)
#
# PURPOSE
#   1) Build endpoint-specific analytic datasets
#   2) Handle missing baseline biomarkers with tailored multiple imputation
#   3) Fit cause-specific Cox models for competing-risk endpoints
#   4) Fit Cox PH models for mortality / survival composites
#   5) Fit Fine-Gray only as a sensitivity analysis for competing-risk endpoints
#   6) Estimate adjusted absolute risks, risk differences, and risk ratios
#   7) Add sex-stratified models, spline models, discrimination metrics, and E-values
#   8) Save bundle objects required by the results pipeline
#
# METHODS
#   - CSC / cause-specific Cox + g-computation for the primary competing-risk analysis
#   - Cox PH for survival endpoints
#   - Fine-Gray retained only as a sensitivity analysis
#   - MICE with biomarker-specific predictor matrices
#   - B = 0 for ATE (influence-function uncertainty; no bootstrap)
#
# Primary Endpoint
# └── ASCVD, HF, PAD
# Secondary Endpoints
# ├── Atrial Fibrillation
# ├── ASCVD + PAD
# └── All-cause Mortality


############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(prodlim)
  library(riskRegression)
  library(splines)
  library(survminer)
  library(lme4)
  library(mice)
  library(lattice)
})

############################################################
# 0) INIT / SAFETY CHECKS / OUTPUT PATHS
############################################################
if (!exists("dt")) {
  stop("Object 'dt' not found. Load your cohort data.table before running this script.")
}
setDT(dt)
set.seed(93277)

base_dir <- "acr_prognosis_final"
paths <- list(
  base = base_dir,
  bundles = file.path(base_dir, "bundles"),
  logs = file.path(base_dir, "logs"),
  figures = file.path(base_dir, "figures"),
  tables = file.path(base_dir, "tables")
)
lapply(paths, dir.create, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(paths$logs, "analysis_log.txt")
log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%F %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
  cat(msg, "
", file = log_file, append = TRUE)
}

log_msg("START: modular analysis pipeline initialised")
log_msg("Cohort size: ", nrow(dt))

############################################################
# 1) GLOBAL CONFIGURATION
############################################################
# Exposure
acr_cat  <- "acr.level"
acr_cont <- "acr"
acr_lev  <- c("A1", "A2", "A3")

# Time horizons (kept in sync with the results pipeline)
peak_times  <- c(1, 3, 5)
curve_times <- seq(0.05, 5, 0.05)

# Baseline biomarkers
# LDL has ~30% missingness; SBP/DBP are ~2% missing and can use simpler predictors.
biomarkers <- c("b.eGFR", "b.SBP", "b.DBP", "b.A1c", "b.LDL", "b.HDL", "b.bmi")

# Full adjustment set
cov_full <- c(
  "age", "sex", "iyear",
  "ane", "lipid", "ast", "copd", "cancer", "dem",
  "ckd", "lupus", "gout",
  biomarkers,
  "m.Sulfony", "m.Thiazide", "m.DPP4", "m.Metformin",
  "m.ARBs", "m.ACEi", "m.CCBlocker", "m.BetaB",
  "m.MRAs", "m.OAC", "m.Diuretics", "m.Antiplatelet",
  "m.Insulins", "m.Statins"
)
analysis_cov <- setdiff(cov_full, biomarkers)
cov_demo <- c("age", "sex")
cov_no_sex <- setdiff(cov_full, "sex")
 
# Endpoint map (single source of truth)
# CVD+P is treated as a competing-risk endpoint here.
eps <- list(
  ascvd    = list(t = "time_ascvd",    s = "status_ascvd",    type = "competing", cause = 1L, label = "ASCVD"),
  hf     = list(t = "time_hf",     s = "status_hf",     type = "competing", cause = 1L, label = "HF"),
  pad    = list(t = "time_pad",    s = "status_pad",    type = "competing", cause = 1L, label = "PAD"),
  acs    = list(t = "time_acs",    s = "status_acs",    type = "competing", cause = 1L, label = "ACS"),
  af     = list(t = "time_af",     s = "status_af",     type = "competing", cause = 1L, label = "AF"),
  stroke = list(t = "time_stroke", s = "status_stroke", type = "competing", cause = 1L, label = "Stroke"),
  CVDp  = list(t = "time_cvdp",  s = "status_cvdp",  type = "competing", cause = 1L, label = "CVD+P"),
  death  = list(t = "time_death",  s = "status_death",  type = "survival",  cause = NA_integer_, label = "Mortality")
)

# QC 1: source-data completeness
req_cols <- unique(c(
  unlist(lapply(eps, `[[`, "t")),
  unlist(lapply(eps, `[[`, "s")),
  acr_cat, acr_cont, cov_full
))
missing_cols <- setdiff(req_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}
log_msg("QC1 OK: all required columns are present")

############################################################
# 2) MODULAR HELPERS
############################################################
`%||%` <- function(x, y) if (is.null(x)) y else x

pretty_endpoint <- function(x) {
  labs <- c(
    cvd = "CVD", acs = "ACS", hf = "HF", af = "AF",
    stroke = "Stroke", pad = "PAD", CVDp = "CVD+P", death = "Mortality"
  )
  ifelse(x %in% names(labs), unname(labs[x]), x)
}

missingness_summary <- function(dat) {
  pct <- sort(sapply(dat, function(x) mean(is.na(x))) * 100, decreasing = TRUE)
  data.table(variable = names(pct), missing_pct = as.numeric(pct))
}

build_hist_formula <- function(time_var, status_var, rhs_str) {
  as.formula(paste0("Hist(", time_var, ", ", status_var, ") ~ ", rhs_str))
}

build_surv_formula <- function(time_var, status_var, rhs_str) {
  as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", rhs_str))
}

# Primary cause-specific Cox model used for the competing-risk analysis.
fit_csc <- function(d, time_var, status_var, rhs_str, cause) {
  fml <- build_hist_formula(time_var, status_var, rhs_str)
  CSC(fml, data = d, cause = cause)
}

# Cox model used for survival endpoints.
fit_cox <- function(d, time_var, status_var, rhs_str) {
  fml <- build_surv_formula(time_var, status_var, rhs_str)
  coxph(fml, data = d, x = TRUE, y = TRUE, model = TRUE)
}

# Fine-Gray sensitivity model (competing-risk endpoints only).
fit_fgr <- function(d, time_var, status_var, rhs_str, cause) {
  fml <- build_hist_formula(time_var, status_var, rhs_str)
  FGR(fml, data = d, cause = cause, y = TRUE)
}

# ATE with B = 0 (influence-function uncertainty; no bootstrap)
run_ate <- function(fit, d, treatment_var, times, cause = NA_integer_) {
  args <- list(
    event     = fit,
    treatment = treatment_var,
    data      = d,
    times     = times,
    se        = TRUE,
    band      = FALSE,
    B         = 0,
    estimator = "GFORMULA"
  )
  if (!is.na(cause)) args$cause <- cause
  a <- do.call(ate, args)
  list(
    obj   = a,
    mean  = model.tables(a, times = peak_times, type = "meanRisk"),
    diff  = model.tables(a, times = peak_times, type = "diffRisk"),
    ratio = model.tables(a, times = peak_times, type = "ratioRisk"),
    curve = model.tables(a, times = curve_times, type = "meanRisk"),
    ci    = confint(a, p.value = TRUE, band = FALSE)
  )
}

run_cuminc <- function(d, time_var, status_var) {
  ci <- cuminc(ftime = d[[time_var]], fstatus = d[[status_var]], group = d[[acr_cat]])
  list(obj = ci, tp = timepoints(ci, times = curve_times))
}

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
  p <- if ("Pr(>|z|)" %in% names(tab)) tab[["Pr(>|z|)"]] else rep(NA_real_, length(beta))

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

# Rubin pooling for scalar estimates.
pool_scalar_rubin <- function(est, se) {
  m <- length(est)
  qbar <- mean(est, na.rm = TRUE)
  ubar <- mean(se^2, na.rm = TRUE)
  b <- stats::var(est, na.rm = TRUE)
  tvar <- ubar + (1 + 1/m) * b
  list(
    est = qbar,
    se = sqrt(tvar),
    lcl = qbar - 1.96 * sqrt(tvar),
    ucl = qbar + 1.96 * sqrt(tvar)
  )
}

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

# E-value calculation for hazard ratios.
calc_evalue_hr <- function(hr, lcl) {
  hr_use <- ifelse(hr < 1, 1 / hr, hr)
  lcl_use <- ifelse(lcl < 1, 1 / lcl, lcl)
  point <- hr_use + sqrt(hr_use * (hr_use - 1))
  ci_e <- ifelse(lcl_use > 1, lcl_use + sqrt(lcl_use * (lcl_use - 1)), 1)
  data.table(evalue = point, evalue_ci = ci_e)
}

# Discrimination / calibration for the fully adjusted model.
# For MI, we compute the score on every imputation, then summarize mean/range
# and Rubin-pooled estimates.
run_performance <- function(fit, d, sp, label = "full") {
  if (sp$type == "competing") {
    sc <- tryCatch(
      Score(
        list(full = fit),
        formula  = as.formula("Hist(t, s) ~ 1"),
        data     = d,
        cause    = sp$cause,
        times    = peak_times,
        plots    = "Calibration",
        metrics  = c("AUC", "Brier"),
        conf.int = TRUE,
        se.fit   = TRUE
      ),
      error = function(e) NULL
    )
    list(score = sc, cindex = NA_real_, auc = sc, label = label)
  } else {
    sc <- tryCatch(
      Score(
        list(full = fit),
        formula  = as.formula("Surv(t, s) ~ 1"),
        data     = d,
        times    = peak_times,
        plots    = "Calibration",
        metrics  = c("AUC", "Brier"),
        conf.int = TRUE,
        se.fit   = TRUE
      ),
      error = function(e) NULL
    )
    cidx <- tryCatch(concordance(fit)$concordance, error = function(e) NA_real_)
    list(score = sc, cindex = cidx, auc = sc, label = label)
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

############################################################
# 3) IMPUTATION MODULES
############################################################
# Missingness summary for reporting and QC.
missingness_table <- function(dat) {
  miss <- missingness_summary(dat)
  miss[, rank := .I]
  miss
}

# Build auxiliary variables for imputation.
build_imputation_data <- function(d) {
  d <- copy(d)
  d[, any_event := as.integer(
    (status_cvd == 1L) | (status_acs == 1L) | (status_hf == 1L) |
      (status_af == 1L) | (status_stroke == 1L) | (status_pad == 1L) |
      (status_cvdp == 1L) | (status_death == 1L)
  )]
  d[, any_time := pmin(
    time_cvd, time_acs, time_hf, time_af,
    time_stroke, time_pad, time_cvdp, time_death,
    na.rm = TRUE
  )]
  d[, primary_event := as.integer(status_cvd == 1L)]
  d[, primary_time := time_cvd]
  d
}

get_existing <- function(x, pool) intersect(x, pool)

# Biomarker-specific predictor matrices.
# LDL gets the richest predictor set because missingness is highest.
# SBP/DBP use narrower, clinically sensible predictor sets because missingness is low.
build_mice_spec <- function(imp_dat) {
  base_preds <- c(
    "age", "sex", "iyear", acr_cat,
    "any_event", "any_time", "primary_event", "primary_time",
    "ane", "lipid", "ast", "copd", "cancer", "dem",
    "ckd", "lupus", "gout"
  )

  biomarker_pred_map <- list(
    "b.LDL" = c(base_preds,
      "b.HDL", "b.A1c", "b.bmi", "b.eGFR", "b.SBP", "b.DBP",
      "m.Statins", "m.Sulfony", "m.Metformin", "m.DPP4", "m.Insulins",
      "m.ARBs", "m.ACEi", "m.Diuretics"
    ),
    "b.HDL" = c(base_preds,
      "b.LDL", "b.A1c", "b.bmi", "b.eGFR", "b.SBP", "b.DBP",
      "m.Statins", "m.Sulfony", "m.Metformin", "m.Insulins"
    ),
    "b.A1c" = c(base_preds,
      "b.bmi", "b.eGFR", "b.SBP", "b.DBP",
      "m.Sulfony", "m.Metformin", "m.DPP4", "m.Insulins"
    ),
    "b.eGFR" = c(base_preds,
      "ht", "b.SBP", "b.DBP", "b.A1c", "b.bmi",
      "m.ACEi", "m.ARBs", "m.Diuretics", "m.BetaB", "m.CCBlocker", "m.MRAs",
      "m.Statins"
    ),
    "b.SBP" = c(base_preds,
      "ht", "ckd", "b.DBP", "b.bmi", "b.eGFR",
      "m.ACEi", "m.ARBs", "m.Diuretics", "m.BetaB", "m.CCBlocker", "m.MRAs"
    ),
    "b.DBP" = c(base_preds,
      "ht", "ckd", "b.SBP", "b.bmi", "b.eGFR",
      "m.ACEi", "m.ARBs", "m.Diuretics", "m.BetaB", "m.CCBlocker", "m.MRAs"
    ),
    "b.bmi" = c(base_preds,
      "b.A1c", "b.SBP", "b.DBP", "b.eGFR",
      "m.Metformin", "m.Sulfony", "m.DPP4", "m.Insulins", "m.Statins"
    )
  )

  meth <- make.method(imp_dat)
  for (v in biomarkers) meth[v] <- "pmm"

  structural <- c(acr_cat, acr_cont, "sex", "iyear", "any_event", "any_time", "primary_event", "primary_time")
  meth[structural] <- ""

  outcome_fields <- unique(c(unlist(lapply(eps, `[[`, "t")), unlist(lapply(eps, `[[`, "s"))))
  meth[get_existing(outcome_fields, names(imp_dat))] <- ""

  pred <- matrix(0, nrow = ncol(imp_dat), ncol = ncol(imp_dat),
                 dimnames = list(names(imp_dat), names(imp_dat)))

  for (v in biomarkers) {
    preds <- get_existing(biomarker_pred_map[[v]] %||% base_preds, names(imp_dat))
    pred[v, preds] <- 1
    pred[v, v] <- 0
  }

  pred[structural, ] <- 0
  pred[get_existing(outcome_fields, names(imp_dat)), ] <- 0

  miss <- sort(sapply(imp_dat[, ..biomarkers], function(x) mean(is.na(x))), decreasing = FALSE)
  visitSequence <- names(miss)

  list(method = meth, predictorMatrix = pred, visitSequence = visitSequence)
}

run_imputation_qc_outputs <- function(imp, imp_dat, mice_spec, out_base, biomarkers) {
  log_msg("IMPUTATION QC: writing diagnostics")

  save(imp, file = file.path(out_base, "imputation_model.RData"))
  fwrite(missingness_table(imp_dat), file.path(out_base, "tables", "Table_Missingness.csv"))
  fwrite(as.data.table(mice_spec$predictorMatrix, keep.rownames = "variable"),
         file.path(out_base, "tables", "Imputation_PredictorMatrix.csv"))
  writeLines(mice_spec$visitSequence, con = file.path(out_base, "tables", "Imputation_VisitSequence.txt"))

  pdf(file.path(out_base, "figures", "Imputation_TracePlots.pdf"), width = 10, height = 8)
  plot(imp)
  dev.off()

  pdf(file.path(out_base, "figures", "Imputation_DensityPlots.pdf"), width = 10, height = 8)
  for (v in biomarkers) {
    try(print(densityplot(imp, as.formula(paste0("~", v)), main = paste("Observed vs imputed:", v))), silent = TRUE)
  }
  dev.off()

  pdf(file.path(out_base, "figures", "Imputation_StripPlots.pdf"), width = 10, height = 8)
  for (v in biomarkers) {
    try(print(stripplot(imp, as.formula(paste0(v, " ~ .imp")), pch = 20, cex = 0.4,
                         main = paste("Imputed values by draw:", v))), silent = TRUE)
  }
  dev.off()

  log_msg("IMPUTATION QC COMPLETE: trace, density, strip plots and tables saved")
}

############################################################
# 4) MODEL / POOLING MODULES
############################################################
make_spline_knots <- function(d) {
  obs_log_acr <- log(d[[acr_cont]] + 0.01)
  k <- quantile(obs_log_acr, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
  b <- quantile(obs_log_acr, probs = c(0.05, 0.95), na.rm = TRUE)
  sprintf(
    "splines::ns(log_acr, knots = c(%s), Boundary.knots = c(%s))",
    paste(format(k, 8), collapse = ","),
    paste(format(b, 8), collapse = ",")
  )
}

build_rhs <- function(spl_term) {
  list(
    un  = "acr.level",
    de  = paste("acr.level", paste(cov_demo, collapse = " + "), sep = " + "),
    fu  = paste("acr.level", paste(cov_full, collapse = " + "), sep = " + "),
    int = paste("acr.level * sex", paste(cov_no_sex, collapse = " + "), sep = " + "),
    spl = paste(spl_term, paste(cov_full, collapse = " + "), sep = " + "),
    spl_sex = paste(sprintf("%s * sex", spl_term), paste(cov_no_sex, collapse = " + "), sep = " + ")
  )
}

fit_primary_family <- function(d, sp, rhs) {
  if (sp$type == "competing") {
    list(
      un  = fit_csc(d, "t", "s", rhs$un,  sp$cause),
      de  = fit_csc(d, "t", "s", rhs$de,  sp$cause),
      fu  = fit_csc(d, "t", "s", rhs$fu,  sp$cause),
      int = fit_csc(d, "t", "s", rhs$int, sp$cause),
      spl = fit_csc(d, "t", "s", rhs$spl, sp$cause)
    )
  } else {
    list(
      un  = fit_cox(d, "t", "s", rhs$un),
      de  = fit_cox(d, "t", "s", rhs$de),
      fu  = fit_cox(d, "t", "s", rhs$fu),
      int = fit_cox(d, "t", "s", rhs$int),
      spl = fit_cox(d, "t", "s", rhs$spl)
    )
  }
}

fit_sex_family <- function(d, sp, rhs) {
  d_m <- d[sex == "Male"]
  d_f <- d[sex == "Female"]
  if (sp$type == "competing") {
    list(
      strat_models = list(
        male   = fit_csc(d_m, "t", "s", rhs$fu, sp$cause),
        female = fit_csc(d_f, "t", "s", rhs$fu, sp$cause)
      ),
      spline_model = fit_csc(d, "t", "s", rhs$spl_sex, sp$cause)
    )
  } else {
    list(
      strat_models = list(
        male   = fit_cox(d_m, "t", "s", rhs$fu),
        female = fit_cox(d_f, "t", "s", rhs$fu)
      ),
      spline_model = fit_cox(d, "t", "s", rhs$spl_sex)
    )
  }
}

fit_fgr_sensitivity <- function(d, sp, rhs) {
  if (sp$type != "competing") return(NULL)
  d_m <- d[sex == "Male"]
  d_f <- d[sex == "Female"]
  list(
    fu = fit_fgr(d,   "t", "s", rhs$fu, sp$cause),
    male = fit_fgr(d_m, "t", "s", rhs$fu, sp$cause),
    female = fit_fgr(d_f, "t", "s", rhs$fu, sp$cause)
  )
}

run_ate_family <- function(models, d, sp) {
  if (sp$type == "competing") {
    list(
      un = run_ate(models$un, d, acr_cat, peak_times, sp$cause),
      de = run_ate(models$de, d, acr_cat, peak_times, sp$cause),
      fu = run_ate(models$fu, d, acr_cat, peak_times, sp$cause)
    )
  } else {
    list(
      un = run_ate(models$un, d, acr_cat, peak_times, NA_integer_),
      de = run_ate(models$de, d, acr_cat, peak_times, NA_integer_),
      fu = run_ate(models$fu, d, acr_cat, peak_times, NA_integer_)
    )
  }
}

run_sex_ate_family <- function(sex_models, d, sp) {
  if (sp$type == "competing") {
    list(
      male   = run_ate(sex_models$strat_models$male,   d[sex == "Male"],   acr_cat, peak_times, sp$cause),
      female = run_ate(sex_models$strat_models$female, d[sex == "Female"], acr_cat, peak_times, sp$cause)
    )
  } else {
    list(
      male   = run_ate(sex_models$strat_models$male,   d[sex == "Male"],   acr_cat, peak_times, NA_integer_),
      female = run_ate(sex_models$strat_models$female, d[sex == "Female"], acr_cat, peak_times, NA_integer_)
    )
  }
}

pool_ate_component <- function(ate_list, component) {
  rows <- lapply(seq_along(ate_list), function(i) {
    x <- standardize_ate_tbl(ate_list[[i]][[component]])
    if (nrow(x) == 0) return(NULL)
    x[, imp := i]
    x
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.table())
  dat <- rbindlist(rows, fill = TRUE)
  if (!all(c("estimate", "se") %in% names(dat))) return(dat)

  grp <- intersect(names(dat), c("time", "treatment", "reference", "comparison"))
  if (length(grp) == 0) grp <- "time"

  dat[, `:=`(
    estimate = as.numeric(estimate),
    se = as.numeric(se)
  )]

  dat[, {
    pooled <- pool_scalar_rubin(estimate, ifelse(is.na(se), 0, se))
    .(estimate = pooled$est, se = pooled$se, lcl = pooled$lcl, ucl = pooled$ucl)
  }, by = grp]
}

pool_ate_family <- function(ate_list) {
  list(
    mean  = pool_ate_component(ate_list, "mean"),
    diff  = pool_ate_component(ate_list, "diff"),
    ratio = pool_ate_component(ate_list, "ratio"),
    curve = pool_ate_component(ate_list, "curve")
  )
}

pool_coef_list <- function(fits, slot = NULL) {
  rows <- lapply(seq_along(fits), function(i) {
    mod <- if (is.null(slot)) fits[[i]] else fits[[i]][[slot]]
    if (is.null(mod)) return(NULL)
    x <- coef_table_from_fit(mod)
    if (nrow(x) == 0) return(NULL)
    x[, imp := i]
    x
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(NULL)
  dat <- rbindlist(rows, fill = TRUE)
  dat[, {
    pooled <- pool_scalar_rubin(beta, se)
    z <- pooled$est / pooled$se
    p <- 2 * pnorm(-abs(z))
    .(beta = pooled$est, se = pooled$se, HR = exp(pooled$est), LCL = exp(pooled$lcl), UCL = exp(pooled$ucl), p = p)
  }, by = Variable]
}

pool_interaction_table <- function(fit_list) {
  # Rubin-pooled interaction coefficients from the interaction model.
  tbl <- pool_coef_list(fit_list, slot = "int")
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(list(table = NULL, global_p = NA_real_))
  }

  interaction_terms <- grepl(":sex", tbl$Variable, fixed = TRUE) | grepl("sex", tbl$Variable)
  interaction_terms <- interaction_terms & grepl("acr\.level", tbl$Variable)
  if (!any(interaction_terms)) {
    return(list(table = tbl, global_p = NA_real_))
  }

  wald <- sum((tbl$beta[interaction_terms] / tbl$se[interaction_terms])^2, na.rm = TRUE)
  global_p <- pchisq(wald, df = sum(interaction_terms), lower.tail = FALSE)
  list(table = tbl, global_p = global_p)
}

make_spline_grid <- function(d) {
  seq(
    quantile(d[[acr_cont]], 0.02, na.rm = TRUE),
    quantile(d[[acr_cont]], 0.98, na.rm = TRUE),
    length.out = 60
  )
}

calc_evalues_table <- function(coef_tbl) {
  if (is.null(coef_tbl) || nrow(coef_tbl) == 0) return(NULL)
  rbindlist(lapply(seq_len(nrow(coef_tbl)), function(i) {
    x <- coef_tbl[i]
    e <- calc_evalue_hr(x$HR, x$LCL)
    data.table(Variable = x$Variable, evalue = e$evalue, evalue_ci = e$evalue_ci)
  }))
}

############################################################
# 5) IMPUTATION INPUTS
############################################################
imp_source <- build_imputation_data(copy(dt))

imp_vars <- unique(c(
  acr_cat, acr_cont, "sex", "iyear",
  biomarkers,
  cov_full,
  "any_event", "any_time", "primary_event", "primary_time",
  unlist(lapply(eps, `[[`, "t")),
  unlist(lapply(eps, `[[`, "s"))
))
imp_dat <- copy(imp_source[, ..imp_vars])
miss_tbl <- missingness_summary(imp_dat)
log_msg("Top missingness variables: ", paste(head(miss_tbl$variable, 8), collapse = ", "))

m_imp <- 30L
maxit_imp <- 20L
mice_spec <- build_mice_spec(imp_dat)

############################################################
# 6) IMPUTATION
############################################################
log_msg("START: MICE imputation with m=", m_imp, ", maxit=", maxit_imp)
imp <- mice(
  imp_dat,
  m = m_imp,
  maxit = maxit_imp,
  method = mice_spec$method,
  predictorMatrix = mice_spec$predictorMatrix,
  visitSequence = mice_spec$visitSequence,
  seed = 93277,
  printFlag = FALSE
)
log_msg("Imputation complete")

# Imputation diagnostics and audit trail.
run_imputation_qc_outputs(imp, imp_dat, mice_spec, base_dir, biomarkers)

############################################################
# 7) PER-ENDPOINT PROCESSOR
############################################################
spl_term <- make_spline_knots(d0)
rhs <- build_rhs(spl_term)
process_endpoint <- function(nm, sp, imp) {
  log_msg("PROCESSING endpoint: ", nm, " (", sp$type, ")")

  analysis_cols <- unique(c(
    sp$t, sp$s, acr_cat, acr_cont, "sex", "iyear",
    "any_event", "any_time", "primary_event", "primary_time",
    biomarkers, analysis_cov
  ))

  # First completed imputation used as the compatibility dataset for the results pipeline.
  d0 <- as.data.table(complete(imp, action = 1))[, ..analysis_cols]
  setnames(d0, analysis_cols)
  d0 <- d0[!is.na(t) & !is.na(s) & !is.na(get(acr_cat))]
  d0[, sex := factor(sex)]
  d0[, (acr_cat) := factor(get(acr_cat), levels = acr_lev)]
  d0[, (acr_cat) := relevel(get(acr_cat), ref = "A1")]
  d0[, iyear := factor(iyear)]
  d0[, log_acr := log(get(acr_cont) + 0.01)]

  qc <- list(
    n = nrow(d0),
    status_values = sort(unique(na.omit(d0$s))),
    positive_time = sum(d0$t > 0, na.rm = TRUE)
  )
  if (qc$n == 0L) stop("No analysis rows left after structural QC for endpoint: ", nm)
  log_msg("QC2 OK: n=", qc$n, ", status codes=", paste(qc$status_values, collapse = ","))

  # Build completed datasets for all imputations.
  comp_list <- lapply(seq_len(imp$m), function(m) {
    dd <- as.data.table(complete(imp, action = m))[, ..analysis_cols]
    setnames(dd, analysis_cols)
    dd <- dd[!is.na(t) & !is.na(s) & !is.na(get(acr_cat))]
    dd[, sex := factor(sex)]
    dd[, (acr_cat) := factor(get(acr_cat), levels = acr_lev)]
    dd[, (acr_cat) := relevel(get(acr_cat), ref = "A1")]
    dd[, iyear := factor(iyear)]
    dd[, log_acr := log(get(acr_cont) + 0.01)]
    dd
  })

  if (any(vapply(comp_list, function(x) anyNA(x[[acr_cat]]) || anyNA(x$sex) || anyNA(x$iyear), logical(1)))) {
    stop("Structural NA introduced after imputation for endpoint: ", nm)
  }

  # Fit model families and compute ATEs on each imputation.
  fit_list <- vector("list", length(comp_list))
  sex_fit_list <- vector("list", length(comp_list))
  fgr_sens_list <- vector("list", length(comp_list))
  ate_list <- vector("list", length(comp_list))
  sex_ate_list <- vector("list", length(comp_list))
  cuminc_list <- vector("list", length(comp_list))
  perf_list <- vector("list", length(comp_list))

  for (m in seq_along(comp_list)) {
    d <- comp_list[[m]]
    fit_list[[m]] <- fit_primary_family(d, sp, rhs)
    sex_fit_list[[m]] <- fit_sex_family(d, sp, rhs)
    ate_list[[m]] <- run_ate_family(fit_list[[m]], d, sp)
    sex_ate_list[[m]] <- run_sex_ate_family(sex_fit_list[[m]], d, sp)

    perf_raw <- run_performance(fit_list[[m]]$fu, d, sp, label = paste0(nm, "_imp", m))
    perf_list[[m]] <- extract_performance_summary(perf_raw$score, model_name = "full")

    if (sp$type == "competing") {
      cuminc_list[[m]] <- run_cuminc(d, "t", "s")
      fgr_sens_list[[m]] <- fit_fgr_sensitivity(d, sp, rhs)
    }

    if (m %% 5 == 0) log_msg("  completed imputation ", m, "/", length(comp_list), " for endpoint ", nm)
  }

  # --- REQUIRED CHANGE: Pool each ATE adjustment level separately ---
  pooled_ate_un <- pool_ate_family(lapply(ate_list, `[[`, "un"))
  pooled_ate_de <- pool_ate_family(lapply(ate_list, `[[`, "de"))
  pooled_ate_fu <- pool_ate_family(lapply(ate_list, `[[`, "fu"))
  
  pooled_coef <- pool_coef_list(fit_list, slot = "fu")
  pooled_interaction <- pool_interaction_table(fit_list)
  pooled_global_p_int <- pooled_interaction$global_p
  pooled_evalues <- calc_evalues_table(pooled_coef)
  pooled_sex_ate <- list(
    male   = pool_ate_family(lapply(sex_ate_list, `[[`, "male")),
    female = pool_ate_family(lapply(sex_ate_list, `[[`, "female"))
  )
  performance_summary <- pool_performance_metrics(perf_list)

  # First-imputation compatibility objects for the results pipeline.
  first_fit <- fit_list[[1]]
  first_sex <- sex_fit_list[[1]]
  cuminc_1 <- if (sp$type == "competing") cuminc_list[[1]] else NULL
  spline_grid <- make_spline_grid(comp_list[[1]])

  # Fine-Gray sensitivity (competing-risk endpoints only).
  fgr_pooled <- if (sp$type == "competing") list(
    fu = pool_coef_list(fgr_sens_list, slot = "fu"),
    male = pool_coef_list(fgr_sens_list, slot = "male"),
    female = pool_coef_list(fgr_sens_list, slot = "female")
  ) else NULL

  # Final bundle
  bundle <- list(
    meta = list(
      ep = nm, label = sp$label, t = "t", s = "s", type = sp$type, cause = sp$cause,
      acr_cat = acr_cat, acr_cont = acr_cont, acr_levels = acr_lev,
      peak_times = peak_times, curve_times = curve_times,
      t_times = peak_times, c_times = curve_times, qc = qc,
      imputation = list(m = imp$m, maxit = maxit_imp, method = "MICE-PMM")
    ),
    data = comp_list[[1]],
    mi = list(
      miss_pct = miss_tbl,
      pooled_coef = pooled_coef,
      pooled_interaction = pooled_interaction$table,
      pooled_global_p_int = pooled_global_p_int,
      pooled_mean = pooled_ate_fu$mean,
      pooled_diff = pooled_ate_fu$diff,
      pooled_ratio = pooled_ate_fu$ratio,
      pooled_curve = pooled_ate_fu$curve,
      pooled_sex_ate = pooled_sex_ate,
      pooled_evalues = pooled_evalues,
      performance_raw = perf_list,
      performance_summary = performance_summary,
      sensitivity_fgr = fgr_pooled
    ),
    models = if (sp$type == "competing") {
      list(csc = first_fit, fgr = fgr_sens_list[[1]], strat_models = first_sex$strat_models, spline = first_sex$spline_model)
    } else {
      list(cox = first_fit, strat_models = first_sex$strat_models, spline = first_sex$spline_model)
    },
    # --- REQUIRED CHANGE: Replace single-imputation list with pooled results ---
    ate = list(
      un = pooled_ate_un,
      de = pooled_ate_de,
      fu = pooled_ate_fu
    ),
    p_int = pooled_global_p_int,
    cuminc = cuminc_1,
    spline = list(model = first_sex$spline_model, grid = spline_grid),
    sex_analysis = list(
      p_int = pooled_global_p_int, strat_models = first_sex$strat_models,
      spline_model = first_sex$spline_model, ate = pooled_sex_ate,
      pooled_interaction = pooled_interaction$table, pooled_global_p_int = pooled_global_p_int,
      pooled_evalues = pooled_evalues, sensitivity_fgr = if (sp$type == "competing") fgr_pooled else NULL
    ),
    sensitivity = list(fgr = if (sp$type == "competing") fgr_sens_list[[1]] else NULL),
    rhs = rhs,
    pooled = list(
      coef = pooled_coef, interaction = pooled_interaction$table, global_p_int = pooled_global_p_int,
      mean = pooled_ate_fu$mean, diff = pooled_ate_fu$diff, ratio = pooled_ate_fu$ratio,
      curve = pooled_ate_fu$curve, evalues = pooled_evalues, sex_ate = pooled_sex_ate,
      performance = performance_summary, performance_raw = perf_list, sensitivity_fgr = fgr_pooled
    )
  )

  # Save endpoint bundle.
  fpath <- file.path(paths$bundles, paste0(nm, ".RData"))
  save(bundle, file = fpath)
  log_msg("Saved bundle: ", basename(fpath))

  list(ep = nm, file = fpath, n = nrow(bundle$data), type = sp$type)
}

############################################################
# 8) RUN ALL ENDPOINTS
############################################################
bundle_paths <- lapply(names(eps), function(nm) {
  process_endpoint(nm, eps[[nm]], imp)
})

manifest <- rbindlist(bundle_paths)
setcolorder(manifest, c("ep", "file", "n", "type"))
save(manifest, file = file.path(paths$base, "manifest.RData"))

save(sessionInfo(), file = file.path(paths$base, "session_info.RData"))
log_msg("END: analysis pipeline complete")
log_msg("Bundles saved to: ", paths$bundles)


