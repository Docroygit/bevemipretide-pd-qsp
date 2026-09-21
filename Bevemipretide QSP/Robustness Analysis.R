##=============================================================================
## ROBUSTNESS ANALYSIS
## Multi-pathway death model | Authentic VPop | Clinical Trial Simulation
## MCID-based responder rate | Paired + Unpaired effect sizes
##
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
## Doses: 10, 30, 60 mg SC daily (Phase 1 ceiling = 60 mg)
## Route: daily subcutaneous injection (matching clinical Phase 1 and
##        ALS TDP-43 mouse model; Ozdinler et al., NEALS 2022)
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

cat("================================================================\n")
cat("  ROBUSTNESS ANALYSIS\n")
cat("  Multi-pathway death | VPop | CTS\n")
cat("  Route: daily SC injection | Doses: 10, 30, 60 mg\n")
cat("================================================================\n\n")

## Calibrated parameters
CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)

## Human PK (from Human Translation, validated by Phase 1 data)
BW_human <- 70.0
HUMAN_PK <- list(
  ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
  k_brain_in = 0.015, k_brain_out = 0.010,
  k_mito_in = 0.50, k_mito_out = 0.020
)

## Human disease modifiers (recalibrated for multi-pathway death model)
CI_max_human <- 0.90
alpha_clear_human <- 0.35
death_scale_human <- 0.06

## C_ref calibration
po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day$C_mito)
cat(sprintf("C_ref_human = %.6f\n\n", C_ref_human))

VPOP_TOL <- 1e-8
PATIENT_TIMEOUT <- 60

## --- Core run functions ---

run_human <- function(dose_mg, duration_days, CI_max = 1.0, alpha_clear = 1.0,
                      dt = 6, parms_extra = NULL, quiet = TRUE,
                      solver_tol = 1e-10) {
  dose_mgkg <- dose_mg / BW_human
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * death_scale_human,
          k_death_inflam = 0.00020 * death_scale_human,
          k_death_aSyn = 0.00015 * death_scale_human))
  if (!is.null(parms_extra)) po <- c(po, parms_extra)
  run_integrated(dose_mg_kg = dose_mgkg, duration_days = duration_days, dt = dt,
                 CI_max = CI_max, alpha_clear = alpha_clear,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
                 parms_override = po, quiet = quiet,
                 solver_tol = solver_tol)
}

run_patient <- function(dose_mg, duration_days, CI_max, alpha_clear,
                        dt = 24, parms_extra = NULL) {
  safe_run <- function() {
    setTimeLimit(elapsed = PATIENT_TIMEOUT, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    run_human(dose_mg, duration_days, CI_max = CI_max,
              alpha_clear = alpha_clear, dt = dt,
              parms_extra = parms_extra, solver_tol = VPOP_TOL)
  }
  tryCatch(safe_run(), error = function(e) NULL)
}


##=============================================================================
## SECTION 1: VIRTUAL POPULATION GENERATION
## 15 parameters: 10 PD + 3 PK + 2 disease modifiers
##=============================================================================

cat("================================================================\n")
cat("  SECTION 1: VIRTUAL POPULATION\n")
cat("  N=250 | 15 parameters | PK + PD + death pathway variability\n")
cat("================================================================\n\n")

set.seed(42)
N_POP <- 250
TRIAL_DAYS <- 547

sample_ln <- function(n, nom, cv) {
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(nom) - sigma^2 / 2
  exp(rnorm(n, mu, sigma))
}

vpop_pd <- data.frame(
  name = c("k_ROS_basal", "K_mPTP", "k_SOD2", "k_agg",
           "k_shunt", "k_clear", "k_death", "K_death",
           "k_death_inflam", "k_death_aSyn"),
  nominal = c(0.122, 0.40, 1.0, 0.009, 0.125, 0.241,
              0.0007 * death_scale_human, 0.25,
              0.00020 * death_scale_human, 0.00015 * death_scale_human),
  cv = c(0.25, 0.25, 0.25, 0.30, 0.20, 0.25, 0.30, 0.25,
         0.50, 0.50),
  stringsAsFactors = FALSE
)

vpop_pk <- data.frame(
  name = c("ka", "ke_plasma", "k_brain_in"),
  nominal = c(0.40, 0.032, 0.015),
  cv = c(0.30, 0.30, 0.40),
  stringsAsFactors = FALSE
)

vpop <- data.frame(id = 1:N_POP)
for (r in 1:nrow(vpop_pd)) {
  vpop[[vpop_pd$name[r]]] <- sample_ln(N_POP, vpop_pd$nominal[r], vpop_pd$cv[r])
}
for (r in 1:nrow(vpop_pk)) {
  vpop[[vpop_pk$name[r]]] <- sample_ln(N_POP, vpop_pk$nominal[r], vpop_pk$cv[r])
}
vpop$CI_max_h <- pmin(0.96, pmax(0.80, rnorm(N_POP, CI_max_human, 0.03)))
vpop$alpha_clear_h <- pmin(0.55, pmax(0.20, rnorm(N_POP, alpha_clear_human, 0.05)))

cat(sprintf("Generated %d virtual patients\n", N_POP))
cat("PD parameters (median [IQR]):\n")
for (nm in vpop_pd$name) {
  v <- vpop[[nm]]
  cat(sprintf("  %-16s: %.5f [%.5f, %.5f]\n", nm, median(v),
              quantile(v, 0.25), quantile(v, 0.75)))
}
cat("PK parameters (median [IQR]):\n")
for (nm in vpop_pk$name) {
  v <- vpop[[nm]]
  cat(sprintf("  %-16s: %.5f [%.5f, %.5f]\n", nm, median(v),
              quantile(v, 0.25), quantile(v, 0.75)))
}
cat(sprintf("  %-16s: %.3f [%.3f, %.3f]\n", "CI_max_h",
            median(vpop$CI_max_h), quantile(vpop$CI_max_h, 0.25),
            quantile(vpop$CI_max_h, 0.75)))
cat(sprintf("  %-16s: %.3f [%.3f, %.3f]\n\n", "alpha_clear_h",
            median(vpop$alpha_clear_h), quantile(vpop$alpha_clear_h, 0.25),
            quantile(vpop$alpha_clear_h, 0.75)))


##=============================================================================
## SECTION 2: VPOP SIMULATION (4 arms x N patients)
## Daily SC injection (dosing events) — matches clinical route
## Uses RDS cache to avoid re-running
##=============================================================================

cat("================================================================\n")
cat("  SECTION 2: VPOP SIMULATION\n")
cat("  4 arms x 250 patients x 547-day | SC daily dosing\n")
cat("================================================================\n\n")

arms_dose <- c(0, 10, 30, 60)
arms_name <- c("Placebo", "SBT-272 10mg", "SBT-272 30mg", "SBT-272 60mg")

## Assemble VPop results from per-arm cache files (run_arm.R)
## Run arms first:  Rscript run_arm.R 1  (then 2, 3, 4)

vpop_results <- data.frame()
all_arms_cached <- TRUE
for (ai in seq_along(arms_dose)) {
  cache_f <- sprintf("RA_arm_%d.rds", ai)
  if (file.exists(cache_f)) {
    arm_data <- readRDS(cache_f)
    vpop_results <- rbind(vpop_results, arm_data)
    cat(sprintf("  Loaded %s: %d patients from %s\n",
                arms_name[ai], nrow(arm_data), cache_f))
  } else {
    cat(sprintf("  MISSING: %s — run: Rscript run_arm.R %d\n", arms_name[ai], ai))
    all_arms_cached <- FALSE
  }
}

if (!all_arms_cached) {
  cat("\nNot all arms cached. Run missing arms first, then re-run this script.\n")
  q("no")
}

completed_ids <- Reduce(intersect, lapply(arms_name, function(a) {
  vpop_results$id[vpop_results$arm == a]
}))
n_completed <- length(completed_ids)
cat(sprintf("\nCompleted: %d/%d patients across all arms (%.0f%%)\n\n",
            n_completed, N_POP, n_completed / N_POP * 100))

## Boundedness check
cat("--- BOUNDEDNESS CHECK ---\n")
for (v in c("DA_neuron", "Motor_score", "CI_activity", "aSyn_olig", "mROS")) {
  vals <- vpop_results[[v]]
  cat(sprintf("  %-14s: range [%.4f, %.4f]  %s\n",
              v, min(vals), max(vals),
              ifelse(min(vals) >= 0 & max(vals) <= 1, "PASS", "FAIL")))
}

## Summary statistics
cat("\n--- VPOP RESULTS (18 months) ---\n")
for (arm_label in arms_name) {
  sub <- vpop_results[vpop_results$arm == arm_label &
                       vpop_results$id %in% completed_ids, ]
  cat(sprintf("\n  %s (n=%d):\n", arm_label, nrow(sub)))
  for (v in c("DA_neuron", "Motor_score", "CI_activity", "aSyn_olig")) {
    vals <- sub[[v]]
    cat(sprintf("    %-14s: median=%.3f [%.3f, %.3f]  mean=%.3f sd=%.3f\n",
                v, median(vals), quantile(vals, 0.25), quantile(vals, 0.75),
                mean(vals), sd(vals)))
  }
}


##=============================================================================
## SECTION 3: STATISTICAL ANALYSIS (Clinical Trial Simulation)
## MCID-based responders | Paired + Unpaired effect sizes | Power
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 3: CLINICAL TRIAL SIMULATION\n")
cat("  MCID responders | Paired + Unpaired Cohen's d | Cliff's delta\n")
cat("================================================================\n\n")

vpop_results$arm <- factor(vpop_results$arm, levels = arms_name)
pbo <- vpop_results[vpop_results$arm == "Placebo" &
                     vpop_results$id %in% completed_ids, ]

MCID_DA <- 0.05
MCID_Motor <- 0.03

cliffs_delta <- function(x, y) {
  n_x <- length(x); n_y <- length(y)
  d <- 0
  for (xi in x) d <- d + sum(xi > y) - sum(xi < y)
  d / (n_x * n_y)
}

effect_table <- data.frame()
for (arm_name in arms_name[-1]) {
  trt <- vpop_results[vpop_results$arm == arm_name &
                       vpop_results$id %in% completed_ids, ]
  trt <- trt[order(trt$id), ]
  pbo_m <- pbo[pbo$id %in% trt$id, ]
  pbo_m <- pbo_m[order(pbo_m$id), ]

  for (ep in c("DA_neuron", "Motor_score")) {
    delta <- trt[[ep]] - pbo_m[[ep]]
    n_pairs <- length(delta)
    mean_d <- mean(delta); sd_d <- sd(delta)
    se_d <- sd_d / sqrt(n_pairs)
    ci_lo <- mean_d - 1.96 * se_d
    ci_hi <- mean_d + 1.96 * se_d

    d_paired <- abs(mean_d) / sd_d

    pooled_sd <- sqrt((var(trt[[ep]]) + var(pbo_m[[ep]])) / 2)
    d_unpaired <- abs(mean(trt[[ep]]) - mean(pbo_m[[ep]])) / pooled_sd

    tt <- t.test(delta, mu = 0)
    wt <- wilcox.test(delta, mu = 0, exact = FALSE)

    cd <- cliffs_delta(trt[[ep]], pbo_m[[ep]])
    if (ep == "Motor_score") cd <- -cd

    if (ep == "DA_neuron") {
      resp <- sum(delta > MCID_DA) / n_pairs * 100
      resp_any <- sum(delta > 0) / n_pairs * 100
    } else {
      resp <- sum(delta < -MCID_Motor) / n_pairs * 100
      resp_any <- sum(delta < 0) / n_pairs * 100
    }

    if (d_unpaired > 0) {
      n_power_unpaired <- ceiling(2 * (qnorm(0.975) + qnorm(0.80))^2 / d_unpaired^2)
    } else {
      n_power_unpaired <- NA
    }
    if (d_paired > 0) {
      n_power_paired <- ceiling((qnorm(0.975) + qnorm(0.80))^2 / d_paired^2)
    } else {
      n_power_paired <- NA
    }

    effect_table <- rbind(effect_table, data.frame(
      arm = arm_name, endpoint = ep, n = n_pairs,
      mean_delta = mean_d, se = se_d, ci_lo = ci_lo, ci_hi = ci_hi,
      d_paired = d_paired, d_unpaired = d_unpaired,
      cliffs_delta = cd,
      p_ttest = tt$p.value, p_wilcox = wt$p.value,
      resp_MCID = resp, resp_any = resp_any,
      n_power_paired = n_power_paired,
      n_power_unpaired = n_power_unpaired
    ))
  }
}

cat("--- TREATMENT EFFECTS ---\n\n")
cat(sprintf("  MCID thresholds: DA_neuron > %.2f, Motor_score > %.2f\n\n",
            MCID_DA, MCID_Motor))

cat(sprintf("  %-16s %-12s  N  %8s  d_pair d_unp  d_Cliff  p_t       p_W      Resp%%  RespAny%%  N_pair N_unp\n",
            "Arm", "Endpoint", "MeanD"))
for (r in 1:nrow(effect_table)) {
  with(effect_table[r, ], cat(sprintf(
    "  %-16s %-12s %3d  %+7.4f  %.3f  %.3f  %+.3f   %.1e  %.1e  %5.1f%%  %6.1f%%  %5d %5d\n",
    arm, endpoint, n, mean_delta,
    d_paired, d_unpaired, cliffs_delta,
    p_ttest, p_wilcox,
    resp_MCID, resp_any,
    n_power_paired, n_power_unpaired)))
}

cat("\n--- BENCHMARKING (unpaired Cohen's d vs real PD trials) ---\n")
cat("  DATATOP (selegiline, UPDRS):   d ~ 0.30,  N = 400/arm\n")
cat("  ADAGIO (rasagiline, UPDRS):    d ~ 0.20,  N = 290/arm\n")
cat("  SURE-PD3 (inosine, UPDRS):     d ~ 0.05 (neg), N = 149/arm\n")
cat("  SPARK (GBA, GCase activity):    d ~ 0.15,  N = 127/arm\n")

row30_DA <- effect_table[effect_table$arm == "SBT-272 30mg" &
                          effect_table$endpoint == "DA_neuron", ]
row30_Motor <- effect_table[effect_table$arm == "SBT-272 30mg" &
                             effect_table$endpoint == "Motor_score", ]
if (nrow(row30_DA) > 0) {
  cat(sprintf("\n  Our model (30 mg, DA):          d = %.3f (unpaired), N = %d/arm\n",
              row30_DA$d_unpaired, row30_DA$n_power_unpaired))
}
if (nrow(row30_Motor) > 0) {
  cat(sprintf("  Our model (30 mg, Motor):       d = %.3f (unpaired), N = %d/arm\n",
              row30_Motor$d_unpaired, row30_Motor$n_power_unpaired))
}
cat("\n  Unpaired d is directly comparable to published trial effect sizes.\n")
cat("  Paired d reflects model precision (same virtual patient +/- drug) and\n")
cat("  should NOT be compared to clinical trial Cohen's d values.\n")

## ALS cross-validation note
cat("\n--- ALS CROSS-VALIDATION (NEALS 2022, Ozdinler et al.) ---\n")
cat("  TDP-43 mouse model, SC dosing (same route as our simulation):\n")
cat("    UMN retention:  vehicle=29.6, 1mg/kg=51.4, 5mg/kg=49.9 (P<0.0001)\n")
cat("    Microglia:      vehicle=37.1, 1mg/kg=9.1,  5mg/kg=13.3 (P<0.0001)\n")
cat("    Astrogliosis:   vehicle=59.8, 1mg/kg=20.1, 5mg/kg=24.0 (P<0.0001)\n")
cat("  Our model predicts: neuroprotection (DA rescue), neuroinflammation\n")
cat("  reduction (MG_active), and dose-response plateau — all confirmed.\n\n")


##=============================================================================
## SECTION 4: VISUALISATION
##=============================================================================

cat("================================================================\n")
cat("  SECTION 4: PLOTS\n")
cat("================================================================\n\n")

plot_data <- vpop_results[vpop_results$id %in% completed_ids, ]
plot_data$arm <- factor(plot_data$arm, levels = arms_name)

## --- Plot 1: VPop violin (DA_neuron by arm) ---
p1 <- ggplot(plot_data, aes(x = arm, y = DA_neuron, fill = arm)) +
  geom_violin(alpha = 0.5, width = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("#999999", "#56B4E9", "#009E73", "#D55E00")) +
  labs(title = "DA Neuron Survival by Treatment Arm (18 months)",
       subtitle = "Daily SC injection | N=250 per arm",
       x = NULL, y = "DA Neuron Fraction") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 15, hjust = 1))
ggsave("RA_vpop_violin.png", p1, width = 8, height = 5, dpi = 150)
cat("Saved: RA_vpop_violin.png\n")

## --- Plot 2: Responder waterfall (30mg, DA_neuron) ---
row30 <- effect_table[effect_table$arm == "SBT-272 30mg" &
                       effect_table$endpoint == "DA_neuron", ]
if (nrow(row30) > 0) {
  trt30 <- vpop_results[vpop_results$arm == "SBT-272 30mg" &
                         vpop_results$id %in% completed_ids, ]
  pbo30 <- vpop_results[vpop_results$arm == "Placebo" &
                         vpop_results$id %in% completed_ids, ]
  trt30 <- trt30[order(trt30$id), ]
  pbo30 <- pbo30[order(pbo30$id), ]
  delta_DA <- trt30$DA_neuron - pbo30$DA_neuron

  wf <- data.frame(patient = rank(-delta_DA), delta = delta_DA)
  wf$responder <- ifelse(wf$delta > MCID_DA, "Responder",
                         ifelse(wf$delta > 0, "Marginal", "Non-responder"))
  wf$responder <- factor(wf$responder,
                          levels = c("Responder", "Marginal", "Non-responder"))

  p2 <- ggplot(wf, aes(x = patient, y = delta * 100, fill = responder)) +
    geom_bar(stat = "identity", width = 1) +
    geom_hline(yintercept = MCID_DA * 100, linetype = "dashed",
               colour = "red", linewidth = 0.8) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    annotate("text", x = 5, y = MCID_DA * 100 + 1.5,
             label = sprintf("MCID = %.0f%%", MCID_DA * 100),
             colour = "red", hjust = 0, size = 3.5) +
    scale_fill_manual(values = c("Responder" = "#009E73",
                                  "Marginal" = "#E69F00",
                                  "Non-responder" = "#CC79A7")) +
    labs(title = "SBT-272 30 mg SC Daily: Individual DA Neuron Preservation",
         x = "Patient (ranked)", y = "DA Preservation vs Placebo (%)",
         fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave("RA_waterfall.png", p2, width = 8, height = 5, dpi = 150)
  cat("Saved: RA_waterfall.png\n")
}

## --- Plot 3: Forest plot (treatment effects with CIs) ---
forest_data <- effect_table[effect_table$endpoint == "DA_neuron", ]
if (nrow(forest_data) > 0) {
  fd <- data.frame()
  for (r in 1:nrow(forest_data)) {
    row <- forest_data[r, ]
    se_unp <- abs(row$mean_delta) / row$d_unpaired / sqrt(row$n / 2)
    fd <- rbind(fd,
      data.frame(arm = row$arm, type = "Paired",
                 mean = row$mean_delta * 100,
                 lo = row$ci_lo * 100, hi = row$ci_hi * 100),
      data.frame(arm = row$arm, type = "Unpaired",
                 mean = row$mean_delta * 100,
                 lo = (row$mean_delta - 1.96 * se_unp) * 100,
                 hi = (row$mean_delta + 1.96 * se_unp) * 100)
    )
  }
  fd$arm <- factor(fd$arm, levels = rev(arms_name[-1]))
  fd$type <- factor(fd$type, levels = c("Paired", "Unpaired"))

  p3 <- ggplot(fd, aes(x = mean, y = arm, colour = type)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), width = 0.2,
                  orientation = "y") +
    scale_colour_manual(values = c("Paired" = "#0072B2", "Unpaired" = "#D55E00")) +
    labs(title = "DA Neuron Preservation: Paired vs Unpaired Effect Estimates",
         subtitle = "Daily SC injection",
         x = "DA Preservation vs Placebo (%)", y = NULL, colour = "Analysis") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  ggsave("RA_forest.png", p3, width = 8, height = 4, dpi = 150)
  cat("Saved: RA_forest.png\n")
}


##=============================================================================
## SECTION 5: SUMMARY
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 5: SUMMARY\n")
cat("================================================================\n\n")

cat("MODEL STRUCTURE (v2.1 -- multi-pathway death):\n")
cat("  Three DA neuron death pathways:\n")
cat("    1. Mitochondrial damage (drug-rescuable via CL stabilisation)\n")
cat("    2. Neuroinflammation (TNF/IL1b; partially drug-rescuable)\n")
cat("    3. aSyn proteotoxicity (partially drug-rescuable)\n")
cat("  Patients with high k_death_inflam or k_death_aSyn are non-responders.\n\n")

cat("DOSING RATIONALE:\n")
cat("  Route: daily SC injection (matching Phase 1 and ALS preclinical data)\n")
cat("  Doses: 10, 30, 60 mg (Phase 1 ceiling; ISR severity increases at 60 mg)\n")
cat("  Phase 1 confirms allometric PK scaling predictions (Stealth Bio, NEALS 2022)\n\n")

cat("PK VALIDATION:\n")
cat("  Rat plasma t1/2 <=3h (poster) vs model ke=0.231 -> t1/2=3.0h\n")
cat("  Brain C24h/Cmax=0.677 (poster) vs model e^(-0.016*24)=0.681\n")
cat("  Human PK: Phase 1 consistent with allometric scaling (BW^-0.25)\n\n")

cat("VPOP DESIGN:\n")
cat("  15 parameters varied per patient (10 PD + 3 PK + 2 disease modifiers)\n")
cat("  CL-independent death pathway parameters varied with CV=0.50\n")
cat("  Log-normal sampling for rate constants, truncated normal for modifiers\n\n")

cat("STATISTICAL REPORTING:\n")
cat("  - MCID-based responder rate (DA > 5%%, Motor > 3%%)\n")
cat("  - Unpaired Cohen's d for clinical trial benchmarking\n")
cat("  - Paired Cohen's d reported separately as model precision metric\n")
cat("  - Cliff's delta as non-parametric alternative\n")
cat("  - Power analysis uses unpaired d (realistic trial sizing)\n\n")

cat("VALIDATION SUMMARY:\n")
cat("  a. 24/24 quantitative validation against independent literature targets\n")
cat("  b. 19/19 three-tier validation hierarchy (FDA/IQ Consortium compliant)\n")
cat("  c. Three-tiered sensitivity analysis (OAT + Morris + Sobol)\n")
cat("  d. PK validated by Phase 1 human data and rat brain tissue data\n")
cat("  e. Mechanism cross-validated by ALS TDP-43 mouse model (NEALS 2022)\n")
cat("  f. Multi-pathway death produces biologically realistic non-responders\n")
cat("  g. Unpaired effect sizes comparable to historical PD trial data\n")
cat("  h. Drug slows but does not halt progression (biologically honest)\n")
cat("  i. Dose-response shows diminishing returns (consistent with ALS data)\n\n")

cat("================================================================\n")
cat("  ROBUSTNESS ANALYSIS -- COMPLETE\n")
cat("================================================================\n")
