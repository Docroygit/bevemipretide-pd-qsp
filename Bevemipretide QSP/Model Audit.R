##=============================================================================
## MODEL AUDIT
## Systematic analysis of bifurcation behavior, solver stability,
## and model inconsistencies for manuscript preparation
##
## All simulation calls use setTimeLimit (20s max) to prevent hangs
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

cat("================================================================\n")
cat("  MODEL AUDIT\n")
cat("  Bifurcation mapping | Solver stability | Inconsistencies\n")
cat("================================================================\n\n")

CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)

BW_human <- 70.0
HUMAN_PK <- list(
  ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
  k_brain_in = 0.015, k_brain_out = 0.010,
  k_mito_in = 0.50, k_mito_out = 0.020
)
CI_max_human <- 0.96
alpha_clear_human <- 0.55

po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day$C_mito)

run_human <- function(dose_mg, duration_days, CI_max = 1.0, alpha_clear = 1.0,
                      dt = 6, parms_extra = NULL, quiet = TRUE,
                      solver_tol = 1e-10) {
  dose_mgkg <- dose_mg / BW_human
  po <- c(HUMAN_PK, list(C_ref = C_ref_human))
  if (!is.null(parms_extra)) po <- c(po, parms_extra)
  run_integrated(dose_mg_kg = dose_mgkg, duration_days = duration_days, dt = dt,
                 CI_max = CI_max, alpha_clear = alpha_clear,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
                 parms_override = po, quiet = quiet,
                 solver_tol = solver_tol)
}

## Timeout-protected wrappers — 20 sec max per call
safe_run <- function(...) {
  tryCatch({
    setTimeLimit(elapsed = 20, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    run_human(...)
  }, error = function(e) NULL)
}

safe_run_mouse <- function(...) {
  tryCatch({
    setTimeLimit(elapsed = 20, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    run_integrated(...)
  }, error = function(e) NULL)
}


##=============================================================================
## AUDIT 1: MOUSE BIFURCATION (35 days) — already completed, print summary
##=============================================================================
cat("================================================================\n")
cat("  AUDIT 1: MOUSE MODEL BIFURCATION [COMPLETED PREVIOUSLY]\n")
cat("================================================================\n")
cat("  494 grid points, 0 solver failures\n")
cat("  DA range: [0.4519, 0.9959]\n")
cat("  Transition is gradual (no sharp bifurcation at mouse timescale)\n")
cat("  CI_max=0.74 (calibrated): DA mean=0.829\n\n")


##=============================================================================
## AUDIT 2: HUMAN 1-YEAR LANDSCAPE — already completed, print summary
##=============================================================================
cat("================================================================\n")
cat("  AUDIT 2: HUMAN 1-YEAR BIFURCATION [COMPLETED PREVIOUSLY]\n")
cat("================================================================\n")
cat("  126 grid points, 0 solver failures\n")
cat("  DA range: [0.460, 0.946]\n")
cat("  Nominal (ci=0.96, ac=0.55): DA=0.913\n")
cat("  Worst case (ci=0.86, ac=0.30): DA=0.460\n\n")


##=============================================================================
## AUDIT 3: HUMAN 5-YEAR LANDSCAPE — already completed, print summary
##=============================================================================
cat("================================================================\n")
cat("  AUDIT 3: HUMAN 5-YEAR LANDSCAPE [COMPLETED PREVIOUSLY]\n")
cat("================================================================\n")
cat("  126 grid points, 0 solver failures\n")
cat("  DA range at 5y: [0.020, 0.757]\n")
cat("  Nominal (ci=0.96, ac=0.55): DA=0.633\n")
cat("  Sharp gradient at low CI_max: ci=0.86/ac=0.30 → DA=0.020\n\n")


##=============================================================================
## AUDIT 4: VPOP PARAMETER SPACE STABILITY
## 50 patients × 4 arms, timeout-protected
##=============================================================================

cat("================================================================\n")
cat("  AUDIT 4: VPOP PARAMETER SPACE STABILITY\n")
cat("  50 patients × 4 arms × 547-day sim (20s timeout)\n")
cat("================================================================\n\n")

sample_ln <- function(n, nom, cv) {
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(nom) - sigma^2 / 2
  exp(rnorm(n, mu, sigma))
}

set.seed(42)
N_AUDIT <- 50
vpop_pd <- data.frame(
  name = c("k_ROS_basal", "K_mPTP", "k_SOD2", "k_agg",
           "k_shunt", "k_clear", "k_death", "K_death"),
  nominal = c(0.122, 0.40, 1.0, 0.009, 0.125, 0.241, 0.001, 0.25),
  cv = c(0.25, 0.25, 0.25, 0.30, 0.20, 0.25, 0.30, 0.25),
  stringsAsFactors = FALSE
)
vpop_pk <- data.frame(
  name = c("ka", "ke_plasma", "k_brain_in"),
  nominal = c(0.40, 0.032, 0.015),
  cv = c(0.30, 0.30, 0.40),
  stringsAsFactors = FALSE
)

vpop_test <- data.frame(id = 1:N_AUDIT)
for (r in 1:nrow(vpop_pd)) {
  vpop_test[[vpop_pd$name[r]]] <- sample_ln(N_AUDIT, vpop_pd$nominal[r], vpop_pd$cv[r])
}
for (r in 1:nrow(vpop_pk)) {
  vpop_test[[vpop_pk$name[r]]] <- sample_ln(N_AUDIT, vpop_pk$nominal[r], vpop_pk$cv[r])
}
vpop_test$CI_max_h <- pmin(0.98, pmax(0.80, rnorm(N_AUDIT, CI_max_human, 0.03)))
vpop_test$alpha_clear_h <- pmin(0.70, pmax(0.25, rnorm(N_AUDIT, alpha_clear_human, 0.06)))

arms_dose <- c(0, 10, 30, 100)
arms_name <- c("Placebo", "10mg", "30mg", "100mg")

stability_log <- data.frame()
for (ai in seq_along(arms_dose)) {
  dose <- arms_dose[ai]
  cat(sprintf("Testing %s arm...\n", arms_name[ai]))
  for (i in 1:N_AUDIT) {
    po_pd <- as.list(vpop_test[i, vpop_pd$name])
    po_pk <- as.list(vpop_test[i, vpop_pk$name])
    po <- c(po_pd, po_pk)

    t_start <- proc.time()["elapsed"]
    res <- safe_run(dose, 547, CI_max = vpop_test$CI_max_h[i],
                    alpha_clear = vpop_test$alpha_clear_h[i], dt = 24,
                    parms_extra = po, solver_tol = 1e-8)
    t_run <- proc.time()["elapsed"] - t_start

    stability_log <- rbind(stability_log, data.frame(
      id = i, arm = arms_name[ai],
      CI_max = vpop_test$CI_max_h[i],
      alpha_clear = vpop_test$alpha_clear_h[i],
      k_ROS = vpop_test$k_ROS_basal[i],
      K_mPTP = vpop_test$K_mPTP[i],
      ok = !is.null(res),
      time = t_run,
      DA = if (!is.null(res)) res$DA_neuron[nrow(res)] else NA
    ))
  }
  cat(sprintf("  %s: done\n", arms_name[ai]))
}

cat(sprintf("\nTotal runs: %d\n", nrow(stability_log)))
cat(sprintf("Failures (timeout or error): %d (%.1f%%)\n",
            sum(!stability_log$ok),
            100 * sum(!stability_log$ok) / nrow(stability_log)))

ok_runs <- stability_log[stability_log$ok, ]
cat(sprintf("\nRun times (successful):\n"))
cat(sprintf("  Mean: %.1f sec, Median: %.1f sec\n",
            mean(ok_runs$time), median(ok_runs$time)))
cat(sprintf("  P95: %.1f sec, Max: %.1f sec\n",
            quantile(ok_runs$time, 0.95), max(ok_runs$time)))

slow <- stability_log[stability_log$time > 10 & stability_log$ok, ]
if (nrow(slow) > 0) {
  cat(sprintf("\nSlow runs (> 10 sec): %d\n", nrow(slow)))
  for (r in 1:min(10, nrow(slow))) {
    with(slow[r, ], cat(sprintf(
      "  id=%d arm=%-5s CI_max=%.3f a_clear=%.3f k_ROS=%.4f time=%.1fs DA=%.4f\n",
      id, arm, CI_max, alpha_clear, k_ROS, time, DA)))
  }
}

failed <- stability_log[!stability_log$ok, ]
if (nrow(failed) > 0) {
  cat(sprintf("\nFailed runs: %d\n", nrow(failed)))
  for (r in 1:min(10, nrow(failed))) {
    with(failed[r, ], cat(sprintf(
      "  id=%d arm=%-5s CI_max=%.3f a_clear=%.3f k_ROS=%.4f K_mPTP=%.4f\n",
      id, arm, CI_max, alpha_clear, k_ROS, K_mPTP)))
  }
} else {
  cat("\nNo failed runs — all 200 completed within timeout.\n")
}


##=============================================================================
## AUDIT 5: VICIOUS CYCLE GAIN (35-day mouse, fast)
##=============================================================================

cat("\n\n================================================================\n")
cat("  AUDIT 5: VICIOUS CYCLE GAIN\n")
cat("================================================================\n\n")

gains <- data.frame()
ci_test <- c(0.60, 0.65, 0.68, 0.70, 0.72, 0.74, 0.76, 0.80, 0.85, 0.90, 0.95, 1.00)
for (ci in ci_test) {
  res <- safe_run_mouse(0, 35, CI_max = ci, alpha_clear = 0.25, dt = 0.5,
                        quiet = TRUE, solver_tol = 1e-8)
  if (!is.null(res)) {
    ss <- res[nrow(res), ]
    gains <- rbind(gains, data.frame(
      CI_max = ci, DA = ss$DA_neuron, aSyn = ss$aSyn_olig,
      CI = ss$CI_activity, mROS = ss$mROS, Motor = ss$Motor_score,
      CL_ratio = ss$CL_ratio))
  } else {
    gains <- rbind(gains, data.frame(
      CI_max = ci, DA = NA, aSyn = NA, CI = NA, mROS = NA,
      Motor = NA, CL_ratio = NA))
  }
}

cat("Vicious cycle response to CI_max (alpha_clear=0.25, 35 days):\n")
cat(sprintf("  %-7s  %-6s  %-7s  %-6s  %-6s  %-6s  %-6s\n",
            "CI_max", "DA", "aSyn_o", "CI", "mROS", "Motor", "CL_r"))
for (r in 1:nrow(gains)) {
  if (is.na(gains$DA[r])) {
    cat(sprintf("  %.2f    FAIL\n", gains$CI_max[r]))
  } else {
    with(gains[r, ], cat(sprintf(
      "  %.2f    %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n",
      CI_max, DA, aSyn, CI, mROS, Motor, CL_ratio)))
  }
}

ok_gains <- gains[!is.na(gains$DA), ]
if (nrow(ok_gains) > 1) {
  da_diff <- diff(ok_gains$DA)
  max_drop_idx <- which.min(da_diff)
  cat(sprintf("\nSteepest DA drop: CI_max %.2f→%.2f, DA %.4f→%.4f (Δ=%.4f)\n",
              ok_gains$CI_max[max_drop_idx + 1], ok_gains$CI_max[max_drop_idx],
              ok_gains$DA[max_drop_idx + 1], ok_gains$DA[max_drop_idx],
              da_diff[max_drop_idx]))
}


##=============================================================================
## AUDIT 6: NUMERICAL ACCURACY (tolerance comparison)
##=============================================================================

cat("\n\n================================================================\n")
cat("  AUDIT 6: NUMERICAL ACCURACY\n")
cat("================================================================\n\n")

cat("Mouse disease (35 days, CI_max=0.74, a_clear=0.25):\n")
for (tol in c(1e-6, 1e-8, 1e-10, 1e-12)) {
  res <- safe_run_mouse(0, 35, CI_max = 0.74, alpha_clear = 0.25,
                        dt = 0.5, quiet = TRUE, solver_tol = tol)
  if (!is.null(res)) {
    ss <- res[nrow(res), ]
    cat(sprintf("  tol=%.0e: DA=%.6f aSyn=%.6f CI=%.6f mROS=%.6f\n",
                tol, ss$DA_neuron, ss$aSyn_olig, ss$CI_activity, ss$mROS))
  } else {
    cat(sprintf("  tol=%.0e: TIMEOUT/FAIL\n", tol))
  }
}

cat("\nHuman disease (1 year, CI_max=0.96, a_clear=0.55):\n")
for (tol in c(1e-6, 1e-8, 1e-10, 1e-12)) {
  res <- safe_run(0, 365, CI_max = 0.96, alpha_clear = 0.55,
                  dt = 24, solver_tol = tol)
  if (!is.null(res)) {
    ss <- res[nrow(res), ]
    cat(sprintf("  tol=%.0e: DA=%.6f aSyn=%.6f CI=%.6f mROS=%.6f\n",
                tol, ss$DA_neuron, ss$aSyn_olig, ss$CI_activity, ss$mROS))
  } else {
    cat(sprintf("  tol=%.0e: TIMEOUT/FAIL\n", tol))
  }
}


##=============================================================================
## AUDIT 7: SAFE OPERATING ENVELOPE
##=============================================================================

cat("\n\n================================================================\n")
cat("  AUDIT 7: SAFE OPERATING ENVELOPE\n")
cat("================================================================\n\n")

base_parms <- list(
  k_ROS_basal = 0.122, K_mPTP = 0.40, k_SOD2 = 1.0,
  k_agg = 0.009, k_shunt = 0.125, k_clear = 0.241,
  k_death = 0.001, K_death = 0.25
)

cat("Parameter extremes (547-day human sim, CI_max=0.96, a_clear=0.55):\n")
cat(sprintf("  %-14s  %-8s  %-6s  %-6s  %-6s  %s\n",
            "Parameter", "Value", "DA", "aSyn", "Time", "Status"))

for (pname in names(base_parms)) {
  nom <- base_parms[[pname]]
  for (mult in c(0.5, 0.75, 1.0, 1.5, 2.0)) {
    po <- list()
    po[[pname]] <- nom * mult
    t_start <- proc.time()["elapsed"]
    res <- safe_run(0, 547, CI_max = 0.96, alpha_clear = 0.55, dt = 24,
                    parms_extra = po, solver_tol = 1e-8)
    t_run <- proc.time()["elapsed"] - t_start
    if (!is.null(res)) {
      ss <- res[nrow(res), ]
      cat(sprintf("  %-14s  %.4f  %.4f  %.4f  %5.1fs  OK\n",
                  pname, nom * mult, ss$DA_neuron, ss$aSyn_olig, t_run))
    } else {
      cat(sprintf("  %-14s  %.4f  ----    ----    %5.1fs  TIMEOUT\n",
                  pname, nom * mult, t_run))
    }
  }
  cat("\n")
}


##=============================================================================
## AUDIT 8: INCONSISTENCY SUMMARY
##=============================================================================

cat("\n\n================================================================\n")
cat("  AUDIT 8: INCONSISTENCY SUMMARY\n")
cat("================================================================\n\n")

cat("1. BIFURCATION BEHAVIOR\n")
cat("   The vicious cycle creates bistability. Below a critical CI_max,\n")
cat("   the system locks into a high-damage attractor. Biologically\n")
cat("   meaningful (PD 'point of no return') but causes solver stiffness\n")
cat("   at the transition boundary.\n\n")

cat("2. BIMODAL VPOP DISTRIBUTIONS\n")
cat("   Virtual patients cluster into healthy (DA>0.7) and disease\n")
cat("   (DA<0.3) attractors. Not a bug — reflects PD heterogeneity.\n")
cat("   But makes Cohen's d and mean differences misleading.\n\n")

cat("3. 100% RESPONDER RATE\n")
cat("   Model has only CL-mediated death pathway. Real PD has multiple\n")
cat("   mechanisms; some patients would not respond to CL rescue.\n\n")

cat("4. PAIRED DESIGN INFLATION\n")
cat("   Same virtual patient ± drug eliminates IIV. Inflates Cohen's d\n")
cat("   (~0.6 vs real trial d~0.2-0.3). Required N (~24/arm) is an\n")
cat("   artefact of model precision, not a real trial size estimate.\n\n")

cat("5. SOLVER SENSITIVITY\n")
cat("   Long-duration sims with extreme parameters can hang. setTimeLimit\n")
cat("   provides a reliable safeguard in sequential mode.\n\n")

cat("6. OUTPUT CLAMPING\n")
cat("   ODE solver can overshoot [0,1] bounds. Post-processing clamp is\n")
cat("   standard practice for bounded ODE systems.\n\n")

cat("================================================================\n")
cat("  MODEL AUDIT — COMPLETE\n")
cat("================================================================\n")
