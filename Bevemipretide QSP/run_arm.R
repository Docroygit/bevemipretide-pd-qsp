## Run a single VPop arm and save to partial cache
## Usage: Rscript run_arm.R <arm_index>  (1=Placebo, 2=10mg, 3=30mg, 4=60mg)

args <- commandArgs(trailingOnly = TRUE)
arm_idx <- as.integer(args[1])

library(deSolve)
suppressPackageStartupMessages(library(dplyr))

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)
BW_human <- 70.0
HUMAN_PK <- list(ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
                 k_brain_in = 0.015, k_brain_out = 0.010,
                 k_mito_in = 0.50, k_mito_out = 0.020)
CI_max_human <- 0.90; alpha_clear_human <- 0.35; death_scale_human <- 0.06

po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day$C_mito)

VPOP_TOL <- 1e-8; PATIENT_TIMEOUT <- 60; N_POP <- 250; TRIAL_DAYS <- 547

run_human <- function(dose_mg, duration_days, CI_max = 1.0, alpha_clear = 1.0,
                      dt = 6, parms_extra = NULL, quiet = TRUE, solver_tol = 1e-10) {
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
                 parms_override = po, quiet = quiet, solver_tol = solver_tol)
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

sample_ln <- function(n, nom, cv) {
  sigma <- sqrt(log(1 + cv^2)); mu <- log(nom) - sigma^2 / 2
  exp(rnorm(n, mu, sigma))
}

set.seed(42)
vpop_pd <- data.frame(
  name = c("k_ROS_basal", "K_mPTP", "k_SOD2", "k_agg", "k_shunt", "k_clear",
           "k_death", "K_death", "k_death_inflam", "k_death_aSyn"),
  nominal = c(0.122, 0.40, 1.0, 0.009, 0.125, 0.241,
              0.0007 * death_scale_human, 0.25,
              0.00020 * death_scale_human, 0.00015 * death_scale_human),
  cv = c(0.25, 0.25, 0.25, 0.30, 0.20, 0.25, 0.30, 0.25, 0.50, 0.50),
  stringsAsFactors = FALSE)
vpop_pk <- data.frame(
  name = c("ka", "ke_plasma", "k_brain_in"),
  nominal = c(0.40, 0.032, 0.015), cv = c(0.30, 0.30, 0.40),
  stringsAsFactors = FALSE)

vpop <- data.frame(id = 1:N_POP)
for (r in 1:nrow(vpop_pd)) vpop[[vpop_pd$name[r]]] <- sample_ln(N_POP, vpop_pd$nominal[r], vpop_pd$cv[r])
for (r in 1:nrow(vpop_pk)) vpop[[vpop_pk$name[r]]] <- sample_ln(N_POP, vpop_pk$nominal[r], vpop_pk$cv[r])
vpop$CI_max_h <- pmin(0.96, pmax(0.80, rnorm(N_POP, CI_max_human, 0.03)))
vpop$alpha_clear_h <- pmin(0.55, pmax(0.20, rnorm(N_POP, alpha_clear_human, 0.05)))

arms_dose <- c(0, 10, 30, 60)
arms_name <- c("Placebo", "SBT-272 10mg", "SBT-272 30mg", "SBT-272 60mg")

dose <- arms_dose[arm_idx]
arm_label <- arms_name[arm_idx]
CACHE_ARM <- sprintf("RA_arm_%d.rds", arm_idx)

if (file.exists(CACHE_ARM)) {
  cat(sprintf("Arm %s already cached (%s). Skipping.\n", arm_label, CACHE_ARM))
  q("no")
}

cat(sprintf("Running %s arm (%d patients, SC daily)...\n", arm_label, N_POP))
arm_results <- data.frame()
n_timeout <- 0
t_start <- proc.time()["elapsed"]

for (i in 1:N_POP) {
  if (i %% 50 == 0) cat(sprintf("  %s: %d/%d\n", arm_label, i, N_POP))
  po_pd <- as.list(vpop[i, vpop_pd$name])
  po_pk <- as.list(vpop[i, vpop_pk$name])
  po <- c(po_pd, po_pk)
  res <- run_patient(dose, TRIAL_DAYS, CI_max = vpop$CI_max_h[i],
                     alpha_clear = vpop$alpha_clear_h[i], parms_extra = po)
  if (is.null(res)) { n_timeout <- n_timeout + 1; next }
  ss <- res[nrow(res), ]
  arm_results <- rbind(arm_results, data.frame(
    id = i, arm = arm_label,
    DA_neuron = ss$DA_neuron, Motor_score = ss$Motor_score,
    CI_activity = ss$CI_activity, aSyn_olig = ss$aSyn_olig,
    mROS = ss$mROS))
}

t_elapsed <- proc.time()["elapsed"] - t_start
cat(sprintf("  %s: done (%d/%d completed, %d timeout, %.0f sec)\n",
            arm_label, nrow(arm_results), N_POP, n_timeout, t_elapsed))
saveRDS(arm_results, CACHE_ARM)
cat(sprintf("  Saved: %s\n", CACHE_ARM))
