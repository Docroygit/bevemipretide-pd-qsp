## ============================================================
##  UNCERTAINTY ANALYSIS
##  Bootstrap CIs on VPop effect sizes + parameter uncertainty
##  propagation using existing simulation data
## ============================================================

library(deSolve)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(tidyr)

cat("\n================================================================\n")
cat("  UNCERTAINTY ANALYSIS\n")
cat("  Bootstrap CIs + Parameter Uncertainty Propagation\n")
cat("================================================================\n\n")

## ---- Load VPop arm results ----
arm1 <- readRDS("RA_arm_1.rds")
arm2 <- readRDS("RA_arm_2.rds")
arm3 <- readRDS("RA_arm_3.rds")
arm4 <- readRDS("RA_arm_4.rds")

## Match patients across arms (same VPop, same id)
ids_common <- Reduce(intersect, list(arm1$id, arm2$id, arm3$id, arm4$id))
cat(sprintf("Matched patients across all 4 arms: %d\n\n", length(ids_common)))

pbo <- arm1[arm1$id %in% ids_common, ] %>% arrange(id)
d10 <- arm2[arm2$id %in% ids_common, ] %>% arrange(id)
d30 <- arm3[arm3$id %in% ids_common, ] %>% arrange(id)
d60 <- arm4[arm4$id %in% ids_common, ] %>% arrange(id)

## ============================================================
##  SECTION 1: BOOTSTRAP CONFIDENCE INTERVALS
##  Nonparametric bootstrap on matched patient pairs
## ============================================================

cat("================================================================\n")
cat("  SECTION 1: BOOTSTRAP CONFIDENCE INTERVALS\n")
cat("  B=10000 resamples | Percentile method | Matched pairs\n")
cat("================================================================\n\n")

B <- 10000
set.seed(123)
n <- nrow(pbo)

## Pre-allocate bootstrap storage
boot_stats <- list()

compute_stats <- function(pbo_b, drug_b, endpoint) {
  if (endpoint == "DA_neuron") {
    delta <- drug_b[[endpoint]] - pbo_b[[endpoint]]
    mcid_thresh <- 0.05
  } else {
    delta <- pbo_b[[endpoint]] - drug_b[[endpoint]]  # motor: lower is better
    mcid_thresh <- 0.03
  }
  mean_delta <- mean(delta)
  sd_delta   <- sd(delta)
  d_paired   <- mean_delta / sd_delta

  pooled_sd  <- sqrt((var(pbo_b[[endpoint]]) + var(drug_b[[endpoint]])) / 2)
  d_unpaired <- abs(mean(drug_b[[endpoint]]) - mean(pbo_b[[endpoint]])) / pooled_sd

  resp_rate  <- mean(abs(delta) > mcid_thresh &
                       ((endpoint == "DA_neuron" & delta > 0) |
                        (endpoint == "Motor_score" & delta > 0)))

  ## Cliff's delta
  n_b <- length(delta)
  pbo_vals <- pbo_b[[endpoint]]
  drg_vals <- drug_b[[endpoint]]
  if (endpoint == "DA_neuron") {
    cliff <- mean(sign(outer(drg_vals, pbo_vals, "-")))
  } else {
    cliff <- mean(sign(outer(pbo_vals, drg_vals, "-")))
  }

  c(mean_delta = mean_delta, d_paired = d_paired, d_unpaired = d_unpaired,
    cliff = cliff, resp_rate = resp_rate)
}

## Run bootstrap for each dose × endpoint
doses <- list(list(data = d10, name = "10mg"),
              list(data = d30, name = "30mg"),
              list(data = d60, name = "60mg"))
endpoints <- c("DA_neuron", "Motor_score")

cat("Running bootstrap (B=10000)...\n")

results_ci <- data.frame()

for (dose_info in doses) {
  for (ep in endpoints) {
    ## Observed statistics
    obs <- compute_stats(pbo, dose_info$data, ep)

    ## Bootstrap
    boot_mat <- matrix(NA, nrow = B, ncol = length(obs))
    colnames(boot_mat) <- names(obs)

    for (b in 1:B) {
      idx <- sample(1:n, n, replace = TRUE)
      boot_mat[b, ] <- compute_stats(pbo[idx, ], dose_info$data[idx, ], ep)
    }

    ## Percentile CIs (2.5%, 97.5%)
    ci_lo <- apply(boot_mat, 2, quantile, probs = 0.025)
    ci_hi <- apply(boot_mat, 2, quantile, probs = 0.975)

    ## Bias-corrected estimate
    boot_mean <- colMeans(boot_mat)
    boot_se   <- apply(boot_mat, 2, sd)

    for (stat_name in names(obs)) {
      results_ci <- rbind(results_ci, data.frame(
        Dose = dose_info$name, Endpoint = ep, Statistic = stat_name,
        Observed = obs[stat_name],
        Boot_Mean = boot_mean[stat_name],
        Boot_SE = boot_se[stat_name],
        CI_lo = ci_lo[stat_name],
        CI_hi = ci_hi[stat_name],
        stringsAsFactors = FALSE
      ))
    }

    ## Store bootstrap distributions for plotting
    boot_stats[[paste(dose_info$name, ep, sep = "_")]] <- as.data.frame(boot_mat)
  }
  cat(sprintf("  %s complete\n", dose_info$name))
}

cat("\n--- BOOTSTRAP 95%% CONFIDENCE INTERVALS ---\n\n")

## Format and print results
for (ep in endpoints) {
  cat(sprintf("  Endpoint: %s\n", ep))
  cat(sprintf("  %-8s  %-12s  %10s  [%10s, %10s]  %10s\n",
              "Dose", "Statistic", "Observed", "CI_lo", "CI_hi", "Boot_SE"))
  cat(paste(rep("-", 72), collapse = ""), "\n")

  sub <- results_ci[results_ci$Endpoint == ep, ]
  for (i in 1:nrow(sub)) {
    r <- sub[i, ]
    cat(sprintf("  %-8s  %-12s  %10.4f  [%10.4f, %10.4f]  %10.4f\n",
                r$Dose, r$Statistic, r$Observed, r$CI_lo, r$CI_hi, r$Boot_SE))
  }
  cat("\n")
}

## ============================================================
##  SECTION 2: PREDICTION INTERVALS ON KEY OUTCOMES
##  What range of treatment effects should we expect?
## ============================================================

cat("================================================================\n")
cat("  SECTION 2: PREDICTION INTERVALS\n")
cat("  Population-level uncertainty on key predictions\n")
cat("================================================================\n\n")

## DA neuron preservation: individual patient prediction interval
cat("--- INDIVIDUAL PATIENT PREDICTION INTERVALS (18 months) ---\n\n")

for (dose_info in doses) {
  delta_DA <- dose_info$data$DA_neuron[dose_info$data$id %in% ids_common] -
              pbo$DA_neuron
  delta_DA <- sort(delta_DA)
  n_d <- length(delta_DA)

  cat(sprintf("  %s DA_neuron improvement:\n", dose_info$name))
  cat(sprintf("    Mean:   %+.4f\n", mean(delta_DA)))
  cat(sprintf("    Median: %+.4f\n", median(delta_DA)))
  cat(sprintf("    80%% PI: [%+.4f, %+.4f]\n",
              quantile(delta_DA, 0.10), quantile(delta_DA, 0.90)))
  cat(sprintf("    95%% PI: [%+.4f, %+.4f]\n",
              quantile(delta_DA, 0.025), quantile(delta_DA, 0.975)))
  cat(sprintf("    Range:  [%+.4f, %+.4f]\n", min(delta_DA), max(delta_DA)))
  cat(sprintf("    %% with any benefit (>0): %.1f%%\n",
              100 * mean(delta_DA > 0)))
  cat(sprintf("    %% with MCID benefit (>5%%): %.1f%%\n\n",
              100 * mean(delta_DA > 0.05)))
}

## ============================================================
##  SECTION 3: PARAMETER UNCERTAINTY PROPAGATION
##  How sensitive are effect sizes to calibration parameter choices?
## ============================================================

cat("================================================================\n")
cat("  SECTION 3: PARAMETER UNCERTAINTY PROPAGATION\n")
cat("  Effect size sensitivity to death_scale, CI_max, alpha_clear\n")
cat("================================================================\n\n")

## Load integrated model for re-running scenarios
SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)
BW_human <- 70.0
HUMAN_PK <- list(ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
                 k_brain_in = 0.015, k_brain_out = 0.010,
                 k_mito_in = 0.50, k_mito_out = 0.020)

## Compute C_ref
po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day$C_mito)

run_scenario <- function(CI_max_h, alpha_clear_h, death_scale_h, dose_mg,
                          duration_days = 547, dt = 24) {
  dose_mgkg <- dose_mg / BW_human
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * death_scale_h,
          k_death_inflam = 0.00020 * death_scale_h,
          k_death_aSyn = 0.00015 * death_scale_h))
  run_integrated(dose_mg_kg = dose_mgkg, duration_days = duration_days, dt = dt,
                 CI_max = CI_max_h, alpha_clear = alpha_clear_h,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
                 parms_override = po, quiet = TRUE)
}

## Sweep key translational parameters around their nominal values
## Nominal: CI_max=0.90, alpha_clear=0.35, death_scale=0.06
cat("Sweeping translational parameters (±25% around nominal)...\n\n")

CI_max_vals      <- c(0.86, 0.88, 0.90, 0.92, 0.94)
alpha_clear_vals <- c(0.28, 0.32, 0.35, 0.38, 0.42)
death_scale_vals <- c(0.04, 0.05, 0.06, 0.07, 0.08)

## Effect of death_scale on predictions (most uncertain parameter)
cat("--- DEATH_SCALE SENSITIVITY (CI_max=0.90, alpha_clear=0.35) ---\n\n")
cat(sprintf("  %-12s  %10s  %10s  %10s  %10s  %10s\n",
            "death_scale", "DA_dis", "DA_drug", "DA_saved", "Motor_dis", "Motor_drug"))
cat(paste(rep("-", 72), collapse = ""), "\n")

ds_results <- data.frame()
for (ds in death_scale_vals) {
  res_dis  <- run_scenario(0.90, 0.35, ds, 0,  547)
  res_drug <- run_scenario(0.90, 0.35, ds, 30, 547)

  da_dis  <- res_dis[nrow(res_dis), "DA_neuron"]
  da_drug <- res_drug[nrow(res_drug), "DA_neuron"]
  mo_dis  <- res_dis[nrow(res_dis), "Motor_score"]
  mo_drug <- res_drug[nrow(res_drug), "Motor_score"]

  cat(sprintf("  %-12.2f  %10.4f  %10.4f  %+10.4f  %10.4f  %10.4f\n",
              ds, da_dis, da_drug, da_drug - da_dis, mo_dis, mo_drug))

  ds_results <- rbind(ds_results, data.frame(
    death_scale = ds, DA_disease = da_dis, DA_drug = da_drug,
    DA_saved = da_drug - da_dis, Motor_disease = mo_dis, Motor_drug = mo_drug))
}

cat(sprintf("\n  DA_saved range across death_scale [0.04-0.08]: [%.4f, %.4f]\n",
            min(ds_results$DA_saved), max(ds_results$DA_saved)))
cat(sprintf("  Relative uncertainty: ±%.0f%% around nominal\n\n",
            100 * (max(ds_results$DA_saved) - min(ds_results$DA_saved)) /
            (2 * ds_results$DA_saved[ds_results$death_scale == 0.06])))

## Effect of CI_max on predictions
cat("--- CI_MAX SENSITIVITY (alpha_clear=0.35, death_scale=0.06) ---\n\n")
cat(sprintf("  %-12s  %10s  %10s  %10s\n", "CI_max", "DA_dis", "DA_drug", "DA_saved"))
cat(paste(rep("-", 50), collapse = ""), "\n")

ci_results <- data.frame()
for (ci in CI_max_vals) {
  res_dis  <- run_scenario(ci, 0.35, 0.06, 0,  547)
  res_drug <- run_scenario(ci, 0.35, 0.06, 30, 547)

  da_dis  <- res_dis[nrow(res_dis), "DA_neuron"]
  da_drug <- res_drug[nrow(res_drug), "DA_neuron"]

  cat(sprintf("  %-12.2f  %10.4f  %10.4f  %+10.4f\n", ci, da_dis, da_drug, da_drug - da_dis))

  ci_results <- rbind(ci_results, data.frame(
    CI_max = ci, DA_disease = da_dis, DA_drug = da_drug, DA_saved = da_drug - da_dis))
}

cat(sprintf("\n  DA_saved range across CI_max [0.86-0.94]: [%.4f, %.4f]\n\n",
            min(ci_results$DA_saved), max(ci_results$DA_saved)))

## Effect of alpha_clear on predictions
cat("--- ALPHA_CLEAR SENSITIVITY (CI_max=0.90, death_scale=0.06) ---\n\n")
cat(sprintf("  %-12s  %10s  %10s  %10s\n", "alpha_clear", "DA_dis", "DA_drug", "DA_saved"))
cat(paste(rep("-", 50), collapse = ""), "\n")

ac_results <- data.frame()
for (ac in alpha_clear_vals) {
  res_dis  <- run_scenario(0.90, ac, 0.06, 0,  547)
  res_drug <- run_scenario(0.90, ac, 0.06, 30, 547)

  da_dis  <- res_dis[nrow(res_dis), "DA_neuron"]
  da_drug <- res_drug[nrow(res_drug), "DA_neuron"]

  cat(sprintf("  %-12.2f  %10.4f  %10.4f  %+10.4f\n", ac, da_dis, da_drug, da_drug - da_dis))

  ac_results <- rbind(ac_results, data.frame(
    alpha_clear = ac, DA_disease = da_dis, DA_drug = da_drug, DA_saved = da_drug - da_dis))
}

cat(sprintf("\n  DA_saved range across alpha_clear [0.28-0.42]: [%.4f, %.4f]\n\n",
            min(ac_results$DA_saved), max(ac_results$DA_saved)))

## ============================================================
##  SECTION 4: COMBINED PREDICTION UNCERTAINTY ENVELOPE
##  Joint sweep of all 3 translational parameters
## ============================================================

cat("================================================================\n")
cat("  SECTION 4: JOINT PARAMETER UNCERTAINTY ENVELOPE\n")
cat("  CI_max × alpha_clear × death_scale (125 combinations)\n")
cat("================================================================\n\n")

joint_results <- data.frame()
combos <- expand.grid(CI_max = CI_max_vals,
                      alpha_clear = alpha_clear_vals,
                      death_scale = death_scale_vals)
n_combos <- nrow(combos)
cat(sprintf("Running %d scenarios...\n", n_combos))

for (i in 1:n_combos) {
  if (i %% 25 == 0) cat(sprintf("  %d/%d\n", i, n_combos))

  ci <- combos$CI_max[i]
  ac <- combos$alpha_clear[i]
  ds <- combos$death_scale[i]

  res_dis  <- run_scenario(ci, ac, ds, 0,  547)
  res_drug <- run_scenario(ci, ac, ds, 30, 547)

  da_dis  <- res_dis[nrow(res_dis), "DA_neuron"]
  da_drug <- res_drug[nrow(res_drug), "DA_neuron"]
  mo_dis  <- res_dis[nrow(res_dis), "Motor_score"]
  mo_drug <- res_drug[nrow(res_drug), "Motor_score"]

  joint_results <- rbind(joint_results, data.frame(
    CI_max = ci, alpha_clear = ac, death_scale = ds,
    DA_disease = da_dis, DA_drug = da_drug, DA_saved = da_drug - da_dis,
    Motor_disease = mo_dis, Motor_drug = mo_drug,
    Motor_saved = mo_dis - mo_drug))
}

cat("\n--- PREDICTION UNCERTAINTY ENVELOPE (30 mg, 18 months) ---\n\n")

cat(sprintf("  DA neuron preservation (DA_saved):\n"))
cat(sprintf("    Nominal (0.90, 0.35, 0.06):  %+.4f\n",
            joint_results$DA_saved[joint_results$CI_max == 0.90 &
                                    joint_results$alpha_clear == 0.35 &
                                    joint_results$death_scale == 0.06]))
cat(sprintf("    Minimum across 125 combos:   %+.4f\n", min(joint_results$DA_saved)))
cat(sprintf("    Maximum across 125 combos:   %+.4f\n", max(joint_results$DA_saved)))
cat(sprintf("    Median:                      %+.4f\n", median(joint_results$DA_saved)))
cat(sprintf("    IQR:                         [%+.4f, %+.4f]\n",
            quantile(joint_results$DA_saved, 0.25),
            quantile(joint_results$DA_saved, 0.75)))
cat(sprintf("    95%% envelope:                [%+.4f, %+.4f]\n\n",
            quantile(joint_results$DA_saved, 0.025),
            quantile(joint_results$DA_saved, 0.975)))

cat(sprintf("  Motor score improvement (Motor_saved):\n"))
cat(sprintf("    Nominal:    %+.4f\n",
            joint_results$Motor_saved[joint_results$CI_max == 0.90 &
                                      joint_results$alpha_clear == 0.35 &
                                      joint_results$death_scale == 0.06]))
cat(sprintf("    95%% envelope: [%+.4f, %+.4f]\n\n",
            quantile(joint_results$Motor_saved, 0.025),
            quantile(joint_results$Motor_saved, 0.975)))

## Approximate unpaired d range
## Use observed VPop SD (pooled) as denominator
pooled_sd_DA <- sqrt((var(pbo$DA_neuron) + var(d30$DA_neuron)) / 2)
pooled_sd_Mo <- sqrt((var(pbo$Motor_score) + var(d30$Motor_score)) / 2)

cat(sprintf("  Approximate unpaired Cohen's d range (using observed pooled SD):\n"))
cat(sprintf("    DA_neuron d:   [%.3f, %.3f]  (nominal %.3f)\n",
            min(joint_results$DA_saved) / pooled_sd_DA,
            max(joint_results$DA_saved) / pooled_sd_DA,
            joint_results$DA_saved[joint_results$CI_max == 0.90 &
                                    joint_results$alpha_clear == 0.35 &
                                    joint_results$death_scale == 0.06] / pooled_sd_DA))
cat(sprintf("    Motor_score d: [%.3f, %.3f]  (nominal %.3f)\n\n",
            min(joint_results$Motor_saved) / pooled_sd_Mo,
            max(joint_results$Motor_saved) / pooled_sd_Mo,
            joint_results$Motor_saved[joint_results$CI_max == 0.90 &
                                      joint_results$alpha_clear == 0.35 &
                                      joint_results$death_scale == 0.06] / pooled_sd_Mo))

## ============================================================
##  SECTION 5: PLOTS
## ============================================================

cat("================================================================\n")
cat("  SECTION 5: PLOTS\n")
cat("================================================================\n\n")

theme_pub <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

## Plot 1: Bootstrap distributions of unpaired d
boot_d_data <- data.frame()
for (dose_info in doses) {
  for (ep in endpoints) {
    key <- paste(dose_info$name, ep, sep = "_")
    bd <- boot_stats[[key]]
    boot_d_data <- rbind(boot_d_data, data.frame(
      Dose = dose_info$name, Endpoint = ep, d_unpaired = bd$d_unpaired))
  }
}

p1 <- ggplot(boot_d_data, aes(x = d_unpaired, fill = Dose)) +
  geom_density(alpha = 0.5, colour = "grey30") +
  facet_wrap(~Endpoint, scales = "free_x",
             labeller = labeller(Endpoint = c(DA_neuron = "DA Neuron Survival",
                                               Motor_score = "Motor Score"))) +
  scale_fill_manual(values = c("10mg" = "#66c2a5", "30mg" = "#fc8d62", "60mg" = "#8da0cb")) +
  geom_vline(xintercept = 0.15, linetype = "dashed", colour = "red", alpha = 0.7) +
  geom_vline(xintercept = 0.20, linetype = "dashed", colour = "blue", alpha = 0.7) +
  annotate("text", x = 0.15, y = Inf, label = "SPARK", vjust = 2, hjust = 1.1,
           colour = "red", size = 3) +
  annotate("text", x = 0.20, y = Inf, label = "ADAGIO", vjust = 2, hjust = -0.1,
           colour = "blue", size = 3) +
  labs(title = "Bootstrap Distribution of Unpaired Cohen's d (B=10,000)",
       x = "Unpaired Cohen's d", y = "Density",
       fill = "Dose") +
  theme_pub

ggsave("UA_bootstrap_d.png", p1, width = 10, height = 5, dpi = 300)
cat("Saved: UA_bootstrap_d.png\n")

## Plot 2: Death scale sensitivity tornado
ds_long <- ds_results %>%
  mutate(DA_saved_pct = 100 * (DA_saved - DA_saved[death_scale == 0.06]) /
           DA_saved[death_scale == 0.06])

p2 <- ggplot(ds_results, aes(x = factor(death_scale), y = DA_saved)) +
  geom_col(fill = "#3288bd", width = 0.6) +
  geom_hline(yintercept = ds_results$DA_saved[ds_results$death_scale == 0.06],
             linetype = "dashed", colour = "red") +
  annotate("text", x = 5.4, y = ds_results$DA_saved[ds_results$death_scale == 0.06],
           label = "Nominal", vjust = -0.5, colour = "red", size = 3.5) +
  labs(title = "DA Neuron Preservation vs Death Rate Scaling",
       subtitle = "30 mg, 18 months | CI_max=0.90, alpha_clear=0.35",
       x = "death_scale_human", y = "DA neurons preserved (drug - placebo)") +
  theme_pub

ggsave("UA_death_scale_sensitivity.png", p2, width = 7, height = 5, dpi = 300)
cat("Saved: UA_death_scale_sensitivity.png\n")

## Plot 3: Joint uncertainty envelope
p3 <- ggplot(joint_results, aes(x = DA_saved)) +
  geom_histogram(bins = 25, fill = "#3288bd", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = joint_results$DA_saved[joint_results$CI_max == 0.90 &
                                                   joint_results$alpha_clear == 0.35 &
                                                   joint_results$death_scale == 0.06],
             colour = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(joint_results$DA_saved, c(0.025, 0.975)),
             colour = "orange", linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = joint_results$DA_saved[joint_results$CI_max == 0.90 &
                                                joint_results$alpha_clear == 0.35 &
                                                joint_results$death_scale == 0.06],
           y = Inf, label = "Nominal", vjust = 2, colour = "red", size = 3.5) +
  labs(title = "Prediction Uncertainty Envelope: DA Neuron Preservation",
       subtitle = "125 combinations of CI_max x alpha_clear x death_scale | 30 mg, 18 months",
       x = "DA neurons preserved (drug - placebo)", y = "Count") +
  theme_pub

ggsave("UA_prediction_envelope.png", p3, width = 8, height = 5, dpi = 300)
cat("Saved: UA_prediction_envelope.png\n")

## Plot 4: Prediction interval per dose (individual patient)
pi_data <- data.frame()
for (dose_info in doses) {
  delta <- dose_info$data$DA_neuron[dose_info$data$id %in% ids_common] - pbo$DA_neuron
  pi_data <- rbind(pi_data, data.frame(Dose = dose_info$name, delta_DA = delta))
}

p4 <- ggplot(pi_data, aes(x = Dose, y = delta_DA, fill = Dose)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0.05, linetype = "dotted", colour = "red") +
  annotate("text", x = 3.4, y = 0.05, label = "MCID (5%)", vjust = -0.5,
           colour = "red", size = 3) +
  scale_fill_manual(values = c("10mg" = "#66c2a5", "30mg" = "#fc8d62", "60mg" = "#8da0cb")) +
  labs(title = "Individual Patient DA Neuron Preservation",
       subtitle = "Matched pairs (drug - placebo) | 18 months",
       x = "Dose", y = "DA neuron change (drug - placebo)") +
  theme_pub +
  theme(legend.position = "none")

ggsave("UA_patient_prediction_interval.png", p4, width = 6, height = 5, dpi = 300)
cat("Saved: UA_patient_prediction_interval.png\n")

## ============================================================
##  SECTION 6: SUMMARY TABLE
## ============================================================

cat("\n================================================================\n")
cat("  SECTION 6: SUMMARY — CONFIDENCE INTERVALS ON ALL KEY PREDICTIONS\n")
cat("================================================================\n\n")

## Extract key CIs from bootstrap
get_ci <- function(dose_name, ep, stat) {
  r <- results_ci[results_ci$Dose == dose_name &
                    results_ci$Endpoint == ep &
                    results_ci$Statistic == stat, ]
  sprintf("%.3f [%.3f, %.3f]", r$Observed, r$CI_lo, r$CI_hi)
}

cat("--- EFFECT SIZE 95% CONFIDENCE INTERVALS (Bootstrap, B=10000) ---\n\n")
cat(sprintf("  %-8s  %-12s  %-28s  %-28s\n", "Dose", "Statistic", "DA_neuron", "Motor_score"))
cat(paste(rep("-", 84), collapse = ""), "\n")

for (stat in c("mean_delta", "d_paired", "d_unpaired", "cliff", "resp_rate")) {
  for (dose_name in c("10mg", "30mg", "60mg")) {
    da_ci <- get_ci(dose_name, "DA_neuron", stat)
    mo_ci <- get_ci(dose_name, "Motor_score", stat)
    cat(sprintf("  %-8s  %-12s  %-28s  %-28s\n", dose_name, stat, da_ci, mo_ci))
  }
  cat("\n")
}

cat("--- PARAMETER UNCERTAINTY ON PREDICTIONS (Joint 125-combo sweep) ---\n\n")

nom_da <- joint_results$DA_saved[joint_results$CI_max == 0.90 &
                                   joint_results$alpha_clear == 0.35 &
                                   joint_results$death_scale == 0.06]
nom_mo <- joint_results$Motor_saved[joint_results$CI_max == 0.90 &
                                     joint_results$alpha_clear == 0.35 &
                                     joint_results$death_scale == 0.06]

cat(sprintf("  DA_saved (30mg, 18mo):     %.4f  [%.4f, %.4f] (95%% envelope)\n",
            nom_da, quantile(joint_results$DA_saved, 0.025),
            quantile(joint_results$DA_saved, 0.975)))
cat(sprintf("  Motor_saved (30mg, 18mo):  %.4f  [%.4f, %.4f] (95%% envelope)\n",
            nom_mo, quantile(joint_results$Motor_saved, 0.025),
            quantile(joint_results$Motor_saved, 0.975)))
cat(sprintf("  Approx d_unpaired range:   %.3f  [%.3f, %.3f]\n\n",
            nom_da / pooled_sd_DA,
            min(joint_results$DA_saved) / pooled_sd_DA,
            max(joint_results$DA_saved) / pooled_sd_DA))

cat("================================================================\n")
cat("  UNCERTAINTY ANALYSIS COMPLETE\n")
cat("================================================================\n")
