##=============================================================================
## HUMAN TRANSLATION
## Allometric PK Scaling | Human Progression | Virtual Population | Clinical Trial
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
##
## Sections:
##   1. Allometric PK scaling (mouse → human)
##   2. Human disease timescale calibration (20-year trajectory)
##   3. Virtual population analysis (SA-informed, N=200)
##   4. Clinical trial simulation (placebo-controlled, 3 dose levels)
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

cat("================================================================\n")
cat("  HUMAN TRANSLATION\n")
cat("  Allometric PK | Human Timescale | VPop | Clinical Trial\n")
cat("================================================================\n\n")

## Calibrated mouse parameters (from Integrated Model)
CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)

## Human death rate scaling (recalibrated for multi-pathway death model)
death_scale_human <- 0.06


##=============================================================================
## SECTION 1: ALLOMETRIC PK SCALING (Mouse → Human)
##=============================================================================

cat("================================================================\n")
cat("  SECTION 1: ALLOMETRIC PK SCALING\n")
cat("================================================================\n\n")

BW_mouse <- 0.025
BW_human <- 70.0
scale_025 <- (BW_human / BW_mouse)^(-0.25)

cat(sprintf("BW ratio = %.0f, rate scaling BW^(-0.25) = %.4f\n",
            BW_human / BW_mouse, scale_025))
cat(sprintf("Predicted t½ scaling BW^(+0.25) = %.1fx\n\n", 1 / scale_025))

## Human PK rate constants
## Elimination/distribution: allometric BW^(-0.25)
## BBB: physiology-based (monkey CSF:plasma = 0.23 vs rat 0.55, ~2x reduction)
## Intracellular (mito): not scaled (biochemistry same across species)
HUMAN_PK <- list(
  ka          = 0.40,     # SC absorption (human tmax ~1.5–2h for peptides)
  ke_plasma   = 0.032,    # t½ ≈ 22h (mouse 3h × 7.3 = 22h)
  k_12        = 0.041,    # Plasma→peripheral (allometric)
  k_21        = 0.027,    # Peripheral→plasma (allometric)
  k_brain_in  = 0.015,    # BBB penetration (reduced; monkey:rat ~0.42)
  k_brain_out = 0.010,    # Brain efflux (longer human brain retention)
  k_mito_in   = 0.50,     # IMM binding (intracellular, same)
  k_mito_out  = 0.020     # IMM release (intracellular, same)
)

mouse_pk_vals <- c(1.5, 0.231, 0.30, 0.20, 0.040, 0.016, 0.50, 0.020)
pk_names <- names(HUMAN_PK)
cat("Parameter scaling:\n")
cat(sprintf("  %-14s  %8s  %8s  %10s  %10s\n",
            "Parameter", "Mouse", "Human", "t½_mouse", "t½_human"))
for (i in seq_along(pk_names)) {
  cat(sprintf("  %-14s  %8.3f  %8.3f  %10.1fh  %10.1fh\n",
              pk_names[i], mouse_pk_vals[i], HUMAN_PK[[pk_names[i]]],
              log(2) / mouse_pk_vals[i], log(2) / HUMAN_PK[[pk_names[i]]]))
}

## FDA Km dose conversion (guidance for industry: HED = animal × Km_animal/Km_human)
Km_mouse <- 3; Km_human <- 37
cat(sprintf("\nFDA Km dose conversion factor: %.4f\n", Km_mouse / Km_human))
cat(sprintf("  Mouse 0.5 mg/kg → Human %.3f mg/kg = %.1f mg (70 kg)\n",
            0.5 * Km_mouse / Km_human, 0.5 * Km_mouse / Km_human * BW_human))
cat(sprintf("  Mouse 5.0 mg/kg → Human %.3f mg/kg = %.1f mg (70 kg)\n\n",
            5.0 * Km_mouse / Km_human, 5.0 * Km_mouse / Km_human * BW_human))

## Clinical dose levels (SC daily)
DOSES_MG <- c(10, 30, 100)
DOSES_MGKG <- DOSES_MG / BW_human
cat("Planned clinical doses:\n")
for (i in seq_along(DOSES_MG)) {
  cat(sprintf("  %3d mg SC daily = %.4f mg/kg (Dose_norm = %.5f)\n",
              DOSES_MG[i], DOSES_MGKG[i], DOSES_MGKG[i] / 5.0))
}

## Find human C_ref: trough C_mito at 30 mg SS (HED of mouse 5 mg/kg)
cat("\nCalibrating human C_ref at 30 mg/day for 120 days...\n")
po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = DOSES_MGKG[2], duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day_ref <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day_ref$C_mito)
cat(sprintf("  C_ref_human = %.6f (trough C_mito at 30 mg SS)\n", C_ref_human))
cat(sprintf("  cf. C_ref_mouse = 4.43 (trough at 5 mg/kg SS)\n"))
cat(sprintf("  Ratio: %.1fx lower (reflects lower BBB penetration at human doses)\n\n",
            4.43 / C_ref_human))

## Helper: run model with human PK
run_human <- function(dose_mg, duration_days, CI_max = 1.0, alpha_clear = 1.0,
                      dt = 6, parms_extra = NULL, quiet = TRUE) {
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
                 parms_override = po, quiet = quiet)
}

## Compare mouse vs human PK (7-day profiles at reference doses)
cat("Comparing single-dose PK profiles...\n")
res_mouse_pk <- run_integrated(dose_mg_kg = 5.0, duration_days = 7, dt = 0.25,
                               CI_max = 1.0, alpha_clear = 1.0, quiet = TRUE)
res_human_pk <- run_human(dose_mg = 30, duration_days = 7, dt = 0.25)

pk_m <- res_mouse_pk %>%
  select(time, C_plasma, C_brain, C_mito) %>%
  mutate(species = "Mouse (5 mg/kg IP)") %>%
  pivot_longer(c(C_plasma, C_brain, C_mito), names_to = "comp", values_to = "conc")
pk_h <- res_human_pk %>%
  select(time, C_plasma, C_brain, C_mito) %>%
  mutate(species = "Human (30 mg SC)") %>%
  pivot_longer(c(C_plasma, C_brain, C_mito), names_to = "comp", values_to = "conc")
pk_all <- rbind(pk_m, pk_h)
pk_all$comp <- factor(pk_all$comp, levels = c("C_plasma", "C_brain", "C_mito"),
                      labels = c("Plasma", "Brain ISF", "Mitochondrial"))

p_pk <- ggplot(pk_all, aes(x = time, y = conc, color = species)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ comp, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Mouse (5 mg/kg IP)" = "#d62728",
                                "Human (30 mg SC)" = "#1f77b4")) +
  labs(title = "Allometric PK Scaling: Mouse vs Human",
       subtitle = "Single dose at species-equivalent reference levels (FDA Km)",
       x = "Time (hours)", y = "Concentration (model units)", color = "Species") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave("HT_PK_comparison.png", p_pk, width = 14, height = 5, dpi = 300, bg = "white")
cat("Saved: HT_PK_comparison.png\n\n")


##=============================================================================
## SECTION 2: HUMAN DISEASE TIMESCALE
##=============================================================================

cat("================================================================\n")
cat("  SECTION 2: HUMAN DISEASE TIMESCALE CALIBRATION\n")
cat("  Target: prodromal 5–10 yr, Bernheimer threshold ~10 yr\n")
cat("================================================================\n\n")

## Human PD is slower than mouse models:
## Mouse: CI_max=0.74, alpha_clear=0.25 → 33% DA loss in 5 weeks
## Human: prodromal ~5–10 yr, symptom onset at ~60–70% DA loss
## Strategy: milder disease modifiers → slower vicious cycle engagement

CI_max_grid   <- seq(0.86, 0.96, by = 0.01)
alpha_clr_grid <- seq(0.25, 0.50, by = 0.05)

cat(sprintf("Sweeping %d CI_max × %d alpha_clear = %d combos (5-year screen)...\n",
            length(CI_max_grid), length(alpha_clr_grid),
            length(CI_max_grid) * length(alpha_clr_grid)))

human_sweep <- data.frame()
for (cm in CI_max_grid) {
  for (ac in alpha_clr_grid) {
    res <- run_human(0, 1825, CI_max = cm, alpha_clear = ac, dt = 24)
    ss <- res[nrow(res), ]

    idx_2y <- which.min(abs(res$time - 730 * 24))

    human_sweep <- rbind(human_sweep, data.frame(
      CI_max = cm, alpha_clear = ac,
      DA_2y = res$DA_neuron[idx_2y],
      DA_5y = ss$DA_neuron,
      DA_loss_5y = (1 - ss$DA_neuron) * 100,
      Motor_5y = ss$Motor_score,
      CI_5y = ss$CI_activity,
      aSyn_5y = ss$aSyn_olig
    ))
  }
}

## Target: DA loss 20–35% at 5 years (on track for 50–65% at 10 years)
human_sweep$err <- abs(human_sweep$DA_loss_5y - 28)
good_h <- human_sweep[human_sweep$DA_loss_5y > 15 &
                       human_sweep$DA_loss_5y < 45 &
                       human_sweep$CI_5y > 0.30, ]

if (nrow(good_h) > 0) {
  cat(sprintf("\n%d candidates found. Top 10 by proximity to 28%% DA loss at 5y:\n", nrow(good_h)))
  top_h <- head(good_h[order(good_h$err), ], 10)
  print(top_h[, c("CI_max", "alpha_clear", "DA_loss_5y", "Motor_5y",
                   "CI_5y", "aSyn_5y")], row.names = FALSE, digits = 3)
  best_h <- good_h[which.min(good_h$err), ]
} else {
  cat("\nNo ideal candidates — using closest match:\n")
  best_h <- human_sweep[which.min(human_sweep$err), ]
}

CI_max_human <- best_h$CI_max
alpha_clear_human <- best_h$alpha_clear
cat(sprintf("\n** Selected: CI_max_human=%.2f, alpha_clear_human=%.2f **\n",
            CI_max_human, alpha_clear_human))
cat(sprintf("   DA loss 5y = %.1f%%, Motor 5y = %.3f, CI 5y = %.3f\n\n",
            best_h$DA_loss_5y, best_h$Motor_5y, best_h$CI_5y))

## Full 20-year trajectories (healthy, disease, drug 30 mg)
cat("Running 20-year trajectories (3 scenarios)...\n")
res_h20 <- run_human(0, 7300, CI_max = 1.0, alpha_clear = 1.0, dt = 24)
res_d20 <- run_human(0, 7300, CI_max = CI_max_human,
                     alpha_clear = alpha_clear_human, dt = 24)
res_rx20 <- run_human(30, 7300, CI_max = CI_max_human,
                      alpha_clear = alpha_clear_human, dt = 24)

## Milestone table
cat("\nHuman PD progression milestones:\n")
cat(sprintf("  %-6s  %-8s  %-8s  %-8s  %-8s  %-10s\n",
            "Year", "DA_dis", "DA_drug", "Motor_d", "Motor_rx", "DA_saved"))
for (yr in c(1, 2, 5, 10, 15, 20)) {
  t_h <- yr * 365 * 24
  id <- which.min(abs(res_d20$time - t_h))
  ir <- which.min(abs(res_rx20$time - t_h))
  cat(sprintf("  %-6d  %.4f    %.4f    %.4f    %.4f    %+.1f%%\n",
              yr, res_d20$DA_neuron[id], res_rx20$DA_neuron[ir],
              res_d20$Motor_score[id], res_rx20$Motor_score[ir],
              (res_rx20$DA_neuron[ir] / res_d20$DA_neuron[id] - 1) * 100))
}

## Find Bernheimer threshold crossing
bern_idx <- which(res_d20$DA_neuron < 0.40)[1]
if (!is.na(bern_idx)) {
  bern_years <- res_d20$time[bern_idx] / 24 / 365
  cat(sprintf("\nBernheimer threshold (60%% DA loss) crossed at %.1f years\n", bern_years))
} else {
  cat("\nBernheimer threshold not crossed in 20 years\n")
}

bern_idx_rx <- which(res_rx20$DA_neuron < 0.40)[1]
if (!is.na(bern_idx_rx)) {
  bern_years_rx <- res_rx20$time[bern_idx_rx] / 24 / 365
  cat(sprintf("With 30 mg drug: threshold crossing delayed to %.1f years\n", bern_years_rx))
  cat(sprintf("Drug extends time to motor onset by ~%.1f years\n\n",
              bern_years_rx - bern_years))
} else {
  cat("With 30 mg drug: threshold not crossed in 20 years\n\n")
}

## Plot: 20-year trajectory
res_h20$scenario  <- "Healthy"
res_d20$scenario  <- "Disease (no drug)"
res_rx20$scenario <- "Disease + SBT-272 (30 mg)"
all20 <- rbind(res_h20, res_d20, res_rx20)
all20$time_years <- all20$time / 24 / 365
all20$scenario <- factor(all20$scenario,
  levels = c("Healthy", "Disease (no drug)", "Disease + SBT-272 (30 mg)"))

colors3 <- c("Healthy" = "#2ca02c", "Disease (no drug)" = "#d62728",
             "Disease + SBT-272 (30 mg)" = "#1f77b4")

prog_vars <- c("DA_neuron", "Motor_score", "CI_activity", "aSyn_olig")
long20 <- all20 %>%
  select(time_years, scenario, all_of(prog_vars)) %>%
  pivot_longer(cols = all_of(prog_vars), names_to = "variable", values_to = "value")
long20$variable <- factor(long20$variable, levels = prog_vars,
  labels = c("DA Neuron Survival", "Motor Score",
             "Complex I Activity", "\u03b1-Synuclein Oligomers"))

p_prog <- ggplot(long20, aes(x = time_years, y = value, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(data = data.frame(
    variable = factor("DA Neuron Survival", levels = levels(long20$variable)),
    yint = 0.40), aes(yintercept = yint), linetype = "dashed", color = "gray50") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = colors3) +
  labs(title = "Human PD Progression: 20-Year Trajectory",
       subtitle = sprintf("Human disease modifiers: CI_max=%.2f, alpha_clear=%.2f | Allometric PK",
                          CI_max_human, alpha_clear_human),
       x = "Time (years)", y = "Normalised level", color = "Scenario",
       caption = "Dashed line: Bernheimer threshold (60% DA loss \u2192 motor symptom onset)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave("HT_human_progression.png", p_prog, width = 12, height = 8, dpi = 300, bg = "white")
cat("Saved: HT_human_progression.png\n\n")


##=============================================================================
## SECTION 3: VIRTUAL POPULATION
##=============================================================================

cat("================================================================\n")
cat("  SECTION 3: VIRTUAL POPULATION ANALYSIS\n")
cat("  N=200 patients | SA-informed parameter variability\n")
cat("================================================================\n\n")

set.seed(42)
N_POP <- 200
TRIAL_DAYS <- 547   # 18 months

## Log-normal sampler
sample_ln <- function(n, nom, cv) {
  sigma <- sqrt(log(1 + cv^2))
  mu <- log(nom) - sigma^2 / 2
  exp(rnorm(n, mu, sigma))
}

## Parameter distributions (based on Sobol ST rankings)
vpop_defs <- data.frame(
  name = c("k_ROS_basal", "K_mPTP", "k_SOD2", "k_agg",
           "k_shunt", "k_clear", "k_death", "K_death"),
  nominal = c(0.122, 0.40, 1.0, 0.009, 0.125, 0.241, 0.001, 0.25),
  cv = c(0.25, 0.25, 0.25, 0.30, 0.20, 0.25, 0.30, 0.25),
  stringsAsFactors = FALSE
)

## Generate virtual patients
vpop <- data.frame(id = 1:N_POP)
for (r in 1:nrow(vpop_defs)) {
  vpop[[vpop_defs$name[r]]] <- sample_ln(N_POP, vpop_defs$nominal[r], vpop_defs$cv[r])
}
vpop$CI_max_h <- pmin(0.96, pmax(0.80, rnorm(N_POP, CI_max_human, 0.03)))
vpop$alpha_clear_h <- pmin(0.55, pmax(0.20, rnorm(N_POP, alpha_clear_human, 0.05)))

cat(sprintf("Generated %d virtual patients\n", N_POP))
cat("Parameter distributions (median [IQR]):\n")
for (nm in vpop_defs$name) {
  v <- vpop[[nm]]
  cat(sprintf("  %-14s: %.4f [%.4f, %.4f]\n", nm, median(v),
              quantile(v, 0.25), quantile(v, 0.75)))
}
cat(sprintf("  %-14s: %.3f [%.3f, %.3f]\n", "CI_max_h",
            median(vpop$CI_max_h), quantile(vpop$CI_max_h, 0.25),
            quantile(vpop$CI_max_h, 0.75)))
cat(sprintf("  %-14s: %.3f [%.3f, %.3f]\n\n", "alpha_clear_h",
            median(vpop$alpha_clear_h), quantile(vpop$alpha_clear_h, 0.25),
            quantile(vpop$alpha_clear_h, 0.75)))

## Run: Placebo + 30 mg
cat(sprintf("Running %d patients × 2 arms (placebo, 30 mg) for %d days...\n",
            N_POP, TRIAL_DAYS))

vpop_results <- data.frame()
for (i in 1:N_POP) {
  if (i %% 50 == 0) cat(sprintf("  VPop: %d/%d patients done\n", i, N_POP))

  po <- as.list(vpop[i, vpop_defs$name])

  res_pbo <- tryCatch(
    run_human(0, TRIAL_DAYS, CI_max = vpop$CI_max_h[i],
              alpha_clear = vpop$alpha_clear_h[i], dt = 24, parms_extra = po),
    error = function(e) NULL)
  res_drg <- tryCatch(
    run_human(30, TRIAL_DAYS, CI_max = vpop$CI_max_h[i],
              alpha_clear = vpop$alpha_clear_h[i], dt = 24, parms_extra = po),
    error = function(e) NULL)

  if (is.null(res_pbo) || is.null(res_drg)) next

  sp <- res_pbo[nrow(res_pbo), ]
  sd <- res_drg[nrow(res_drg), ]

  vpop_results <- rbind(vpop_results, data.frame(
    id = rep(i, 2), arm = c("Placebo", "SBT-272 30mg"),
    DA_neuron = c(sp$DA_neuron, sd$DA_neuron),
    Motor_score = c(sp$Motor_score, sd$Motor_score),
    CI_activity = c(sp$CI_activity, sd$CI_activity),
    aSyn_olig = c(sp$aSyn_olig, sd$aSyn_olig),
    mROS = c(sp$mROS, sd$mROS)
  ))
}
cat(sprintf("  VPop: %d/%d patients done\n\n", N_POP, N_POP))

n_completed <- length(unique(vpop_results$id))
cat(sprintf("Completed: %d/%d patients (%.0f%%)\n\n", n_completed, N_POP,
            n_completed / N_POP * 100))

## Summary
cat("Virtual Population Results (18 months):\n")
for (arm_name in c("Placebo", "SBT-272 30mg")) {
  sub <- vpop_results[vpop_results$arm == arm_name, ]
  cat(sprintf("\n  %s (n=%d):\n", arm_name, nrow(sub)))
  for (v in c("DA_neuron", "Motor_score", "CI_activity", "aSyn_olig")) {
    vals <- sub[[v]]
    cat(sprintf("    %-14s: mean=%.3f sd=%.3f median=%.3f [%.3f, %.3f]\n",
                v, mean(vals), sd(vals), median(vals),
                quantile(vals, 0.25), quantile(vals, 0.75)))
  }
}

## Plot: VPop endpoint distributions
vpop_ep <- vpop_results %>%
  select(id, arm, DA_neuron, Motor_score) %>%
  pivot_longer(c(DA_neuron, Motor_score), names_to = "endpoint", values_to = "value")
vpop_ep$endpoint <- factor(vpop_ep$endpoint,
  levels = c("DA_neuron", "Motor_score"),
  labels = c("DA Neuron Survival", "Motor Score"))
vpop_ep$arm <- factor(vpop_ep$arm, levels = c("Placebo", "SBT-272 30mg"))

p_vpop <- ggplot(vpop_ep, aes(x = value, fill = arm)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  facet_wrap(~ endpoint, scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Placebo" = "#d62728", "SBT-272 30mg" = "#1f77b4")) +
  labs(title = "Virtual Population: Endpoint Distributions at 18 Months",
       subtitle = sprintf("N=%d patients per arm | SA-informed variability", n_completed),
       x = "Endpoint Value", y = "Density", fill = "Treatment Arm") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave("HT_vpop_distribution.png", p_vpop, width = 12, height = 5, dpi = 300, bg = "white")
cat("\nSaved: HT_vpop_distribution.png\n\n")


##=============================================================================
## SECTION 4: CLINICAL TRIAL SIMULATION
##=============================================================================

cat("================================================================\n")
cat("  SECTION 4: CLINICAL TRIAL SIMULATION\n")
cat("  Parallel-group | Placebo + 3 doses | 18 months\n")
cat("================================================================\n\n")

## Run additional dose arms (10 mg, 100 mg) for same virtual patients
completed_ids <- unique(vpop_results$id)

for (dose in c(10, 100)) {
  dose_label <- sprintf("SBT-272 %dmg", dose)
  cat(sprintf("Running %s arm (%d patients)...\n", dose_label, length(completed_ids)))

  for (i in completed_ids) {
    if (i %% 50 == 0) cat(sprintf("  %s: %d/%d\n", dose_label, i, max(completed_ids)))

    po <- as.list(vpop[i, vpop_defs$name])

    res <- tryCatch(
      run_human(dose, TRIAL_DAYS, CI_max = vpop$CI_max_h[i],
                alpha_clear = vpop$alpha_clear_h[i], dt = 24, parms_extra = po),
      error = function(e) NULL)
    if (is.null(res)) next

    ss <- res[nrow(res), ]
    vpop_results <- rbind(vpop_results, data.frame(
      id = i, arm = dose_label,
      DA_neuron = ss$DA_neuron, Motor_score = ss$Motor_score,
      CI_activity = ss$CI_activity, aSyn_olig = ss$aSyn_olig,
      mROS = ss$mROS
    ))
  }
  cat(sprintf("  %s: done\n", dose_label))
}

## Treatment effect analysis
cat("\n--- TREATMENT EFFECTS (vs Placebo, paired) ---\n\n")

arms <- c("Placebo", "SBT-272 10mg", "SBT-272 30mg", "SBT-272 100mg")
vpop_results$arm <- factor(vpop_results$arm, levels = arms)
pbo <- vpop_results[vpop_results$arm == "Placebo", ]

effect_table <- data.frame()
for (arm_name in arms[-1]) {
  trt <- vpop_results[vpop_results$arm == arm_name, ]
  trt <- trt[trt$id %in% pbo$id, ]
  pbo_matched <- pbo[pbo$id %in% trt$id, ]

  trt <- trt[order(trt$id), ]
  pbo_matched <- pbo_matched[order(pbo_matched$id), ]

  for (ep in c("DA_neuron", "Motor_score")) {
    delta <- trt[[ep]] - pbo_matched[[ep]]
    n_pairs <- length(delta)
    mean_d <- mean(delta); sd_d <- sd(delta)
    se_d <- sd_d / sqrt(n_pairs)
    ci_lo <- mean_d - 1.96 * se_d
    ci_hi <- mean_d + 1.96 * se_d

    d_paired <- abs(mean_d) / sd_d
    tt <- t.test(delta, mu = 0)

    if (ep == "DA_neuron") {
      resp <- sum(delta > 0) / n_pairs * 100
    } else {
      resp <- sum(delta < 0) / n_pairs * 100
    }

    effect_table <- rbind(effect_table, data.frame(
      arm = arm_name, endpoint = ep, n = n_pairs,
      mean_delta = mean_d, se = se_d, ci_lo = ci_lo, ci_hi = ci_hi,
      cohens_d = d_paired, p_value = tt$p.value, responder_pct = resp
    ))
  }
}

cat(sprintf("  %-18s %-14s  N   %8s %8s  [%7s, %7s]  d=%5s  %8s  Resp%%\n",
            "Arm", "Endpoint", "Mean\u0394", "SE", "CI_lo", "CI_hi", "Cohen", "p-value"))
for (r in 1:nrow(effect_table)) {
  with(effect_table[r, ], cat(sprintf(
    "  %-18s %-14s %3d  %+8.4f %8.4f  [%+7.4f, %+7.4f]  d=%.2f  %.1e  %.0f%%\n",
    arm, endpoint, n, mean_delta, se, ci_lo, ci_hi, cohens_d, p_value, responder_pct)))
}

## Power analysis
cat("\n--- POWER ANALYSIS (80%% power, alpha=0.05, paired t-test) ---\n\n")
for (r in 1:nrow(effect_table)) {
  d <- effect_table$cohens_d[r]
  if (d > 0.01) {
    n_needed <- ceiling((qnorm(0.975) + qnorm(0.80))^2 / d^2)
  } else {
    n_needed <- NA
  }
  cat(sprintf("  %-18s %-14s: d=%.3f → N=%s per arm\n",
              effect_table$arm[r], effect_table$endpoint[r], d,
              ifelse(!is.na(n_needed) & n_needed < 10000,
                     as.character(n_needed), ">10,000")))
}

## Plot 1: Forest plot (Motor score)
forest_data <- effect_table[effect_table$endpoint == "Motor_score", ]
forest_data$arm <- factor(forest_data$arm, levels = rev(arms[-1]))

p_forest <- ggplot(forest_data, aes(x = mean_delta, y = arm)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.8) +
  geom_point(size = 3, color = "#1f77b4") +
  labs(title = "Clinical Trial: Motor Score Treatment Effect (18 Months)",
       subtitle = sprintf("Paired difference vs placebo (N=%d per arm, 95%% CI)", n_completed),
       x = "Change in Motor Score vs Placebo (negative = improvement)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())
ggsave("HT_trial_forest.png", p_forest, width = 10, height = 4, dpi = 300, bg = "white")
cat("\nSaved: HT_trial_forest.png\n")

## Plot 2: Waterfall (DA neuron, 30 mg arm)
trt30 <- vpop_results[vpop_results$arm == "SBT-272 30mg", ]
trt30 <- trt30[order(trt30$id), ]
pbo_w <- pbo[pbo$id %in% trt30$id, ]
pbo_w <- pbo_w[order(pbo_w$id), ]

wf <- data.frame(
  id = trt30$id,
  DA_benefit = (trt30$DA_neuron - pbo_w$DA_neuron) * 100
)
wf <- wf[order(wf$DA_benefit), ]
wf$rank <- 1:nrow(wf)
wf$resp <- ifelse(wf$DA_benefit > 0, "Responder", "Non-responder")

p_wf <- ggplot(wf, aes(x = rank, y = DA_benefit, fill = resp)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("Responder" = "#1f77b4", "Non-responder" = "#d62728")) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(title = "Individual Patient Response: SBT-272 30 mg vs Placebo",
       subtitle = sprintf("DA neuron preservation at 18 months (N=%d)", nrow(wf)),
       x = "Patient (ranked)", y = "DA Neuron Change vs Placebo (%)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("HT_trial_waterfall.png", p_wf, width = 12, height = 5, dpi = 300, bg = "white")
cat("Saved: HT_trial_waterfall.png\n")

## Plot 3: Dose-response
dr <- vpop_results %>%
  group_by(arm) %>%
  summarise(DA_m = mean(DA_neuron), DA_se = sd(DA_neuron) / sqrt(n()),
            Mo_m = mean(Motor_score), Mo_se = sd(Motor_score) / sqrt(n()),
            .groups = "drop")
dr$dose_mg <- c(0, 10, 30, 100)

p_dr <- ggplot(dr, aes(x = dose_mg)) +
  geom_ribbon(aes(ymin = DA_m - 1.96 * DA_se, ymax = DA_m + 1.96 * DA_se),
              fill = "#1f77b4", alpha = 0.2) +
  geom_line(aes(y = DA_m), color = "#1f77b4", linewidth = 1) +
  geom_point(aes(y = DA_m), color = "#1f77b4", size = 3) +
  geom_errorbar(aes(ymin = DA_m - 1.96 * DA_se, ymax = DA_m + 1.96 * DA_se),
                width = 3, color = "#1f77b4") +
  scale_x_continuous(breaks = c(0, 10, 30, 100)) +
  labs(title = "Dose-Response: DA Neuron Survival at 18 Months",
       subtitle = sprintf("Mean \u00b1 95%% CI | N=%d per arm", n_completed),
       x = "Daily Dose (mg SC)", y = "DA Neuron Survival (fraction)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())
ggsave("HT_trial_dose_response.png", p_dr, width = 8, height = 5, dpi = 300, bg = "white")
cat("Saved: HT_trial_dose_response.png\n\n")


##=============================================================================
## SUMMARY
##=============================================================================

cat("================================================================\n")
cat("  HUMAN TRANSLATION — COMPLETE\n")
cat("================================================================\n\n")

cat("SECTION 1: Allometric PK Scaling\n")
cat(sprintf("  Plasma t\u00bd: mouse %.1fh \u2192 human %.1fh\n",
            log(2) / 0.231, log(2) / HUMAN_PK$ke_plasma))
cat(sprintf("  C_ref: mouse %.2f \u2192 human %.6f\n", 4.43, C_ref_human))
cat(sprintf("  Human reference dose: 30 mg SC daily (HED of mouse 5 mg/kg)\n\n"))

cat("SECTION 2: Human Disease Timescale\n")
cat(sprintf("  CI_max_human = %.2f (mouse: %.2f)\n", CI_max_human, CAL$CI_max))
cat(sprintf("  alpha_clear_human = %.2f (mouse: %.2f)\n", alpha_clear_human, CAL$alpha_clear))
if (!is.na(bern_idx)) {
  cat(sprintf("  Bernheimer threshold: %.1f years (disease), ", bern_years))
  if (!is.na(bern_idx_rx)) {
    cat(sprintf("%.1f years (drug)\n", bern_years_rx))
  } else {
    cat("not reached (drug)\n")
  }
}
cat("\n")

cat(sprintf("SECTION 3: Virtual Population (N=%d)\n", n_completed))
cat("  SA-informed log-normal variability | 8 PD + 2 disease parameters\n\n")

cat("SECTION 4: Clinical Trial Simulation\n")
cat("  18-month, placebo-controlled, 3 dose levels (10, 30, 100 mg SC)\n")
cat(sprintf("  Best effect (100 mg, Motor): d=%.2f, p=%.1e\n",
            effect_table$cohens_d[effect_table$arm == "SBT-272 100mg" &
                                  effect_table$endpoint == "Motor_score"],
            effect_table$p_value[effect_table$arm == "SBT-272 100mg" &
                                 effect_table$endpoint == "Motor_score"]))
cat("\n")

cat("Plots generated:\n")
cat("  HT_PK_comparison.png       \u2014 Mouse vs human PK profiles\n")
cat("  HT_human_progression.png   \u2014 20-year DA/Motor trajectories\n")
cat("  HT_vpop_distribution.png   \u2014 Virtual population endpoint distributions\n")
cat("  HT_trial_forest.png        \u2014 Forest plot of treatment effects\n")
cat("  HT_trial_waterfall.png     \u2014 Individual patient DA response\n")
cat("  HT_trial_dose_response.png \u2014 Dose-response curve\n")
cat("================================================================\n")
