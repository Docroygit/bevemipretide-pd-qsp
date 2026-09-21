##=============================================================================
## PHASED SIMULATION
## Clinically Realistic Disease Progression → Diagnosis → Treatment
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
##
## Rationale:
##   PD is a chronic progressive disease. Motor symptoms appear after
##   substantial DA neuron loss (~30-50% SNpc loss; Fearnley & Lees 1991,
##   Bezard et al. 2001). Diagnosis follows with a median delay of ~12
##   months (Breen et al. 2013). Treatment therefore begins at a disease
##   state far from healthy baseline.
##
##   This script implements phased simulation:
##     Phase 1: Prodromal disease (no drug) — run until symptom threshold
##     Phase 2: Diagnostic delay (no drug) — continue for delay period
##     Phase 3: Treatment — drug from the disease state at diagnosis
##
## Key literature:
##   Fearnley & Lees 1991 (Brain): 48% SNpc loss at symptom onset
##   Bezard et al. 2001 (J Neurosci): 43% TH+ neuron loss at symptom onset
##   Breen et al. 2013 (J Neurol): median 12-month diagnostic delay
##   Bernheimer et al. 1973: 60% striatal DA loss for established PD
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

cat("================================================================\n")
cat("  PHASED SIMULATION\n")
cat("  Prodromal → Diagnosis → Treatment\n")
cat("  Clinically realistic disease-treatment timeline\n")
cat("================================================================\n\n")


##=============================================================================
## SETUP: Human parameters (from Human Translation / Robustness Analysis)
##=============================================================================

CAL <- list(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
            K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)

BW_human <- 70.0
HUMAN_PK <- list(
  ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
  k_brain_in = 0.015, k_brain_out = 0.010,
  k_mito_in = 0.50, k_mito_out = 0.020
)

CI_max_human      <- 0.90
alpha_clear_human <- 0.35
death_scale_human <- 0.06

## Calibrate C_ref_human
po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE)
last_day <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_day$C_mito)
cat(sprintf("C_ref_human = %.6f\n\n", C_ref_human))


##=============================================================================
## CORE FUNCTIONS
##=============================================================================

## Helper: run human model with optional custom initial conditions
run_human_phased <- function(dose_mg, duration_days, CI_max = 1.0,
                             alpha_clear = 1.0, dt = 24,
                             parms_extra = NULL, init_state = NULL,
                             quiet = TRUE) {
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
                 init_override = init_state)
}


## Extract state vector from a simulation result row
extract_state <- function(sim_result, row_idx = nrow(sim_result)) {
  state_names <- c("Depot", "C_plasma", "C_periph", "C_brain", "C_mito",
                   "CL_n", "CL_ox", "CL_ext", "ALCAT1", "TAZ",
                   "aSyn_mono", "aSyn_olig", "aSyn_ext",
                   "CI_activity", "SC_integrity", "delta_psi_m", "ATP", "mPTP_open",
                   "mROS", "SOD2_act", "GPx4_act",
                   "PINK1_act", "mito_damage", "DA_neuron",
                   "MG_active", "TNF", "IL1b",
                   "Motor_score")
  state <- as.numeric(sim_result[row_idx, state_names])
  names(state) <- state_names
  ## PK states must be zero when starting a new phase without drug
  state[c("Depot", "C_plasma", "C_periph", "C_brain", "C_mito")] <- 0
  return(state)
}


## Find the time (in years) when DA_neuron crosses a threshold
find_threshold_crossing <- function(sim_result, threshold, variable = "DA_neuron") {
  idx <- which(sim_result[[variable]] < threshold)
  if (length(idx) == 0) return(NA)
  return(sim_result$time[idx[1]] / 24 / 365)
}


##=============================================================================
## SECTION 1: PRODROMAL PHASE — How long until symptoms appear?
##=============================================================================

cat("================================================================\n")
cat("  SECTION 1: PRODROMAL PHASE CHARACTERISATION\n")
cat("  Disease-only trajectory to determine symptom onset\n")
cat("================================================================\n\n")

## Symptom onset thresholds (from literature)
THRESHOLD_EARLY_SYMPTOMS <- 0.70   # 30% SNpc loss — first motor signs
THRESHOLD_CLINICAL_PD    <- 0.50   # 50% SNpc loss — clinical PD diagnosis
THRESHOLD_BERNHEIMER     <- 0.40   # 60% loss — established PD

## Run 20-year disease-only trajectory
cat("Running 20-year disease-only trajectory...\n")
res_disease_20y <- run_human_phased(dose_mg = 0, duration_days = 7300,
                                    CI_max = CI_max_human,
                                    alpha_clear = alpha_clear_human, dt = 24)

## Find threshold crossings
t_early    <- find_threshold_crossing(res_disease_20y, THRESHOLD_EARLY_SYMPTOMS)
t_clinical <- find_threshold_crossing(res_disease_20y, THRESHOLD_CLINICAL_PD)
t_bernheim <- find_threshold_crossing(res_disease_20y, THRESHOLD_BERNHEIMER)

cat("\n--- PRODROMAL PHASE MILESTONES ---\n")
cat(sprintf("  First motor signs (DA < %.0f%%): %.1f years\n",
            THRESHOLD_EARLY_SYMPTOMS * 100, t_early))
cat(sprintf("  Clinical PD diagnosis grade (DA < %.0f%%): %.1f years\n",
            THRESHOLD_CLINICAL_PD * 100, t_clinical))
cat(sprintf("  Bernheimer threshold (DA < %.0f%%): %.1f years\n",
            THRESHOLD_BERNHEIMER * 100, t_bernheim))

## DA neuron level at each milestone
idx_early <- which.min(abs(res_disease_20y$time - t_early * 365 * 24))
idx_clin  <- which.min(abs(res_disease_20y$time - t_clinical * 365 * 24))

cat(sprintf("\n  At first motor signs (%.1f yr):\n", t_early))
cat(sprintf("    DA_neuron = %.4f, Motor_score = %.4f\n",
            res_disease_20y$DA_neuron[idx_early],
            res_disease_20y$Motor_score[idx_early]))
cat(sprintf("    CI = %.3f, aSyn_olig = %.3f, mROS = %.3f\n",
            res_disease_20y$CI_activity[idx_early],
            res_disease_20y$aSyn_olig[idx_early],
            res_disease_20y$mROS[idx_early]))

cat(sprintf("\n  At clinical PD grade (%.1f yr):\n", t_clinical))
cat(sprintf("    DA_neuron = %.4f, Motor_score = %.4f\n",
            res_disease_20y$DA_neuron[idx_clin],
            res_disease_20y$Motor_score[idx_clin]))
cat(sprintf("    CI = %.3f, aSyn_olig = %.3f, mROS = %.3f\n",
            res_disease_20y$CI_activity[idx_clin],
            res_disease_20y$aSyn_olig[idx_clin],
            res_disease_20y$mROS[idx_clin]))


##=============================================================================
## SECTION 2: DIAGNOSTIC DELAY — How much more is lost?
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 2: DIAGNOSTIC DELAY IMPACT\n")
cat("  What happens during the 12-month lag from symptoms to diagnosis?\n")
cat("================================================================\n\n")

DIAG_DELAY_MONTHS <- 12   # Breen et al. 2013: median 12 months
DIAG_DELAY_DAYS   <- round(DIAG_DELAY_MONTHS * 30.44)

## State at symptom onset
state_at_symptoms <- extract_state(res_disease_20y, idx_early)
cat(sprintf("State at symptom onset (t = %.1f years):\n", t_early))
cat(sprintf("  DA_neuron = %.4f (%.1f%% loss)\n",
            state_at_symptoms["DA_neuron"],
            (1 - state_at_symptoms["DA_neuron"]) * 100))

## Run disease for 12 more months from symptom onset state
res_diag_delay <- run_human_phased(dose_mg = 0, duration_days = DIAG_DELAY_DAYS,
                                   CI_max = CI_max_human,
                                   alpha_clear = alpha_clear_human,
                                   dt = 24, init_state = state_at_symptoms)

state_at_diagnosis <- extract_state(res_diag_delay)
t_diagnosis <- t_early + DIAG_DELAY_MONTHS / 12

cat(sprintf("\nState at diagnosis (t = %.1f years, after %d-month delay):\n",
            t_diagnosis, DIAG_DELAY_MONTHS))
cat(sprintf("  DA_neuron = %.4f (%.1f%% loss)\n",
            state_at_diagnosis["DA_neuron"],
            (1 - state_at_diagnosis["DA_neuron"]) * 100))
cat(sprintf("  Additional DA loss during diagnostic delay: %.1f%%\n",
            (state_at_symptoms["DA_neuron"] - state_at_diagnosis["DA_neuron"]) * 100))
cat(sprintf("  Motor_score = %.4f (vs %.4f at symptom onset)\n",
            state_at_diagnosis["Motor_score"], state_at_symptoms["Motor_score"]))

## Disease state comparison table
cat("\n--- DISEASE STATE: SYMPTOM ONSET vs DIAGNOSIS ---\n")
compare_vars <- c("DA_neuron", "CI_activity", "aSyn_olig", "mROS",
                   "CL_n", "mPTP_open", "MG_active", "Motor_score")
cat(sprintf("  %-16s  %10s  %10s  %10s\n",
            "Variable", "Symptom", "Diagnosis", "Change"))
for (v in compare_vars) {
  cat(sprintf("  %-16s  %10.4f  %10.4f  %+10.4f\n",
              v, state_at_symptoms[v], state_at_diagnosis[v],
              state_at_diagnosis[v] - state_at_symptoms[v]))
}


##=============================================================================
## SECTION 3: PHASED TREATMENT — Drug starts at diagnosis
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 3: PHASED TREATMENT SIMULATION\n")
cat("  Drug starts at diagnosis | 3-year treatment duration\n")
cat("================================================================\n\n")

TREATMENT_DAYS <- 1095   # 3 years of treatment
DOSES_MG <- c(10, 30, 60)

## Run disease-only continuation for 3 years from diagnosis (no treatment control)
res_notrt <- run_human_phased(dose_mg = 0, duration_days = TREATMENT_DAYS,
                              CI_max = CI_max_human,
                              alpha_clear = alpha_clear_human,
                              dt = 24, init_state = state_at_diagnosis)

## Run treatment at each dose from diagnosis state
treatment_results <- list()
treatment_results[["No treatment"]] <- res_notrt

for (d in DOSES_MG) {
  cat(sprintf("  Running %d mg treatment for %d days from diagnosis state...\n",
              d, TREATMENT_DAYS))
  res_trt <- run_human_phased(dose_mg = d, duration_days = TREATMENT_DAYS,
                              CI_max = CI_max_human,
                              alpha_clear = alpha_clear_human,
                              dt = 24, init_state = state_at_diagnosis)
  treatment_results[[sprintf("SBT-272 %d mg", d)]] <- res_trt
}

## Summary table
cat("\n--- TREATMENT OUTCOMES (3 years from diagnosis) ---\n")
cat(sprintf("  Starting DA_neuron at diagnosis: %.4f (%.1f%% loss)\n\n",
            state_at_diagnosis["DA_neuron"],
            (1 - state_at_diagnosis["DA_neuron"]) * 100))

cat(sprintf("  %-20s  %10s  %10s  %10s  %10s  %10s\n",
            "Scenario", "DA_final", "DA_saved", "Motor_f", "Motor_Δ", "CI_final"))
for (nm in names(treatment_results)) {
  res <- treatment_results[[nm]]
  ss <- res[nrow(res), ]
  da_saved <- ss$DA_neuron - treatment_results[["No treatment"]][nrow(treatment_results[["No treatment"]]), "DA_neuron"]
  motor_delta <- ss$Motor_score - treatment_results[["No treatment"]][nrow(treatment_results[["No treatment"]]), "Motor_score"]
  cat(sprintf("  %-20s  %10.4f  %+10.4f  %10.4f  %+10.4f  %10.4f\n",
              nm, ss$DA_neuron, da_saved, ss$Motor_score, motor_delta, ss$CI_activity))
}


##=============================================================================
## SECTION 4: SCENARIO ANALYSIS — Treatment delay after diagnosis
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 4: SCENARIO ANALYSIS\n")
cat("  Treatment delay: at diagnosis, +6m, +1y, +2y\n")
cat("  Treatment duration: 18 months from start\n")
cat("================================================================\n\n")

TREATMENT_DURATION <- 547   # 18-month trial (consistent with RA)
DELAY_MONTHS <- c(0, 6, 12, 24)
REFERENCE_DOSE <- 30   # mg

## For each delay scenario:
## 1. Continue disease from diagnosis for delay period
## 2. Start drug from that state
## 3. Run for 18 months
## Also run disease-only for the full period as control

scenario_results <- data.frame()

for (delay_m in DELAY_MONTHS) {
  delay_days <- round(delay_m * 30.44)

  ## State at treatment start
  if (delay_m == 0) {
    state_trt_start <- state_at_diagnosis
    da_at_start <- state_at_diagnosis["DA_neuron"]
  } else {
    res_delay <- run_human_phased(dose_mg = 0, duration_days = delay_days,
                                  CI_max = CI_max_human,
                                  alpha_clear = alpha_clear_human,
                                  dt = 24, init_state = state_at_diagnosis)
    state_trt_start <- extract_state(res_delay)
    da_at_start <- state_trt_start["DA_neuron"]
  }

  ## Run treatment (drug) for 18 months
  res_drug <- run_human_phased(dose_mg = REFERENCE_DOSE,
                               duration_days = TREATMENT_DURATION,
                               CI_max = CI_max_human,
                               alpha_clear = alpha_clear_human,
                               dt = 24, init_state = state_trt_start)

  ## Run disease-only for same 18 months (matched control)
  res_ctrl <- run_human_phased(dose_mg = 0,
                               duration_days = TREATMENT_DURATION,
                               CI_max = CI_max_human,
                               alpha_clear = alpha_clear_human,
                               dt = 24, init_state = state_trt_start)

  da_drug_end  <- res_drug[nrow(res_drug), "DA_neuron"]
  da_ctrl_end  <- res_ctrl[nrow(res_ctrl), "DA_neuron"]
  da_saved     <- da_drug_end - da_ctrl_end
  motor_drug   <- res_drug[nrow(res_drug), "Motor_score"]
  motor_ctrl   <- res_ctrl[nrow(res_ctrl), "Motor_score"]
  motor_saved  <- motor_ctrl - motor_drug

  scenario_results <- rbind(scenario_results, data.frame(
    delay_months     = delay_m,
    DA_at_start      = as.numeric(da_at_start),
    DA_loss_at_start = (1 - as.numeric(da_at_start)) * 100,
    DA_drug_end      = da_drug_end,
    DA_ctrl_end      = da_ctrl_end,
    DA_preserved     = da_saved,
    Motor_drug_end   = motor_drug,
    Motor_ctrl_end   = motor_ctrl,
    Motor_improvement = motor_saved,
    stringsAsFactors = FALSE
  ))
}

## Calculate efficacy retention (relative to no-delay scenario)
max_benefit <- scenario_results$DA_preserved[1]
scenario_results$efficacy_retained <- scenario_results$DA_preserved / max_benefit * 100

cat("--- SCENARIO ANALYSIS RESULTS (30 mg, 18-month treatment) ---\n\n")
cat(sprintf("  %-8s  %10s  %10s  %10s  %10s  %10s  %10s\n",
            "Delay", "DA@start", "DA_drug", "DA_ctrl", "DA_saved",
            "Motor_Δ", "Efficacy%"))
for (i in 1:nrow(scenario_results)) {
  r <- scenario_results[i, ]
  cat(sprintf("  %-8s  %10.4f  %10.4f  %10.4f  %+10.4f  %+10.4f  %9.1f%%\n",
              sprintf("%dm", r$delay_months),
              r$DA_at_start, r$DA_drug_end, r$DA_ctrl_end,
              r$DA_preserved, r$Motor_improvement, r$efficacy_retained))
}

cat(sprintf("\n  Key finding: Every %d-month delay costs ~%.1f%% of maximum drug benefit\n",
            6, 100 - scenario_results$efficacy_retained[2]))


##=============================================================================
## SECTION 5: MULTI-DOSE SCENARIO ANALYSIS
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 5: MULTI-DOSE × MULTI-DELAY ANALYSIS\n")
cat("  Doses: 10, 30, 60 mg × Delays: 0, 6, 12, 24 months\n")
cat("================================================================\n\n")

full_scenario <- data.frame()

for (dose in DOSES_MG) {
  for (delay_m in DELAY_MONTHS) {
    delay_days <- round(delay_m * 30.44)

    if (delay_m == 0) {
      state_trt_start <- state_at_diagnosis
    } else {
      res_delay <- run_human_phased(dose_mg = 0, duration_days = delay_days,
                                    CI_max = CI_max_human,
                                    alpha_clear = alpha_clear_human,
                                    dt = 24, init_state = state_at_diagnosis)
      state_trt_start <- extract_state(res_delay)
    }

    res_drug <- run_human_phased(dose_mg = dose,
                                 duration_days = TREATMENT_DURATION,
                                 CI_max = CI_max_human,
                                 alpha_clear = alpha_clear_human,
                                 dt = 24, init_state = state_trt_start)
    res_ctrl <- run_human_phased(dose_mg = 0,
                                 duration_days = TREATMENT_DURATION,
                                 CI_max = CI_max_human,
                                 alpha_clear = alpha_clear_human,
                                 dt = 24, init_state = state_trt_start)

    da_saved <- res_drug[nrow(res_drug), "DA_neuron"] - res_ctrl[nrow(res_ctrl), "DA_neuron"]
    motor_saved <- res_ctrl[nrow(res_ctrl), "Motor_score"] - res_drug[nrow(res_drug), "Motor_score"]

    full_scenario <- rbind(full_scenario, data.frame(
      dose_mg = dose,
      delay_months = delay_m,
      DA_at_start = as.numeric(state_trt_start["DA_neuron"]),
      DA_preserved = da_saved,
      Motor_improvement = motor_saved,
      stringsAsFactors = FALSE
    ))
  }
}

## Efficacy retention per dose (relative to 0-delay for that dose)
for (d in DOSES_MG) {
  ref_val <- full_scenario$DA_preserved[full_scenario$dose_mg == d &
                                         full_scenario$delay_months == 0]
  idx <- full_scenario$dose_mg == d
  full_scenario$efficacy_retained[idx] <- full_scenario$DA_preserved[idx] / ref_val * 100
}

cat("--- FULL DOSE × DELAY MATRIX ---\n\n")
cat(sprintf("  %-8s  %-8s  %10s  %10s  %10s  %10s\n",
            "Dose", "Delay", "DA@start", "DA_saved", "Motor_Δ", "Efficacy%"))
for (i in 1:nrow(full_scenario)) {
  r <- full_scenario[i, ]
  cat(sprintf("  %-8s  %-8s  %10.4f  %+10.4f  %+10.4f  %9.1f%%\n",
              sprintf("%dmg", r$dose_mg), sprintf("%dm", r$delay_months),
              r$DA_at_start, r$DA_preserved, r$Motor_improvement,
              r$efficacy_retained))
}


##=============================================================================
## SECTION 6: PLOTS
##=============================================================================

cat("\n\n================================================================\n")
cat("  SECTION 6: PUBLICATION PLOTS\n")
cat("================================================================\n\n")

## --- Plot 1: Full disease trajectory with phase annotations ---

res_disease_20y$time_years <- res_disease_20y$time / 24 / 365
phase_df <- data.frame(
  xmin = c(0, t_early, t_diagnosis),
  xmax = c(t_early, t_diagnosis, 20),
  phase = factor(c("Prodromal\n(asymptomatic)",
                   "Diagnostic\ndelay",
                   "Post-diagnosis"),
                 levels = c("Prodromal\n(asymptomatic)",
                            "Diagnostic\ndelay",
                            "Post-diagnosis"))
)

prog_vars <- c("DA_neuron", "Motor_score", "CI_activity", "aSyn_olig")
long_dis <- res_disease_20y %>%
  select(time_years, all_of(prog_vars)) %>%
  pivot_longer(cols = all_of(prog_vars), names_to = "variable", values_to = "value")
long_dis$variable <- factor(long_dis$variable, levels = prog_vars,
                            labels = c("DA Neuron Survival", "Motor Score",
                                       "Complex I Activity",
                                       "α-Synuclein Oligomers"))

p1 <- ggplot(long_dis, aes(x = time_years, y = value)) +
  geom_rect(data = phase_df, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
            alpha = 0.15) +
  geom_line(linewidth = 0.9, color = "#d62728") +
  geom_vline(xintercept = t_early, linetype = "dashed", color = "orange", linewidth = 0.6) +
  geom_vline(xintercept = t_diagnosis, linetype = "dashed", color = "red", linewidth = 0.6) +
  geom_hline(data = data.frame(
    variable = factor("DA Neuron Survival", levels = levels(long_dis$variable)),
    yint = c(THRESHOLD_EARLY_SYMPTOMS, THRESHOLD_CLINICAL_PD, THRESHOLD_BERNHEIMER)),
    aes(yintercept = yint), linetype = "dotted", color = "gray50") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Prodromal\n(asymptomatic)" = "#2ca02c",
                                "Diagnostic\ndelay" = "#ff7f0e",
                                "Post-diagnosis" = "#d62728")) +
  labs(title = "Human PD Progression: Clinical Phase Annotations",
       subtitle = sprintf("CI_max=%.2f, alpha_clear=%.2f | Symptom onset: %.1f yr | Diagnosis: %.1f yr",
                          CI_max_human, alpha_clear_human, t_early, t_diagnosis),
       x = "Time (years)", y = "Normalised level", fill = "Disease Phase",
       caption = "Dashed lines: symptom onset (orange), diagnosis (red). Dotted: DA thresholds (70%, 50%, 40%)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave("PS_disease_phases.png", p1, width = 12, height = 8, dpi = 300, bg = "white")
cat("Saved: PS_disease_phases.png\n")


## --- Plot 2: Treatment from diagnosis — dose comparison ---

trt_long <- data.frame()
for (nm in names(treatment_results)) {
  df <- treatment_results[[nm]]
  df$scenario <- nm
  df$time_months <- df$time / 24 / 30.44
  trt_long <- rbind(trt_long, df)
}

trt_vars <- c("DA_neuron", "Motor_score")
trt_plot <- trt_long %>%
  select(time_months, scenario, all_of(trt_vars)) %>%
  pivot_longer(cols = all_of(trt_vars), names_to = "variable", values_to = "value")
trt_plot$variable <- factor(trt_plot$variable, levels = trt_vars,
                            labels = c("DA Neuron Survival", "Motor Score"))
trt_plot$scenario <- factor(trt_plot$scenario,
                            levels = c("No treatment", "SBT-272 10 mg",
                                       "SBT-272 30 mg", "SBT-272 60 mg"))

trt_colors <- c("No treatment" = "#d62728", "SBT-272 10 mg" = "#9467bd",
                "SBT-272 30 mg" = "#1f77b4", "SBT-272 60 mg" = "#2ca02c")

p2 <- ggplot(trt_plot, aes(x = time_months, y = value, color = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = trt_colors) +
  labs(title = "Treatment from Diagnosis: Dose Comparison (3-Year)",
       subtitle = sprintf("Treatment starts at diagnosis (%.1f yr from disease onset, DA=%.1f%%)",
                          t_diagnosis, state_at_diagnosis["DA_neuron"] * 100),
       x = "Time from treatment start (months)", y = "Normalised level",
       color = "Scenario") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave("PS_treatment_from_diagnosis.png", p2, width = 12, height = 5, dpi = 300, bg = "white")
cat("Saved: PS_treatment_from_diagnosis.png\n")


## --- Plot 3: Scenario analysis — Efficacy vs delay ---

p3_data <- scenario_results
p3_data$delay_label <- sprintf("%d months", p3_data$delay_months)
p3_data$delay_label <- factor(p3_data$delay_label,
                               levels = sprintf("%d months", DELAY_MONTHS))

p3 <- ggplot(p3_data, aes(x = delay_months, y = efficacy_retained)) +
  geom_line(linewidth = 1.0, color = "#1f77b4") +
  geom_point(size = 3, color = "#1f77b4") +
  geom_text(aes(label = sprintf("%.1f%%", efficacy_retained)),
            vjust = -1.2, size = 3.5, fontface = "bold") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(DELAY_MONTHS), y = 52,
           label = "50% efficacy threshold", hjust = 1, color = "gray50", size = 3) +
  scale_x_continuous(breaks = DELAY_MONTHS) +
  scale_y_continuous(limits = c(0, 110)) +
  labs(title = "Treatment Delay vs Efficacy Retention (SBT-272 30 mg)",
       subtitle = sprintf("18-month treatment | Diagnosis at %.1f yr (DA=%.1f%%)",
                          t_diagnosis, state_at_diagnosis["DA_neuron"] * 100),
       x = "Treatment delay after diagnosis (months)",
       y = "Efficacy retained (%)",
       caption = "Efficacy = DA neurons preserved (drug vs no-drug), relative to treatment at diagnosis") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        panel.grid.minor = element_blank())
ggsave("PS_efficacy_vs_delay.png", p3, width = 8, height = 6, dpi = 300, bg = "white")
cat("Saved: PS_efficacy_vs_delay.png\n")


## --- Plot 4: Multi-dose × multi-delay heatmap ---

full_scenario$dose_label <- sprintf("%d mg", full_scenario$dose_mg)
full_scenario$delay_label <- sprintf("%d months", full_scenario$delay_months)
full_scenario$dose_label <- factor(full_scenario$dose_label,
                                    levels = sprintf("%d mg", DOSES_MG))
full_scenario$delay_label <- factor(full_scenario$delay_label,
                                     levels = sprintf("%d months", DELAY_MONTHS))

p4 <- ggplot(full_scenario, aes(x = delay_label, y = dose_label,
                                 fill = DA_preserved)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("+%.3f\n(%.0f%%)", DA_preserved, efficacy_retained)),
            size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#d62728", mid = "#ffffbf", high = "#1a9850",
                       midpoint = median(full_scenario$DA_preserved),
                       name = "DA neurons\npreserved") +
  labs(title = "DA Neuron Preservation: Dose × Treatment Delay",
       subtitle = "18-month treatment from each delay point | Values: absolute DA preserved (% efficacy retained)",
       x = "Treatment delay after diagnosis", y = "Dose") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 11))
ggsave("PS_dose_delay_heatmap.png", p4, width = 10, height = 5, dpi = 300, bg = "white")
cat("Saved: PS_dose_delay_heatmap.png\n")


## --- Plot 5: Composite trajectory — all phases stitched together ---

## Stitch together: prodromal (to symptoms) → delay → treatment (30 mg vs none)
## Time axis: absolute years from disease onset

# Phase 1: prodromal (disease only, 0 → symptom onset)
p1_data <- res_disease_20y[res_disease_20y$time_years <= t_early, ]
p1_data$phase <- "Prodromal"
p1_data$abs_time <- p1_data$time_years

# Phase 2: diagnostic delay (symptom onset → diagnosis)
p2_data <- res_diag_delay
p2_data$time_years_rel <- p2_data$time / 24 / 365
p2_data$abs_time <- t_early + p2_data$time_years_rel
p2_data$phase <- "Diagnostic delay"

# Phase 3a: no treatment (18 months from diagnosis)
p3_notrt <- treatment_results[["No treatment"]]
p3_notrt$time_years_rel <- p3_notrt$time / 24 / 365
p3_notrt$abs_time <- t_diagnosis + p3_notrt$time_years_rel
p3_notrt$phase <- "No treatment"

# Phase 3b: 30 mg treatment (18 months from diagnosis)
p3_drug <- treatment_results[["SBT-272 30 mg"]]
p3_drug$time_years_rel <- p3_drug$time / 24 / 365
p3_drug$abs_time <- t_diagnosis + p3_drug$time_years_rel
p3_drug$phase <- "SBT-272 30 mg"

# Combine for DA_neuron trajectory
composite <- rbind(
  data.frame(abs_time = p1_data$abs_time, DA_neuron = p1_data$DA_neuron,
             Motor_score = p1_data$Motor_score, trajectory = "Disease course"),
  data.frame(abs_time = p2_data$abs_time, DA_neuron = p2_data$DA_neuron,
             Motor_score = p2_data$Motor_score, trajectory = "Disease course"),
  data.frame(abs_time = p3_notrt$abs_time, DA_neuron = p3_notrt$DA_neuron,
             Motor_score = p3_notrt$Motor_score, trajectory = "No treatment"),
  data.frame(abs_time = p3_drug$abs_time, DA_neuron = p3_drug$DA_neuron,
             Motor_score = p3_drug$Motor_score, trajectory = "SBT-272 30 mg")
)

composite_long <- composite %>%
  pivot_longer(cols = c(DA_neuron, Motor_score),
               names_to = "variable", values_to = "value")
composite_long$variable <- factor(composite_long$variable,
                                   levels = c("DA_neuron", "Motor_score"),
                                   labels = c("DA Neuron Survival", "Motor Score"))
composite_long$trajectory <- factor(composite_long$trajectory,
                                     levels = c("Disease course", "No treatment",
                                                "SBT-272 30 mg"))

comp_colors <- c("Disease course" = "#ff7f0e",
                 "No treatment" = "#d62728",
                 "SBT-272 30 mg" = "#1f77b4")

p5 <- ggplot(composite_long, aes(x = abs_time, y = value, color = trajectory)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = t_early, linetype = "dashed", color = "orange",
             linewidth = 0.5) +
  geom_vline(xintercept = t_diagnosis, linetype = "dashed", color = "red",
             linewidth = 0.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = comp_colors) +
  annotate("text", x = t_early, y = Inf, label = "Symptom\nonset",
           vjust = 1.5, hjust = 1.1, size = 3, color = "orange") +
  annotate("text", x = t_diagnosis, y = Inf, label = "Diagnosis",
           vjust = 1.5, hjust = -0.1, size = 3, color = "red") +
  labs(title = "Complete PD Timeline: Prodromal → Diagnosis → Treatment",
       subtitle = sprintf("Prodromal: %.1f yr | Diagnostic delay: %d mo | Treatment: 3 yr",
                          t_early, DIAG_DELAY_MONTHS),
       x = "Time from disease onset (years)", y = "Normalised level",
       color = "Trajectory") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
ggsave("PS_composite_timeline.png", p5, width = 14, height = 6, dpi = 300, bg = "white")
cat("Saved: PS_composite_timeline.png\n")


##=============================================================================
## SUMMARY
##=============================================================================

cat("\n\n================================================================\n")
cat("  SUMMARY\n")
cat("================================================================\n\n")

cat("DISEASE TIMELINE:\n")
cat(sprintf("  Prodromal duration (to 30%% DA loss): %.1f years\n", t_early))
cat(sprintf("  Diagnostic delay: %d months\n", DIAG_DELAY_MONTHS))
cat(sprintf("  DA loss at diagnosis: %.1f%%\n",
            (1 - state_at_diagnosis["DA_neuron"]) * 100))
cat(sprintf("  Time to Bernheimer threshold (untreated): %.1f years\n", t_bernheim))

cat("\nTREATMENT EFFICACY (30 mg, 18 months from diagnosis):\n")
r0 <- scenario_results[1, ]
cat(sprintf("  DA preserved vs no treatment: +%.4f (%.1f%% of remaining)\n",
            r0$DA_preserved, r0$DA_preserved / r0$DA_ctrl_end * 100))
cat(sprintf("  Motor score improvement: %.4f\n", r0$Motor_improvement))

cat("\nEFFICACY DECAY WITH TREATMENT DELAY:\n")
for (i in 1:nrow(scenario_results)) {
  r <- scenario_results[i, ]
  cat(sprintf("  %2d-month delay: %.1f%% efficacy retained\n",
              r$delay_months, r$efficacy_retained))
}

cat("\nPLOTS GENERATED:\n")
cat("  PS_disease_phases.png       - Disease trajectory with phase annotations\n")
cat("  PS_treatment_from_diagnosis.png - Dose comparison from diagnosis state\n")
cat("  PS_efficacy_vs_delay.png    - Efficacy retention vs treatment delay\n")
cat("  PS_dose_delay_heatmap.png   - Dose × delay matrix heatmap\n")
cat("  PS_composite_timeline.png   - Complete stitched timeline\n")

cat("\n================================================================\n")
cat("  PHASED SIMULATION COMPLETE\n")
cat("================================================================\n")
