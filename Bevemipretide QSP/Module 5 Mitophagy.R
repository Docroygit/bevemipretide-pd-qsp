##=============================================================================
## MODULE M5: MITOPHAGY, QUALITY CONTROL, AND DA NEURON SURVIVAL
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
## Version : 1.0
##
## State variables (3):
##   PINK1_act   — PINK1/Parkin pathway activation [0,1]
##   mito_damage — Net mitochondrial damage accumulation [0,1]
##   DA_neuron   — Dopaminergic neuron survival fraction [0,1]
##
## Upstream coupling:
##   delta_psi_m ← M3 (depolarisation activates PINK1)
##   mPTP_open   ← M3 (permeability transition drives damage)
##   ATP         ← M3 (powers biogenesis/repair)
##   CL_ext      ← M1 (externalised CL as "eat me" mitophagy signal)
##
## Key output:
##   DA_neuron → M7 (motor score via Hill function, Bernheimer threshold)
##
## Biology:
##   When delta_psi_m drops, PINK1 stabilises on the OMM, recruits
##   Parkin, and tags damaged mitochondria for autophagy (mitophagy).
##   CL_ext on the OMM is an additional "eat me" signal (Kagan 2015).
##   Mitophagy is PROTECTIVE — it removes damaged organelles.
##
##   Damage accumulates from mPTP opening (cytochrome c release),
##   energy failure (low ATP), and depolarisation. Repair is driven
##   by PINK1-mediated mitophagy and ATP-dependent biogenesis.
##
##   DA neuron death follows a Hill threshold on mito_damage. Below
##   K_death, death is negligible (healthy). Above K_death, death
##   accelerates dramatically. This captures the Bernheimer observation
##   that motor symptoms appear after 50-70% DA neuron loss — a
##   threshold-crossing event.
##
##   DA neurons CANNOT proliferate in adults, so DA_neuron can only
##   decline. Its value at a given timepoint (not steady state) is the
##   meaningful output.
##
## Standalone calibration note:
##   M3 standalone disease values are mild (mPTP=0.003, dpsi=0.935)
##   because the vicious cycle is not running. Standalone M5 therefore
##   shows modest DA_neuron decline (~5% over 5 weeks at disease).
##   In the integrated model, vicious cycle amplification pushes
##   mPTP to ~0.10-0.20 and ATP to ~0.65-0.75, producing
##   physiologically meaningful DA loss (50%+ over weeks, matching
##   Bido mouse model observations).
##
## Data sources:
##   Narendra 2010 (doi:10.1083/jcb.200911141): PINK1/Parkin mechanism
##   Kagan 2015: CL externalisation as mitophagy signal
##   Fearnley & Lees 1991: 5%/year DA neuron loss in PD
##   Bernheimer 1973: 50-70% DA loss threshold for motor symptoms
##   Bido et al. poster: SBT-272 rescues TH+ neurons in PD mouse model
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M5: MITOPHAGY & DA NEURON SURVIVAL\n")
cat("  Version 1.0 | Upstream: M3 (ETC), M1 (CL_ext)\n")
cat("  Output: DA_neuron -> M7 (motor endpoint)\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M5_parameters <- function(delta_psi_m = 0.960, mPTP_open = 0.002,
                          ATP = 0.954, CL_ext = 0.002) {

  p <- c(
    ## --- PINK1/Parkin dynamics ---
    k_PINK1_on  = 2.0,       # PINK1 activation by depolarisation (h^-1)
                               # ASSUMED: PINK1 stabilises on OMM when
                               # delta_psi_m drops; Narendra 2010
    k_PINK1_CL  = 20.0,      # CL_ext enhancement of PINK1 (h^-1)
                               # ASSUMED: CL_ext is "eat me" signal that
                               # promotes mitophagy; Kagan 2015
                               # High coefficient because CL_ext is small
    k_PINK1_off = 0.50,      # PINK1 deactivation by repolarisation (h^-1)
                               # ASSUMED: normal delta_psi_m promotes
                               # PINK1 import and degradation

    ## --- Mitochondrial damage dynamics ---
    k_damage_mPTP   = 10.0,  # mPTP-driven damage (h^-1)
                               # ESTIMATED: mPTP opening releases cytochrome c,
                               # disrupts cristae; high coefficient for
                               # sensitivity to small mPTP changes
    k_damage_energy = 1.0,   # Energy failure damage (h^-1)
                               # ASSUMED: low ATP impairs protein quality
                               # control, ion homeostasis
    k_damage_dpsi   = 2.0,   # Depolarisation-driven damage (h^-1)
                               # ASSUMED: loss of membrane potential disrupts
                               # import, ETC function, Ca2+ buffering
    k_repair        = 1.0,   # PINK1-dependent mitophagy repair (h^-1)
                               # ASSUMED: mitophagy removes damaged
                               # mitochondria; slow relative to biogenesis
    k_biogen        = 3.0,   # ATP-dependent biogenesis/repair (h^-1)
                               # ASSUMED: PGC-1a-driven mitochondrial
                               # biogenesis replaces damaged organelles;
                               # dominant repair mechanism

    ## --- DA neuron death ---
    k_death   = 0.001,       # Base death rate (h^-1)
                               # ESTIMATED: tuned for ~5% DA loss at
                               # standalone disease over 5 weeks (840h);
                               # integrated model produces 30-50% loss
    K_death   = 0.20,        # Hill half-max for death threshold
                               # ESTIMATED: threshold for clinically
                               # significant damage; below this, death
                               # is negligible
    n_death   = 3,           # Hill cooperativity for death threshold
                               # ASSUMED: cooperative transition from
                               # compensated to decompensated state

    ## --- Upstream coupling (scenario parameters) ---
    delta_psi_m = delta_psi_m,
    mPTP_open   = mPTP_open,
    ATP         = ATP,
    CL_ext      = CL_ext
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M5_init <- function() {
  c(
    PINK1_act   = 0.10,    # Low basal PINK1 activation
    mito_damage = 0.05,    # Near healthy SS
    DA_neuron   = 1.00     # Full neuronal complement at start
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M5_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states ---
    PINK1_act   <- max(0, min(1, PINK1_act))
    mito_damage <- max(0, min(1, mito_damage))
    DA_neuron   <- max(0, min(1, DA_neuron))

    ## --- ODE 1: PINK1/Parkin pathway activation ---
    ## Depolarisation (low delta_psi_m) stabilises PINK1 on OMM
    ## CL_ext on OMM enhances recognition and mitophagy
    ## Repolarisation promotes PINK1 import and degradation
    dPINK1 <- (k_PINK1_on * (1 - delta_psi_m) * (1 - PINK1_act)
               + k_PINK1_CL * CL_ext * (1 - PINK1_act)
               - k_PINK1_off * delta_psi_m * PINK1_act)

    ## --- ODE 2: Mitochondrial damage accumulation ---
    ## Damage from mPTP, energy failure, depolarisation
    ## Repair from PINK1-driven mitophagy and ATP-dependent biogenesis
    damage_in <- (k_damage_mPTP * mPTP_open
                  + k_damage_energy * (1 - ATP)
                  + k_damage_dpsi * (1 - delta_psi_m))

    repair_out <- k_repair * PINK1_act + k_biogen * ATP

    dmito_damage <- damage_in * (1 - mito_damage) - repair_out * mito_damage

    ## --- ODE 3: DA neuron survival ---
    ## Hill-threshold death: negligible at low damage, catastrophic at high
    ## DA neurons do not proliferate — only decline
    Hill_death <- mito_damage^n_death /
                  (K_death^n_death + mito_damage^n_death)

    death_rate <- k_death * Hill_death

    dDA_neuron <- -death_rate * DA_neuron

    ## --- Derived quantities ---
    mitophagy_flux <- k_repair * PINK1_act * mito_damage
    biogen_flux    <- k_biogen * ATP * mito_damage
    DA_loss_pct    <- (1 - DA_neuron) * 100

    list(
      c(dPINK1, dmito_damage, dDA_neuron),
      damage_in      = damage_in,
      repair_out     = repair_out,
      Hill_death     = Hill_death,
      death_rate     = death_rate,
      mitophagy_flux = mitophagy_flux,
      biogen_flux    = biogen_flux,
      DA_loss_pct    = DA_loss_pct
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M5 <- function(delta_psi_m = 0.960, mPTP_open = 0.002,
                   ATP = 0.954, CL_ext = 0.002,
                   duration_h = 840, dt = 0.5, y0 = NULL) {

  parms <- M5_parameters(delta_psi_m, mPTP_open, ATP, CL_ext)
  if (is.null(y0)) y0 <- M5_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M5_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$delta_psi_m <- delta_psi_m
  result$mPTP_open   <- mPTP_open
  result$ATP         <- ATP
  result$CL_ext      <- CL_ext

  ss <- result[nrow(result), ]
  cat(sprintf(
    "M5 | dpsi=%.3f mPTP=%.4f ATP=%.3f CL_ext=%.3f | PINK1=%.3f D=%.4f DA=%.4f (%.1f%% loss at %dh)\n",
    delta_psi_m, mPTP_open, ATP, CL_ext,
    ss$PINK1_act, ss$mito_damage, ss$DA_neuron,
    (1 - ss$DA_neuron) * 100, duration_h))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE AND DA DECAY
##-----------------------------------------------------------------------------

M5_analytical_ss <- function(delta_psi_m = 0.960, mPTP_open = 0.002,
                             ATP = 0.954, CL_ext = 0.002) {

  p <- M5_parameters(delta_psi_m, mPTP_open, ATP, CL_ext)

  ## PINK1: closed-form SS
  act <- as.numeric(p["k_PINK1_on"]) * (1 - delta_psi_m) +
         as.numeric(p["k_PINK1_CL"]) * CL_ext
  deact <- as.numeric(p["k_PINK1_off"]) * delta_psi_m

  PINK1_ss <- act / (act + deact)

  ## mito_damage: closed-form SS
  a <- (as.numeric(p["k_damage_mPTP"]) * mPTP_open +
        as.numeric(p["k_damage_energy"]) * (1 - ATP) +
        as.numeric(p["k_damage_dpsi"]) * (1 - delta_psi_m))

  b <- as.numeric(p["k_repair"]) * PINK1_ss +
       as.numeric(p["k_biogen"]) * ATP

  D_ss <- a / (a + b)

  ## Hill death rate at SS
  Hill_ss <- D_ss^as.numeric(p["n_death"]) /
             (as.numeric(p["K_death"])^as.numeric(p["n_death"]) +
              D_ss^as.numeric(p["n_death"]))

  death_rate_ss <- as.numeric(p["k_death"]) * Hill_ss

  return(c(PINK1_act   = PINK1_ss,
           mito_damage = D_ss,
           Hill_death  = Hill_ss,
           death_rate  = death_rate_ss))
}


## Analytical DA neuron survival at time t
## Given constant upstream inputs, PINK1 and mito_damage reach SS quickly;
## DA_neuron then follows exponential decay on the slow manifold
M5_DA_analytical <- function(t, delta_psi_m = 0.960, mPTP_open = 0.002,
                             ATP = 0.954, CL_ext = 0.002, DA_0 = 1.0) {

  ss <- M5_analytical_ss(delta_psi_m, mPTP_open, ATP, CL_ext)
  DA_0 * exp(-ss["death_rate"] * t)
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M5 <- function() {

  cat("=== MODULE M5 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6
  t_bido      <- 840    # 5 weeks in hours (Bido protocol)


  ## --- V1: Healthy steady state ---
  ## PINK1 low (basal mitophagy), damage low, DA stable
  ## Upstream from M3/M1 healthy: dpsi=0.960, mPTP=0.002, ATP=0.954
  cat("V1: Healthy steady state and DA stability...\n")
  res_h <- run_M5(delta_psi_m = 0.960, mPTP_open = 0.002,
                  ATP = 0.954, CL_ext = 0.002,
                  duration_h = t_bido)
  ss_h  <- res_h[nrow(res_h), ]
  ss_an_h <- M5_analytical_ss(0.960, 0.002, 0.954, 0.002)

  v1 <- (ss_h$PINK1_act > 0.10 && ss_h$PINK1_act < 0.30 &&
         ss_h$mito_damage > 0.02 && ss_h$mito_damage < 0.10 &&
         ss_h$DA_neuron > 0.95)
  cat(sprintf("    PINK1_act   = %.4f (expect 0.10-0.30)\n", ss_h$PINK1_act))
  cat(sprintf("    mito_damage = %.4f (expect 0.02-0.10)\n", ss_h$mito_damage))
  cat(sprintf("    DA_neuron   = %.4f at %dh (expect > 0.95)\n",
              ss_h$DA_neuron, t_bido))
  cat(sprintf("    DA loss     = %.2f%%\n", (1 - ss_h$DA_neuron) * 100))
  cat(sprintf("    V1 Healthy SS & DA stability: %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation ---
  ## Upstream from M3/M1 disease: dpsi=0.935, mPTP=0.003, ATP=0.891
  cat("V2: Disease perturbation...\n")
  res_d <- run_M5(delta_psi_m = 0.935, mPTP_open = 0.003,
                  ATP = 0.891, CL_ext = 0.003,
                  duration_h = t_bido)
  ss_d  <- res_d[nrow(res_d), ]

  v2 <- (ss_d$mito_damage > ss_h$mito_damage &&
         ss_d$PINK1_act > ss_h$PINK1_act &&
         ss_d$DA_neuron < ss_h$DA_neuron)
  cat(sprintf("    mito_damage: %.4f -> %.4f (increased: %s)\n",
              ss_h$mito_damage, ss_d$mito_damage,
              ifelse(ss_d$mito_damage > ss_h$mito_damage, "YES", "NO")))
  cat(sprintf("    PINK1_act:   %.4f -> %.4f (stress response: %s)\n",
              ss_h$PINK1_act, ss_d$PINK1_act,
              ifelse(ss_d$PINK1_act > ss_h$PINK1_act, "YES", "NO")))
  cat(sprintf("    DA_neuron:   %.4f -> %.4f (loss: %.2f%% vs %.2f%%)\n",
              ss_h$DA_neuron, ss_d$DA_neuron,
              (1 - ss_d$DA_neuron) * 100,
              (1 - ss_h$DA_neuron) * 100))
  cat(sprintf("    Note: standalone effect is modest; vicious cycle amplifies\n"))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue ---
  ## Drug improves M3 outputs: dpsi 0.935->0.940, ATP 0.891->0.899
  cat("V3: Drug rescue (improved ETC from M1->M3)...\n")
  res_drug <- run_M5(delta_psi_m = 0.940, mPTP_open = 0.003,
                     ATP = 0.899, CL_ext = 0.002,
                     duration_h = t_bido)
  ss_rx <- res_drug[nrow(res_drug), ]

  v3 <- (ss_rx$mito_damage < ss_d$mito_damage &&
         ss_rx$DA_neuron > ss_d$DA_neuron)
  cat(sprintf("    mito_damage: disease=%.4f -> drug=%.4f\n",
              ss_d$mito_damage, ss_rx$mito_damage))
  cat(sprintf("    DA_neuron:   disease=%.4f -> drug=%.4f\n",
              ss_d$DA_neuron, ss_rx$DA_neuron))
  cat(sprintf("    DA loss improvement: %.2f pp\n",
              ((1 - ss_d$DA_neuron) - (1 - ss_rx$DA_neuron)) * 100))
  cat(sprintf("    V3 Drug rescue:             %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero upstream damage → no damage, DA stable ---
  cat("V4: Zero upstream damage -> zero mito_damage...\n")
  res_zero <- run_M5(delta_psi_m = 1.0, mPTP_open = 0.0,
                     ATP = 1.0, CL_ext = 0.0,
                     duration_h = t_bido,
                     y0 = c(PINK1_act = 0.0, mito_damage = 0.0,
                            DA_neuron = 1.0))
  ss_z <- res_zero[nrow(res_zero), ]

  v4 <- (abs(ss_z$PINK1_act) < 1e-10 &&
         abs(ss_z$mito_damage) < 1e-10 &&
         abs(ss_z$DA_neuron - 1.0) < 1e-10)
  cat(sprintf("    PINK1_act   = %.2e (expect 0)\n", ss_z$PINK1_act))
  cat(sprintf("    mito_damage = %.2e (expect 0)\n", ss_z$mito_damage))
  cat(sprintf("    DA_neuron   = %.10f (expect 1.0)\n", ss_z$DA_neuron))
  cat(sprintf("    V4 Zero damage:             %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness [0,1] ---
  cat("V5: Non-negativity and boundedness...\n")
  all_res <- rbind(res_h, res_d, res_drug)
  v5 <- (all(all_res$PINK1_act >= -1e-10) &&
         all(all_res$PINK1_act <= 1 + 1e-10) &&
         all(all_res$mito_damage >= -1e-10) &&
         all(all_res$mito_damage <= 1 + 1e-10) &&
         all(all_res$DA_neuron >= -1e-10) &&
         all(all_res$DA_neuron <= 1 + 1e-10))
  cat(sprintf("    PINK1 range:  [%.4f, %.4f]\n",
              min(all_res$PINK1_act), max(all_res$PINK1_act)))
  cat(sprintf("    Damage range: [%.4f, %.4f]\n",
              min(all_res$mito_damage), max(all_res$mito_damage)))
  cat(sprintf("    DA range:     [%.4f, %.4f]\n",
              min(all_res$DA_neuron), max(all_res$DA_neuron)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical vs analytical match ---
  ## PINK1 and mito_damage: analytical SS match (<1%)
  ## DA_neuron: analytical exponential decay match
  cat("V6: Numerical vs analytical steady state...\n")
  ss_an_d <- M5_analytical_ss(0.935, 0.003, 0.891, 0.003)

  ## PINK1 and damage: check at midpoint (t=100h) when SS is reached
  mid_h <- res_h[abs(res_h$time - 100) < 0.3, ][1, ]
  mid_d <- res_d[abs(res_d$time - 100) < 0.3, ][1, ]

  err_PINK1_h <- abs(mid_h$PINK1_act - ss_an_h["PINK1_act"]) /
                 ss_an_h["PINK1_act"] * 100
  err_PINK1_d <- abs(mid_d$PINK1_act - ss_an_d["PINK1_act"]) /
                 ss_an_d["PINK1_act"] * 100
  err_D_h     <- abs(mid_h$mito_damage - ss_an_h["mito_damage"]) /
                 ss_an_h["mito_damage"] * 100
  err_D_d     <- abs(mid_d$mito_damage - ss_an_d["mito_damage"]) /
                 ss_an_d["mito_damage"] * 100

  ## DA_neuron: compare numerical at t=840 with analytical decay
  DA_an_h <- as.numeric(
    M5_DA_analytical(t_bido, 0.960, 0.002, 0.954, 0.002, DA_0 = 1.0))
  DA_an_d <- as.numeric(
    M5_DA_analytical(t_bido, 0.935, 0.003, 0.891, 0.003, DA_0 = 1.0))

  err_DA_h <- abs(ss_h$DA_neuron - DA_an_h) / DA_an_h * 100
  err_DA_d <- abs(ss_d$DA_neuron - DA_an_d) / DA_an_d * 100

  v6 <- (err_PINK1_h < 1.0 && err_PINK1_d < 1.0 &&
         err_D_h < 1.0 && err_D_d < 1.0 &&
         err_DA_h < 1.0 && err_DA_d < 1.0)
  cat(sprintf("    Healthy PINK1: num=%.5f ana=%.5f (err: %.4f%%)\n",
              mid_h$PINK1_act, ss_an_h["PINK1_act"], err_PINK1_h))
  cat(sprintf("    Disease PINK1: num=%.5f ana=%.5f (err: %.4f%%)\n",
              mid_d$PINK1_act, ss_an_d["PINK1_act"], err_PINK1_d))
  cat(sprintf("    Healthy D:     num=%.5f ana=%.5f (err: %.4f%%)\n",
              mid_h$mito_damage, ss_an_h["mito_damage"], err_D_h))
  cat(sprintf("    Disease D:     num=%.5f ana=%.5f (err: %.4f%%)\n",
              mid_d$mito_damage, ss_an_d["mito_damage"], err_D_d))
  cat(sprintf("    Healthy DA(%dh): num=%.5f ana=%.5f (err: %.4f%%)\n",
              t_bido, ss_h$DA_neuron, DA_an_h, err_DA_h))
  cat(sprintf("    Disease DA(%dh): num=%.5f ana=%.5f (err: %.4f%%)\n",
              t_bido, ss_d$DA_neuron, DA_an_d, err_DA_d))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M5 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M5_scenarios <- function() {

  t_bido <- 840  # 5 weeks

  ## Scenario definitions
  ## Upstream values from M3 standalone + M1 (CL_ext)
  ## Drug effect is INDIRECT via M1->M3: improved CL -> better ETC
  scenarios <- list(
    healthy  = list(dpsi = 0.960, mPTP = 0.002, ATP = 0.954,
                    CL_ext = 0.002, label = "Healthy"),
    disease  = list(dpsi = 0.935, mPTP = 0.003, ATP = 0.891,
                    CL_ext = 0.003, label = "Disease (no drug)"),
    drug_low = list(dpsi = 0.937, mPTP = 0.003, ATP = 0.894,
                    CL_ext = 0.003, label = "Disease + SBT-272 low"),
    drug_hi  = list(dpsi = 0.940, mPTP = 0.003, ATP = 0.899,
                    CL_ext = 0.002, label = "Disease + SBT-272 high")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M5(delta_psi_m = s$dpsi, mPTP_open = s$mPTP,
                  ATP = s$ATP, CL_ext = s$CL_ext,
                  duration_h = t_bido)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))

  ## Convert time to days for readability
  all_data$time_days <- all_data$time / 24


  ## --- Plot 1: Time course of all state variables ---
  plot_vars <- c("PINK1_act", "mito_damage", "DA_neuron")

  long_data <- all_data %>%
    select(time_days, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("PINK1/Parkin activation",
               "Mitochondrial damage",
               "DA neuron survival"))

  p1 <- ggplot(long_data, aes(x = time_days, y = value, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (no drug)" = "#d62728",
                                  "Disease + SBT-272 low" = "#ff7f0e",
                                  "Disease + SBT-272 high" = "#1f77b4")) +
    labs(
      title = "Module M5: Mitophagy & DA Neuron Survival",
      subtitle = paste("5-week simulation (Bido protocol) |",
                        "Upstream: M3 (ETC), M1 (CL_ext)"),
      x = "Time (days)",
      y = "Normalised level",
      color = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave("M5_time_course.png", p1,
         width = 10, height = 10, dpi = 300, bg = "white")
  cat("Saved: M5_time_course.png\n")


  ## --- Plot 2: mPTP sweep → mito_damage and DA at 840h ---
  ## Shows what integration would produce at higher mPTP values
  cat("Running mPTP sweep (projected integration range)...\n")
  mPTP_sweep <- seq(0.0, 0.30, by = 0.005)
  sweep_data <- data.frame()

  for (mp in mPTP_sweep) {
    ss_vals <- M5_analytical_ss(delta_psi_m = 0.935, mPTP_open = mp,
                                ATP = 0.891, CL_ext = 0.003)
    DA_5wk <- as.numeric(
      M5_DA_analytical(t_bido, 0.935, mp, 0.891, 0.003, DA_0 = 1.0))
    sweep_data <- rbind(sweep_data, data.frame(
      mPTP_open   = mp,
      mito_damage = ss_vals["mito_damage"],
      Hill_death  = ss_vals["Hill_death"],
      DA_5wk      = DA_5wk,
      DA_loss_pct = (1 - DA_5wk) * 100
    ))
  }

  p2 <- ggplot(sweep_data) +
    geom_line(aes(x = mPTP_open, y = DA_5wk), linewidth = 1,
              color = "#1f77b4") +
    geom_line(aes(x = mPTP_open, y = mito_damage), linewidth = 1,
              color = "#d62728", linetype = "dashed") +
    geom_hline(yintercept = 0.50, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = c(0.003, 0.10), linetype = "dotted",
               color = "gray50") +
    annotate("text", x = 0.003, y = 0.95,
             label = "Standalone\ndisease", hjust = -0.1, size = 3,
             color = "#d62728") +
    annotate("text", x = 0.10, y = 0.95,
             label = "Projected\nintegrated", hjust = -0.1, size = 3,
             color = "gray40") +
    annotate("text", x = 0.25, y = 0.52,
             label = "Bernheimer\nthreshold", size = 3, color = "gray50") +
    annotate("text", x = 0.20, y = 0.25,
             label = "DA neuron survival", size = 3.5, color = "#1f77b4") +
    annotate("text", x = 0.25, y = 0.85,
             label = "Mito damage (SS)", size = 3.5, color = "#d62728") +
    labs(
      title = "M5: DA Neuron Survival vs mPTP Opening",
      subtitle = "5-week projection | Shows threshold for catastrophic DA loss",
      x = "mPTP_open (from M3)",
      y = "Value at 840 hours",
      caption = paste("At standalone disease (mPTP=0.003), DA loss is modest.",
                       "Integration amplifies mPTP to ~0.10, driving",
                       "clinically meaningful DA loss.")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M5_mPTP_sweep.png", p2,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: M5_mPTP_sweep.png\n")


  ## --- Plot 3: Sensitivity tornado (OAT ±50% on mito_damage) ---
  cat("Running sensitivity analysis...\n")
  base_ss <- M5_analytical_ss(0.935, 0.003, 0.891, 0.003)
  base_D  <- base_ss["mito_damage"]

  p_base <- M5_parameters(0.935, 0.003, 0.891, 0.003)

  param_names <- c("k_PINK1_on", "k_PINK1_CL", "k_PINK1_off",
                    "k_damage_mPTP", "k_damage_energy", "k_damage_dpsi",
                    "k_repair", "k_biogen", "k_death", "K_death")

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      ## Recompute PINK1 SS
      act_mod <- as.numeric(p_mod["k_PINK1_on"]) * (1 - 0.935) +
                 as.numeric(p_mod["k_PINK1_CL"]) * 0.003
      deact_mod <- as.numeric(p_mod["k_PINK1_off"]) * 0.935
      PINK1_mod <- act_mod / (act_mod + deact_mod)

      ## Recompute damage SS
      a_mod <- (as.numeric(p_mod["k_damage_mPTP"]) * 0.003 +
                as.numeric(p_mod["k_damage_energy"]) * 0.109 +
                as.numeric(p_mod["k_damage_dpsi"]) * 0.065)
      b_mod <- as.numeric(p_mod["k_repair"]) * PINK1_mod +
               as.numeric(p_mod["k_biogen"]) * 0.891
      D_mod <- a_mod / (a_mod + b_mod)

      pct_change <- as.numeric((D_mod - base_D) / base_D * 100)

      sens_data <- rbind(sens_data, data.frame(
        parameter  = pname,
        direction  = ifelse(direction > 0, "+50%", "-50%"),
        mito_damage = D_mod,
        pct_change = pct_change
      ))
    }
  }

  tornado <- sens_data %>%
    group_by(parameter) %>%
    summarise(
      low  = min(pct_change, na.rm = TRUE),
      high = max(pct_change, na.rm = TRUE),
      span = abs(max(pct_change, na.rm = TRUE) -
                 min(pct_change, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(span)

  tornado$parameter <- factor(tornado$parameter, levels = tornado$parameter)

  p3 <- ggplot(tornado, aes(y = parameter)) +
    geom_segment(aes(x = low, xend = high, yend = parameter),
                 linewidth = 6, color = "#e377c2", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M5: Parameter Sensitivity on Disease mito_damage",
      subtitle = "OAT +/-50% perturbation (dpsi=0.935, mPTP=0.003, ATP=0.891)",
      x = "Change in mito_damage_ss (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M5_sensitivity.png", p3,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: M5_sensitivity.png\n\n")


  return(list(data = all_data, sweep = sweep_data, sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

## Run verification
all_passed <- verify_M5()

## Run scenarios and generate plots
if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Running scenarios anyway for diagnostic purposes...\n\n")
}

results <- run_all_M5_scenarios()

## Print SS summary and projected integration impact
cat("Analytical SS and DA survival (5 weeks):\n")
t_5wk <- 840
for (label in c("Healthy", "Disease", "Drug_low", "Drug_high")) {
  if (label == "Healthy") {
    ss <- M5_analytical_ss(0.960, 0.002, 0.954, 0.002)
    DA <- M5_DA_analytical(t_5wk, 0.960, 0.002, 0.954, 0.002)
  } else if (label == "Disease") {
    ss <- M5_analytical_ss(0.935, 0.003, 0.891, 0.003)
    DA <- M5_DA_analytical(t_5wk, 0.935, 0.003, 0.891, 0.003)
  } else if (label == "Drug_low") {
    ss <- M5_analytical_ss(0.937, 0.003, 0.894, 0.003)
    DA <- M5_DA_analytical(t_5wk, 0.937, 0.003, 0.894, 0.003)
  } else {
    ss <- M5_analytical_ss(0.940, 0.003, 0.899, 0.002)
    DA <- M5_DA_analytical(t_5wk, 0.940, 0.003, 0.899, 0.002)
  }
  cat(sprintf("  %-10s PINK1=%.3f D=%.4f death_rate=%.2e DA(5wk)=%.4f (%.1f%% loss)\n",
              label, ss["PINK1_act"], ss["mito_damage"],
              ss["death_rate"], DA, (1-DA)*100))
}

## Projected integration impact
cat("\nProjected integration impact (amplified mPTP values):\n")
for (mp in c(0.003, 0.05, 0.10, 0.15, 0.20)) {
  ss <- M5_analytical_ss(0.935, mp, 0.891, 0.003)
  DA <- M5_DA_analytical(t_5wk, 0.935, mp, 0.891, 0.003)
  cat(sprintf("  mPTP=%.3f -> D=%.3f death_rate=%.4f DA(5wk)=%.3f (%.0f%% loss)\n",
              mp, ss["mito_damage"], ss["death_rate"], DA, (1-DA)*100))
}

cat("\nSensitivity analysis results:\n")
print(results$sensitivity[, c("parameter", "direction", "pct_change")],
      row.names = FALSE)

cat("\n========================================================\n")
cat("  MODULE M5 COMPLETE\n")
cat("  Output: DA_neuron -> M7 (motor score)\n")
cat("  Upstream: delta_psi_m, mPTP_open, ATP <- M3\n")
cat("            CL_ext <- M1\n")
cat("  Standalone: modest DA loss (1-5%% over 5 weeks)\n")
cat("  Integrated: projected 30-50%% DA loss via vicious cycle\n")
cat("========================================================\n")
