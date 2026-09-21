##=============================================================================
## MODULE M1: CARDIOLIPIN DYNAMICS (v2 — upstream-coupled)
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## State variables (5):
##   CL_n    — Native (mature) cardiolipin fraction [0,1]
##   CL_ox   — Oxidised/damaged CL fraction [0,1]
##   CL_ext  — Externalised CL on outer mitochondrial membrane [0,1]
##   ALCAT1  — ALCAT1 enzyme active fraction [0,1]
##   TAZ     — Tafazzin active fraction [0,1]
##
## Upstream coupling:
##   mROS        ← M4 (mitochondrial ROS, normalised)
##   C_mito_norm ← M0 (drug mitochondrial conc, normalised by C_ref = 4.43)
##   aSyn_oligo  — scenario parameter until M2 is coded
##
## Key outputs:
##   CL_ratio = CL_n / (CL_n + CL_ox) → M3 (supercomplex stability)
##   CL_ext   → M2 (aSyn seeding), M5 (mitophagy signal)
##
## Calibration targets:
##   Healthy  CL_ratio ≈ 0.94  (CL_n ≈ 0.91)
##   Disease  CL_ratio ≈ 0.79  (CL_n ~22% reduction, Gao 2017)
##   Drug rescue improves CL_ratio toward healthy
##
## Data sources:
##   Gao 2017 (doi:10.1093/hmg/ddx100): CL_n 23% reduction in PD model
##   Song et al. 2019 (doi:10.1111/acel.12941): ALCAT1 upregulation
##   Bolte 2004: CL half-life ≈ 27 days
##   Kagan 2015: CL externalisation as apoptotic/mitophagy signal
##   Schlame & Greenberg 2017: tafazzin remodelling
##   Bayir 2009: aSyn–CL interaction and peroxidation
##   Tung et al. 2025 (doi:10.3390/ijms26030944): bevemipretide mechanism
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M1: CARDIOLIPIN DYNAMICS — Bevemipretide QSP\n")
cat("  Version 2.0 | Upstream-coupled to M0 (PK), M4 (ROS)\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M1_parameters <- function(mROS = 0.10, aSyn_oligo = 0.05, C_mito_norm = 0.0) {

  p <- c(
    ## --- CL synthesis and capacity ---
    k_syn       = 0.007,    # CL synthesis rate (h^-1), logistic growth
                             # ESTIMATED: derived from CL t1/2 ~ 27d (Bolte 2004),
                             # scaled for logistic carrying-capacity model
    K_CL        = 1.0,      # CL carrying capacity (normalised maximum)

    ## --- CL oxidation ---
    k_ox        = 0.020,    # mROS-driven CL oxidation rate (h^-1)
                             # ESTIMATED: tuned so disease mROS=0.176 gives
                             # ~22% CL_n reduction (Gao 2017)
    k_ALCAT     = 5.0,      # ALCAT1 amplification of oxidation (dimensionless)
                             # ASSUMED: Song et al. 2019
                             # MPTP -> 3x ALCAT1 -> amplified CL remodelling
    k_aSyn_ox   = 0.020,    # Direct aSyn-mediated CL damage (h^-1)
                             # ASSUMED: Bayir 2009, aSyn binds CL -> peroxidation

    ## --- CL repair ---
    k_red       = 0.070,    # TAZ-mediated CL repair/remodelling (h^-1)
                             # ASSUMED: tafazzin remodels CL_ox -> CL_n
                             # Schlame & Greenberg 2017
    k_drug      = 0.100,    # Bevemipretide CL stabilisation (h^-1)
                             # ESTIMATED: to be fitted to Bido et al. poster
                             # Tung et al. 2025

    ## --- CL degradation ---
    k_degrad_ox = 0.010,    # CL_ox degradation/mitophagy (h^-1)
                             # ASSUMED: ~70h clearance for damaged CL
                             # Chu 2018 (doi:10.1016/j.neulet.2018.04.004)

    ## --- CL externalisation ---
    k_ext       = 0.001,    # CL flip to OMM (h^-1)
                             # ASSUMED: Kagan 2015, apoptotic CL signal
    k_clear_ext = 0.050,    # CL_ext clearance (h^-1)
                             # ASSUMED: phagocytic/mitophagy clearance

    ## --- ALCAT1 dynamics ---
    k_act_A     = 0.50,     # ALCAT1 activation by mROS (h^-1)
                             # ESTIMATED: Song et al. 2019
    k_deact_A   = 0.45,     # ALCAT1 deactivation (h^-1)
                             # ESTIMATED: sets ALCAT1_healthy ~ 0.10

    ## --- Tafazzin dynamics ---
    k_act_T     = 0.50,     # TAZ constitutive activation (h^-1)
                             # ASSUMED: maintains TAZ ~ 0.90 at healthy mROS
    k_deact_T   = 0.556,    # TAZ inactivation by mROS (h^-1)
                             # ESTIMATED: Falabella et al. 2021

    ## --- Upstream coupling (scenario parameters) ---
    mROS        = mROS,
    aSyn_oligo  = aSyn_oligo,
    C_mito_norm = C_mito_norm
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M1_init <- function() {
  c(
    CL_n   = 0.90,    # Near healthy SS
    CL_ox  = 0.05,    # Low basal oxidation
    CL_ext = 0.001,   # Minimal externalised CL
    ALCAT1 = 0.10,    # Low basal ALCAT1
    TAZ    = 0.90     # High tafazzin activity
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M1_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states to [0, 1] ---
    CL_n   <- max(0, min(1, CL_n))
    CL_ox  <- max(0, min(1, CL_ox))
    CL_ext <- max(0, min(1, CL_ext))
    ALCAT1 <- max(0, min(1, ALCAT1))
    TAZ    <- max(0, min(1, TAZ))

    ## Proxy for mPTP/delta_psi damage signal until M3 is built
    damage_signal <- mROS

    ## --- ODE 1: Native cardiolipin ---
    dCL_n <- (k_syn * (1 - CL_n / K_CL)
              - k_ox * mROS * CL_n * (1 + k_ALCAT * ALCAT1)
              - k_aSyn_ox * aSyn_oligo * CL_n
              + k_red * TAZ * CL_ox
              + k_drug * C_mito_norm * CL_ox
              - k_ext * CL_n * damage_signal)

    ## --- ODE 2: Oxidised cardiolipin ---
    dCL_ox <- (k_ox * mROS * CL_n * (1 + k_ALCAT * ALCAT1)
               + k_aSyn_ox * aSyn_oligo * CL_n
               - k_red * TAZ * CL_ox
               - k_drug * C_mito_norm * CL_ox
               - k_degrad_ox * CL_ox)

    ## --- ODE 3: Externalised cardiolipin ---
    dCL_ext <- (k_ext * CL_n * damage_signal
                - k_clear_ext * CL_ext)

    ## --- ODE 4: ALCAT1 active fraction ---
    dALCAT1 <- (k_act_A * mROS * (1 - ALCAT1)
                - k_deact_A * ALCAT1)

    ## --- ODE 5: Tafazzin active fraction ---
    dTAZ <- (k_act_T * (1 - TAZ)
             - k_deact_T * mROS * TAZ)

    ## --- Derived quantities ---
    CL_total   <- CL_n + CL_ox
    CL_ratio   <- CL_n / max(CL_total, 1e-8)
    drug_flux  <- k_drug * C_mito_norm * CL_ox
    ox_flux    <- k_ox * mROS * CL_n * (1 + k_ALCAT * ALCAT1)
    ext_flux   <- k_ext * CL_n * damage_signal

    list(
      c(dCL_n, dCL_ox, dCL_ext, dALCAT1, dTAZ),
      CL_total  = CL_total,
      CL_ratio  = CL_ratio,
      drug_flux = drug_flux,
      ox_flux   = ox_flux,
      ext_flux  = ext_flux
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M1 <- function(mROS = 0.10, aSyn_oligo = 0.05, C_mito_norm = 0.0,
                   duration_h = 2000, dt = 0.5, y0 = NULL) {

  parms <- M1_parameters(mROS, aSyn_oligo, C_mito_norm)
  if (is.null(y0)) y0 <- M1_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M1_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$mROS        <- mROS
  result$aSyn_oligo  <- aSyn_oligo
  result$C_mito_norm <- C_mito_norm

  ss <- result[nrow(result), ]
  cat(sprintf(
    "M1 | mROS=%.3f aSyn=%.2f C_mito=%.2f | SS: CL_n=%.3f CL_ox=%.3f CL_ext=%.4f ALCAT1=%.3f TAZ=%.3f CL_ratio=%.3f\n",
    mROS, aSyn_oligo, C_mito_norm,
    ss$CL_n, ss$CL_ox, ss$CL_ext, ss$ALCAT1, ss$TAZ, ss$CL_ratio))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE
##-----------------------------------------------------------------------------

M1_analytical_ss <- function(mROS = 0.10, aSyn_oligo = 0.05, C_mito_norm = 0.0) {
  ## Fully analytical: ALCAT1/TAZ decouple, then CL_n/CL_ox is a 2x2 linear system

  p <- M1_parameters(mROS, aSyn_oligo, C_mito_norm)

  ## ALCAT1 and TAZ steady states (independent of CL)
  ALCAT1_ss <- as.numeric(p["k_act_A"] * mROS /
                           (p["k_act_A"] * mROS + p["k_deact_A"]))

  TAZ_ss <- as.numeric(p["k_act_T"] /
                        (p["k_act_T"] + p["k_deact_T"] * mROS))

  ## Oxidation and clearance coefficients
  alpha <- as.numeric(
    p["k_ox"] * mROS * (1 + p["k_ALCAT"] * ALCAT1_ss) +
    p["k_aSyn_ox"] * aSyn_oligo)

  beta <- as.numeric(
    p["k_red"] * TAZ_ss + p["k_drug"] * C_mito_norm + p["k_degrad_ox"])

  ## CL_n from mass balance: dCL_n + dCL_ox = 0 at SS
  ## k_syn*(1 - CL_n/K_CL) = k_degrad_ox*(alpha/beta)*CL_n + k_ext*CL_n*mROS
  denom <- as.numeric(
    p["k_syn"] / p["K_CL"] +
    p["k_degrad_ox"] * alpha / beta +
    p["k_ext"] * mROS)

  CL_n_ss   <- as.numeric(p["k_syn"] / denom)
  CL_ox_ss  <- (alpha / beta) * CL_n_ss
  CL_ext_ss <- as.numeric(p["k_ext"]) * CL_n_ss * mROS /
               as.numeric(p["k_clear_ext"])

  CL_ratio_ss <- CL_n_ss / (CL_n_ss + CL_ox_ss)

  return(c(CL_n    = CL_n_ss,
           CL_ox   = CL_ox_ss,
           CL_ext  = CL_ext_ss,
           ALCAT1  = ALCAT1_ss,
           TAZ     = TAZ_ss,
           CL_ratio = CL_ratio_ss))
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M1 <- function() {

  cat("=== MODULE M1 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6
  t_verify    <- 3000


  ## --- V1: Healthy steady state in expected range ---
  cat("V1: Healthy steady state...\n")
  res_h <- run_M1(mROS = 0.10, aSyn_oligo = 0.05, C_mito_norm = 0.0,
                  duration_h = t_verify)
  ss_h  <- res_h[nrow(res_h), ]

  v1 <- (ss_h$CL_n > 0.85 && ss_h$CL_n < 0.95 &&
         ss_h$CL_ox > 0.03 && ss_h$CL_ox < 0.10 &&
         ss_h$ALCAT1 > 0.05 && ss_h$ALCAT1 < 0.20 &&
         ss_h$TAZ > 0.80 && ss_h$TAZ < 0.95 &&
         ss_h$CL_ratio > 0.90 && ss_h$CL_ratio < 0.98)
  cat(sprintf("    CL_n = %.4f (expect 0.85-0.95)\n", ss_h$CL_n))
  cat(sprintf("    CL_ox = %.4f (expect 0.03-0.10)\n", ss_h$CL_ox))
  cat(sprintf("    CL_ext = %.5f\n", ss_h$CL_ext))
  cat(sprintf("    ALCAT1 = %.4f (expect 0.05-0.20)\n", ss_h$ALCAT1))
  cat(sprintf("    TAZ = %.4f (expect 0.80-0.95)\n", ss_h$TAZ))
  cat(sprintf("    CL_ratio = %.4f (expect 0.90-0.98)\n", ss_h$CL_ratio))
  cat(sprintf("    V1 Healthy SS:              %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation ---
  ## mROS = 0.176 from M4 (CI=0.49), aSyn = 0.50 (pathological)
  cat("V2: Disease perturbation (mROS=0.176, aSyn=0.50)...\n")
  res_d <- run_M1(mROS = 0.176, aSyn_oligo = 0.50, C_mito_norm = 0.0,
                  duration_h = t_verify)
  ss_d  <- res_d[nrow(res_d), ]

  pct_drop <- (1 - ss_d$CL_n / ss_h$CL_n) * 100
  v2 <- (ss_d$CL_n < ss_h$CL_n &&
         ss_d$CL_ox > ss_h$CL_ox &&
         ss_d$ALCAT1 > ss_h$ALCAT1 &&
         ss_d$TAZ < ss_h$TAZ &&
         ss_d$CL_ratio < ss_h$CL_ratio)
  cat(sprintf("    CL_n: %.4f -> %.4f (%.1f%% drop, target ~23%%)\n",
              ss_h$CL_n, ss_d$CL_n, pct_drop))
  cat(sprintf("    CL_ox: %.4f -> %.4f\n", ss_h$CL_ox, ss_d$CL_ox))
  cat(sprintf("    CL_ratio: %.4f -> %.4f\n", ss_h$CL_ratio, ss_d$CL_ratio))
  cat(sprintf("    ALCAT1: %.4f -> %.4f (upregulated: %s)\n",
              ss_h$ALCAT1, ss_d$ALCAT1,
              ifelse(ss_d$ALCAT1 > ss_h$ALCAT1, "YES", "NO")))
  cat(sprintf("    TAZ: %.4f -> %.4f (inhibited: %s)\n",
              ss_h$TAZ, ss_d$TAZ,
              ifelse(ss_d$TAZ < ss_h$TAZ, "YES", "NO")))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue (high dose C_mito_norm = 1.0) ---
  cat("V3: Drug rescue (C_mito_norm = 1.0)...\n")
  res_drug <- run_M1(mROS = 0.176, aSyn_oligo = 0.50, C_mito_norm = 1.0,
                     duration_h = t_verify)
  ss_rx <- res_drug[nrow(res_drug), ]

  v3 <- (ss_rx$CL_n > ss_d$CL_n &&
         ss_rx$CL_ox < ss_d$CL_ox &&
         ss_rx$CL_ratio > ss_d$CL_ratio)
  cat(sprintf("    CL_n: disease=%.4f -> drug=%.4f\n", ss_d$CL_n, ss_rx$CL_n))
  cat(sprintf("    CL_ox: disease=%.4f -> drug=%.4f\n", ss_d$CL_ox, ss_rx$CL_ox))
  cat(sprintf("    CL_ratio: disease=%.4f -> drug=%.4f\n",
              ss_d$CL_ratio, ss_rx$CL_ratio))
  cat(sprintf("    Improvement from disease: %.1f%%\n",
              (ss_rx$CL_ratio / ss_d$CL_ratio - 1) * 100))
  cat(sprintf("    V3 Drug rescue:             %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero drug concentration -> zero drug flux ---
  cat("V4: Zero drug concentration -> zero drug flux...\n")
  v4 <- all(abs(res_h$drug_flux) < 1e-15)
  cat(sprintf("    Max drug_flux (healthy, C_mito=0) = %.2e\n",
              max(abs(res_h$drug_flux))))
  cat(sprintf("    V4 Zero drug flux:          %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness [0,1] ---
  cat("V5: Non-negativity and boundedness...\n")
  all_res <- rbind(res_h, res_d, res_drug)
  v5 <- (all(all_res$CL_n >= -1e-10) && all(all_res$CL_n <= 1 + 1e-10) &&
         all(all_res$CL_ox >= -1e-10) && all(all_res$CL_ox <= 1 + 1e-10) &&
         all(all_res$CL_ext >= -1e-10) && all(all_res$CL_ext <= 1 + 1e-10) &&
         all(all_res$ALCAT1 >= -1e-10) && all(all_res$ALCAT1 <= 1 + 1e-10) &&
         all(all_res$TAZ >= -1e-10) && all(all_res$TAZ <= 1 + 1e-10) &&
         all((all_res$CL_n + all_res$CL_ox) <= 1.0 + 1e-6))
  cat(sprintf("    CL_n range:   [%.4e, %.4f]\n",
              min(all_res$CL_n), max(all_res$CL_n)))
  cat(sprintf("    CL_ox range:  [%.4e, %.4f]\n",
              min(all_res$CL_ox), max(all_res$CL_ox)))
  cat(sprintf("    CL_ext range: [%.4e, %.6f]\n",
              min(all_res$CL_ext), max(all_res$CL_ext)))
  cat(sprintf("    Max CL_total: %.4f\n",
              max(all_res$CL_n + all_res$CL_ox)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical SS matches analytical SS (<1% error) ---
  cat("V6: Numerical vs analytical steady state...\n")
  ss_an_h <- M1_analytical_ss(mROS = 0.10, aSyn_oligo = 0.05, C_mito_norm = 0.0)
  ss_an_d <- M1_analytical_ss(mROS = 0.176, aSyn_oligo = 0.50, C_mito_norm = 0.0)

  err_h_ratio <- abs(ss_h$CL_ratio - ss_an_h["CL_ratio"]) /
                 ss_an_h["CL_ratio"] * 100
  err_d_ratio <- abs(ss_d$CL_ratio - ss_an_d["CL_ratio"]) /
                 ss_an_d["CL_ratio"] * 100

  err_h_n <- abs(ss_h$CL_n - ss_an_h["CL_n"]) / ss_an_h["CL_n"] * 100
  err_d_n <- abs(ss_d$CL_n - ss_an_d["CL_n"]) / ss_an_d["CL_n"] * 100

  v6 <- (err_h_ratio < 1.0 && err_d_ratio < 1.0 &&
         err_h_n < 1.0 && err_d_n < 1.0)
  cat(sprintf("    Healthy CL_ratio: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$CL_ratio, ss_an_h["CL_ratio"], err_h_ratio))
  cat(sprintf("    Disease CL_ratio: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$CL_ratio, ss_an_d["CL_ratio"], err_d_ratio))
  cat(sprintf("    Healthy CL_n:     num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$CL_n, ss_an_h["CL_n"], err_h_n))
  cat(sprintf("    Disease CL_n:     num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$CL_n, ss_an_d["CL_n"], err_d_n))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M1 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M1_scenarios <- function() {

  ## Scenario definitions
  ## Drug effect in standalone M1 is DIRECT (CL stabilisation only)
  ## mROS held at disease level for drug scenarios — upstream M4 feedback
  ## will reduce mROS further in the integrated model
  scenarios <- list(
    healthy  = list(mROS = 0.10,  aSyn = 0.05, C_mito = 0.0,
                    label = "Healthy"),
    disease  = list(mROS = 0.176, aSyn = 0.50, C_mito = 0.0,
                    label = "Disease (no drug)"),
    drug_low = list(mROS = 0.176, aSyn = 0.50, C_mito = 0.10,
                    label = "Disease + SBT-272 low"),
    drug_hi  = list(mROS = 0.176, aSyn = 0.50, C_mito = 1.00,
                    label = "Disease + SBT-272 high")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M1(mROS = s$mROS, aSyn_oligo = s$aSyn, C_mito_norm = s$C_mito,
                  duration_h = 2000)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))


  ## --- Plot 1: Time course (all state variables + CL_ratio) ---
  plot_vars <- c("CL_n", "CL_ox", "CL_ext", "ALCAT1", "TAZ", "CL_ratio")

  long_data <- all_data %>%
    select(time, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("CL native", "CL oxidised", "CL externalised",
               "ALCAT1", "Tafazzin", "CL functional ratio"))

  p1 <- ggplot(long_data, aes(x = time, y = value, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (no drug)" = "#d62728",
                                  "Disease + SBT-272 low" = "#ff7f0e",
                                  "Disease + SBT-272 high" = "#1f77b4")) +
    labs(
      title = "Module M1: Cardiolipin Dynamics",
      subtitle = paste("Upstream: mROS from M4 (0.10 healthy / 0.176 disease),",
                        "C_mito from M0 (normalised by C_ref)"),
      x = "Time (hours)",
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

  ggsave("M1_time_course.png", p1,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Saved: M1_time_course.png\n")


  ## --- Plot 2: Dose-response at steady state (analytical) ---
  cat("Running dose-response sweep (analytical SS)...\n")
  dose_levels <- seq(0, 1.5, by = 0.025)
  ss_data <- data.frame()

  for (dose in dose_levels) {
    ss_vals <- M1_analytical_ss(mROS = 0.176, aSyn_oligo = 0.50,
                                C_mito_norm = dose)
    ss_data <- rbind(ss_data, data.frame(
      C_mito_norm = dose,
      CL_n     = ss_vals["CL_n"],
      CL_ox    = ss_vals["CL_ox"],
      CL_ratio = ss_vals["CL_ratio"],
      CL_ext   = ss_vals["CL_ext"]
    ))
  }

  p2 <- ggplot(ss_data, aes(x = C_mito_norm)) +
    geom_line(aes(y = CL_n, color = "CL native"), linewidth = 1) +
    geom_line(aes(y = CL_ox, color = "CL oxidised"), linewidth = 1) +
    geom_line(aes(y = CL_ratio, color = "CL functional ratio"),
              linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = c(0.10, 1.00), linetype = "dotted",
               color = "gray50", linewidth = 0.5) +
    annotate("text", x = 0.12, y = 0.92, label = "Low dose",
             size = 3, hjust = 0, color = "gray40") +
    annotate("text", x = 1.02, y = 0.92, label = "High dose",
             size = 3, hjust = 0, color = "gray40") +
    scale_color_manual(values = c("CL native" = "#2ca02c",
                                  "CL oxidised" = "#d62728",
                                  "CL functional ratio" = "#1f77b4")) +
    labs(
      title = "M1 Dose-Response: Steady-State CL Restoration",
      subtitle = "Disease background (mROS = 0.176 from M4, aSyn = 0.50)",
      x = "Bevemipretide C_mito (normalised, C_ref = 4.43)",
      y = "Normalised CL level at steady state",
      color = "Variable"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M1_dose_response.png", p2,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: M1_dose_response.png\n")


  ## --- Plot 3: Sensitivity tornado (OAT +/-50% on CL_ratio) ---
  cat("Running one-at-a-time sensitivity analysis...\n")
  base_ss <- M1_analytical_ss(mROS = 0.176, aSyn_oligo = 0.50, C_mito_norm = 1.0)
  base_CL_ratio <- base_ss["CL_ratio"]

  p_base <- M1_parameters(mROS = 0.176, aSyn_oligo = 0.50, C_mito_norm = 1.0)

  param_names <- c("k_syn", "k_ox", "k_ALCAT", "k_aSyn_ox",
                    "k_red", "k_drug", "k_degrad_ox",
                    "k_ext", "k_act_A", "k_deact_A",
                    "k_act_T", "k_deact_T")

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      ## Recompute analytical SS with modified parameter
      mod_ALCAT1 <- p_mod["k_act_A"] * 0.176 /
                    (p_mod["k_act_A"] * 0.176 + p_mod["k_deact_A"])
      mod_TAZ    <- p_mod["k_act_T"] /
                    (p_mod["k_act_T"] + p_mod["k_deact_T"] * 0.176)
      mod_alpha  <- p_mod["k_ox"] * 0.176 *
                    (1 + p_mod["k_ALCAT"] * mod_ALCAT1) +
                    p_mod["k_aSyn_ox"] * 0.50
      mod_beta   <- p_mod["k_red"] * mod_TAZ +
                    p_mod["k_drug"] * 1.0 +
                    p_mod["k_degrad_ox"]
      mod_denom  <- p_mod["k_syn"] / p_mod["K_CL"] +
                    p_mod["k_degrad_ox"] * mod_alpha / mod_beta +
                    p_mod["k_ext"] * 0.176
      mod_CL_n   <- as.numeric(p_mod["k_syn"] / mod_denom)
      mod_CL_ox  <- as.numeric((mod_alpha / mod_beta) * mod_CL_n)
      mod_ratio  <- mod_CL_n / (mod_CL_n + mod_CL_ox)

      pct_change <- as.numeric(
        (mod_ratio - base_CL_ratio) / base_CL_ratio * 100)

      sens_data <- rbind(sens_data, data.frame(
        parameter  = pname,
        direction  = ifelse(direction > 0, "+50%", "-50%"),
        CL_ratio   = mod_ratio,
        pct_change = pct_change
      ))
    }
  }

  ## Build tornado bars
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
                 linewidth = 6, color = "#1f77b4", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M1: Parameter Sensitivity on CL Functional Ratio",
      subtitle = "OAT +/-50% perturbation, disease + high dose (C_mito = 1.0)",
      x = "Change in CL_ratio (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M1_sensitivity.png", p3,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: M1_sensitivity.png\n\n")


  return(list(data = all_data, dose_response = ss_data,
              sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

## Run verification
all_passed <- verify_M1()

## Run scenarios and generate plots
if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Running scenarios anyway for diagnostic purposes...\n\n")
}

results <- run_all_M1_scenarios()

## Print dose-response summary
cat("Dose-response summary (selected points):\n")
dr <- results$dose_response
sel <- dr$C_mito_norm %in% c(0.0, 0.10, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50)
print(dr[sel, c("C_mito_norm", "CL_n", "CL_ox", "CL_ratio")], row.names = FALSE)

cat("\nSensitivity analysis results:\n")
print(results$sensitivity[, c("parameter", "direction", "pct_change")],
      row.names = FALSE)

cat("\n========================================================\n")
cat("  MODULE M1 COMPLETE\n")
cat("  Outputs: CL_ratio -> M3 (supercomplex stability)\n")
cat("           CL_ext   -> M2 (aSyn seeding), M5 (mitophagy)\n")
cat("  Upstream: mROS <- M4, C_mito <- M0, aSyn <- scenario\n")
cat("========================================================\n")
