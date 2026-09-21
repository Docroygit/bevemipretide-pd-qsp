##=============================================================================
## MODULE M3: ETC BIOENERGETICS
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## State variables (5):
##   CI_activity  — Complex I fractional activity [0,1]
##   SC_integrity — Supercomplex (respirasome) integrity [0,1]
##   delta_psi_m  — Mitochondrial membrane potential (normalised) [0,1]
##   ATP          — Mitochondrial ATP level (normalised) [0,1]
##   mPTP_open    — mPTP open fraction [0,1]
##
## Upstream coupling:
##   CL_ratio    <- M1 (supercomplex stability via cardiolipin)
##   mROS        <- M4 (CI damage, mPTP opening)
##   aSyn_oligo  — scenario parameter until M2 is coded
##
## Key outputs:
##   CI_activity  -> M4 (electron leak -> ROS production)
##   SC_integrity -> M4 (electron channeling)
##   delta_psi_m  -> M5 (PINK1 activation)
##   mPTP_open    -> M1 (damage_signal), M5 (cell death)
##   ATP          -> M2 (clearance), M5 (repair capacity)
##
## Calibration targets:
##   Healthy: CI ~ 0.95, SC ~ 0.97, dpsi ~ 0.96, ATP ~ 0.95, mPTP ~ 0.02
##   Disease: CI ~ 0.49 (Gao 2017), dpsi/ATP reduced, mPTP increased
##
## Spec deviation: aSyn damage to CI uses aSyn^2 (cooperative binding)
## instead of linear aSyn, to achieve the ~20x dynamic range needed
## between healthy (aSyn=0.05) and disease (aSyn=0.50) states.
## dpsi and ATP equations include (1-x) saturation for boundedness.
##
## Data sources:
##   Gao 2017 (doi:10.1093/hmg/ddx100): CI = 49% of WT in PD model
##   Devi 2008: aSyn binds CI subunits, inhibits activity
##   Kembro 2013 (Biophys J 104:332): proton leak g_H
##   Cortassa 2004 (Biophys J 84:2734): mitochondrial energetics
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M3: ETC BIOENERGETICS — Bevemipretide QSP\n")
cat("  Version 1.0 | Upstream-coupled to M1 (CL), M4 (ROS)\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M3_parameters <- function(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05) {

  p <- c(
    ## --- Complex I dynamics ---
    k_CI_repair = 0.050,    # CI repair/turnover rate (h^-1)
                             # ASSUMED: slow mitochondrial protein turnover
    k_CI_ROS    = 0.020,    # ROS-mediated CI damage (h^-1)
                             # ESTIMATED: secondary driver of CI loss
    k_CI_aSyn   = 0.200,    # aSyn-mediated CI damage (h^-1, applied to aSyn^2)
                             # ESTIMATED: Devi 2008, aSyn binds CI subunits
                             # Cooperative: requires oligomer accumulation

    ## --- Supercomplex dynamics ---
    k_SC_form   = 0.10,     # CL-dependent SC assembly (h^-1)
                             # ASSUMED: CL is required for respirasome stability
    k_SC_dissoc = 0.05,     # SC disassembly rate (h^-1)
                             # ASSUMED: proportional to (1 - CL_ratio)

    ## --- Membrane potential ---
    k_resp       = 1.0,     # ETC respiration rate (h^-1)
                             # ESTIMATED: overall electron transport capacity
    k_ATPase     = 0.020,   # ATP synthase proton consumption (h^-1)
                             # ASSUMED: Cortassa 2004 framework
    k_leak       = 0.010,   # Proton leak (h^-1)
                             # LITERATURE: Kembro 2013, g_H parameter
    k_mPTP_depol = 0.50,    # mPTP-mediated depolarisation (h^-1)
                             # ASSUMED: mPTP opening collapses membrane potential

    ## --- ATP dynamics ---
    k_ATP_syn     = 1.20,   # ATP synthesis rate (h^-1)
                             # ESTIMATED: tuned for ATP_healthy ~ 0.95
    k_ATP_consume = 0.025,  # Basal ATP consumption (h^-1)
                             # ASSUMED: cellular ATP demand
    k_ATP_mPTP    = 2.0,    # mPTP-mediated ATP depletion (h^-1)
                             # ASSUMED: cytochrome c release + uncoupling

    ## --- mPTP dynamics ---
    k_mPTP_open  = 0.50,    # mPTP opening rate (h^-1)
                             # ASSUMED: phenomenological
    k_mPTP_close = 0.50,    # mPTP closing rate (h^-1)
                             # ASSUMED: dpsi-dependent closure
    K_mPTP       = 0.40,    # Hill half-max for ROS-driven opening
                             # ASSUMED: threshold for mPTP activation
    n_Hill       = 3,       # Hill coefficient (cooperativity)
                             # ASSUMED: cooperative Ca2+/ROS gating

    ## --- Upstream coupling (scenario parameters) ---
    CL_ratio   = CL_ratio,
    mROS       = mROS,
    aSyn_oligo = aSyn_oligo
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M3_init <- function() {
  c(
    CI_activity  = 0.95,    # Near healthy SS
    SC_integrity = 0.97,    # Near healthy SS
    delta_psi_m  = 0.96,    # Near healthy SS
    ATP          = 0.95,    # Near healthy SS
    mPTP_open    = 0.015    # Nearly closed
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M3_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states to [0, 1] ---
    CI_activity  <- max(0, min(1, CI_activity))
    SC_integrity <- max(0, min(1, SC_integrity))
    delta_psi_m  <- max(0, min(1, delta_psi_m))
    ATP          <- max(0, min(1, ATP))
    mPTP_open    <- max(0, min(1, mPTP_open))

    ## --- Hill function for mPTP ROS sensitivity ---
    Hill_ROS <- mROS^n_Hill / (K_mPTP^n_Hill + mROS^n_Hill)

    ## --- ODE 1: Complex I activity ---
    ## Repair toward full activity; damage by ROS and aSyn oligomers
    ## aSyn^2 gives cooperative damage (Devi 2008)
    dCI <- (k_CI_repair * (1 - CI_activity)
            - k_CI_ROS * mROS * CI_activity
            - k_CI_aSyn * aSyn_oligo^2 * CI_activity)

    ## --- ODE 2: Supercomplex integrity ---
    ## CL-dependent assembly; (1-CL_ratio)-dependent disassembly
    dSC <- (k_SC_form * CL_ratio * (1 - SC_integrity)
            - k_SC_dissoc * (1 - CL_ratio) * SC_integrity)

    ## --- ODE 3: Membrane potential ---
    ## Respiration charges; ATPase + leak + mPTP discharge
    ## (1 - dpsi) saturation ensures boundedness
    dPsi <- (k_resp * CI_activity * SC_integrity * (1 - delta_psi_m)
             - k_ATPase * delta_psi_m
             - k_leak * delta_psi_m
             - k_mPTP_depol * mPTP_open * delta_psi_m)

    ## --- ODE 4: ATP level ---
    ## Synthesis driven by pmf * ETC capacity; consumption + mPTP loss
    ## (1 - ATP) saturation ensures boundedness
    dATP <- (k_ATP_syn * delta_psi_m * CI_activity * SC_integrity * (1 - ATP)
             - k_ATP_consume * ATP
             - k_ATP_mPTP * mPTP_open * ATP)

    ## --- ODE 5: mPTP open fraction ---
    ## ROS-driven opening (Hill); dpsi-dependent closure
    dmPTP <- (k_mPTP_open * Hill_ROS * (1 - mPTP_open)
              - k_mPTP_close * mPTP_open * delta_psi_m)

    ## --- Derived quantities ---
    ETC_capacity <- CI_activity * SC_integrity
    resp_flux    <- k_resp * CI_activity * SC_integrity * (1 - delta_psi_m)
    ATP_flux     <- k_ATP_syn * delta_psi_m * CI_activity * SC_integrity * (1 - ATP)

    list(
      c(dCI, dSC, dPsi, dATP, dmPTP),
      ETC_capacity = ETC_capacity,
      Hill_ROS     = Hill_ROS,
      resp_flux    = resp_flux,
      ATP_flux     = ATP_flux
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M3 <- function(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05,
                   duration_h = 200, dt = 0.1, y0 = NULL) {

  parms <- M3_parameters(CL_ratio, mROS, aSyn_oligo)
  if (is.null(y0)) y0 <- M3_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M3_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$CL_ratio   <- CL_ratio
  result$mROS        <- mROS
  result$aSyn_oligo  <- aSyn_oligo

  ss <- result[nrow(result), ]
  cat(sprintf(
    "M3 | CL=%.3f mROS=%.3f aSyn=%.2f | SS: CI=%.3f SC=%.3f dpsi=%.3f ATP=%.3f mPTP=%.4f\n",
    CL_ratio, mROS, aSyn_oligo,
    ss$CI_activity, ss$SC_integrity, ss$delta_psi_m, ss$ATP, ss$mPTP_open))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE
##-----------------------------------------------------------------------------

M3_analytical_ss <- function(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05) {

  p <- M3_parameters(CL_ratio, mROS, aSyn_oligo)

  ## CI and SC decouple from other M3 states
  CI_ss <- as.numeric(
    p["k_CI_repair"] /
    (p["k_CI_repair"] + p["k_CI_ROS"] * mROS + p["k_CI_aSyn"] * aSyn_oligo^2))

  SC_ss <- as.numeric(
    p["k_SC_form"] * CL_ratio /
    (p["k_SC_form"] * CL_ratio + p["k_SC_dissoc"] * (1 - CL_ratio)))

  ## dpsi and mPTP are coupled — solve via uniroot
  A <- as.numeric(p["k_resp"]) * CI_ss * SC_ss
  B <- as.numeric(p["k_ATPase"] + p["k_leak"])
  k_dep <- as.numeric(p["k_mPTP_depol"])
  k_o   <- as.numeric(p["k_mPTP_open"])
  k_c   <- as.numeric(p["k_mPTP_close"])

  H <- mROS^as.numeric(p["n_Hill"]) /
       (as.numeric(p["K_mPTP"])^as.numeric(p["n_Hill"]) +
        mROS^as.numeric(p["n_Hill"]))

  if (H < 1e-15) {
    mPTP_ss   <- 0
    dpsi_ss   <- A / (A + B)
  } else {
    f_mPTP <- function(m) {
      dpsi <- A / (A + B + k_dep * m)
      k_o * H * (1 - m) - k_c * m * dpsi
    }
    sol <- uniroot(f_mPTP, interval = c(0, 1 - 1e-10), tol = 1e-12)
    mPTP_ss <- sol$root
    dpsi_ss <- A / (A + B + k_dep * mPTP_ss)
  }

  ## ATP from saturating SS
  P_atp <- as.numeric(p["k_ATP_syn"]) * dpsi_ss * CI_ss * SC_ss
  Q_atp <- as.numeric(p["k_ATP_consume"]) +
           as.numeric(p["k_ATP_mPTP"]) * mPTP_ss
  ATP_ss <- P_atp / (P_atp + Q_atp)

  return(c(CI_activity  = CI_ss,
           SC_integrity = SC_ss,
           delta_psi_m  = dpsi_ss,
           ATP          = ATP_ss,
           mPTP_open    = mPTP_ss))
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M3 <- function() {

  cat("=== MODULE M3 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6
  t_verify    <- 500


  ## --- V1: Healthy steady state in expected range ---
  cat("V1: Healthy steady state...\n")
  res_h <- run_M3(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05,
                  duration_h = t_verify)
  ss_h  <- res_h[nrow(res_h), ]

  v1 <- (ss_h$CI_activity > 0.90 && ss_h$CI_activity <= 1.0 &&
         ss_h$SC_integrity > 0.95 &&
         ss_h$delta_psi_m > 0.90 &&
         ss_h$ATP > 0.90 &&
         ss_h$mPTP_open < 0.05)
  cat(sprintf("    CI_activity  = %.4f (expect > 0.90)\n", ss_h$CI_activity))
  cat(sprintf("    SC_integrity = %.4f (expect > 0.95)\n", ss_h$SC_integrity))
  cat(sprintf("    delta_psi_m  = %.4f (expect > 0.90)\n", ss_h$delta_psi_m))
  cat(sprintf("    ATP          = %.4f (expect > 0.90)\n", ss_h$ATP))
  cat(sprintf("    mPTP_open    = %.4f (expect < 0.05)\n", ss_h$mPTP_open))
  cat(sprintf("    V1 Healthy SS:              %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation ---
  ## CL_ratio=0.807 from M1, mROS=0.176 from M4, aSyn=0.50
  cat("V2: Disease perturbation...\n")
  res_d <- run_M3(CL_ratio = 0.807, mROS = 0.176, aSyn_oligo = 0.50,
                  duration_h = t_verify)
  ss_d  <- res_d[nrow(res_d), ]

  v2 <- (ss_d$CI_activity < ss_h$CI_activity &&
         ss_d$SC_integrity < ss_h$SC_integrity &&
         ss_d$delta_psi_m < ss_h$delta_psi_m &&
         ss_d$ATP < ss_h$ATP &&
         ss_d$mPTP_open > ss_h$mPTP_open)
  cat(sprintf("    CI:   %.4f -> %.4f (target ~0.49)\n",
              ss_h$CI_activity, ss_d$CI_activity))
  cat(sprintf("    SC:   %.4f -> %.4f\n",
              ss_h$SC_integrity, ss_d$SC_integrity))
  cat(sprintf("    dpsi: %.4f -> %.4f\n",
              ss_h$delta_psi_m, ss_d$delta_psi_m))
  cat(sprintf("    ATP:  %.4f -> %.4f\n",
              ss_h$ATP, ss_d$ATP))
  cat(sprintf("    mPTP: %.4f -> %.4f\n",
              ss_h$mPTP_open, ss_d$mPTP_open))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue (CL_ratio improvement from M1) ---
  ## Drug acts through M1: CL_ratio improves 0.807 -> 0.911
  ## This only affects SC directly; CI unchanged (mROS/aSyn same)
  cat("V3: Drug rescue (CL_ratio 0.807 -> 0.911 from M1)...\n")
  res_drug <- run_M3(CL_ratio = 0.911, mROS = 0.176, aSyn_oligo = 0.50,
                     duration_h = t_verify)
  ss_rx <- res_drug[nrow(res_drug), ]

  v3 <- (ss_rx$SC_integrity > ss_d$SC_integrity &&
         ss_rx$delta_psi_m > ss_d$delta_psi_m &&
         ss_rx$ATP > ss_d$ATP)
  cat(sprintf("    SC:   disease=%.4f -> drug=%.4f\n",
              ss_d$SC_integrity, ss_rx$SC_integrity))
  cat(sprintf("    dpsi: disease=%.4f -> drug=%.4f\n",
              ss_d$delta_psi_m, ss_rx$delta_psi_m))
  cat(sprintf("    ATP:  disease=%.4f -> drug=%.4f\n",
              ss_d$ATP, ss_rx$ATP))
  cat(sprintf("    CI:   disease=%.4f -> drug=%.4f (unchanged, expected)\n",
              ss_d$CI_activity, ss_rx$CI_activity))
  cat(sprintf("    V3 Drug rescue (via SC):    %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero damage -> full ETC function ---
  cat("V4: Zero damage -> full ETC function...\n")
  res_zero <- run_M3(CL_ratio = 1.0, mROS = 0.0, aSyn_oligo = 0.0,
                     duration_h = t_verify,
                     y0 = c(CI_activity = 1.0, SC_integrity = 1.0,
                            delta_psi_m = 0.97, ATP = 0.95, mPTP_open = 0.0))
  ss_z <- res_zero[nrow(res_zero), ]

  v4 <- (abs(ss_z$CI_activity - 1.0) < 1e-6 &&
         abs(ss_z$mPTP_open) < 1e-10)
  cat(sprintf("    CI_activity = %.6f (expect 1.0)\n", ss_z$CI_activity))
  cat(sprintf("    mPTP_open   = %.2e (expect 0)\n", ss_z$mPTP_open))
  cat(sprintf("    SC = %.4f  dpsi = %.4f  ATP = %.4f\n",
              ss_z$SC_integrity, ss_z$delta_psi_m, ss_z$ATP))
  cat(sprintf("    V4 Zero damage:             %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness [0,1] ---
  cat("V5: Non-negativity and boundedness...\n")
  all_res <- rbind(res_h, res_d, res_drug)
  v5 <- (all(all_res$CI_activity >= -1e-10) &&
         all(all_res$CI_activity <= 1 + 1e-10) &&
         all(all_res$SC_integrity >= -1e-10) &&
         all(all_res$SC_integrity <= 1 + 1e-10) &&
         all(all_res$delta_psi_m >= -1e-10) &&
         all(all_res$delta_psi_m <= 1 + 1e-10) &&
         all(all_res$ATP >= -1e-10) &&
         all(all_res$ATP <= 1 + 1e-10) &&
         all(all_res$mPTP_open >= -1e-10) &&
         all(all_res$mPTP_open <= 1 + 1e-10))
  cat(sprintf("    CI range:   [%.4f, %.4f]\n",
              min(all_res$CI_activity), max(all_res$CI_activity)))
  cat(sprintf("    SC range:   [%.4f, %.4f]\n",
              min(all_res$SC_integrity), max(all_res$SC_integrity)))
  cat(sprintf("    dpsi range: [%.4f, %.4f]\n",
              min(all_res$delta_psi_m), max(all_res$delta_psi_m)))
  cat(sprintf("    ATP range:  [%.4f, %.4f]\n",
              min(all_res$ATP), max(all_res$ATP)))
  cat(sprintf("    mPTP range: [%.4e, %.4f]\n",
              min(all_res$mPTP_open), max(all_res$mPTP_open)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical SS matches analytical SS (<1%) ---
  cat("V6: Numerical vs analytical steady state...\n")
  ss_an_h <- M3_analytical_ss(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05)
  ss_an_d <- M3_analytical_ss(CL_ratio = 0.807, mROS = 0.176, aSyn_oligo = 0.50)

  err_CI_h  <- abs(ss_h$CI_activity - ss_an_h["CI_activity"]) /
               ss_an_h["CI_activity"] * 100
  err_CI_d  <- abs(ss_d$CI_activity - ss_an_d["CI_activity"]) /
               ss_an_d["CI_activity"] * 100
  err_psi_h <- abs(ss_h$delta_psi_m - ss_an_h["delta_psi_m"]) /
               ss_an_h["delta_psi_m"] * 100
  err_psi_d <- abs(ss_d$delta_psi_m - ss_an_d["delta_psi_m"]) /
               ss_an_d["delta_psi_m"] * 100

  v6 <- (err_CI_h < 1.0 && err_CI_d < 1.0 &&
         err_psi_h < 1.0 && err_psi_d < 1.0)
  cat(sprintf("    Healthy CI:   num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$CI_activity, ss_an_h["CI_activity"], err_CI_h))
  cat(sprintf("    Disease CI:   num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$CI_activity, ss_an_d["CI_activity"], err_CI_d))
  cat(sprintf("    Healthy dpsi: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$delta_psi_m, ss_an_h["delta_psi_m"], err_psi_h))
  cat(sprintf("    Disease dpsi: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$delta_psi_m, ss_an_d["delta_psi_m"], err_psi_d))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M3 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M3_scenarios <- function() {

  ## Scenario definitions
  ## Drug effect is INDIRECT via M1: CL_ratio improves -> SC improves
  ## mROS and aSyn held at disease level (not yet feedback-coupled)
  scenarios <- list(
    healthy  = list(CL = 0.948, mROS = 0.10,  aSyn = 0.05,
                    label = "Healthy"),
    disease  = list(CL = 0.807, mROS = 0.176, aSyn = 0.50,
                    label = "Disease (no drug)"),
    drug_low = list(CL = 0.827, mROS = 0.176, aSyn = 0.50,
                    label = "Disease + SBT-272 low"),
    drug_hi  = list(CL = 0.911, mROS = 0.176, aSyn = 0.50,
                    label = "Disease + SBT-272 high")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M3(CL_ratio = s$CL, mROS = s$mROS, aSyn_oligo = s$aSyn,
                  duration_h = 200)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))


  ## --- Plot 1: Time course of all states ---
  plot_vars <- c("CI_activity", "SC_integrity", "delta_psi_m",
                 "ATP", "mPTP_open")

  long_data <- all_data %>%
    select(time, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("Complex I activity", "Supercomplex integrity",
               "Membrane potential", "ATP", "mPTP open fraction"))

  p1 <- ggplot(long_data, aes(x = time, y = value, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (no drug)" = "#d62728",
                                  "Disease + SBT-272 low" = "#ff7f0e",
                                  "Disease + SBT-272 high" = "#1f77b4")) +
    labs(
      title = "Module M3: ETC Bioenergetics",
      subtitle = paste("Upstream: CL_ratio from M1, mROS from M4,",
                        "aSyn scenario parameter"),
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

  ggsave("M3_ETC_time_course.png", p1,
         width = 12, height = 10, dpi = 300, bg = "white")
  cat("Saved: M3_ETC_time_course.png\n")


  ## --- Plot 2: mROS sweep -> all SS outputs ---
  cat("Running mROS sweep (analytical SS)...\n")
  mROS_sweep <- seq(0.05, 0.30, by = 0.005)
  sweep_data <- data.frame()

  for (mr in mROS_sweep) {
    ss_vals <- M3_analytical_ss(CL_ratio = 0.807, mROS = mr, aSyn_oligo = 0.50)
    sweep_data <- rbind(sweep_data, data.frame(
      mROS         = mr,
      CI_activity  = ss_vals["CI_activity"],
      SC_integrity = ss_vals["SC_integrity"],
      delta_psi_m  = ss_vals["delta_psi_m"],
      ATP          = ss_vals["ATP"],
      mPTP_open    = ss_vals["mPTP_open"]
    ))
  }

  sweep_long <- sweep_data %>%
    pivot_longer(cols = -mROS, names_to = "variable", values_to = "value")

  sweep_long$variable <- factor(sweep_long$variable,
    levels = c("CI_activity", "SC_integrity", "delta_psi_m", "ATP", "mPTP_open"),
    labels = c("CI activity", "SC integrity", "Membrane potential",
               "ATP", "mPTP open"))

  p2 <- ggplot(sweep_long, aes(x = mROS, y = value)) +
    geom_line(linewidth = 1, color = "#d62728") +
    facet_wrap(~ variable, scales = "free_y", ncol = 2) +
    geom_vline(xintercept = 0.10, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0.176, linetype = "dotted", color = "#d62728") +
    labs(
      title = "M3: Steady-State ETC Response to mROS",
      subtitle = "Disease background (CL_ratio = 0.807, aSyn = 0.50)",
      x = "mROS (normalised, from M4)",
      y = "Value",
      caption = "Dashed: healthy mROS | Dotted: disease mROS"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave("M3_mROS_sweep.png", p2,
         width = 10, height = 8, dpi = 300, bg = "white")
  cat("Saved: M3_mROS_sweep.png\n")


  ## --- Plot 3: Sensitivity tornado on ATP_ss ---
  cat("Running sensitivity analysis...\n")
  base_ss <- M3_analytical_ss(CL_ratio = 0.807, mROS = 0.176, aSyn_oligo = 0.50)
  base_ATP <- base_ss["ATP"]

  p_base <- M3_parameters(CL_ratio = 0.807, mROS = 0.176, aSyn_oligo = 0.50)

  param_names <- c("k_CI_repair", "k_CI_ROS", "k_CI_aSyn",
                    "k_SC_form", "k_SC_dissoc",
                    "k_resp", "k_ATPase", "k_leak", "k_mPTP_depol",
                    "k_ATP_syn", "k_ATP_consume", "k_ATP_mPTP",
                    "k_mPTP_open", "k_mPTP_close", "K_mPTP")

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      ## Recompute SS with modified parameter
      mod_CI <- p_mod["k_CI_repair"] /
                (p_mod["k_CI_repair"] + p_mod["k_CI_ROS"] * 0.176 +
                 p_mod["k_CI_aSyn"] * 0.50^2)
      mod_SC <- p_mod["k_SC_form"] * 0.807 /
                (p_mod["k_SC_form"] * 0.807 + p_mod["k_SC_dissoc"] * 0.193)

      mod_A <- p_mod["k_resp"] * mod_CI * mod_SC
      mod_B <- p_mod["k_ATPase"] + p_mod["k_leak"]
      mod_H <- 0.176^p_mod["n_Hill"] /
               (p_mod["K_mPTP"]^p_mod["n_Hill"] + 0.176^p_mod["n_Hill"])

      if (mod_H < 1e-15) {
        mod_mPTP <- 0
      } else {
        f_mod <- function(m) {
          dpsi <- as.numeric(mod_A / (mod_A + mod_B + p_mod["k_mPTP_depol"]*m))
          as.numeric(p_mod["k_mPTP_open"] * mod_H * (1-m) -
                     p_mod["k_mPTP_close"] * m * dpsi)
        }
        sol <- tryCatch(
          uniroot(f_mod, interval = c(0, 1 - 1e-10), tol = 1e-12),
          error = function(e) list(root = NA))
        mod_mPTP <- sol$root
      }

      if (is.na(mod_mPTP)) {
        pct_change <- NA
      } else {
        mod_dpsi <- as.numeric(mod_A / (mod_A + mod_B +
                    p_mod["k_mPTP_depol"] * mod_mPTP))
        mod_P <- as.numeric(p_mod["k_ATP_syn"] * mod_dpsi * mod_CI * mod_SC)
        mod_Q <- as.numeric(p_mod["k_ATP_consume"] +
                 p_mod["k_ATP_mPTP"] * mod_mPTP)
        mod_ATP <- mod_P / (mod_P + mod_Q)
        pct_change <- as.numeric((mod_ATP - base_ATP) / base_ATP * 100)
      }

      sens_data <- rbind(sens_data, data.frame(
        parameter  = pname,
        direction  = ifelse(direction > 0, "+50%", "-50%"),
        pct_change = pct_change
      ))
    }
  }

  tornado <- sens_data %>%
    filter(!is.na(pct_change)) %>%
    group_by(parameter) %>%
    summarise(
      low  = min(pct_change),
      high = max(pct_change),
      span = abs(max(pct_change) - min(pct_change)),
      .groups = "drop"
    ) %>%
    arrange(span)

  tornado$parameter <- factor(tornado$parameter, levels = tornado$parameter)

  p3 <- ggplot(tornado, aes(y = parameter)) +
    geom_segment(aes(x = low, xend = high, yend = parameter),
                 linewidth = 6, color = "#5cb85c", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M3: Parameter Sensitivity on Disease ATP",
      subtitle = "OAT +/-50% perturbation (CL = 0.807, mROS = 0.176, aSyn = 0.50)",
      x = "Change in ATP_ss (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M3_sensitivity.png", p3,
         width = 10, height = 7, dpi = 300, bg = "white")
  cat("Saved: M3_sensitivity.png\n\n")


  return(list(data = all_data, sweep = sweep_data, sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

all_passed <- verify_M3()

if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Running scenarios anyway for diagnostic purposes...\n\n")
}

results <- run_all_M3_scenarios()

## Summary table
cat("Analytical SS for key scenarios:\n")
for (label in c("Healthy", "Disease", "Drug_high")) {
  if (label == "Healthy") {
    ss <- M3_analytical_ss(CL_ratio = 0.948, mROS = 0.10, aSyn_oligo = 0.05)
  } else if (label == "Disease") {
    ss <- M3_analytical_ss(CL_ratio = 0.807, mROS = 0.176, aSyn_oligo = 0.50)
  } else {
    ss <- M3_analytical_ss(CL_ratio = 0.911, mROS = 0.176, aSyn_oligo = 0.50)
  }
  cat(sprintf("  %-10s CI=%.3f SC=%.3f dpsi=%.3f ATP=%.3f mPTP=%.4f\n",
              label, ss["CI_activity"], ss["SC_integrity"],
              ss["delta_psi_m"], ss["ATP"], ss["mPTP_open"]))
}

cat("\nSensitivity analysis (top parameters):\n")
print(results$sensitivity[order(-abs(results$sensitivity$pct_change)), ],
      row.names = FALSE)

cat("\n========================================================\n")
cat("  MODULE M3 COMPLETE\n")
cat("  Outputs: CI_activity, SC_integrity -> M4 (ROS)\n")
cat("           delta_psi_m, mPTP_open    -> M5 (mitophagy)\n")
cat("           ATP                       -> M2 (clearance)\n")
cat("  Upstream: CL_ratio <- M1, mROS <- M4, aSyn <- scenario\n")
cat("========================================================\n")
