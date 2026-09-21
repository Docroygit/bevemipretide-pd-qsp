##=============================================================================
## MODULE M6: NEUROINFLAMMATION
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
## Version : 1.0
##
## State variables (3):
##   MG_active — Microglial activation (M1 pro-inflammatory) [0,1]
##   TNF       — TNF-alpha level (normalised) [0,1]
##   IL1b      — IL-1-beta level (normalised) [0,1]
##
## Upstream coupling:
##   aSyn_ext ← M2 (extracellular alpha-synuclein triggers TLR2/TLR4)
##
## Key outputs:
##   MG_active → M2 (impairs aSyn clearance via k_impair)
##   TNF, IL1b → M7 (inflammatory contribution to neuronal damage)
##
## Biology:
##   Extracellular aSyn (released from damaged DA neurons) activates
##   microglia through TLR2/TLR4 signaling (Ivanova 2024). Activated
##   microglia release pro-inflammatory cytokines (TNF-alpha, IL-1-beta).
##   These cytokines sustain microglial activation (positive feedback),
##   impair aSyn clearance (M2 coupling), and contribute to neuronal
##   damage (M7 coupling). This creates the third feedback loop:
##   aSyn -> neuroinflammation -> impaired clearance -> more aSyn.
##
##   SBT-272 does not directly target microglia; its effect is INDIRECT
##   via CL rescue -> better ETC -> less ROS -> less aSyn -> less
##   aSyn_ext -> reduced inflammation.
##
## Standalone calibration note:
##   M2 standalone disease produces aSyn_ext ~ 0.005, which drives
##   modest MG activation (~0.18). In the integrated model, the
##   aSyn->inflammation->impaired clearance feedback amplifies both
##   aSyn and MG to higher levels. At projected integrated aSyn_ext
##   ~ 0.03, MG reaches ~ 0.45, producing 22% clearance impairment
##   in M2 — consistent with the Ivanova 2024 TLR2 KO validation
##   target (23-45% aSyn reduction when inflammation is blocked).
##
## Data sources:
##   Ivanova 2024: TLR2 KO = 23-45% aSyn reduction in mouse model
##   Gerhard 2006 (doi:10.1016/j.nbd.2005.08.002): PET evidence of
##     microglial activation in PD (11C-PK11195 binding)
##   Codolo 2013 (doi:10.1371/journal.pone.0055تو07): aSyn activates
##     NLRP3 inflammasome -> IL-1-beta release
##   Tracey 1987: TNF-alpha half-life ~ 1-7h in tissue
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M6: NEUROINFLAMMATION\n")
cat("  Version 1.0 | Upstream: M2 (aSyn_ext)\n")
cat("  Output: MG_active -> M2; TNF, IL1b -> M7\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M6_parameters <- function(aSyn_ext = 0.003) {

  p <- c(
    ## --- Microglial activation ---
    k_act   = 20.0,       # aSyn_ext-driven activation via TLR2/TLR4 (h^-1)
                           # ESTIMATED: tuned for MG_healthy ~ 0.08-0.14;
                           # high coefficient because aSyn_ext is small;
                           # Ivanova 2024 TLR2 signaling
    k_auto  = 0.50,       # Auto-amplification via cytokines (h^-1)
                           # ASSUMED: TNF + IL1b sustain M1 microglial
                           # phenotype; creates positive feedback that
                           # amplifies disease signal during integration;
                           # subcritical (G < k_deact) to prevent bistability
    k_deact = 1.0,        # Deactivation / resolution (h^-1)
                           # ASSUMED: anti-inflammatory pathways (IL-10,
                           # TGF-beta, CX3CR1) resolve activation;
                           # t_resolution ~ 1h timescale

    ## --- TNF-alpha dynamics ---
    k_TNF_rel = 0.20,     # TNF release by activated microglia (h^-1)
                           # ESTIMATED: calibrated for TNF_healthy ~ 0.10-0.25
    k_TNF_deg = 0.10,     # TNF degradation (h^-1)
                           # LITERATURE: TNF-alpha half-life ~ 1-7h in brain
                           # tissue; k_deg = ln(2)/7 ~ 0.10; Tracey 1987

    ## --- IL-1-beta dynamics ---
    k_IL1b_rel = 0.15,    # IL-1-beta release via NLRP3 inflammasome (h^-1)
                           # ESTIMATED: NLRP3 activated by aSyn via MG;
                           # lower than TNF (NLRP3 activation is slower);
                           # Codolo 2013
    k_IL1b_deg = 0.10,    # IL-1-beta degradation (h^-1)
                           # LITERATURE: similar half-life to TNF-alpha

    ## --- Upstream coupling (scenario parameter) ---
    aSyn_ext = aSyn_ext
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M6_init <- function() {
  c(
    MG_active = 0.05,     # Near resting state
    TNF       = 0.10,     # Low basal
    IL1b      = 0.08      # Low basal
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M6_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states ---
    MG_active <- max(0, min(1, MG_active))
    TNF       <- max(0, min(1, TNF))
    IL1b      <- max(0, min(1, IL1b))

    ## --- ODE 1: Microglial activation ---
    ## TLR2/TLR4-driven activation by extracellular aSyn
    ## Auto-amplification: cytokines sustain M1 phenotype
    ## Resolution: return to resting/M2 phenotype
    cytokine_signal <- (TNF + IL1b) / 2

    dMG <- (k_act * aSyn_ext * (1 - MG_active)
            + k_auto * cytokine_signal * (1 - MG_active)
            - k_deact * MG_active)

    ## --- ODE 2: TNF-alpha ---
    ## Released by activated microglia; saturating production
    ## ensures TNF stays bounded [0,1]
    dTNF <- k_TNF_rel * MG_active * (1 - TNF) - k_TNF_deg * TNF

    ## --- ODE 3: IL-1-beta ---
    ## Released via NLRP3 inflammasome in activated microglia;
    ## saturating production for boundedness
    dIL1b <- k_IL1b_rel * MG_active * (1 - IL1b) - k_IL1b_deg * IL1b

    ## --- Derived quantities ---
    clearance_impairment <- 0.50 * MG_active  # k_impair from M2

    list(
      c(dMG, dTNF, dIL1b),
      cytokine_signal       = cytokine_signal,
      clearance_impairment  = clearance_impairment
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M6 <- function(aSyn_ext = 0.003, duration_h = 200, dt = 0.1,
                   y0 = NULL) {

  parms <- M6_parameters(aSyn_ext)
  if (is.null(y0)) y0 <- M6_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M6_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$aSyn_ext <- aSyn_ext

  ss <- result[nrow(result), ]
  cat(sprintf(
    "M6 | aSyn_ext=%.4f | MG=%.4f TNF=%.4f IL1b=%.4f | clearance impair=%.1f%%\n",
    aSyn_ext, ss$MG_active, ss$TNF, ss$IL1b,
    0.50 * ss$MG_active * 100))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE (semi-analytical via uniroot)
##-----------------------------------------------------------------------------

M6_analytical_ss <- function(aSyn_ext = 0.003) {

  p <- M6_parameters(aSyn_ext)

  ## Special case: zero input
  if (aSyn_ext <= 0) {
    return(c(MG_active = 0, TNF = 0, IL1b = 0))
  }

  ## For given MG, TNF and IL1b are determined (bounded production):
  ## TNF(MG) = k_TNF_rel * MG / (k_TNF_rel * MG + k_TNF_deg)
  ## IL1b(MG) = k_IL1b_rel * MG / (k_IL1b_rel * MG + k_IL1b_deg)
  ##
  ## MG equation at SS:
  ## k_act * aSyn_ext * (1-MG) + k_auto * C(MG) * (1-MG) - k_deact * MG = 0
  ## where C(MG) = (TNF(MG) + IL1b(MG)) / 2

  root_fn <- function(MG) {
    TNF_val  <- as.numeric(p["k_TNF_rel"]) * MG /
                (as.numeric(p["k_TNF_rel"]) * MG + as.numeric(p["k_TNF_deg"]))
    IL1b_val <- as.numeric(p["k_IL1b_rel"]) * MG /
                (as.numeric(p["k_IL1b_rel"]) * MG + as.numeric(p["k_IL1b_deg"]))
    C_val    <- (TNF_val + IL1b_val) / 2

    (as.numeric(p["k_act"]) * aSyn_ext * (1 - MG)
     + as.numeric(p["k_auto"]) * C_val * (1 - MG)
     - as.numeric(p["k_deact"]) * MG)
  }

  MG_ss <- uniroot(root_fn, interval = c(1e-12, 1 - 1e-12),
                    tol = 1e-12)$root

  TNF_ss  <- as.numeric(p["k_TNF_rel"]) * MG_ss /
             (as.numeric(p["k_TNF_rel"]) * MG_ss + as.numeric(p["k_TNF_deg"]))
  IL1b_ss <- as.numeric(p["k_IL1b_rel"]) * MG_ss /
             (as.numeric(p["k_IL1b_rel"]) * MG_ss + as.numeric(p["k_IL1b_deg"]))

  return(c(MG_active = MG_ss, TNF = TNF_ss, IL1b = IL1b_ss))
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M6 <- function() {

  cat("=== MODULE M6 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6


  ## --- V1: Healthy steady state ---
  ## aSyn_ext from M2 healthy (~0.003)
  ## MG low, cytokines at basal levels
  cat("V1: Healthy steady state...\n")
  res_h <- run_M6(aSyn_ext = 0.003, duration_h = 200)
  ss_h  <- res_h[nrow(res_h), ]

  v1 <- (ss_h$MG_active > 0.05 && ss_h$MG_active < 0.25 &&
         ss_h$TNF > 0.05 && ss_h$TNF < 0.30 &&
         ss_h$IL1b > 0.05 && ss_h$IL1b < 0.25)
  cat(sprintf("    MG_active = %.4f (expect 0.05-0.25)\n", ss_h$MG_active))
  cat(sprintf("    TNF       = %.4f (expect 0.05-0.30)\n", ss_h$TNF))
  cat(sprintf("    IL1b      = %.4f (expect 0.05-0.25)\n", ss_h$IL1b))
  cat(sprintf("    V1 Healthy SS:              %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation ---
  ## aSyn_ext from M2 disease (~0.005)
  ## Increased MG activation, cytokine elevation
  cat("V2: Disease perturbation...\n")
  res_d <- run_M6(aSyn_ext = 0.005, duration_h = 200)
  ss_d  <- res_d[nrow(res_d), ]

  v2 <- (ss_d$MG_active > ss_h$MG_active &&
         ss_d$TNF > ss_h$TNF &&
         ss_d$IL1b > ss_h$IL1b)
  cat(sprintf("    MG: %.4f -> %.4f (ratio: %.2fx)\n",
              ss_h$MG_active, ss_d$MG_active,
              ss_d$MG_active / ss_h$MG_active))
  cat(sprintf("    TNF: %.4f -> %.4f (ratio: %.2fx)\n",
              ss_h$TNF, ss_d$TNF, ss_d$TNF / ss_h$TNF))
  cat(sprintf("    IL1b: %.4f -> %.4f (ratio: %.2fx)\n",
              ss_h$IL1b, ss_d$IL1b, ss_d$IL1b / ss_h$IL1b))
  cat(sprintf("    Clearance impairment: %.1f%% -> %.1f%%\n",
              0.50 * ss_h$MG_active * 100,
              0.50 * ss_d$MG_active * 100))
  cat(sprintf("    Note: standalone effect is modest; aSyn-inflammation\n"))
  cat(sprintf("          feedback loop amplifies during integration\n"))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue ---
  ## Drug reduces aSyn_ext indirectly via CL -> ETC -> ROS -> aSyn
  cat("V3: Drug rescue (reduced aSyn_ext via upstream CL rescue)...\n")
  res_rx <- run_M6(aSyn_ext = 0.004, duration_h = 200)
  ss_rx  <- res_rx[nrow(res_rx), ]

  v3 <- (ss_rx$MG_active < ss_d$MG_active &&
         ss_rx$TNF < ss_d$TNF &&
         ss_rx$IL1b < ss_d$IL1b)
  cat(sprintf("    MG: disease=%.4f -> drug=%.4f\n",
              ss_d$MG_active, ss_rx$MG_active))
  cat(sprintf("    TNF: disease=%.4f -> drug=%.4f\n",
              ss_d$TNF, ss_rx$TNF))
  cat(sprintf("    IL1b: disease=%.4f -> drug=%.4f\n",
              ss_d$IL1b, ss_rx$IL1b))
  cat(sprintf("    V3 Drug rescue:             %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero aSyn_ext → no activation ---
  cat("V4: Zero aSyn_ext -> zero inflammation...\n")
  res_z <- run_M6(aSyn_ext = 0.0, duration_h = 200,
                  y0 = c(MG_active = 0.0, TNF = 0.0, IL1b = 0.0))
  ss_z  <- res_z[nrow(res_z), ]

  v4 <- (abs(ss_z$MG_active) < 1e-10 &&
         abs(ss_z$TNF) < 1e-10 &&
         abs(ss_z$IL1b) < 1e-10)
  cat(sprintf("    MG_active = %.2e (expect 0)\n", ss_z$MG_active))
  cat(sprintf("    TNF       = %.2e (expect 0)\n", ss_z$TNF))
  cat(sprintf("    IL1b      = %.2e (expect 0)\n", ss_z$IL1b))
  cat(sprintf("    V4 Zero aSyn_ext:           %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness [0,1] ---
  cat("V5: Non-negativity and boundedness...\n")
  all_res <- rbind(res_h, res_d, res_rx)
  v5 <- (all(all_res$MG_active >= -1e-10) &&
         all(all_res$MG_active <= 1 + 1e-10) &&
         all(all_res$TNF >= -1e-10) &&
         all(all_res$TNF <= 1 + 1e-10) &&
         all(all_res$IL1b >= -1e-10) &&
         all(all_res$IL1b <= 1 + 1e-10))
  cat(sprintf("    MG range:   [%.4f, %.4f]\n",
              min(all_res$MG_active), max(all_res$MG_active)))
  cat(sprintf("    TNF range:  [%.4f, %.4f]\n",
              min(all_res$TNF), max(all_res$TNF)))
  cat(sprintf("    IL1b range: [%.4f, %.4f]\n",
              min(all_res$IL1b), max(all_res$IL1b)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical vs analytical match ---
  ## Semi-analytical: uniroot for MG, then closed-form TNF/IL1b
  cat("V6: Numerical vs analytical steady state...\n")
  ss_an_h <- M6_analytical_ss(0.003)
  ss_an_d <- M6_analytical_ss(0.005)

  err_MG_h   <- abs(ss_h$MG_active - ss_an_h["MG_active"]) /
                ss_an_h["MG_active"] * 100
  err_MG_d   <- abs(ss_d$MG_active - ss_an_d["MG_active"]) /
                ss_an_d["MG_active"] * 100
  err_TNF_h  <- abs(ss_h$TNF - ss_an_h["TNF"]) / ss_an_h["TNF"] * 100
  err_TNF_d  <- abs(ss_d$TNF - ss_an_d["TNF"]) / ss_an_d["TNF"] * 100
  err_IL1b_h <- abs(ss_h$IL1b - ss_an_h["IL1b"]) / ss_an_h["IL1b"] * 100
  err_IL1b_d <- abs(ss_d$IL1b - ss_an_d["IL1b"]) / ss_an_d["IL1b"] * 100

  v6 <- (err_MG_h < 1.0 && err_MG_d < 1.0 &&
         err_TNF_h < 1.0 && err_TNF_d < 1.0 &&
         err_IL1b_h < 1.0 && err_IL1b_d < 1.0)
  cat(sprintf("    Healthy MG:   num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$MG_active, ss_an_h["MG_active"], err_MG_h))
  cat(sprintf("    Disease MG:   num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$MG_active, ss_an_d["MG_active"], err_MG_d))
  cat(sprintf("    Healthy TNF:  num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$TNF, ss_an_h["TNF"], err_TNF_h))
  cat(sprintf("    Disease TNF:  num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$TNF, ss_an_d["TNF"], err_TNF_d))
  cat(sprintf("    Healthy IL1b: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$IL1b, ss_an_h["IL1b"], err_IL1b_h))
  cat(sprintf("    Disease IL1b: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$IL1b, ss_an_d["IL1b"], err_IL1b_d))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M6 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M6_scenarios <- function() {

  ## Scenario definitions
  ## aSyn_ext from M2 analytical SS at each condition
  ## aSyn_ext = 0.05 * aSyn_olig (k_secrete/k_clear_ext)
  scenarios <- list(
    healthy  = list(aSyn_ext = 0.003,  label = "Healthy"),
    disease  = list(aSyn_ext = 0.005,  label = "Disease (no drug)"),
    drug_low = list(aSyn_ext = 0.0045, label = "Disease + SBT-272 low"),
    drug_hi  = list(aSyn_ext = 0.004,  label = "Disease + SBT-272 high")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M6(aSyn_ext = s$aSyn_ext, duration_h = 200)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))


  ## --- Plot 1: Time course ---
  plot_vars <- c("MG_active", "TNF", "IL1b")

  long_data <- all_data %>%
    select(time, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("Microglial activation",
               "TNF-α",
               "IL-1β"))

  p1 <- ggplot(long_data, aes(x = time, y = value, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (no drug)" = "#d62728",
                                  "Disease + SBT-272 low" = "#ff7f0e",
                                  "Disease + SBT-272 high" = "#1f77b4")) +
    labs(
      title = "Module M6: Neuroinflammation",
      subtitle = "Upstream: aSyn_ext from M2 | Auto-amplification via TNF + IL-1β",
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

  ggsave("M6_time_course.png", p1,
         width = 10, height = 10, dpi = 300, bg = "white")
  cat("Saved: M6_time_course.png\n")


  ## --- Plot 2: aSyn_ext sweep (projected integration range) ---
  cat("Running aSyn_ext sweep (projected integration range)...\n")
  aSyn_sweep <- seq(0.001, 0.05, by = 0.001)
  sweep_data <- data.frame()

  for (ae in aSyn_sweep) {
    ss_vals <- M6_analytical_ss(ae)
    sweep_data <- rbind(sweep_data, data.frame(
      aSyn_ext  = ae,
      MG_active = ss_vals["MG_active"],
      TNF       = ss_vals["TNF"],
      IL1b      = ss_vals["IL1b"],
      clearance_impair_pct = 0.50 * ss_vals["MG_active"] * 100
    ))
  }

  p2 <- ggplot(sweep_data) +
    geom_line(aes(x = aSyn_ext, y = MG_active), linewidth = 1,
              color = "#d62728") +
    geom_line(aes(x = aSyn_ext, y = TNF), linewidth = 1,
              color = "#ff7f0e", linetype = "dashed") +
    geom_line(aes(x = aSyn_ext, y = IL1b), linewidth = 1,
              color = "#9467bd", linetype = "dotdash") +
    geom_vline(xintercept = c(0.003, 0.005), linetype = "dotted",
               color = "gray50") +
    geom_vline(xintercept = 0.03, linetype = "dotted",
               color = "gray50") +
    annotate("text", x = 0.003, y = 0.65, label = "Healthy",
             hjust = -0.1, size = 3, color = "#2ca02c") +
    annotate("text", x = 0.005, y = 0.65, label = "Disease\n(standalone)",
             hjust = -0.1, size = 3, color = "#d62728") +
    annotate("text", x = 0.03, y = 0.65, label = "Projected\nintegrated",
             hjust = -0.1, size = 3, color = "gray40") +
    annotate("text", x = 0.045, y = 0.56,
             label = "MG_active", size = 3.5, color = "#d62728") +
    annotate("text", x = 0.045, y = 0.48,
             label = "TNF", size = 3.5, color = "#ff7f0e") +
    annotate("text", x = 0.045, y = 0.38,
             label = "IL-1β", size = 3.5, color = "#9467bd") +
    labs(
      title = "M6: Neuroinflammation vs Extracellular α-Synuclein",
      subtitle = paste("Shows how vicious cycle amplification of aSyn_ext",
                        "drives stronger neuroinflammation"),
      x = "aSyn_ext (from M2)",
      y = "Steady-state level",
      caption = paste("At standalone disease (aSyn_ext=0.005), MG activation",
                       "is modest (~0.18).",
                       "Projected integration (aSyn_ext~0.03) drives MG~0.45,",
                       "impairing aSyn clearance by ~23%.")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M6_aSyn_sweep.png", p2,
         width = 10, height = 6, dpi = 300, bg = "white")
  cat("Saved: M6_aSyn_sweep.png\n")


  ## --- Plot 3: Sensitivity tornado on MG_active (disease) ---
  cat("Running sensitivity analysis...\n")
  ss_base <- M6_analytical_ss(0.005)
  base_MG <- ss_base["MG_active"]

  param_names <- c("k_act", "k_auto", "k_deact",
                    "k_TNF_rel", "k_TNF_deg", "k_IL1b_rel", "k_IL1b_deg")

  p_base <- M6_parameters(0.005)

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      ## Solve modified SS via uniroot
      root_fn <- function(MG) {
        TNF_v  <- as.numeric(p_mod["k_TNF_rel"]) * MG /
                  (as.numeric(p_mod["k_TNF_rel"]) * MG +
                   as.numeric(p_mod["k_TNF_deg"]))
        IL1b_v <- as.numeric(p_mod["k_IL1b_rel"]) * MG /
                  (as.numeric(p_mod["k_IL1b_rel"]) * MG +
                   as.numeric(p_mod["k_IL1b_deg"]))
        C_v    <- (TNF_v + IL1b_v) / 2

        (as.numeric(p_mod["k_act"]) * 0.005 * (1 - MG)
         + as.numeric(p_mod["k_auto"]) * C_v * (1 - MG)
         - as.numeric(p_mod["k_deact"]) * MG)
      }

      MG_mod <- tryCatch(
        uniroot(root_fn, interval = c(1e-12, 1 - 1e-12),
                tol = 1e-12)$root,
        error = function(e) NA
      )

      if (!is.na(MG_mod)) {
        pct_change <- as.numeric((MG_mod - base_MG) / base_MG * 100)
        sens_data <- rbind(sens_data, data.frame(
          parameter  = pname,
          direction  = ifelse(direction > 0, "+50%", "-50%"),
          MG_active  = MG_mod,
          pct_change = pct_change
        ))
      }
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
                 linewidth = 6, color = "#ff9896", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M6: Parameter Sensitivity on Disease MG_active",
      subtitle = "OAT +/-50% perturbation (aSyn_ext = 0.005)",
      x = "Change in MG_active_ss (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M6_sensitivity.png", p3,
         width = 10, height = 5, dpi = 300, bg = "white")
  cat("Saved: M6_sensitivity.png\n\n")


  return(list(data = all_data, sweep = sweep_data, sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

## Run verification
all_passed <- verify_M6()

## Run scenarios and plots
if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Running scenarios anyway for diagnostics...\n\n")
}

results <- run_all_M6_scenarios()

## Print analytical SS summary
cat("Analytical SS summary:\n")
for (ae in c(0.003, 0.005, 0.0045, 0.004)) {
  ss <- M6_analytical_ss(ae)
  label <- ifelse(ae == 0.003, "Healthy",
           ifelse(ae == 0.005, "Disease",
           ifelse(ae == 0.0045, "Drug_low", "Drug_high")))
  cat(sprintf("  %-10s (aSyn_ext=%.4f): MG=%.4f TNF=%.4f IL1b=%.4f | impair=%.1f%%\n",
              label, ae, ss["MG_active"], ss["TNF"], ss["IL1b"],
              0.50 * ss["MG_active"] * 100))
}

## Projected integration impact
cat("\nProjected integration impact (amplified aSyn_ext):\n")
for (ae in c(0.005, 0.010, 0.020, 0.030, 0.050)) {
  ss <- M6_analytical_ss(ae)
  cat(sprintf("  aSyn_ext=%.3f -> MG=%.3f TNF=%.3f IL1b=%.3f | impair=%.0f%%\n",
              ae, ss["MG_active"], ss["TNF"], ss["IL1b"],
              0.50 * ss["MG_active"] * 100))
}

## Feedback loop analysis
cat("\nFeedback loop #3 analysis (aSyn->inflammation->impaired clearance):\n")
cat("  In M2, clearance = k_clear * ATP * (1 - 0.50 * MG_active)\n")
cat("  At standalone disease MG=0.18: clearance impaired by 9%\n")
cat("  At projected integrated MG=0.45: clearance impaired by 23%\n")
cat("  Ivanova 2024 target: TLR2 KO reduces aSyn by 23-45%\n")
cat("  -> Blocking MG_active (TLR2 KO) in integrated model should match\n")

cat("\nSensitivity analysis results:\n")
print(results$sensitivity[, c("parameter", "direction", "pct_change")],
      row.names = FALSE)

cat("\n========================================================\n")
cat("  MODULE M6 COMPLETE\n")
cat("  Output: MG_active -> M2 (clearance impairment)\n")
cat("          TNF, IL1b -> M7 (inflammatory neuronal damage)\n")
cat("  Feedback loop #3: aSyn -> neuroinflammation -> impaired\n")
cat("                    clearance -> more aSyn (CLOSED)\n")
cat("========================================================\n")
