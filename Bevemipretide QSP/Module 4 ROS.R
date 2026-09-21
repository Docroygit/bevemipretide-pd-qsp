##=============================================================================
## MODULE M4: MITOCHONDRIAL ROS
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
## Version : 1.0
##
## State variables (3):
##   mROS     — Mitochondrial ROS level [0,1] (normalised)
##   SOD2_act — MnSOD (SOD2) active fraction [0,1]
##   GPx4_act — GPx4 active fraction [0,1]
##
## Key output to downstream:
##   mROS → M1 (CL oxidation), M2 (aSyn aggregation),
##          M3 (mPTP opening), M7 (direct neurotoxicity)
##
## Upstream inputs (fixed as scenario parameters until M3 is built):
##   CI_activity  — Complex I activity [0,1], from M3
##   SC_integrity — Supercomplex integrity [0,1], from M3
##
## Biology:
##   Superoxide production from ETC electron leak (CI damage), buffered
##   by SOD2 (matrix O₂⁻ → H₂O₂) and GPx4 (H₂O₂ → H₂O).
##   Antioxidant capacity is exhaustible under sustained oxidative stress.
##   mROS represents the integrated oxidative burden, not a single species.
##
## Validation targets:
##   Healthy: mROS ≈ 0.10  (basal)
##   Disease: mROS ≈ 0.177 (177% of basal, Choi 2022)
##   Kembro 2013: [O₂⁻]m ≈ 6.4 nM, [H₂O₂]m ≈ 82 nM at baseline
##
## Data sources:
##   Cortassa 2004 (Biophys J 84:2734): mitochondrial energetics ODE model,
##     SOD kinetics k1 = 2.4×10⁶ mM⁻¹s⁻¹, shunt fraction ~5% of V_O2
##   Kembro/Cortassa 2013 (Biophys J 104:332): two-compartment ME-R model,
##     SOD2 E_T = 0.3 μM, GPx E_T = 0.1 μM, Φ₁ = 5×10⁻³
##   Choi 2022: H₂O₂ = 177% of basal in PD mitochondria
##   Gao 2017: Complex I activity = 49% of WT in PD model
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M4: MITOCHONDRIAL ROS — Bevemipretide QSP\n")
cat("  AIIMS Bhubaneswar | August 2026\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M4_parameters <- function(CI_activity = 1.0, SC_integrity = 1.0) {

  p <- c(
    ## --- ROS production ---
    k_ROS_basal = 0.122,     # Basal mROS production (h⁻¹)
                              # ESTIMATED: tuned so healthy mROS_ss ≈ 0.10
                              # Kembro: 0.15–2% of total O₂ flux
    k_shunt     = 0.125,     # ETC electron leak from CI damage (h⁻¹)
                              # ESTIMATED: Cortassa shunt concept, normalised
                              # to hit disease mROS ≈ 177% basal (Choi 2022)
                              # with CI_activity = 0.49 (Gao 2017)

    ## --- Antioxidant clearance ---
    k_SOD2      = 1.0,       # SOD2 catalytic rate constant (h⁻¹)
                              # LITERATURE: Cortassa k1=2.4×10⁶ mM⁻¹s⁻¹;
                              # Kembro E_T=0.3 μM; lumped & normalised
    k_GPx4      = 0.5,       # GPx4 catalytic rate constant (h⁻¹)
                              # LITERATURE: Kembro E_T=0.1 μM, Φ₁=5×10⁻³;
                              # lumped & normalised
    k_ROS_other = 0.05,      # Other antioxidant clearance (h⁻¹)
                              # ASSUMED: catalase, thioredoxin, Prx3, etc.

    ## --- SOD2 dynamics ---
    k_SOD2_on   = 0.50,      # SOD2 recovery/reactivation rate (h⁻¹)
                              # ASSUMED: enzyme resynthesis + refolding
    k_SOD2_off  = 1.0,       # SOD2 inactivation by ROS (h⁻¹)
                              # ESTIMATED: H₂O₂ product inhibition
                              # Cortassa: Ki ≈ 0.5 mM

    ## --- GPx4 dynamics ---
    k_GPx4_on   = 0.30,      # GPx4 recovery rate (h⁻¹)
                              # ASSUMED: GSH-dependent recycling
    k_GPx4_off  = 1.5,       # GPx4 inactivation by ROS (h⁻¹)
                              # ESTIMATED: GPx4 more vulnerable than SOD2
                              # (selenocysteine oxidation, ferroptosis link)

    ## --- Upstream coupling (scenario parameters) ---
    CI_activity  = CI_activity,   # From M3: Complex I fractional activity
    SC_integrity = SC_integrity   # From M3: supercomplex assembly [0,1]
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M4_init <- function() {
  c(
    mROS     = 0.10,   # Healthy basal ROS level
    SOD2_act = 0.833,  # Near full activity at healthy mROS
    GPx4_act = 0.667   # Active fraction at healthy mROS
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M4_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states to [0, 1] for numerical safety ---
    mROS     <- max(0, min(1, mROS))
    SOD2_act <- max(0, min(1, SOD2_act))
    GPx4_act <- max(0, min(1, GPx4_act))

    ## --- ODE 1: Mitochondrial ROS ---
    ## Production: basal + electron leak from damaged CI
    ## Clearance: SOD2, GPx4, other antioxidants
    ROS_production <- k_ROS_basal +
                      k_shunt * (1 - CI_activity) * SC_integrity

    ROS_clearance  <- (k_SOD2 * SOD2_act +
                       k_GPx4 * GPx4_act +
                       k_ROS_other) * mROS

    dmROS <- ROS_production - ROS_clearance

    ## --- ODE 2: SOD2 active fraction ---
    ## Recovery towards full activity, inactivation by ROS
    dSOD2_act <- k_SOD2_on * (1 - SOD2_act) -
                 k_SOD2_off * mROS * SOD2_act

    ## --- ODE 3: GPx4 active fraction ---
    ## GSH-dependent recovery, ROS-driven inactivation
    dGPx4_act <- k_GPx4_on * (1 - GPx4_act) -
                 k_GPx4_off * mROS * GPx4_act

    ## --- Derived quantities ---
    total_clearance <- k_SOD2 * SOD2_act + k_GPx4 * GPx4_act + k_ROS_other
    mROS_pct_basal  <- mROS / 0.10 * 100   # percentage of healthy basal

    list(
      c(dmROS, dSOD2_act, dGPx4_act),
      ROS_production  = ROS_production,
      ROS_clearance   = ROS_clearance,
      total_clearance = total_clearance,
      mROS_pct_basal  = mROS_pct_basal
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M4 <- function(CI_activity = 1.0, SC_integrity = 1.0,
                   duration_h = 72, dt = 0.1, y0 = NULL) {

  parms <- M4_parameters(CI_activity, SC_integrity)
  if (is.null(y0)) y0 <- M4_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M4_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$CI_activity  <- CI_activity
  result$SC_integrity <- SC_integrity

  ## Report SS values (last time point)
  ss <- result[nrow(result), ]
  cat(sprintf("M4 | CI=%.2f SC=%.2f | SS: mROS=%.4f (%.0f%% basal) SOD2=%.3f GPx4=%.3f\n",
              CI_activity, SC_integrity,
              ss$mROS, ss$mROS / 0.10 * 100,
              ss$SOD2_act, ss$GPx4_act))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE
##-----------------------------------------------------------------------------

M4_analytical_ss <- function(CI_activity = 1.0, SC_integrity = 1.0) {
  ## Solve the nonlinear SS equation numerically
  p <- M4_parameters(CI_activity, SC_integrity)

  P <- p["k_ROS_basal"] + p["k_shunt"] * (1 - CI_activity) * SC_integrity

  f <- function(m) {
    s2 <- p["k_SOD2_on"] / (p["k_SOD2_on"] + p["k_SOD2_off"] * m)
    g4 <- p["k_GPx4_on"] / (p["k_GPx4_on"] + p["k_GPx4_off"] * m)
    clearance <- (p["k_SOD2"] * s2 + p["k_GPx4"] * g4 + p["k_ROS_other"]) * m
    as.numeric(P - clearance)
  }

  sol <- uniroot(f, interval = c(1e-8, 1.0), tol = 1e-12)

  mROS_ss <- sol$root
  SOD2_ss <- as.numeric(p["k_SOD2_on"] / (p["k_SOD2_on"] + p["k_SOD2_off"] * mROS_ss))
  GPx4_ss <- as.numeric(p["k_GPx4_on"] / (p["k_GPx4_on"] + p["k_GPx4_off"] * mROS_ss))

  return(c(mROS = mROS_ss, SOD2_act = SOD2_ss, GPx4_act = GPx4_ss))
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M4 <- function() {

  cat("=== MODULE M4 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6


  ## --- V1: Healthy steady state in expected range ---
  cat("V1: Healthy steady state (mROS ≈ 0.10)...\n")
  res_h <- run_M4(CI_activity = 1.0, SC_integrity = 1.0, duration_h = 72)
  ss_h  <- res_h[nrow(res_h), ]

  v1 <- (ss_h$mROS > 0.05 && ss_h$mROS < 0.15 &&
         ss_h$SOD2_act > 0.70 && ss_h$GPx4_act > 0.50)
  cat(sprintf("    mROS = %.4f (target: 0.05–0.15)\n", ss_h$mROS))
  cat(sprintf("    SOD2_act = %.3f (expect > 0.70)\n", ss_h$SOD2_act))
  cat(sprintf("    GPx4_act = %.3f (expect > 0.50)\n", ss_h$GPx4_act))
  cat(sprintf("    V1 Healthy SS:              %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation increases mROS ---
  ## Gao 2017: CI = 49% WT; Choi 2022: H₂O₂ = 177% basal
  cat("V2: Disease perturbation (CI = 0.49)...\n")
  res_d <- run_M4(CI_activity = 0.49, SC_integrity = 1.0, duration_h = 72)
  ss_d  <- res_d[nrow(res_d), ]

  fold_increase <- ss_d$mROS / ss_h$mROS
  v2 <- (fold_increase > 1.5 && fold_increase < 2.5 &&
         ss_d$SOD2_act < ss_h$SOD2_act &&
         ss_d$GPx4_act < ss_h$GPx4_act)
  cat(sprintf("    Disease mROS = %.4f (%.0f%% of basal)\n",
              ss_d$mROS, fold_increase * 100))
  cat(sprintf("    Target: 150–250%% of basal (Choi 2022: 177%%)\n"))
  cat(sprintf("    SOD2_act = %.3f (was %.3f, depleted: %s)\n",
              ss_d$SOD2_act, ss_h$SOD2_act,
              ifelse(ss_d$SOD2_act < ss_h$SOD2_act, "YES", "NO")))
  cat(sprintf("    GPx4_act = %.3f (was %.3f, depleted: %s)\n",
              ss_d$GPx4_act, ss_h$GPx4_act,
              ifelse(ss_d$GPx4_act < ss_h$GPx4_act, "YES", "NO")))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue reduces mROS below disease ---
  ## Drug acts indirectly: CL stabilisation → better CI/SC → less leak
  cat("V3: Drug rescue (CI improved to 0.80)...\n")
  res_drug <- run_M4(CI_activity = 0.80, SC_integrity = 1.0, duration_h = 72)
  ss_dr <- res_drug[nrow(res_drug), ]

  v3 <- (ss_dr$mROS < ss_d$mROS && ss_dr$mROS > ss_h$mROS * 0.8)
  cat(sprintf("    Drug mROS = %.4f (disease = %.4f, healthy = %.4f)\n",
              ss_dr$mROS, ss_d$mROS, ss_h$mROS))
  cat(sprintf("    Reduction from disease: %.1f%%\n",
              (1 - ss_dr$mROS / ss_d$mROS) * 100))
  cat(sprintf("    V3 Drug rescue:             %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero production → zero mROS ---
  cat("V4: Zero production → zero mROS...\n")
  p_zero <- M4_parameters(CI_activity = 1.0, SC_integrity = 1.0)
  p_zero["k_ROS_basal"] <- 0
  y0_zero <- c(mROS = 0.0, SOD2_act = 1.0, GPx4_act = 1.0)
  times_z <- seq(0, 72, by = 0.5)

  out_z <- ode(y = y0_zero, times = times_z, func = M4_odes, parms = p_zero,
               method = "lsoda", atol = 1e-10, rtol = 1e-10)
  df_z <- as.data.frame(out_z)

  v4 <- (max(df_z$mROS) < 1e-10)
  cat(sprintf("    Max mROS with zero production = %.2e\n", max(df_z$mROS)))
  cat(sprintf("    V4 Zero production:         %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness ---
  cat("V5: Non-negativity and boundedness [0,1]...\n")
  all_res <- rbind(res_h, res_d, res_drug)
  v5 <- (all(all_res$mROS >= -1e-10) && all(all_res$mROS <= 1 + 1e-10) &&
         all(all_res$SOD2_act >= -1e-10) && all(all_res$SOD2_act <= 1 + 1e-10) &&
         all(all_res$GPx4_act >= -1e-10) && all(all_res$GPx4_act <= 1 + 1e-10))
  cat(sprintf("    mROS range:     [%.4e, %.4f]\n",
              min(all_res$mROS), max(all_res$mROS)))
  cat(sprintf("    SOD2_act range: [%.4f, %.4f]\n",
              min(all_res$SOD2_act), max(all_res$SOD2_act)))
  cat(sprintf("    GPx4_act range: [%.4f, %.4f]\n",
              min(all_res$GPx4_act), max(all_res$GPx4_act)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical SS matches analytical SS ---
  cat("V6: Numerical vs analytical steady state...\n")
  ss_analytical_h <- M4_analytical_ss(CI_activity = 1.0, SC_integrity = 1.0)
  ss_analytical_d <- M4_analytical_ss(CI_activity = 0.49, SC_integrity = 1.0)

  err_h <- abs(ss_h$mROS - ss_analytical_h["mROS"]) / ss_analytical_h["mROS"] * 100
  err_d <- abs(ss_d$mROS - ss_analytical_d["mROS"]) / ss_analytical_d["mROS"] * 100

  v6 <- (err_h < 1.0 && err_d < 1.0)
  cat(sprintf("    Healthy: numerical = %.6f, analytical = %.6f (err: %.3f%%)\n",
              ss_h$mROS, ss_analytical_h["mROS"], err_h))
  cat(sprintf("    Disease: numerical = %.6f, analytical = %.6f (err: %.3f%%)\n",
              ss_d$mROS, ss_analytical_d["mROS"], err_d))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M4 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M4_scenarios <- function() {

  ## Scenario definitions
  ## Drug effect is INDIRECT: CL stabilisation → improved CI/SC → less leak
  scenarios <- list(
    healthy          = list(CI = 1.00, SC = 1.0, label = "Healthy"),
    disease          = list(CI = 0.49, SC = 1.0, label = "Disease (CI = 49%)"),
    disease_drug_low = list(CI = 0.60, SC = 1.0, label = "Disease + Low dose"),
    disease_drug_hi  = list(CI = 0.80, SC = 1.0, label = "Disease + High dose")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M4(CI_activity = s$CI, SC_integrity = s$SC, duration_h = 72)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))


  ## --- Plot 1: mROS time course ---
  p1 <- ggplot(all_data, aes(x = time, y = mROS, color = scenario)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (CI = 49%)" = "#d62728",
                                  "Disease + Low dose" = "#ff7f0e",
                                  "Disease + High dose" = "#1f77b4")) +
    labs(
      title = "Module M4: Mitochondrial ROS Dynamics",
      subtitle = "Response to CI damage ± indirect drug rescue",
      x = "Time (hours)",
      y = "mROS (normalised)",
      color = "Scenario",
      caption = paste("AIIMS Bhubaneswar | Roy & Padhy | August 2026\n",
                      "ROS kinetics: Cortassa 2004 / Kembro 2013")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M4_ROS_time_course.png", p1, width = 10, height = 6, dpi = 300)
  cat("Saved: M4_ROS_time_course.png\n")


  ## --- Plot 2: All state variables faceted ---
  plot_vars <- c("mROS", "SOD2_act", "GPx4_act")
  long_data <- all_data %>%
    select(time, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("mROS (oxidative burden)",
               "SOD2 active fraction",
               "GPx4 active fraction"))

  p2 <- ggplot(long_data, aes(x = time, y = value, color = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (CI = 49%)" = "#d62728",
                                  "Disease + Low dose" = "#ff7f0e",
                                  "Disease + High dose" = "#1f77b4")) +
    labs(
      title = "M4: State Variables Across Scenarios",
      subtitle = "Antioxidant depletion under oxidative stress",
      x = "Time (hours)",
      y = "Value [0, 1]",
      color = "Scenario"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave("M4_ROS_states.png", p2, width = 10, height = 10, dpi = 300)
  cat("Saved: M4_ROS_states.png\n")


  ## --- Plot 3: CI_activity sweep → mROS steady state ---
  ci_sweep <- seq(0.30, 1.00, by = 0.01)
  sweep_data <- data.frame()

  for (ci in ci_sweep) {
    ss_vals <- M4_analytical_ss(CI_activity = ci, SC_integrity = 1.0)
    sweep_data <- rbind(sweep_data, data.frame(
      CI_activity = ci,
      mROS     = ss_vals["mROS"],
      SOD2_act = ss_vals["SOD2_act"],
      GPx4_act = ss_vals["GPx4_act"],
      mROS_pct = ss_vals["mROS"] / 0.1 * 100
    ))
  }

  p3 <- ggplot(sweep_data, aes(x = CI_activity, y = mROS)) +
    geom_line(linewidth = 1, color = "#d62728") +
    geom_hline(yintercept = 0.10, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = 0.177, linetype = "dotted", color = "#d62728",
               linewidth = 0.5) +
    geom_vline(xintercept = 0.49, linetype = "dotted", color = "gray50") +
    annotate("text", x = 0.49, y = max(sweep_data$mROS) * 0.95,
             label = "Gao: CI = 49%", hjust = -0.1, size = 3, color = "gray40") +
    annotate("text", x = 0.95, y = 0.105,
             label = "Healthy basal", size = 3, color = "gray50") +
    annotate("text", x = 0.32, y = 0.182,
             label = "Choi: 177%", size = 3, color = "#d62728") +
    labs(
      title = "M4: Steady-State mROS vs Complex I Activity",
      subtitle = "Nonlinear amplification via antioxidant depletion",
      x = "CI activity (fraction of WT)",
      y = "mROS (normalised)",
      caption = "SC_integrity = 1.0 (standalone M4)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M4_CI_sweep.png", p3, width = 8, height = 5, dpi = 300)
  cat("Saved: M4_CI_sweep.png\n")


  ## --- Plot 4: Tornado sensitivity (±50% on mROS_ss) ---
  base_mROS <- M4_analytical_ss(CI_activity = 0.49, SC_integrity = 1.0)["mROS"]
  p_base    <- M4_parameters(CI_activity = 0.49, SC_integrity = 1.0)

  param_names <- c("k_ROS_basal", "k_shunt", "k_SOD2", "k_GPx4",
                    "k_ROS_other", "k_SOD2_on", "k_SOD2_off",
                    "k_GPx4_on", "k_GPx4_off")

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      P_mod <- p_mod["k_ROS_basal"] +
               p_mod["k_shunt"] * (1 - 0.49) * 1.0

      f_mod <- function(m) {
        s2 <- p_mod["k_SOD2_on"] / (p_mod["k_SOD2_on"] + p_mod["k_SOD2_off"] * m)
        g4 <- p_mod["k_GPx4_on"] / (p_mod["k_GPx4_on"] + p_mod["k_GPx4_off"] * m)
        clr <- (p_mod["k_SOD2"] * s2 + p_mod["k_GPx4"] * g4 + p_mod["k_ROS_other"]) * m
        as.numeric(P_mod - clr)
      }

      sol <- tryCatch(
        uniroot(f_mod, interval = c(1e-8, 1.0), tol = 1e-12),
        error = function(e) list(root = NA)
      )

      pct_change <- (sol$root - base_mROS) / base_mROS * 100

      sens_data <- rbind(sens_data, data.frame(
        parameter = pname,
        direction = ifelse(direction > 0, "+50%", "-50%"),
        mROS_ss   = sol$root,
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
      span = abs(max(pct_change, na.rm = TRUE) - min(pct_change, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(span)

  tornado$parameter <- factor(tornado$parameter, levels = tornado$parameter)

  p4 <- ggplot(tornado, aes(y = parameter)) +
    geom_segment(aes(x = low, xend = high, yend = parameter),
                 linewidth = 6, color = "#d62728", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M4: Parameter Sensitivity on Disease mROS",
      subtitle = "OAT ±50% perturbation at CI = 0.49",
      x = "Change in mROS_ss (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M4_sensitivity.png", p4, width = 8, height = 5, dpi = 300)
  cat("Saved: M4_sensitivity.png\n\n")


  return(list(data = all_data, sweep = sweep_data, sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

## Run verification
all_passed <- verify_M4()

## Run scenarios and generate plots
if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
  results <- run_all_M4_scenarios()
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Review output above before running scenarios.\n")
  cat("Running scenarios anyway for diagnostic purposes...\n\n")
  results <- run_all_M4_scenarios()
}
