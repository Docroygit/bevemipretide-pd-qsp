##=============================================================================
## MODULE M2: α-SYNUCLEIN AGGREGATION AND CLEARANCE
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Authors : Dr. Sayantan Shankar Roy, Prof. Biswa Mohan Padhy
## Institute: AIIMS Bhubaneswar, Department of Pharmacology
## Date    : August 2026
## Version : 1.0
##
## State variables (3):
##   aSyn_mono — Monomeric (native) α-synuclein [0,1]
##   aSyn_olig — Oligomeric (pathological) α-synuclein [0,1]
##   aSyn_ext  — Extracellular α-synuclein [0,1]
##
## Upstream coupling:
##   mROS       ← M4 (oxidative stress promotes aggregation)
##   CL_ratio   ← M1 (CL dysfunction removes chaperone protection)
##   CL_ext     ← M1 (externalised CL seeds aggregation on OMM surface)
##   ATP        ← M3 (powers autophagy/UPS clearance)
##   MG_active  — scenario parameter until M6 is coded
##
## Key outputs:
##   aSyn_olig → M1 (CL oxidation via k_aSyn_ox), M3 (CI damage via k_CI_aSyn)
##   aSyn_ext  → M6 (microglial activation trigger)
##
## Biology:
##   α-Synuclein aggregation is the central proteinopathy in PD. Native
##   monomeric aSyn misfolds into toxic oligomers under oxidative stress
##   (mROS, from M4) and loss of cardiolipin-mediated membrane stability
##   (CL_ratio, from M1). CL on the inner mitochondrial membrane
##   interacts with aSyn to stabilise its native fold; loss of functional
##   CL removes this chaperone-like protection. Externalised CL on the
##   outer membrane provides a templating surface for aSyn nucleation.
##
##   Oligomers are cleared by ATP-dependent autophagy (macroautophagy,
##   CMA) and proteasomal degradation. Neuroinflammation (MG_active,
##   from M6) impairs this clearance. Secreted aSyn_olig triggers
##   microglial activation via TLR2.
##
##   Drug effect is INDIRECT: bevemipretide stabilises CL (M1), which:
##   (1) improves CL_ratio → less aggregation, more disaggregation
##   (2) improves ETC (M3) → better ATP → more clearance
##   (3) lowers mROS (M4) → less oxidative misfolding
##
## Calibration note:
##   Standalone M2 produces ~1.8x increase in aSyn_olig (healthy→disease).
##   This is amplified to ~5-10x by the vicious cycle M2→M1→M3→M4→M2
##   during integration. Individual modules showing modest effects is
##   expected; the cascade amplification IS the disease mechanism.
##
## Data sources:
##   Bhatt 2009 (doi:10.1523/JNEUROSCI.5515-08.2009): aSyn t½ ≈ 50 h
##   Iljina 2016 (doi:10.1073/pnas.1524128113): aSyn aggregation kinetics
##   Bayir 2009: aSyn-CL interaction, CL-dependent refolding/peroxidation
##   Webb 2003 / Cuervo 2004: aSyn clearance via autophagy/CMA/UPS
##   Ivanova 2024 (doi:10.1002/psp4.13223): QSP aSyn/neuroinflammation
##   Bido et al. poster: SBT-272 reduces aSyn in PD mouse model
##   Kagan 2015: CL externalisation as apoptotic/mitophagy signal
##=============================================================================


library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)


cat("========================================================\n")
cat("  MODULE M2: a-SYNUCLEIN DYNAMICS — Bevemipretide QSP\n")
cat("  Version 1.0 | Upstream-coupled to M1 (CL), M3 (ATP),\n")
cat("                 M4 (ROS); output to M1, M3, M6\n")
cat("========================================================\n\n")


##-----------------------------------------------------------------------------
## 1. PARAMETERS
##-----------------------------------------------------------------------------

M2_parameters <- function(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                          ATP = 0.95, MG_active = 0.0) {

  p <- c(
    ## --- aSyn synthesis and turnover ---
    k_syn       = 0.50,      # aSyn synthesis rate (h^-1), self-limiting
                              # ESTIMATED: derived from aSyn t½ ≈ 50h
                              # (Bhatt 2009) and target SS monomer ≈ 0.95
    k_turn      = 0.014,     # Monomeric aSyn turnover (h^-1)
                              # LITERATURE: aSyn t½ ≈ 50h in mouse brain
                              # → ln(2)/50 ≈ 0.014 (Bhatt 2009)

    ## --- Aggregation ---
    k_agg       = 0.009,     # Basal aggregation rate constant (h^-1)
                              # ESTIMATED: tuned for aSyn_olig_healthy ≈ 0.05
                              # Iljina 2016: nucleation rate constants
                              # (in vitro μM; normalised for [0,1] model)
    k_ROS_agg   = 2.0,       # ROS amplification of aggregation (dimensionless)
                              # ASSUMED: oxidative stress promotes misfolding
                              # ROS modifies Met residues → aggregation-prone
    k_CL_loss   = 5.0,       # CL dysfunction amplifies aggregation (dimensionless)
                              # ASSUMED: Bayir 2009, CL on IMM stabilises native
                              # aSyn fold; loss of CL_ratio removes this protection
    k_seed_CL   = 1.0,       # CL_ext seeding of aggregation (dimensionless)
                              # ASSUMED: Kagan 2015, externalised CL on OMM
                              # provides templating surface for aSyn nucleation

    ## --- Disaggregation ---
    k_refold    = 0.010,     # CL-dependent disaggregation rate (h^-1)
                              # ASSUMED: slow, CL-mediated refolding of oligomers
                              # Bayir 2009: CL can promote aSyn refolding
                              # at physiological conditions

    ## --- Clearance ---
    k_clear     = 0.241,     # ATP-dependent oligomer clearance (h^-1)
                              # ESTIMATED: autophagy (macro/CMA) and UPS
                              # Webb 2003, Cuervo 2004; tuned for mass balance
    k_impair    = 0.50,      # Inflammation impairment factor (dimensionless)
                              # ASSUMED: TNF/pro-inflammatory cytokines
                              # suppress autophagy flux; MG_active scales [0,1]

    ## --- Secretion and extracellular clearance ---
    k_secrete   = 0.005,     # Oligomer secretion rate (h^-1)
                              # ASSUMED: exosomal/non-conventional secretion,
                              # slow minor pathway
    k_clear_ext = 0.10,      # Extracellular aSyn clearance (h^-1)
                              # ASSUMED: ISF drainage, glymphatic clearance

    ## --- Upstream coupling (scenario parameters) ---
    mROS       = mROS,
    CL_ratio   = CL_ratio,
    CL_ext     = CL_ext,
    ATP        = ATP,
    MG_active  = MG_active
  )

  return(p)
}


##-----------------------------------------------------------------------------
## 2. INITIAL CONDITIONS
##-----------------------------------------------------------------------------

M2_init <- function() {
  c(
    aSyn_mono = 0.95,    # Near healthy SS (mostly monomeric)
    aSyn_olig = 0.05,    # Low basal oligomer level
    aSyn_ext  = 0.003    # Minimal extracellular aSyn
  )
}


##-----------------------------------------------------------------------------
## 3. ODE SYSTEM
##-----------------------------------------------------------------------------

M2_odes <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    ## --- Clamp states to [0, 1] ---
    aSyn_mono <- max(0, min(1, aSyn_mono))
    aSyn_olig <- max(0, min(1, aSyn_olig))
    aSyn_ext  <- max(0, min(1, aSyn_ext))

    ## --- Aggregation stress function ---
    ## Combines oxidative stress, CL dysfunction, and CL_ext seeding
    ## Basal rate of 1.0 always present (constitutive misfolding)
    f_stress <- (1
                 + k_ROS_agg * mROS
                 + k_CL_loss * (1 - CL_ratio)
                 + k_seed_CL * CL_ext)

    ## Effective aggregation rate
    agg_rate <- k_agg * f_stress

    ## Effective clearance (ATP-dependent, inflammation-impaired)
    clear_eff <- k_clear * ATP * (1 - k_impair * MG_active)

    ## --- ODE 1: Monomeric α-synuclein ---
    ## Synthesis: self-limiting, approaches carrying capacity
    ## Loss: aggregation (stress-driven), normal turnover
    ## Gain: disaggregation (CL-dependent refolding of oligomers)
    daSyn_mono <- (k_syn * (1 - aSyn_mono)
                   - agg_rate * aSyn_mono
                   + k_refold * CL_ratio * aSyn_olig
                   - k_turn * aSyn_mono)

    ## --- ODE 2: Oligomeric α-synuclein ---
    ## Gain: aggregation from monomers
    ## Loss: disaggregation, autophagy/UPS clearance, secretion
    daSyn_olig <- (agg_rate * aSyn_mono
                   - k_refold * CL_ratio * aSyn_olig
                   - clear_eff * aSyn_olig
                   - k_secrete * aSyn_olig)

    ## --- ODE 3: Extracellular α-synuclein ---
    ## From oligomer secretion; cleared by ISF drainage/glymphatic
    daSyn_ext <- (k_secrete * aSyn_olig
                  - k_clear_ext * aSyn_ext)

    ## --- Derived quantities ---
    olig_fraction  <- aSyn_olig / max(aSyn_mono + aSyn_olig, 1e-8)
    agg_flux       <- agg_rate * aSyn_mono
    clear_flux     <- clear_eff * aSyn_olig
    refold_flux    <- k_refold * CL_ratio * aSyn_olig
    secrete_flux   <- k_secrete * aSyn_olig

    list(
      c(daSyn_mono, daSyn_olig, daSyn_ext),
      f_stress       = f_stress,
      agg_flux       = agg_flux,
      clear_flux     = clear_flux,
      refold_flux    = refold_flux,
      secrete_flux   = secrete_flux,
      olig_fraction  = olig_fraction,
      clear_eff      = clear_eff
    )
  })
}


##-----------------------------------------------------------------------------
## 4. SOLVER
##-----------------------------------------------------------------------------

run_M2 <- function(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                   ATP = 0.95, MG_active = 0.0,
                   duration_h = 200, dt = 0.1, y0 = NULL) {

  parms <- M2_parameters(mROS, CL_ratio, CL_ext, ATP, MG_active)
  if (is.null(y0)) y0 <- M2_init()
  times <- seq(0, duration_h, by = dt)

  out <- ode(y = y0, times = times, func = M2_odes, parms = parms,
             method = "lsoda", atol = 1e-10, rtol = 1e-10)

  result <- as.data.frame(out)
  result$mROS      <- mROS
  result$CL_ratio  <- CL_ratio
  result$CL_ext    <- CL_ext
  result$ATP       <- ATP
  result$MG_active <- MG_active

  ss <- result[nrow(result), ]
  cat(sprintf(
    "M2 | mROS=%.3f CL=%.3f ATP=%.2f MG=%.1f | SS: mono=%.3f olig=%.4f ext=%.4f\n",
    mROS, CL_ratio, ATP, MG_active,
    ss$aSyn_mono, ss$aSyn_olig, ss$aSyn_ext))

  return(result)
}


##-----------------------------------------------------------------------------
## 5. ANALYTICAL STEADY STATE
##-----------------------------------------------------------------------------

M2_analytical_ss <- function(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                             ATP = 0.95, MG_active = 0.0) {

  p <- M2_parameters(mROS, CL_ratio, CL_ext, ATP, MG_active)

  ## Stress-dependent aggregation rate
  f_stress <- as.numeric(
    1 + p["k_ROS_agg"] * mROS +
    p["k_CL_loss"] * (1 - CL_ratio) +
    p["k_seed_CL"] * CL_ext)

  alpha <- as.numeric(p["k_agg"]) * f_stress  # effective aggregation rate

  ## Effective clearance
  clear_eff <- as.numeric(p["k_clear"]) * ATP *
               (1 - as.numeric(p["k_impair"]) * MG_active)

  ## Total removal rate from oligomer pool
  C <- clear_eff + as.numeric(p["k_secrete"])

  ## Total removal including disaggregation
  beta <- as.numeric(p["k_refold"]) * CL_ratio + C

  ## Synthesis + turnover
  S <- as.numeric(p["k_syn"]) + as.numeric(p["k_turn"])

  ## Analytical solution (linear system, fully closed-form)
  ## From simultaneous daSyn_mono = 0 and daSyn_olig = 0:
  aSyn_olig_ss <- alpha * as.numeric(p["k_syn"]) / (S * beta + alpha * C)

  aSyn_mono_ss <- (as.numeric(p["k_syn"]) - C * aSyn_olig_ss) / S

  aSyn_ext_ss <- as.numeric(p["k_secrete"]) * aSyn_olig_ss /
                 as.numeric(p["k_clear_ext"])

  return(c(aSyn_mono  = aSyn_mono_ss,
           aSyn_olig  = aSyn_olig_ss,
           aSyn_ext   = aSyn_ext_ss,
           f_stress   = f_stress))
}


##-----------------------------------------------------------------------------
## 6. VERIFICATION TESTS (FDA QSP Guidance Section IV)
##-----------------------------------------------------------------------------

verify_M2 <- function() {

  cat("=== MODULE M2 VERIFICATION TESTS ===\n\n")
  pass_count  <- 0
  total_tests <- 6
  t_verify    <- 200


  ## --- V1: Healthy steady state in expected range ---
  ## Upstream from M1/M3/M4: mROS=0.10, CL_ratio=0.948, ATP=0.95
  cat("V1: Healthy steady state...\n")
  res_h <- run_M2(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                  ATP = 0.95, MG_active = 0.0, duration_h = t_verify)
  ss_h  <- res_h[nrow(res_h), ]

  v1 <- (ss_h$aSyn_mono > 0.90 && ss_h$aSyn_mono < 0.98 &&
         ss_h$aSyn_olig > 0.03 && ss_h$aSyn_olig < 0.08 &&
         ss_h$aSyn_ext > 0.001 && ss_h$aSyn_ext < 0.010)
  cat(sprintf("    aSyn_mono = %.4f (expect 0.90-0.98)\n", ss_h$aSyn_mono))
  cat(sprintf("    aSyn_olig = %.4f (expect 0.03-0.08)\n", ss_h$aSyn_olig))
  cat(sprintf("    aSyn_ext  = %.5f (expect 0.001-0.010)\n", ss_h$aSyn_ext))
  cat(sprintf("    f_stress  = %.3f\n", ss_h$f_stress))
  cat(sprintf("    V1 Healthy SS:              %s\n\n",
              ifelse(v1, "PASS", "FAIL")))
  if (v1) pass_count <- pass_count + 1


  ## --- V2: Disease perturbation increases aSyn_olig ---
  ## Upstream from M1/M3/M4 at disease:
  ##   mROS = 0.176 (Choi 2022, 177% basal)
  ##   CL_ratio = 0.807 (Gao 2017, ~22% CL_n reduction)
  ##   ATP ≈ 0.82 (M3 disease SS estimate)
  cat("V2: Disease perturbation (mROS=0.176, CL_ratio=0.807)...\n")
  res_d <- run_M2(mROS = 0.176, CL_ratio = 0.807, CL_ext = 0.003,
                  ATP = 0.82, MG_active = 0.0, duration_h = t_verify)
  ss_d  <- res_d[nrow(res_d), ]

  fold_increase <- ss_d$aSyn_olig / ss_h$aSyn_olig
  v2 <- (ss_d$aSyn_olig > ss_h$aSyn_olig &&
         ss_d$aSyn_mono < ss_h$aSyn_mono &&
         ss_d$aSyn_ext > ss_h$aSyn_ext &&
         fold_increase > 1.3)
  cat(sprintf("    aSyn_olig: %.4f -> %.4f (%.1fx increase)\n",
              ss_h$aSyn_olig, ss_d$aSyn_olig, fold_increase))
  cat(sprintf("    aSyn_mono: %.4f -> %.4f (depleted: %s)\n",
              ss_h$aSyn_mono, ss_d$aSyn_mono,
              ifelse(ss_d$aSyn_mono < ss_h$aSyn_mono, "YES", "NO")))
  cat(sprintf("    aSyn_ext:  %.5f -> %.5f\n", ss_h$aSyn_ext, ss_d$aSyn_ext))
  cat(sprintf("    f_stress:  %.3f -> %.3f (%.1fx)\n",
              ss_h$f_stress, ss_d$f_stress,
              ss_d$f_stress / ss_h$f_stress))
  cat(sprintf("    Note: integrated model amplifies via vicious cycle\n"))
  cat(sprintf("    V2 Disease perturbation:    %s\n\n",
              ifelse(v2, "PASS", "FAIL")))
  if (v2) pass_count <- pass_count + 1


  ## --- V3: Drug rescue reduces aSyn_olig ---
  ## Drug acts via M1: CL_ratio improves 0.807 -> 0.911 (high dose)
  ## Also improves ATP slightly (better ETC), mROS stays at disease
  cat("V3: Drug rescue (CL_ratio 0.807 -> 0.911 from M1)...\n")
  res_drug <- run_M2(mROS = 0.176, CL_ratio = 0.911, CL_ext = 0.002,
                     ATP = 0.85, MG_active = 0.0, duration_h = t_verify)
  ss_rx <- res_drug[nrow(res_drug), ]

  v3 <- (ss_rx$aSyn_olig < ss_d$aSyn_olig &&
         ss_rx$aSyn_olig > ss_h$aSyn_olig * 0.8)
  cat(sprintf("    aSyn_olig: disease=%.4f -> drug=%.4f\n",
              ss_d$aSyn_olig, ss_rx$aSyn_olig))
  cat(sprintf("    Reduction from disease: %.1f%%\n",
              (1 - ss_rx$aSyn_olig / ss_d$aSyn_olig) * 100))
  cat(sprintf("    V3 Drug rescue:             %s\n\n",
              ifelse(v3, "PASS", "FAIL")))
  if (v3) pass_count <- pass_count + 1


  ## --- V4: Zero aggregation → zero oligomers ---
  ## With k_agg = 0, no oligomers should form from zero IC
  cat("V4: Zero aggregation -> zero oligomers...\n")
  p_zero <- M2_parameters(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                          ATP = 0.95, MG_active = 0.0)
  p_zero["k_agg"] <- 0
  y0_zero <- c(aSyn_mono = 0.95, aSyn_olig = 0.0, aSyn_ext = 0.0)
  times_z <- seq(0, 200, by = 0.5)

  out_z <- ode(y = y0_zero, times = times_z, func = M2_odes, parms = p_zero,
               method = "lsoda", atol = 1e-10, rtol = 1e-10)
  df_z <- as.data.frame(out_z)

  v4 <- (max(df_z$aSyn_olig) < 1e-10 && max(df_z$aSyn_ext) < 1e-10)
  cat(sprintf("    Max aSyn_olig with k_agg=0: %.2e\n", max(df_z$aSyn_olig)))
  cat(sprintf("    Max aSyn_ext  with k_agg=0: %.2e\n", max(df_z$aSyn_ext)))
  cat(sprintf("    V4 Zero aggregation:        %s\n\n",
              ifelse(v4, "PASS", "FAIL")))
  if (v4) pass_count <- pass_count + 1


  ## --- V5: Non-negativity and boundedness [0,1] ---
  cat("V5: Non-negativity and boundedness...\n")
  all_res <- rbind(res_h, res_d, res_drug)
  v5 <- (all(all_res$aSyn_mono >= -1e-10) &&
         all(all_res$aSyn_mono <= 1 + 1e-10) &&
         all(all_res$aSyn_olig >= -1e-10) &&
         all(all_res$aSyn_olig <= 1 + 1e-10) &&
         all(all_res$aSyn_ext >= -1e-10) &&
         all(all_res$aSyn_ext <= 1 + 1e-10))
  cat(sprintf("    aSyn_mono range: [%.4f, %.4f]\n",
              min(all_res$aSyn_mono), max(all_res$aSyn_mono)))
  cat(sprintf("    aSyn_olig range: [%.4e, %.4f]\n",
              min(all_res$aSyn_olig), max(all_res$aSyn_olig)))
  cat(sprintf("    aSyn_ext range:  [%.4e, %.6f]\n",
              min(all_res$aSyn_ext), max(all_res$aSyn_ext)))
  cat(sprintf("    V5 Boundedness:             %s\n\n",
              ifelse(v5, "PASS", "FAIL")))
  if (v5) pass_count <- pass_count + 1


  ## --- V6: Numerical SS matches analytical SS (<1%) ---
  cat("V6: Numerical vs analytical steady state...\n")
  ss_an_h <- M2_analytical_ss(mROS = 0.10, CL_ratio = 0.948, CL_ext = 0.002,
                              ATP = 0.95, MG_active = 0.0)
  ss_an_d <- M2_analytical_ss(mROS = 0.176, CL_ratio = 0.807, CL_ext = 0.003,
                              ATP = 0.82, MG_active = 0.0)

  err_olig_h <- abs(ss_h$aSyn_olig - ss_an_h["aSyn_olig"]) /
                ss_an_h["aSyn_olig"] * 100
  err_olig_d <- abs(ss_d$aSyn_olig - ss_an_d["aSyn_olig"]) /
                ss_an_d["aSyn_olig"] * 100
  err_mono_h <- abs(ss_h$aSyn_mono - ss_an_h["aSyn_mono"]) /
                ss_an_h["aSyn_mono"] * 100
  err_mono_d <- abs(ss_d$aSyn_mono - ss_an_d["aSyn_mono"]) /
                ss_an_d["aSyn_mono"] * 100

  v6 <- (err_olig_h < 1.0 && err_olig_d < 1.0 &&
         err_mono_h < 1.0 && err_mono_d < 1.0)
  cat(sprintf("    Healthy olig: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$aSyn_olig, ss_an_h["aSyn_olig"], err_olig_h))
  cat(sprintf("    Disease olig: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$aSyn_olig, ss_an_d["aSyn_olig"], err_olig_d))
  cat(sprintf("    Healthy mono: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_h$aSyn_mono, ss_an_h["aSyn_mono"], err_mono_h))
  cat(sprintf("    Disease mono: num=%.5f ana=%.5f (err: %.4f%%)\n",
              ss_d$aSyn_mono, ss_an_d["aSyn_mono"], err_mono_d))
  cat(sprintf("    V6 Analytical match (<1%%): %s\n\n",
              ifelse(v6, "PASS", "FAIL")))
  if (v6) pass_count <- pass_count + 1


  ## --- Summary ---
  cat(sprintf("=== M2 VERIFICATION RESULT: %d/%d tests passed ===\n\n",
              pass_count, total_tests))

  return(pass_count == total_tests)
}


##-----------------------------------------------------------------------------
## 7. MULTI-SCENARIO SIMULATION AND PLOTTING
##-----------------------------------------------------------------------------

run_all_M2_scenarios <- function() {

  ## Scenario definitions
  ## Upstream values from M1/M3/M4 standalone steady states
  ## Drug effect is INDIRECT via M1: CL_ratio improves -> less aggregation
  ## mROS and ATP held at disease level for drug scenarios
  ## (in integrated model, drug also improves mROS via M4 and ATP via M3)
  scenarios <- list(
    healthy  = list(mROS = 0.10,  CL = 0.948, CL_ext = 0.002,
                    ATP = 0.95, MG = 0.0,
                    label = "Healthy"),
    disease  = list(mROS = 0.176, CL = 0.807, CL_ext = 0.003,
                    ATP = 0.82, MG = 0.0,
                    label = "Disease (no drug)"),
    drug_low = list(mROS = 0.176, CL = 0.827, CL_ext = 0.003,
                    ATP = 0.83, MG = 0.0,
                    label = "Disease + SBT-272 low"),
    drug_hi  = list(mROS = 0.176, CL = 0.911, CL_ext = 0.002,
                    ATP = 0.85, MG = 0.0,
                    label = "Disease + SBT-272 high")
  )

  ## Run all scenarios
  all_data <- data.frame()
  for (name in names(scenarios)) {
    s <- scenarios[[name]]
    res <- run_M2(mROS = s$mROS, CL_ratio = s$CL, CL_ext = s$CL_ext,
                  ATP = s$ATP, MG_active = s$MG, duration_h = 200)
    res$scenario <- s$label
    all_data <- rbind(all_data, res)
  }

  all_data$scenario <- factor(all_data$scenario,
    levels = sapply(scenarios, function(s) s$label))


  ## --- Plot 1: Time course of all state variables ---
  plot_vars <- c("aSyn_mono", "aSyn_olig", "aSyn_ext")

  long_data <- all_data %>%
    select(time, scenario, all_of(plot_vars)) %>%
    pivot_longer(cols = all_of(plot_vars),
                 names_to = "variable", values_to = "value")

  long_data$variable <- factor(long_data$variable,
    levels = plot_vars,
    labels = c("aSyn monomeric (native)",
               "aSyn oligomeric (pathological)",
               "aSyn extracellular"))

  p1 <- ggplot(long_data, aes(x = time, y = value, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("Healthy" = "#2ca02c",
                                  "Disease (no drug)" = "#d62728",
                                  "Disease + SBT-272 low" = "#ff7f0e",
                                  "Disease + SBT-272 high" = "#1f77b4")) +
    labs(
      title = "Module M2: a-Synuclein Dynamics",
      subtitle = paste("Upstream: mROS from M4, CL_ratio from M1,",
                        "ATP from M3"),
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

  ggsave("M2_aSyn_time_course.png", p1,
         width = 10, height = 10, dpi = 300, bg = "white")
  cat("Saved: M2_aSyn_time_course.png\n")


  ## --- Plot 2: CL_ratio sweep → aSyn_olig SS ---
  ## Shows how CL stabilisation by drug reduces aSyn_olig
  cat("Running CL_ratio sweep (analytical SS)...\n")
  cl_sweep <- seq(0.60, 1.00, by = 0.005)
  sweep_CL <- data.frame()

  for (cl in cl_sweep) {
    ss_vals <- M2_analytical_ss(mROS = 0.176, CL_ratio = cl, CL_ext = 0.003,
                                ATP = 0.82, MG_active = 0.0)
    sweep_CL <- rbind(sweep_CL, data.frame(
      CL_ratio   = cl,
      aSyn_olig  = ss_vals["aSyn_olig"],
      aSyn_mono  = ss_vals["aSyn_mono"],
      aSyn_ext   = ss_vals["aSyn_ext"],
      f_stress   = ss_vals["f_stress"]
    ))
  }

  p2 <- ggplot(sweep_CL, aes(x = CL_ratio, y = aSyn_olig)) +
    geom_line(linewidth = 1, color = "#d62728") +
    geom_vline(xintercept = c(0.807, 0.948), linetype = "dashed",
               color = "gray50") +
    geom_vline(xintercept = 0.911, linetype = "dotted",
               color = "#1f77b4") +
    annotate("text", x = 0.807, y = max(sweep_CL$aSyn_olig) * 0.95,
             label = "Disease", hjust = 1.1, size = 3, color = "#d62728") +
    annotate("text", x = 0.948, y = max(sweep_CL$aSyn_olig) * 0.95,
             label = "Healthy", hjust = -0.1, size = 3, color = "#2ca02c") +
    annotate("text", x = 0.911, y = max(sweep_CL$aSyn_olig) * 0.80,
             label = "SBT-272\nhigh dose", hjust = -0.1, size = 3,
             color = "#1f77b4") +
    labs(
      title = "M2: aSyn Oligomer vs CL Functional Ratio",
      subtitle = "Disease background (mROS = 0.176, ATP = 0.82)",
      x = "CL_ratio (from M1)",
      y = "aSyn_olig at steady state",
      caption = "Drug improves CL_ratio -> reduces aSyn_olig"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M2_CL_sweep.png", p2,
         width = 8, height = 5, dpi = 300, bg = "white")
  cat("Saved: M2_CL_sweep.png\n")


  ## --- Plot 3: mROS sweep → aSyn_olig SS ---
  cat("Running mROS sweep (analytical SS)...\n")
  mros_sweep <- seq(0.05, 0.30, by = 0.005)
  sweep_ROS <- data.frame()

  for (mr in mros_sweep) {
    ss_vals <- M2_analytical_ss(mROS = mr, CL_ratio = 0.807,
                                CL_ext = 0.003, ATP = 0.82, MG_active = 0.0)
    sweep_ROS <- rbind(sweep_ROS, data.frame(
      mROS       = mr,
      aSyn_olig  = ss_vals["aSyn_olig"],
      aSyn_mono  = ss_vals["aSyn_mono"],
      f_stress   = ss_vals["f_stress"]
    ))
  }

  p3 <- ggplot(sweep_ROS, aes(x = mROS, y = aSyn_olig)) +
    geom_line(linewidth = 1, color = "#9467bd") +
    geom_vline(xintercept = c(0.10, 0.176), linetype = "dashed",
               color = "gray50") +
    annotate("text", x = 0.10, y = max(sweep_ROS$aSyn_olig) * 0.95,
             label = "Healthy\nmROS", hjust = 1.1, size = 3,
             color = "#2ca02c") +
    annotate("text", x = 0.176, y = max(sweep_ROS$aSyn_olig) * 0.80,
             label = "Disease\nmROS", hjust = -0.1, size = 3,
             color = "#d62728") +
    labs(
      title = "M2: aSyn Oligomer vs Mitochondrial ROS",
      subtitle = "Disease CL background (CL_ratio = 0.807, ATP = 0.82)",
      x = "mROS (from M4)",
      y = "aSyn_olig at steady state",
      caption = "Drug reduces mROS indirectly via CL -> ETC -> ROS"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40")
    )

  ggsave("M2_mROS_sweep.png", p3,
         width = 8, height = 5, dpi = 300, bg = "white")
  cat("Saved: M2_mROS_sweep.png\n")


  ## --- Plot 4: Sensitivity tornado (OAT ±50% on aSyn_olig_ss) ---
  cat("Running one-at-a-time sensitivity analysis...\n")
  base_ss <- M2_analytical_ss(mROS = 0.176, CL_ratio = 0.807,
                               CL_ext = 0.003, ATP = 0.82, MG_active = 0.0)
  base_olig <- base_ss["aSyn_olig"]

  p_base <- M2_parameters(mROS = 0.176, CL_ratio = 0.807,
                           CL_ext = 0.003, ATP = 0.82, MG_active = 0.0)

  param_names <- c("k_syn", "k_turn", "k_agg", "k_ROS_agg", "k_CL_loss",
                    "k_seed_CL", "k_refold", "k_clear", "k_secrete")

  sens_data <- data.frame()
  for (pname in param_names) {
    for (direction in c(-0.5, 0.5)) {
      p_mod <- p_base
      p_mod[pname] <- p_base[pname] * (1 + direction)

      ## Recompute analytical SS with modified parameter
      f_mod <- as.numeric(
        1 + p_mod["k_ROS_agg"] * 0.176 +
        p_mod["k_CL_loss"] * (1 - 0.807) +
        p_mod["k_seed_CL"] * 0.003)

      alpha_mod <- as.numeric(p_mod["k_agg"]) * f_mod
      clear_mod <- as.numeric(p_mod["k_clear"]) * 0.82
      C_mod     <- clear_mod + as.numeric(p_mod["k_secrete"])
      beta_mod  <- as.numeric(p_mod["k_refold"]) * 0.807 + C_mod
      S_mod     <- as.numeric(p_mod["k_syn"]) + as.numeric(p_mod["k_turn"])

      olig_mod  <- alpha_mod * as.numeric(p_mod["k_syn"]) /
                   (S_mod * beta_mod + alpha_mod * C_mod)

      pct_change <- as.numeric(
        (olig_mod - base_olig) / base_olig * 100)

      sens_data <- rbind(sens_data, data.frame(
        parameter  = pname,
        direction  = ifelse(direction > 0, "+50%", "-50%"),
        aSyn_olig  = olig_mod,
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

  p4 <- ggplot(tornado, aes(y = parameter)) +
    geom_segment(aes(x = low, xend = high, yend = parameter),
                 linewidth = 6, color = "#9467bd", alpha = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title = "M2: Parameter Sensitivity on Disease aSyn_olig",
      subtitle = "OAT +/-50% perturbation (mROS=0.176, CL=0.807, ATP=0.82)",
      x = "Change in aSyn_olig_ss (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "gray40"),
      panel.grid.minor = element_blank()
    )

  ggsave("M2_sensitivity.png", p4,
         width = 10, height = 5, dpi = 300, bg = "white")
  cat("Saved: M2_sensitivity.png\n\n")


  return(list(data = all_data, sweep_CL = sweep_CL,
              sweep_ROS = sweep_ROS, sensitivity = sens_data))
}


##-----------------------------------------------------------------------------
## 8. EXECUTE
##-----------------------------------------------------------------------------

## Run verification
all_passed <- verify_M2()

## Run scenarios and generate plots
if (all_passed) {
  cat("All verification tests passed. Running scenarios...\n\n")
} else {
  cat("WARNING: Not all verification tests passed.\n")
  cat("Running scenarios anyway for diagnostic purposes...\n\n")
}

results <- run_all_M2_scenarios()

## Print SS summary table
cat("Analytical SS for key scenarios:\n")
for (label in c("Healthy", "Disease", "Drug_low", "Drug_high")) {
  if (label == "Healthy") {
    ss <- M2_analytical_ss(mROS=0.10, CL_ratio=0.948, CL_ext=0.002,
                           ATP=0.95, MG_active=0.0)
  } else if (label == "Disease") {
    ss <- M2_analytical_ss(mROS=0.176, CL_ratio=0.807, CL_ext=0.003,
                           ATP=0.82, MG_active=0.0)
  } else if (label == "Drug_low") {
    ss <- M2_analytical_ss(mROS=0.176, CL_ratio=0.827, CL_ext=0.003,
                           ATP=0.83, MG_active=0.0)
  } else {
    ss <- M2_analytical_ss(mROS=0.176, CL_ratio=0.911, CL_ext=0.002,
                           ATP=0.85, MG_active=0.0)
  }
  cat(sprintf("  %-10s mono=%.3f olig=%.4f ext=%.4f f_stress=%.3f\n",
              label, ss["aSyn_mono"], ss["aSyn_olig"],
              ss["aSyn_ext"], ss["f_stress"]))
}

## Inflammation sensitivity preview (M6 coupling)
cat("\nInflammation sensitivity (MG_active sweep, disease background):\n")
for (mg in c(0.0, 0.25, 0.50, 0.75)) {
  ss_mg <- M2_analytical_ss(mROS=0.176, CL_ratio=0.807, CL_ext=0.003,
                             ATP=0.82, MG_active=mg)
  cat(sprintf("  MG_active=%.2f -> aSyn_olig=%.4f (%.0f%% vs MG=0)\n",
              mg, ss_mg["aSyn_olig"],
              ss_mg["aSyn_olig"] /
              M2_analytical_ss(mROS=0.176, CL_ratio=0.807, CL_ext=0.003,
                               ATP=0.82, MG_active=0.0)["aSyn_olig"] * 100))
}

cat("\nSensitivity analysis results:\n")
print(results$sensitivity[, c("parameter", "direction", "pct_change")],
      row.names = FALSE)

cat("\n========================================================\n")
cat("  MODULE M2 COMPLETE\n")
cat("  Outputs: aSyn_olig -> M1 (CL oxidation), M3 (CI damage)\n")
cat("           aSyn_ext  -> M6 (microglial activation)\n")
cat("  Upstream: mROS <- M4, CL_ratio/CL_ext <- M1,\n")
cat("            ATP <- M3, MG_active <- M6 (scenario)\n")
cat("  Note: standalone ~1.8x increase; vicious cycle amplifies\n")
cat("========================================================\n")
