##=============================================================================
## SENSITIVITY ANALYSIS: Bevemipretide PD QSP Integrated Model
##
## Pipeline:
##   1. OAT (One-At-a-Time) ±50% → Tornado diagram
##   2. Morris screening (r=20 trajectories) → μ* vs σ plot
##   3. Sobol (top 10 parameters from Morris) → First-order & total indices
##
## Outputs tracked: DA_neuron, Motor_score, CI_activity, aSyn_olig, mROS
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

cat("================================================================\n")
cat("  SENSITIVITY ANALYSIS: Bevemipretide PD QSP\n")
cat("  OAT → Morris → Sobol pipeline\n")
cat("================================================================\n\n")

## Source integrated model functions (without running calibration)
SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

## Calibrated parameter values
CAL <- list(
  K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4,
  CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35
)


##-----------------------------------------------------------------------------
## PARAMETER SPACE DEFINITION
##-----------------------------------------------------------------------------

param_space <- data.frame(
  name = c(
    ## M1: Cardiolipin
    "k_ox", "k_red", "k_drug", "k_aSyn_ox", "k_ext", "k_syn_CL",
    ## M2: alpha-Synuclein
    "k_agg", "k_ROS_agg", "k_CL_loss", "k_seed_CL", "k_clear", "k_syn_aSyn",
    ## M3: ETC
    "k_CI_repair", "k_CI_ROS", "k_CI_aSyn", "k_mPTP_open", "K_mPTP", "k_resp",
    ## M4: ROS
    "k_ROS_basal", "k_shunt", "k_SOD2",
    ## M5: Mitophagy
    "k_death", "K_death", "k_damage_mPTP", "k_repair", "k_biogen",
    ## M6: Neuroinflammation
    "k_act_MG", "k_auto", "k_deact_MG",
    ## M7: Motor
    "K_motor", "k_inflam",
    ## Disease modifiers
    "CI_max", "alpha_clear", "k_impair"
  ),
  module = c(
    rep("M1:CL", 6), rep("M2:aSyn", 6), rep("M3:ETC", 6),
    rep("M4:ROS", 3), rep("M5:Death", 5), rep("M6:Inflam", 3),
    rep("M7:Motor", 2), rep("Disease", 3)
  ),
  stringsAsFactors = FALSE
)

## Get baseline values from calibrated parameter set
base_parms <- integrated_parameters(
  dose_mg_kg = 0, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear,
  k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
  k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death
)

param_space$baseline <- sapply(param_space$name, function(nm) {
  as.numeric(base_parms[nm])
})

n_params <- nrow(param_space)
cat(sprintf("Parameters to analyse: %d\n\n", n_params))

## Key outputs to track
OUTPUT_NAMES <- c("DA_neuron", "Motor_score", "CI_activity",
                  "aSyn_olig", "mROS")


##-----------------------------------------------------------------------------
## HELPER: Run model with parameter override
##-----------------------------------------------------------------------------

run_sa <- function(parms_vec) {
  override <- as.list(parms_vec)
  names(override) <- param_space$name

  res <- tryCatch({
    run_integrated(
      dose_mg_kg = 0, duration_days = 35,
      CI_max = override$CI_max, alpha_clear = override$alpha_clear,
      k_impair_val = override$k_impair,
      K_death_val = override$K_death,
      k_damage_mPTP_val = override$k_damage_mPTP,
      n_death_val = CAL$n_death,
      parms_override = override,
      quiet = TRUE
    )
  }, error = function(e) NULL)

  if (is.null(res)) return(rep(NA, length(OUTPUT_NAMES)))

  ss <- res[nrow(res), ]
  sapply(OUTPUT_NAMES, function(v) as.numeric(ss[[v]]))
}


##-----------------------------------------------------------------------------
## 1. OAT (One-At-a-Time) ±50% TORNADO
##-----------------------------------------------------------------------------

cat("================================================================\n")
cat("  STAGE 1: OAT ±50%% TORNADO ANALYSIS\n")
cat(sprintf("  %d parameters × 2 perturbations = %d runs\n",
            n_params, n_params * 2))
cat("================================================================\n\n")

## Baseline run
baseline_vals <- run_sa(param_space$baseline)
cat("Baseline values:\n")
for (i in seq_along(OUTPUT_NAMES)) {
  cat(sprintf("  %-14s = %.4f\n", OUTPUT_NAMES[i], baseline_vals[i]))
}
cat("\n")

## OAT sweep
oat_results <- data.frame()

for (i in seq_len(n_params)) {
  for (direction in c(-0.5, 0.5)) {
    pv <- param_space$baseline
    pv[i] <- pv[i] * (1 + direction)

    ## Ensure non-negative and n_death stays integer-ish
    pv[i] <- max(pv[i], 1e-6)

    vals <- run_sa(pv)

    for (j in seq_along(OUTPUT_NAMES)) {
      delta <- vals[j] - baseline_vals[j]
      pct   <- ifelse(baseline_vals[j] != 0,
                       delta / baseline_vals[j] * 100, 0)

      oat_results <- rbind(oat_results, data.frame(
        param    = param_space$name[i],
        module   = param_space$module[i],
        direction = ifelse(direction > 0, "+50%", "-50%"),
        output   = OUTPUT_NAMES[j],
        value    = vals[j],
        delta    = delta,
        pct_change = pct,
        stringsAsFactors = FALSE
      ))
    }
  }
  if (i %% 10 == 0) cat(sprintf("  OAT: %d/%d parameters done\n", i, n_params))
}

cat(sprintf("  OAT: %d/%d parameters done\n\n", n_params, n_params))

## Compute max |%change| per param-output pair for ranking
oat_summary <- oat_results %>%
  group_by(param, module, output) %>%
  summarise(max_abs_pct = max(abs(pct_change), na.rm = TRUE),
            .groups = "drop") %>%
  arrange(desc(max_abs_pct))

## Print top 10 most influential for each output
cat("--- TOP 10 MOST SENSITIVE PARAMETERS (per output) ---\n\n")
for (out in OUTPUT_NAMES) {
  top10 <- oat_summary %>%
    filter(output == out) %>%
    head(10)
  cat(sprintf("  %s:\n", out))
  for (r in seq_len(nrow(top10))) {
    cat(sprintf("    %-18s (%s)  max |Δ| = %.1f%%\n",
                top10$param[r], top10$module[r], top10$max_abs_pct[r]))
  }
  cat("\n")
}


## --- TORNADO PLOT (DA_neuron) ---
tornado_da <- oat_results %>%
  filter(output == "DA_neuron") %>%
  group_by(param, module) %>%
  summarise(
    low  = pct_change[direction == "-50%"],
    high = pct_change[direction == "+50%"],
    span = abs(high - low),
    .groups = "drop"
  ) %>%
  arrange(desc(span)) %>%
  head(20)

tornado_da$param <- factor(tornado_da$param,
                           levels = rev(tornado_da$param))

tornado_long <- tornado_da %>%
  select(param, module, low, high) %>%
  pivot_longer(cols = c("low", "high"),
               names_to = "direction", values_to = "pct_change")

p_tornado <- ggplot(tornado_long,
                    aes(x = pct_change, y = param, fill = direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("low" = "#d62728", "high" = "#1f77b4"),
                    labels = c("low" = "-50%", "high" = "+50%")) +
  labs(title = "OAT Tornado: DA Neuron Survival (Day 35)",
       subtitle = "Top 20 parameters by sensitivity | ±50% perturbation",
       x = "% Change in DA_neuron", y = NULL, fill = "Perturbation") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_tornado_DA.png", p_tornado,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: SA_tornado_DA.png\n")


## --- TORNADO PLOT (Motor_score) ---
tornado_motor <- oat_results %>%
  filter(output == "Motor_score") %>%
  group_by(param, module) %>%
  summarise(
    low  = pct_change[direction == "-50%"],
    high = pct_change[direction == "+50%"],
    span = abs(high - low),
    .groups = "drop"
  ) %>%
  arrange(desc(span)) %>%
  head(20)

tornado_motor$param <- factor(tornado_motor$param,
                              levels = rev(tornado_motor$param))

tornado_m_long <- tornado_motor %>%
  select(param, module, low, high) %>%
  pivot_longer(cols = c("low", "high"),
               names_to = "direction", values_to = "pct_change")

p_tornado_m <- ggplot(tornado_m_long,
                      aes(x = pct_change, y = param, fill = direction)) +
  geom_col(position = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("low" = "#d62728", "high" = "#1f77b4"),
                    labels = c("low" = "-50%", "high" = "+50%")) +
  labs(title = "OAT Tornado: Motor Score (Day 35)",
       subtitle = "Top 20 parameters by sensitivity | ±50% perturbation",
       x = "% Change in Motor_score", y = NULL, fill = "Perturbation") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_tornado_Motor.png", p_tornado_m,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("Saved: SA_tornado_Motor.png\n\n")


##-----------------------------------------------------------------------------
## 2. MORRIS SCREENING
##-----------------------------------------------------------------------------

cat("================================================================\n")
cat("  STAGE 2: MORRIS ELEMENTARY EFFECTS SCREENING\n")
cat("================================================================\n\n")

## Morris method: sample r trajectories across p-level grid
## Each trajectory has (k+1) points, varying one param at a time
## Total evaluations: r × (k+1)

morris_r <- 20      # trajectories
morris_p <- 4       # levels
morris_delta <- morris_p / (2 * (morris_p - 1))  # step size in [0,1]

## Parameter bounds: [0.5×baseline, 1.5×baseline]
p_lower <- param_space$baseline * 0.5
p_upper <- param_space$baseline * 1.5

## Generate Morris trajectories
set.seed(42)
morris_ee <- array(NA, dim = c(morris_r, n_params, length(OUTPUT_NAMES)))

cat(sprintf("  Morris: %d trajectories × %d params = %d model evaluations\n",
            morris_r, n_params, morris_r * (n_params + 1)))

for (traj in seq_len(morris_r)) {

  ## Random starting point on the grid
  x_start <- sapply(seq_len(n_params), function(i) {
    grid_vals <- seq(0, 1 - morris_delta, length.out = morris_p)
    sample(grid_vals, 1)
  })

  ## Map to parameter space
  to_real <- function(x01) p_lower + x01 * (p_upper - p_lower)

  ## Evaluate at starting point
  pv0 <- to_real(x_start)
  y0  <- run_sa(pv0)

  ## Random permutation of parameter order
  perm <- sample(n_params)

  x_curr <- x_start
  y_curr <- y0

  for (step in seq_along(perm)) {
    i <- perm[step]

    ## Perturb parameter i by ±delta
    x_new <- x_curr
    if (x_curr[i] + morris_delta <= 1) {
      x_new[i] <- x_curr[i] + morris_delta
    } else {
      x_new[i] <- x_curr[i] - morris_delta
    }

    pv_new <- to_real(x_new)
    y_new  <- run_sa(pv_new)

    ## Elementary effect for parameter i
    dx <- x_new[i] - x_curr[i]
    for (j in seq_along(OUTPUT_NAMES)) {
      if (!is.na(y_new[j]) && !is.na(y_curr[j]) && dx != 0) {
        morris_ee[traj, i, j] <- (y_new[j] - y_curr[j]) / dx
      }
    }

    x_curr <- x_new
    y_curr <- y_new
  }

  if (traj %% 5 == 0) {
    cat(sprintf("  Morris: %d/%d trajectories done\n", traj, morris_r))
  }
}

## Compute Morris statistics: μ* (mean of |EE|) and σ (sd of EE)
morris_stats <- data.frame()
for (i in seq_len(n_params)) {
  for (j in seq_along(OUTPUT_NAMES)) {
    ee_vals <- morris_ee[, i, j]
    ee_vals <- ee_vals[!is.na(ee_vals)]
    if (length(ee_vals) > 1) {
      morris_stats <- rbind(morris_stats, data.frame(
        param  = param_space$name[i],
        module = param_space$module[i],
        output = OUTPUT_NAMES[j],
        mu_star = mean(abs(ee_vals)),
        sigma   = sd(ee_vals),
        stringsAsFactors = FALSE
      ))
    }
  }
}

## Print top 10 for DA_neuron
cat("\n--- MORRIS: TOP 10 PARAMETERS FOR DA_NEURON ---\n")
cat("  (μ* = overall importance, σ = nonlinearity/interactions)\n\n")
morris_da <- morris_stats %>%
  filter(output == "DA_neuron") %>%
  arrange(desc(mu_star)) %>%
  head(10)
cat(sprintf("  %-18s %-10s %8s %8s\n", "Parameter", "Module", "μ*", "σ"))
cat(sprintf("  %-18s %-10s %8s %8s\n", "---------", "------", "--", "-"))
for (r in seq_len(nrow(morris_da))) {
  cat(sprintf("  %-18s %-10s %8.4f %8.4f\n",
              morris_da$param[r], morris_da$module[r],
              morris_da$mu_star[r], morris_da$sigma[r]))
}

cat("\n--- MORRIS: TOP 10 PARAMETERS FOR MOTOR_SCORE ---\n\n")
morris_motor <- morris_stats %>%
  filter(output == "Motor_score") %>%
  arrange(desc(mu_star)) %>%
  head(10)
cat(sprintf("  %-18s %-10s %8s %8s\n", "Parameter", "Module", "μ*", "σ"))
cat(sprintf("  %-18s %-10s %8s %8s\n", "---------", "------", "--", "-"))
for (r in seq_len(nrow(morris_motor))) {
  cat(sprintf("  %-18s %-10s %8.4f %8.4f\n",
              morris_motor$param[r], morris_motor$module[r],
              morris_motor$mu_star[r], morris_motor$sigma[r]))
}

## --- Morris μ* vs σ scatter plot ---
morris_da_all <- morris_stats %>% filter(output == "DA_neuron")

p_morris <- ggplot(morris_da_all, aes(x = mu_star, y = sigma)) +
  geom_point(aes(color = module), size = 3) +
  geom_text(aes(label = param), size = 2.5, vjust = -0.8, hjust = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Morris Screening: DA Neuron Survival",
       subtitle = "μ* = importance | σ = nonlinearity/interactions | above diagonal = interactive",
       x = "μ* (mean of |elementary effects|)",
       y = "σ (std dev of elementary effects)",
       color = "Module") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_morris_DA.png", p_morris,
       width = 11, height = 7, dpi = 300, bg = "white")
cat("\nSaved: SA_morris_DA.png\n")

p_morris_m <- ggplot(
  morris_stats %>% filter(output == "Motor_score"),
  aes(x = mu_star, y = sigma)) +
  geom_point(aes(color = module), size = 3) +
  geom_text(aes(label = param), size = 2.5, vjust = -0.8, hjust = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Morris Screening: Motor Score",
       subtitle = "μ* = importance | σ = nonlinearity/interactions",
       x = "μ* (mean of |elementary effects|)",
       y = "σ (std dev of elementary effects)",
       color = "Module") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_morris_Motor.png", p_morris_m,
       width = 11, height = 7, dpi = 300, bg = "white")
cat("Saved: SA_morris_Motor.png\n\n")


##-----------------------------------------------------------------------------
## 3. SOBOL INDICES (Top 10 parameters from Morris)
##-----------------------------------------------------------------------------

cat("================================================================\n")
cat("  STAGE 3: SOBOL VARIANCE-BASED ANALYSIS (Top 10 params)\n")
cat("================================================================\n\n")

## Select top 10 parameters across all outputs by max μ*
top_params <- morris_stats %>%
  group_by(param, module) %>%
  summarise(max_mu = max(mu_star, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_mu)) %>%
  head(10)

cat("Top 10 parameters for Sobol analysis:\n")
for (r in seq_len(nrow(top_params))) {
  cat(sprintf("  %2d. %-18s (%s)  max μ* = %.4f\n",
              r, top_params$param[r], top_params$module[r], top_params$max_mu[r]))
}

sobol_params <- top_params$param
sobol_idx    <- match(sobol_params, param_space$name)
n_sobol      <- length(sobol_params)

## Sobol via Saltelli sampling (2 matrices × N × (k+2) evaluations)
## N=500 gives 500×(10+2)=6000 runs — feasible at ~3s each = 5 hours
## Use N=200 for tractability: 200×12=2400 runs ≈ 2 hours
## Actually, let me use a simpler bootstrap approach: Latin Hypercube

sobol_N <- 256  # base sample size
total_evals <- sobol_N * (n_sobol + 2)
cat(sprintf("\nSobol: N=%d, k=%d, total evaluations=%d\n", sobol_N, n_sobol, total_evals))
cat("Estimated time: ~", round(total_evals * 3 / 60, 0), "minutes\n\n")

## Generate two independent uniform samples (Saltelli method)
set.seed(123)
A <- matrix(runif(sobol_N * n_sobol), nrow = sobol_N, ncol = n_sobol)
B <- matrix(runif(sobol_N * n_sobol), nrow = sobol_N, ncol = n_sobol)

## Scale to parameter bounds
scale_sample <- function(u, i_sobol) {
  idx <- sobol_idx[i_sobol]
  p_lower[idx] + u * (p_upper[idx] - p_lower[idx])
}

## Build full parameter vector from Sobol subset
make_parms <- function(sobol_vals) {
  pv <- param_space$baseline
  for (s in seq_len(n_sobol)) {
    pv[sobol_idx[s]] <- sobol_vals[s]
  }
  pv
}

## Evaluate model for a sample matrix
eval_matrix <- function(mat, label = "") {
  n <- nrow(mat)
  results <- matrix(NA, nrow = n, ncol = length(OUTPUT_NAMES))
  for (row in seq_len(n)) {
    sobol_vals <- sapply(seq_len(n_sobol), function(s) {
      scale_sample(mat[row, s], s)
    })
    pv <- make_parms(sobol_vals)
    results[row, ] <- run_sa(pv)
    if (row %% 100 == 0) {
      cat(sprintf("  Sobol %s: %d/%d evaluations\n", label, row, n))
    }
  }
  results
}

## Evaluate A and B matrices
cat("  Evaluating matrix A...\n")
yA <- eval_matrix(A, "A")
cat("  Evaluating matrix B...\n")
yB <- eval_matrix(B, "B")

## Evaluate AB_i matrices (A with column i from B)
yABi <- list()
for (i in seq_len(n_sobol)) {
  ABi <- A
  ABi[, i] <- B[, i]
  cat(sprintf("  Evaluating AB_%d (%s)...\n", i, sobol_params[i]))
  yABi[[i]] <- eval_matrix(ABi, sprintf("AB_%d", i))
}

## Compute Sobol indices (Jansen estimator)
sobol_results <- data.frame()

for (j in seq_along(OUTPUT_NAMES)) {
  fA  <- yA[, j]
  fB  <- yB[, j]

  ## Remove NAs
  valid_AB <- !is.na(fA) & !is.na(fB)
  fA_v <- fA[valid_AB]
  fB_v <- fB[valid_AB]
  n_valid <- sum(valid_AB)

  if (n_valid < 10) next

  f0  <- mean(c(fA_v, fB_v))
  VY  <- var(c(fA_v, fB_v))

  if (VY < 1e-12) next

  for (i in seq_len(n_sobol)) {
    fABi <- yABi[[i]][, j]
    valid_i <- valid_AB & !is.na(fABi)
    fA_i   <- fA[valid_i]
    fABi_i <- fABi[valid_i]
    fB_i   <- fB[valid_i]
    n_i    <- sum(valid_i)

    if (n_i < 10) next

    ## First-order (Saltelli 2010)
    Si <- mean(fB_i * (fABi_i - fA_i)) / VY

    ## Total-order (Jansen 1999)
    STi <- mean((fA_i - fABi_i)^2) / (2 * VY)

    sobol_results <- rbind(sobol_results, data.frame(
      param  = sobol_params[i],
      module = param_space$module[sobol_idx[i]],
      output = OUTPUT_NAMES[j],
      S1     = max(0, Si),
      ST     = max(0, STi),
      interaction = max(0, STi - max(0, Si)),
      stringsAsFactors = FALSE
    ))
  }
}

## Print Sobol results for DA_neuron
cat("\n--- SOBOL INDICES: DA_NEURON ---\n")
cat("  S1 = first-order (main effect), ST = total (including interactions)\n\n")
sobol_da <- sobol_results %>%
  filter(output == "DA_neuron") %>%
  arrange(desc(ST))
cat(sprintf("  %-18s %-10s %6s %6s %6s\n",
            "Parameter", "Module", "S1", "ST", "Inter"))
cat(sprintf("  %-18s %-10s %6s %6s %6s\n",
            "---------", "------", "--", "--", "-----"))
for (r in seq_len(nrow(sobol_da))) {
  cat(sprintf("  %-18s %-10s %6.3f %6.3f %6.3f\n",
              sobol_da$param[r], sobol_da$module[r],
              sobol_da$S1[r], sobol_da$ST[r], sobol_da$interaction[r]))
}

cat("\n--- SOBOL INDICES: MOTOR_SCORE ---\n\n")
sobol_motor <- sobol_results %>%
  filter(output == "Motor_score") %>%
  arrange(desc(ST))
cat(sprintf("  %-18s %-10s %6s %6s %6s\n",
            "Parameter", "Module", "S1", "ST", "Inter"))
cat(sprintf("  %-18s %-10s %6s %6s %6s\n",
            "---------", "------", "--", "--", "-----"))
for (r in seq_len(nrow(sobol_motor))) {
  cat(sprintf("  %-18s %-10s %6.3f %6.3f %6.3f\n",
              sobol_motor$param[r], sobol_motor$module[r],
              sobol_motor$S1[r], sobol_motor$ST[r], sobol_motor$interaction[r]))
}


## --- Sobol bar plot ---
sobol_da_plot <- sobol_results %>%
  filter(output == "DA_neuron") %>%
  arrange(desc(ST)) %>%
  select(param, S1, interaction) %>%
  pivot_longer(cols = c("S1", "interaction"),
               names_to = "type", values_to = "index")
sobol_da_plot$param <- factor(sobol_da_plot$param,
                              levels = rev(unique(sobol_da_plot$param)))
sobol_da_plot$type <- factor(sobol_da_plot$type,
                             levels = c("interaction", "S1"),
                             labels = c("Interactions", "First-order"))

p_sobol <- ggplot(sobol_da_plot, aes(x = index, y = param, fill = type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("First-order" = "#1f77b4",
                               "Interactions" = "#ff7f0e")) +
  labs(title = "Sobol Indices: DA Neuron Survival",
       subtitle = "First-order (S1) + interaction effects = total (ST)",
       x = "Sobol Index", y = NULL, fill = "Contribution") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_sobol_DA.png", p_sobol,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("\nSaved: SA_sobol_DA.png\n")

sobol_motor_plot <- sobol_results %>%
  filter(output == "Motor_score") %>%
  arrange(desc(ST)) %>%
  select(param, S1, interaction) %>%
  pivot_longer(cols = c("S1", "interaction"),
               names_to = "type", values_to = "index")
sobol_motor_plot$param <- factor(sobol_motor_plot$param,
                                 levels = rev(unique(sobol_motor_plot$param)))
sobol_motor_plot$type <- factor(sobol_motor_plot$type,
                                levels = c("interaction", "S1"),
                                labels = c("Interactions", "First-order"))

p_sobol_m <- ggplot(sobol_motor_plot, aes(x = index, y = param, fill = type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("First-order" = "#1f77b4",
                               "Interactions" = "#ff7f0e")) +
  labs(title = "Sobol Indices: Motor Score",
       subtitle = "First-order (S1) + interaction effects = total (ST)",
       x = "Sobol Index", y = NULL, fill = "Contribution") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"))
ggsave("SA_sobol_Motor.png", p_sobol_m,
       width = 10, height = 6, dpi = 300, bg = "white")
cat("Saved: SA_sobol_Motor.png\n\n")


##-----------------------------------------------------------------------------
## SUMMARY
##-----------------------------------------------------------------------------

cat("================================================================\n")
cat("  SENSITIVITY ANALYSIS COMPLETE\n")
cat("================================================================\n\n")
cat("  Plots generated:\n")
cat("    SA_tornado_DA.png      — OAT tornado for DA neuron\n")
cat("    SA_tornado_Motor.png   — OAT tornado for Motor score\n")
cat("    SA_morris_DA.png       — Morris μ* vs σ for DA neuron\n")
cat("    SA_morris_Motor.png    — Morris μ* vs σ for Motor score\n")
cat("    SA_sobol_DA.png        — Sobol indices for DA neuron\n")
cat("    SA_sobol_Motor.png     — Sobol indices for Motor score\n")
cat("\n")

cat("  Key findings:\n")
cat("  Top parameters by Morris μ* for DA_neuron:\n")
top5_da <- morris_stats %>%
  filter(output == "DA_neuron") %>%
  arrange(desc(mu_star)) %>%
  head(5)
for (r in seq_len(nrow(top5_da))) {
  cat(sprintf("    %d. %s (%s)\n", r, top5_da$param[r], top5_da$module[r]))
}

cat("\n  Top parameters by Morris μ* for Motor_score:\n")
top5_m <- morris_stats %>%
  filter(output == "Motor_score") %>%
  arrange(desc(mu_star)) %>%
  head(5)
for (r in seq_len(nrow(top5_m))) {
  cat(sprintf("    %d. %s (%s)\n", r, top5_m$param[r], top5_m$module[r]))
}

cat("\n================================================================\n")
