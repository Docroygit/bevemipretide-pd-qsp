##=============================================================================
## FAST PROFILE LIKELIHOOD — FigS5 generation
## Optimised re-implementation of Section 9 from Identifiability Analysis.R:
##   1. Nelder-Mead (derivative-free) instead of L-BFGS-B — eliminates the
##      5-parameter numerical-gradient overhead (~6 extra ODE solves/step)
##   2. Pre-compute healthy scenario once (CI_max=1, alpha_clear=1 — invariant
##      to calibrated parameter perturbation) and skip drug rescue scenario
##      (not used in SSE targets) → 2 ODE solves per eval instead of 4
##   3. 7 grid points (−35% to +35%) instead of 11
##   Combined: ~6-10× faster than the L-BFGS-B original.
##=============================================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(parallel)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

out_dir <- "../New publication plots"

## ---- Calibrated-parameter nominals (same as main script) ----
nom6 <- c(CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35,
          K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4.0)
cal_names <- names(nom6)
DAYS <- 35

## ---- Pre-compute healthy baseline (invariant to calibrated params) ----
cat("Pre-computing healthy baseline (once)...\n")
t_start <- Sys.time()
res_h <- run_integrated(dose_mg_kg = 0, duration_days = DAYS, dt = 1.0,
                         CI_max = 1.0, alpha_clear = 1.0, quiet = TRUE)
h_end  <- res_h[nrow(res_h), ]
h_CL_n <- h_end$CL_n
h_mROS <- h_end$mROS
cat(sprintf("  Healthy: CL_n=%.4f, mROS=%.4f\n", h_CL_n, h_mROS))

## ---- Compute nominal SSE targets from calibrated fit ----
cat("Computing nominal disease + KO outputs...\n")
res_d <- run_integrated(dose_mg_kg = 0, duration_days = DAYS, dt = 1.0,
                         CI_max = 0.74, alpha_clear = 0.25,
                         k_impair_val = 0.35, K_death_val = 0.25,
                         k_damage_mPTP_val = 4.0, n_death_val = 4.0, quiet = TRUE)
res_ko <- run_integrated(dose_mg_kg = 0, duration_days = DAYS, dt = 1.0,
                          CI_max = 0.74, alpha_clear = 0.25,
                          k_impair_val = 0.35, K_death_val = 0.25,
                          k_damage_mPTP_val = 4.0, n_death_val = 4.0,
                          parms_override = list(k_act_MG = 0), quiet = TRUE)
d_end  <- res_d[nrow(res_d), ]
ko_end <- res_ko[nrow(res_ko), ]

target_vals <- c(
  CI_act   = d_end$CI_activity,
  CL_drop  = (1 - d_end$CL_n / h_CL_n) * 100,
  mROS_rat = d_end$mROS / h_mROS * 100,
  DA_loss  = (1 - d_end$DA_neuron) * 100,
  aSyn_red = (1 - ko_end$aSyn_olig / d_end$aSyn_olig) * 100
)
target_sigma <- c(
  CI_act   = (0.58 - 0.40) / 2,
  CL_drop  = (35 - 15)     / 2,
  mROS_rat = (250 - 150)   / 2,
  DA_loss  = (45 - 25)     / 2,
  aSyn_red = (50 - 15)     / 2
)
cat("Nominal SSE targets:\n"); print(round(target_vals, 4))

## ---- Grid: 7 points, −35% to +35% ----
profile_frac <- seq(0.65, 1.35, length.out = 7)
tasks <- expand.grid(pname = cal_names, gi = seq_along(profile_frac),
                      stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("\nProfile likelihood: %d params x %d grid = %d Nelder-Mead optimisations\n",
            length(cal_names), length(profile_frac), n_tasks))

## ---- Parallel execution ----
n_cores <- max(1, detectCores() - 2)
cat(sprintf("Launching parallel cluster: %d workers\n", n_cores))

cl <- makeCluster(n_cores, rscript_args = "--vanilla")
invisible(clusterEvalQ(cl, {
  library(deSolve)
  SOURCED_FOR_FUNCTIONS <- TRUE
  source("Integrated Model.R", local = FALSE)
}))
clusterExport(cl, c("nom6", "cal_names", "target_vals", "target_sigma",
                      "h_CL_n", "h_mROS", "profile_frac", "tasks", "DAYS"))

t_pl <- Sys.time()
results <- parLapply(cl, seq_len(n_tasks), function(k) {
  pname <- tasks$pname[k]
  gi    <- tasks$gi[k]
  gval  <- nom6[[pname]] * profile_frac[gi]
  free_names <- setdiff(cal_names, pname)

  obj <- function(free_vec) {
    full <- nom6
    full[pname]      <- gval
    full[free_names] <- free_vec

    res_d2 <- tryCatch(
      run_integrated(dose_mg_kg = 0, duration_days = DAYS, dt = 1.0,
                      CI_max = full[["CI_max"]],
                      alpha_clear = full[["alpha_clear"]],
                      k_impair_val = full[["k_impair"]],
                      K_death_val = full[["K_death"]],
                      k_damage_mPTP_val = full[["k_damage_mPTP"]],
                      n_death_val = full[["n_death"]], quiet = TRUE),
      error = function(e) NULL)
    if (is.null(res_d2)) return(1e6)

    res_ko2 <- tryCatch(
      run_integrated(dose_mg_kg = 0, duration_days = DAYS, dt = 1.0,
                      CI_max = full[["CI_max"]],
                      alpha_clear = full[["alpha_clear"]],
                      k_impair_val = full[["k_impair"]],
                      K_death_val = full[["K_death"]],
                      k_damage_mPTP_val = full[["k_damage_mPTP"]],
                      n_death_val = full[["n_death"]],
                      parms_override = list(k_act_MG = 0), quiet = TRUE),
      error = function(e) NULL)
    if (is.null(res_ko2)) return(1e6)

    dd <- res_d2[nrow(res_d2), ]
    ko <- res_ko2[nrow(res_ko2), ]
    e <- c(CI_act   = dd$CI_activity,
            CL_drop  = (1 - dd$CL_n / h_CL_n) * 100,
            mROS_rat = dd$mROS / h_mROS * 100,
            DA_loss  = (1 - dd$DA_neuron) * 100,
            aSyn_red = (1 - ko$aSyn_olig / dd$aSyn_olig) * 100)
    r <- (e[names(target_vals)] - target_vals) / target_sigma
    sum(r^2)
  }

  opt <- optim(nom6[free_names], obj, method = "Nelder-Mead",
               control = list(maxit = 500, reltol = 1e-6))
  data.frame(param = pname, value = gval, frac = profile_frac[gi],
             sse = opt$value, convergence = opt$convergence)
})
stopCluster(cl)

elapsed <- as.numeric(difftime(Sys.time(), t_pl, units = "secs"))
cat(sprintf("Profile likelihood done: %.1f s (%.1f min)\n", elapsed, elapsed / 60))
cat(sprintf("Total wall time: %.1f s\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

## ---- Assemble results ----
profile_df <- do.call(rbind, results)
profile_df <- profile_df[order(profile_df$param, profile_df$frac), ]
profile_df <- profile_df %>%
  group_by(param) %>%
  mutate(delta_sse = sse - min(sse)) %>%
  ungroup()

CHI2 <- 3.84
profile_df$above_threshold <- profile_df$delta_sse > CHI2

verdict <- profile_df %>%
  group_by(param) %>%
  summarise(max_delta_sse = max(delta_sse),
            practically_identifiable = any(above_threshold), .groups = "drop")
cat("\nProfile likelihood verdict:\n")
print(as.data.frame(verdict))

n_converged <- sum(profile_df$convergence == 0)
cat(sprintf("Convergence: %d / %d optimisations converged\n", n_converged, n_tasks))

saveRDS(profile_df, "identifiability_profile_likelihood.rds")

## ---- Figure S5 ----
nom_df <- data.frame(param = cal_names, nom = as.numeric(nom6))
psS5 <- ggplot(profile_df, aes(x = value, y = delta_sse)) +
  geom_line(colour = "#1565C0", linewidth = 0.8) +
  geom_point(colour = "#1565C0", size = 1.5) +
  geom_hline(yintercept = CHI2, colour = "#B71C1C", linetype = "dashed") +
  geom_vline(data = nom_df, aes(xintercept = nom),
             colour = "grey40", linetype = "dotted") +
  facet_wrap(~ param, scales = "free_x", nrow = 2) +
  labs(title = "Profile likelihood: practical identifiability of calibrated parameters",
       subtitle = sprintf(
         "Dashed = 95%% chi-squared threshold (%.2f, 1 d.f.); dotted = calibrated nominal",
         CHI2),
       x = "Parameter value",
       y = expression(Delta * "SSE (re-optimised)")) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "FigS5_profile_likelihood.png"), psS5,
       width = 12, height = 7, dpi = 300, bg = "white")
cat("Saved FigS5_profile_likelihood.png\n")
