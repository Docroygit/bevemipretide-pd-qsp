##=============================================================================
## STRUCTURAL & PRACTICAL IDENTIFIABILITY ANALYSIS (R re-implementation)
## Bevemipretide PD QSP -- replaces a lost Python pipeline (Sep 2026 rebuild).
## Produces FigS4 (structural: FIM rank + eigenvalue spectrum + sensitivity
## heatmap), FigS5 (practical: profile likelihood of the 6 calibrated
## parameters), and FigS6 (collinearity + sensitivity-correlation
## diagnostics).
##
## Parameter set and classification counts are read LIVE from Table ST1
## (exported to ST1_classification.json by export_st1_classification.py) so
## this script can never silently drift from the manuscript's own parameter
## table again -- a prior hardcoded version (37/32/6/15 literature/estimated/
## calibrated/assumed) went stale after ST1 was corrected (Sep 2026: 4 PK
## parameters reclassified Literature->Assumed, 1 reclassified
## Literature->Estimated), producing a 90-parameter total that no longer
## matched its own classification bar chart, and an identifiability sweep
## that silently omitted k_brain_in after it was reclassified as Estimated.
##
## Methodology (Methods 2.7):
##  - Structural: central finite differences (Delta = 1%) on the calibrated +
##    estimated parameters (from live ST1), evaluated against 7 calibration
##    outputs (5 quantitative targets, the healthy DA baseline, and the drug
##    rescue endpoint). FIM = S'S; rank via eigenvalues.
##  - Practical: profile likelihood on the 6 calibrated parameters. Each is
##    swept +/-35% around its calibrated value while the other 5 are
##    re-optimised (L-BFGS-B) to minimise SSE against the model's own nominal
##    fit, weighted by each output's literature-derived acceptable half-range
##    (Table 3). Delta-SSE relative to the profile minimum is compared to the
##    95% chi-squared threshold (3.84, 1 d.f.).
##  - Collinearity index per Brun et al. (2001), gamma_K = 1/sqrt(min
##    eigenvalue of the parameter subset's normalised sensitivity correlation
##    matrix).
##=============================================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(jsonlite)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

out_dir <- "../New publication plots"

CAL <- list(K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4,
            CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35)
DAYS_CAL <- 35   # matches Table 3 / Results 3.1 calibration checkpoint

## -----------------------------------------------------------------------
## 1. Parameter set, read live from Table ST1's Classification column.
##    R-name mapping: ST1 shows "k_a"/"α_clear" (Greek stripped); the
##    internal parms vector uses "ka"/"alpha_clear".
## -----------------------------------------------------------------------
st1 <- jsonlite::fromJSON("ST1_classification.json")
name_map <- c("k_a" = "ka", "α_clear" = "alpha_clear")
remap <- function(v) unname(ifelse(v %in% names(name_map), name_map[v], v))

calibrated_params <- remap(st1$params_by_class$Calibrated)
estimated_params  <- remap(st1$params_by_class$Estimated)
all_id_params <- c(calibrated_params, estimated_params)
stopifnot(length(all_id_params) == st1$counts$Calibrated + st1$counts$Estimated)
cat(sprintf("Identifiability parameter set (live from ST1): %d calibrated + %d estimated = %d total\n",
            length(calibrated_params), length(estimated_params), length(all_id_params)))

nominal_p0 <- integrated_parameters(0, CAL$CI_max, CAL$alpha_clear, CAL$k_impair,
                                     CAL$K_death, CAL$k_damage_mPTP, CAL$n_death)
get_nominal <- function(pname) {
  if (pname %in% names(CAL)) return(CAL[[pname]])
  as.numeric(nominal_p0[pname])
}
nominal_vals <- setNames(sapply(all_id_params, get_nominal), all_id_params)

## -----------------------------------------------------------------------
## 2. Run the 4 scenarios needed for the 7 calibration outputs, given a
##    named override list applied to the *shared* parms vector, and given
##    (possibly perturbed) values for the 6 calibrated arguments.
## -----------------------------------------------------------------------
TC_DAYS <- c(7, 14, 21, 28, 35)
TC_VARS <- c("CI_activity", "CL_ratio", "mROS", "aSyn_olig", "DA_neuron", "Motor_score")

run_scenarios <- function(ovr = list(), cal = CAL, days = DAYS_CAL, want_tc = FALSE) {
  res_h <- run_integrated(dose_mg_kg = 0, duration_days = days, dt = 1.0,
                           CI_max = 1.0, alpha_clear = 1.0,
                           parms_override = ovr, quiet = TRUE)
  res_d <- run_integrated(dose_mg_kg = 0, duration_days = days, dt = 1.0,
                           CI_max = cal$CI_max, alpha_clear = cal$alpha_clear,
                           k_impair_val = cal$k_impair, K_death_val = cal$K_death,
                           k_damage_mPTP_val = cal$k_damage_mPTP,
                           n_death_val = cal$n_death,
                           parms_override = ovr, quiet = TRUE)
  ovr_ko <- ovr; ovr_ko[["k_act_MG"]] <- 0
  res_ko <- run_integrated(dose_mg_kg = 0, duration_days = days, dt = 1.0,
                            CI_max = cal$CI_max, alpha_clear = cal$alpha_clear,
                            k_impair_val = cal$k_impair, K_death_val = cal$K_death,
                            k_damage_mPTP_val = cal$k_damage_mPTP,
                            n_death_val = cal$n_death,
                            parms_override = ovr_ko, quiet = TRUE)
  res_rx <- run_integrated(dose_mg_kg = 5.0, duration_days = days, dt = 1.0,
                            CI_max = cal$CI_max, alpha_clear = cal$alpha_clear,
                            k_impair_val = cal$k_impair, K_death_val = cal$K_death,
                            k_damage_mPTP_val = cal$k_damage_mPTP,
                            n_death_val = cal$n_death,
                            parms_override = ovr, quiet = TRUE)
  h <- res_h[nrow(res_h), ]; dd <- res_d[nrow(res_d), ]
  ko <- res_ko[nrow(res_ko), ]; rx <- res_rx[nrow(res_rx), ]
  endpoint <- c(
    CI_act    = dd$CI_activity,
    mROS_rat  = dd$mROS / h$mROS * 100,
    CL_drop   = (1 - dd$CL_n / h$CL_n) * 100,
    aSyn_red  = (1 - ko$aSyn_olig / dd$aSyn_olig) * 100,
    DA_loss   = (1 - dd$DA_neuron) * 100,
    DA_healthy = h$DA_neuron,
    DA_rescue = rx$DA_neuron - dd$DA_neuron
  )
  if (!want_tc) return(endpoint)
  tc <- numeric(0)
  for (dday in TC_DAYS) {
    ridx <- which.min(abs(res_d$time_days - dday))
    row <- res_d[ridx, ]
    vals <- setNames(as.numeric(row[TC_VARS]), sprintf("%s_d%d", TC_VARS, dday))
    tc <- c(tc, vals)
  }
  list(endpoint = endpoint, tc = tc)
}

cat("Computing nominal outputs...\n")
y0_full <- run_scenarios(want_tc = TRUE)
y0 <- y0_full$endpoint
y0_tc <- y0_full$tc
print(round(y0, 4))
outputs <- names(y0)
tc_names <- names(y0_tc)

## -----------------------------------------------------------------------
## 3. Central finite-difference sensitivity matrix S (7 outputs x N params)
##    normalised: S_ij = d(ln y_i) / d(ln p_j)
## -----------------------------------------------------------------------
delta <- 0.01
S <- matrix(0, nrow = length(outputs), ncol = length(all_id_params),
            dimnames = list(outputs, all_id_params))

S_tc <- matrix(0, nrow = length(tc_names), ncol = length(all_id_params),
               dimnames = list(tc_names, all_id_params))

cat(sprintf("Running finite-difference perturbations (%d params x 2 directions, parallel)...\n",
            length(all_id_params)))
t_fd <- Sys.time()

n_cores <- max(1, parallel::detectCores() - 2)
cl <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl, {
  library(deSolve); SOURCED_FOR_FUNCTIONS <- TRUE
  source("Integrated Model.R", local = FALSE)
}))
parallel::clusterExport(cl, c("CAL", "run_scenarios", "TC_DAYS", "TC_VARS", "delta",
                               "all_id_params", "nominal_vals", "DAYS_CAL"))
cat("Parallel cluster:", n_cores, "workers\n")

fd_results <- parallel::parLapply(cl, seq_along(all_id_params), function(j) {
  pname <- all_id_params[j]
  p0 <- unname(nominal_vals[pname])
  cal_up <- CAL; cal_dn <- CAL
  ovr_up <- list(); ovr_dn <- list()
  if (pname %in% names(CAL)) {
    cal_up[[pname]] <- p0 * (1 + delta)
    cal_dn[[pname]] <- p0 * (1 - delta)
  } else {
    ovr_up[[pname]] <- p0 * (1 + delta)
    ovr_dn[[pname]] <- p0 * (1 - delta)
  }
  r_up <- run_scenarios(ovr = ovr_up, cal = cal_up, want_tc = TRUE)
  r_dn <- run_scenarios(ovr = ovr_dn, cal = cal_dn, want_tc = TRUE)
  list(endpoint_diff = r_up$endpoint - r_dn$endpoint, tc_diff = r_up$tc - r_dn$tc)
})
parallel::stopCluster(cl)

for (j in seq_along(all_id_params)) {
  for (i in seq_along(outputs)) {
    y_ref <- if (abs(y0[i]) > 1e-8) y0[i] else 1e-8
    S[i, j] <- (fd_results[[j]]$endpoint_diff[i] / (2 * delta)) / y_ref
  }
  for (i in seq_along(tc_names)) {
    y_ref <- if (abs(y0_tc[i]) > 1e-8) y0_tc[i] else 1e-8
    S_tc[i, j] <- (fd_results[[j]]$tc_diff[i] / (2 * delta)) / y_ref
  }
}
cat(sprintf("Finite-difference perturbations done (%.1f s elapsed)\n",
            as.numeric(difftime(Sys.time(), t_fd, units = "secs"))))
S[!is.finite(S)] <- 0
S_tc[!is.finite(S_tc)] <- 0
saveRDS(S, "identifiability_sensitivity_matrix.rds")
saveRDS(S_tc, "identifiability_sensitivity_matrix_tc.rds")
cat("Sensitivity matrices saved.\n")

## -----------------------------------------------------------------------
## 4. Structural identifiability: FIM = S'S, rank via eigen decomposition
## -----------------------------------------------------------------------
FIM <- t(S) %*% S
eig <- eigen(FIM, symmetric = TRUE)
ev <- sort(Re(eig$values), decreasing = TRUE)
ev[ev < 0] <- 0
rank_endpoint <- sum(ev > 1e-8 * max(ev))
cat(sprintf("\nEndpoint FIM rank: %d (of %d params, %d outputs)\n",
            rank_endpoint, length(all_id_params), length(outputs)))
cat("Top eigenvalues:", paste(sprintf("%.4g", head(ev, 10)), collapse = ", "), "\n")

FIM_tc <- t(S_tc) %*% S_tc
eig_tc <- eigen(FIM_tc, symmetric = TRUE)
ev_tc <- sort(Re(eig_tc$values), decreasing = TRUE)
ev_tc[ev_tc < 0] <- 0
rank_tc <- sum(ev_tc > 1e-8 * max(ev_tc))
cat(sprintf("Time-course FIM rank: %d (of %d params, %d time-sampled outputs)\n",
            rank_tc, length(all_id_params), length(tc_names)))

cond_number <- ev[1] / ev[rank_endpoint]
pct_leading <- 100 * ev[1] / sum(ev[1:rank_endpoint])
cat(sprintf("Condition number (within rank): %.3g\n", cond_number))
cat(sprintf("Leading eigenvalue as %% of retained information: %.1f%%\n", pct_leading))

## -----------------------------------------------------------------------
## 5. Parameter classification bar chart (Panel A) -- read live from Table
##    ST1 via ST1_classification.json (Section 1), so it can never drift
##    from the manuscript's own parameter table again.
## -----------------------------------------------------------------------
n_literature <- st1$counts[["Literature-derived"]]
n_calibrated <- st1$counts[["Calibrated"]]
n_estimated  <- st1$counts[["Estimated"]]
n_assumed    <- st1$counts[["Assumed"]]
n_total      <- st1$counts[["Total"]]
stopifnot(n_total == n_literature + n_calibrated + n_estimated + n_assumed)

class_df <- data.frame(
  category = factor(c("Literature\n(fixed)", "Estimated\n(calibrated)",
                       "Estimated\n(structural)", "Assumed\n(defaults)"),
                     levels = rev(c("Literature\n(fixed)", "Estimated\n(calibrated)",
                                    "Estimated\n(structural)", "Assumed\n(defaults)"))),
  n = c(n_literature, n_calibrated, n_estimated, n_assumed))

psA <- ggplot(class_df, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n), hjust = -0.3, fontface = "bold", size = 5) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c("Literature\n(fixed)" = "#2E7D32",
                                "Estimated\n(calibrated)" = "#C62828",
                                "Estimated\n(structural)" = "#1565C0",
                                "Assumed\n(defaults)" = "#F57C00")) +
  scale_y_continuous(limits = c(0, max(class_df$n) * 1.18)) +
  labs(title = sprintf("(A) Parameter classification (N = %d)", n_total),
       x = NULL, y = "Number of parameters") +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12))

## -----------------------------------------------------------------------
## 6. FIM eigenvalue spectrum plot (Panel B)
## -----------------------------------------------------------------------
spec_df <- data.frame(idx = seq_along(ev), eigenvalue = pmax(ev, 1e-13),
                       grp = ifelse(seq_along(ev) <= rank_endpoint,
                                    ifelse(seq_along(ev) <= 2, "top2", "rank"), "null"))
psB <- ggplot(spec_df, aes(x = idx, y = eigenvalue, fill = grp)) +
  geom_col(width = 0.7) +
  scale_y_log10() +
  scale_fill_manual(values = c("top2" = "#1565C0", "rank" = "#F9A825", "null" = "#B71C1C"),
                     guide = "none") +
  geom_hline(yintercept = 0.01 * ev[1], colour = "#B71C1C", linetype = "dashed") +
  annotate("text", x = length(ev) * 0.75, y = 0.01 * ev[1] * 3, label = "1% of max",
           colour = "#B71C1C", size = 3) +
  annotate("label", x = length(ev) * 0.6, y = ev[1],
           label = sprintf("Endpoint: rank = %d\nTime-course: rank = %d", rank_endpoint, rank_tc),
           size = 3.2, fill = "#FFFDE7") +
  labs(title = sprintf("(B) FIM eigenvalue spectrum (rank = %d)", rank_endpoint),
       x = "Eigenvalue index", y = "FIM eigenvalue (log scale)")

## -----------------------------------------------------------------------
## 7. Normalised sensitivity heatmap (Panel C, top 20 by mean |S|)
## -----------------------------------------------------------------------
mean_abs_S <- colMeans(abs(S))
top20 <- names(sort(mean_abs_S, decreasing = TRUE))[1:20]
S_top <- S[, top20]
heat_df <- as.data.frame(S_top) %>%
  mutate(output = rownames(S_top)) %>%
  pivot_longer(-output, names_to = "param", values_to = "sens")
heat_df$param <- factor(heat_df$param, levels = top20)
heat_df$output <- factor(heat_df$output, levels = rev(outputs))
param_labels <- ifelse(top20 %in% calibrated_params, paste0(top20, "*"), top20)
names(param_labels) <- top20

psC <- ggplot(heat_df, aes(x = param, y = output, fill = sens)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#1565C0", mid = "grey95", high = "#B71C1C", midpoint = 0,
                        name = "Normalised\nsensitivity") +
  scale_x_discrete(labels = param_labels) +
  labs(title = "(C) Normalised sensitivity matrix (top 20, * = calibrated)",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9))

figS4_top <- arrangeGrob(psA, psB, ncol = 2)
figS4 <- arrangeGrob(figS4_top, psC, nrow = 2, heights = c(1, 1.1))
ggsave(file.path(out_dir, "FigS4_structural_identifiability.png"), figS4,
       width = 13, height = 11, dpi = 300, bg = "white")
cat("Saved FigS4_structural_identifiability.png\n")

## -----------------------------------------------------------------------
## 8. Collinearity index (Brun et al. 2001) for all pairs among the 12 most
##    sensitive parameters (endpoint S), plus the 6 calibrated parameters --
##    matches the parameter set actually discussed in the manuscript text.
##    gamma_K = 1 / sqrt(min eigenvalue of the column-normalised correlation
##    matrix restricted to subset K).
## -----------------------------------------------------------------------
col_norm <- function(M) {
  scale(M, center = FALSE, scale = sqrt(colSums(M^2)))
}
Sn <- col_norm(S)   # normalise each parameter's sensitivity column to unit length

collin_set <- union(calibrated_params, top20[1:12])
collin_set <- intersect(collin_set, colnames(Sn))  # preserve order, dedupe

gamma_pair <- function(p1, p2) {
  sub <- Sn[, c(p1, p2), drop = FALSE]
  Ksub <- t(sub) %*% sub
  ev_sub <- eigen(Ksub, symmetric = TRUE, only.values = TRUE)$values
  1 / sqrt(max(min(Re(ev_sub)), 1e-12))
}

pairs <- t(combn(collin_set, 2))
collin_df <- data.frame(p1 = pairs[, 1], p2 = pairs[, 2],
                         gamma = apply(pairs, 1, function(pr) gamma_pair(pr[1], pr[2])))
collin_df <- collin_df[order(-collin_df$gamma), ]
cat("\nTop 10 collinearity pairs (of parameters actually plotted):\n")
print(head(collin_df, 10))

## count over ALL pairs among the full identifiability parameter set (for the
## "gamma>20 out of N pairs" text stat)
all_pairs <- t(combn(all_id_params, 2))
gamma_all <- apply(all_pairs, 1, function(pr) gamma_pair(pr[1], pr[2]))
n_high <- sum(gamma_all > 20)
cat(sprintf("\nPairs with gamma > 20 (all %d params, %d pairs): %d\n",
            length(all_id_params), nrow(all_pairs), n_high))
top_all_idx <- order(-gamma_all)[1:5]
cat("Top 5 pairs overall:\n")
print(data.frame(p1 = all_pairs[top_all_idx, 1], p2 = all_pairs[top_all_idx, 2],
                  gamma = gamma_all[top_all_idx]))

## Panel A: pairwise collinearity heatmap (only the plotted subset)
mat_df <- expand.grid(p1 = collin_set, p2 = collin_set, stringsAsFactors = FALSE)
mat_df$gamma <- mapply(function(a, b) if (a == b) NA_real_ else gamma_pair(a, b),
                        mat_df$p1, mat_df$p2)
mat_df$p1 <- factor(mat_df$p1, levels = collin_set)
mat_df$p2 <- factor(mat_df$p2, levels = rev(collin_set))
bold_labels <- setNames(ifelse(collin_set %in% calibrated_params,
                                paste0("**", collin_set, "**"), collin_set), collin_set)

psS6a <- ggplot(mat_df, aes(x = p1, y = p2, fill = pmin(gamma, 50))) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "#FFF9C4", high = "#B71C1C", na.value = "white",
                       limits = c(0, 50), name = "Collinearity\nindex") +
  scale_x_discrete(labels = collin_set) +
  scale_y_discrete(labels = rev(collin_set)) +
  labs(title = "(A) Pairwise collinearity (bold axis labels = calibrated)", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
          face = ifelse(collin_set %in% calibrated_params, "bold", "plain")),
        axis.text.y = element_text(
          face = ifelse(rev(collin_set) %in% calibrated_params, "bold", "plain")))

## Panel B: sensitivity correlation matrix (cosine similarity of columns)
corr_mat <- cor(S[, collin_set])
corr_df <- as.data.frame(as.table(corr_mat))
names(corr_df) <- c("p1", "p2", "cosine")
corr_df$p1 <- factor(corr_df$p1, levels = collin_set)
corr_df$p2 <- factor(corr_df$p2, levels = rev(collin_set))

psS6b <- ggplot(corr_df, aes(x = p1, y = p2, fill = cosine)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#B71C1C", midpoint = 0,
                        limits = c(-1, 1), name = "Cosine\nsimilarity") +
  labs(title = "(B) Sensitivity correlation matrix", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figS6 <- arrangeGrob(psS6a, psS6b, ncol = 2)
ggsave(file.path(out_dir, "FigS6_collinearity_diagnostics.png"), figS6,
       width = 16, height = 7.5, dpi = 300, bg = "white")
cat("Saved FigS6_collinearity_diagnostics.png\n")

## -----------------------------------------------------------------------
## 9. Practical identifiability: profile likelihood (Figure S5)
##    For each of the 6 calibrated parameters, fix its value across a grid
##    spanning +/-35% of the nominal calibration, re-optimise the remaining
##    5 via L-BFGS-B to minimise SSE against the model's own nominal fit
##    (weighted by each output's literature-derived acceptable half-range,
##    Table 3), and report Delta-SSE relative to the profile minimum against
##    the 95% chi-squared threshold (Delta-chi^2 = 3.84, 1 d.f.).
## -----------------------------------------------------------------------
nom6 <- c(CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear,
          k_impair = CAL$k_impair, K_death = CAL$K_death,
          k_damage_mPTP = CAL$k_damage_mPTP, n_death = CAL$n_death)
cal_names <- names(nom6)

target_vals  <- c(CI_act = unname(y0[["CI_act"]]), CL_drop = unname(y0[["CL_drop"]]),
                   mROS_rat = unname(y0[["mROS_rat"]]), DA_loss = unname(y0[["DA_loss"]]),
                   aSyn_red = unname(y0[["aSyn_red"]]))
target_sigma <- c(CI_act = (0.58 - 0.40) / 2, CL_drop = (35 - 15) / 2,
                   mROS_rat = (250 - 150) / 2, DA_loss = (45 - 25) / 2,
                   aSyn_red = (50 - 15) / 2)

sse_of <- function(cal6) {
  cal <- list(CI_max = cal6[["CI_max"]], alpha_clear = cal6[["alpha_clear"]],
              k_impair = cal6[["k_impair"]], K_death = cal6[["K_death"]],
              k_damage_mPTP = cal6[["k_damage_mPTP"]], n_death = cal6[["n_death"]])
  e <- tryCatch(run_scenarios(cal = cal, want_tc = FALSE), error = function(e) NULL)
  if (is.null(e)) return(1e6)
  r <- (e[names(target_vals)] - target_vals) / target_sigma
  sum(r^2)
}

profile_frac  <- seq(0.65, 1.35, length.out = 11)  # -35% to +35%, nominal exact at index 6
profile_lower <- nom6 * 0.65
profile_upper <- nom6 * 1.35

profile_tasks <- expand.grid(pname = cal_names, gi = seq_along(profile_frac),
                              stringsAsFactors = FALSE)

cat(sprintf("\nRunning profile likelihood (%d params x %d grid points, parallel)...\n",
            length(cal_names), length(profile_frac)))
t_pl <- Sys.time()
cl2 <- parallel::makeCluster(n_cores)
invisible(parallel::clusterEvalQ(cl2, {
  library(deSolve); SOURCED_FOR_FUNCTIONS <- TRUE
  source("Integrated Model.R", local = FALSE)
}))
parallel::clusterExport(cl2, c("nom6", "cal_names", "target_vals", "target_sigma",
                                "sse_of", "run_scenarios", "profile_frac",
                                "profile_lower", "profile_upper", "DAYS_CAL",
                                "profile_tasks"))

pl_results <- parallel::parLapply(cl2, seq_len(nrow(profile_tasks)), function(k) {
  pname <- profile_tasks$pname[k]; gi <- profile_tasks$gi[k]
  gval <- nom6[[pname]] * profile_frac[gi]
  free_names <- setdiff(cal_names, pname)
  obj <- function(free_vec) {
    names(free_vec) <- free_names
    full <- nom6; full[pname] <- gval; full[free_names] <- free_vec
    sse_of(full)
  }
  opt <- optim(nom6[free_names], obj, method = "L-BFGS-B",
               lower = profile_lower[free_names], upper = profile_upper[free_names],
               control = list(maxit = 300))
  data.frame(param = pname, value = gval, frac = profile_frac[gi], sse = opt$value)
})
parallel::stopCluster(cl2)
cat(sprintf("Profile likelihood done (%.1f s elapsed)\n",
            as.numeric(difftime(Sys.time(), t_pl, units = "secs"))))

profile_df <- do.call(rbind, pl_results)
profile_df <- profile_df[order(profile_df$param, profile_df$frac), ]
profile_df <- profile_df %>%
  group_by(param) %>%
  mutate(delta_sse = sse - min(sse)) %>%
  ungroup()

CHI2_THRESH <- 3.84  # 95% chi-squared, 1 d.f.
profile_df$above_threshold <- profile_df$delta_sse > CHI2_THRESH

verdict <- profile_df %>%
  group_by(param) %>%
  summarise(max_delta_sse = max(delta_sse),
            practically_identifiable_in_range = any(above_threshold), .groups = "drop")
cat("\nProfile likelihood verdict (practical identifiability within +/-35% range):\n")
print(as.data.frame(verdict))

saveRDS(profile_df, "identifiability_profile_likelihood.rds")

nom_df <- data.frame(param = cal_names, nom = as.numeric(nom6))
psS5 <- ggplot(profile_df, aes(x = value, y = delta_sse)) +
  geom_line(colour = "#1565C0", linewidth = 0.8) +
  geom_point(colour = "#1565C0", size = 1.5) +
  geom_hline(yintercept = CHI2_THRESH, colour = "#B71C1C", linetype = "dashed") +
  geom_vline(data = nom_df, aes(xintercept = nom), colour = "grey40", linetype = "dotted") +
  facet_wrap(~ param, scales = "free_x", nrow = 2) +
  labs(title = "Profile likelihood: practical identifiability of calibrated parameters",
       subtitle = sprintf("Dashed = 95%% chi-squared threshold (%.2f, 1 d.f.); dotted = calibrated nominal", CHI2_THRESH),
       x = "Parameter value", y = expression(Delta*"SSE (re-optimised)")) +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(out_dir, "FigS5_profile_likelihood.png"), psS5,
       width = 12, height = 7, dpi = 300, bg = "white")
cat("Saved FigS5_profile_likelihood.png\n")

cat("\n=== SUMMARY FOR MANUSCRIPT TEXT UPDATE ===\n")
cat(sprintf("Parameter set: %d calibrated + %d estimated = %d (identifiability sweep)\n",
            n_calibrated, n_estimated, length(all_id_params)))
cat(sprintf("Classification (live from ST1): %d literature (%.1f%%), %d calibrated, %d estimated, %d assumed, %d total\n",
            n_literature, 100 * n_literature / n_total, n_calibrated, n_estimated, n_assumed, n_total))
cat(sprintf("Endpoint FIM rank = %d\n", rank_endpoint))
cat(sprintf("Time-course FIM rank = %d\n", rank_tc))
cat(sprintf("Condition number = %.3g\n", cond_number))
cat(sprintf("Leading eigenvalue %% of info = %.1f%%\n", pct_leading))
cat(sprintf("Pairs gamma>20 (of %d total): %d\n", nrow(all_pairs), n_high))
cat("Top overall collinearity pair:", all_pairs[top_all_idx[1], 1], "-", all_pairs[top_all_idx[1], 2],
    sprintf(" (gamma=%.2f)\n", gamma_all[top_all_idx[1]]))
cat("Profile likelihood verdict per calibrated parameter:\n")
print(as.data.frame(verdict))
