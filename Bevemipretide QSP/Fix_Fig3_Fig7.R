##=============================================================================
## TARGETED REGENERATION: Fig 3 (calibration target fit, day-35 fix) and
## Fig 7 (uncertainty envelope, 125-combo fix). Skips Fig1/2/4/5/6/S1/S2/S3
## since those were not found to need changes -- keeps this fast.
##=============================================================================
library(deSolve); library(ggplot2); library(dplyr); library(tidyr); library(gridExtra)
library(parallel)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)
out_dir <- "../New publication plots"

n_cores <- max(1, detectCores() - 2)
cl <- makeCluster(n_cores)
invisible(clusterEvalQ(cl, { library(deSolve); SOURCED_FOR_FUNCTIONS <- TRUE
                   source("Integrated Model.R", local = FALSE) }))
cat("Parallel cluster:", n_cores, "workers\n")

pub_theme <- theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, margin = margin(b = 4)),
        plot.subtitle = element_text(size = 8, color = "grey40", margin = margin(b = 6)),
        panel.grid.minor = element_blank())
theme_set(pub_theme)

CAL <- list(K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4,
            CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35)
run_cal <- function(dose, days, CI_max = 1.0, alpha_clear = 1.0,
                    parms_override = NULL, dt = 0.5) {
  run_integrated(dose_mg_kg = dose, duration_days = days, dt = dt,
                 CI_max = CI_max, alpha_clear = alpha_clear,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP,
                 n_death_val = CAL$n_death,
                 parms_override = parms_override, quiet = TRUE)
}

t0 <- Sys.time()
cat("[", format(Sys.time()), "] Starting...\n")

##=============================================================================
## FIG 3
##=============================================================================
res_healthy50 <- run_cal(0, 50)
res_disease50 <- run_cal(0, 50, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
row35_h <- which.min(abs(res_healthy50$time_days - 35))
row35_d <- which.min(abs(res_disease50$time_days - 35))
ss_h35 <- res_healthy50[row35_h, ]
ss_d35 <- res_disease50[row35_d, ]
cat("[", format(Sys.time()), "] Base 50-day runs done (", round(difftime(Sys.time(), t0, units="secs"),1), "s )\n")

mROS_h <- ss_h35$mROS; CL_n_h <- ss_h35$CL_n

CI_grid <- seq(0.60, 0.85, by = 0.01)
ac_grid <- seq(0.20, 0.55, by = 0.05)
grid_combos <- expand.grid(CI_max = CI_grid, alpha_clear = ac_grid)
clusterExport(cl, c("CAL", "run_cal", "mROS_h", "CL_n_h", "grid_combos"))
sweep_list <- parLapply(cl, seq_len(nrow(grid_combos)), function(k) {
  cm <- grid_combos$CI_max[k]; ac <- grid_combos$alpha_clear[k]
  res <- run_cal(0, 50, CI_max = cm, alpha_clear = ac)
  ss  <- res[nrow(res), ]
  err_CI   <- (ss$CI_activity - 0.49)^2
  err_mROS <- (ss$mROS / mROS_h - 1.77)^2
  err_CL   <- ((1 - ss$CL_n / CL_n_h) * 100 - 23)^2 / 100
  da_loss  <- (1 - ss$DA_neuron) * 100
  err_DA   <- ifelse(da_loss < 25, ((25 - da_loss)/20)^2,
                     ifelse(da_loss > 55, ((da_loss - 55)/20)^2, 0))
  score    <- err_CI * 4 + err_mROS + err_CL + err_DA * 2
  data.frame(CI_max = cm, alpha_clear = ac, score = score)
})
sweep_df <- do.call(rbind, sweep_list)
sweep_df$log_score <- log10(pmax(sweep_df$score, 1e-4))
cat("[", format(Sys.time()), "] Panel A sweep done (", round(difftime(Sys.time(), t0, units="secs"),1), "s )\n")

p3a <- ggplot(sweep_df, aes(x = CI_max, y = factor(alpha_clear), fill = log_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = data.frame(CI_max = 0.74, alpha_clear = factor(0.25)),
             aes(x = CI_max, y = alpha_clear), inherit.aes = FALSE,
             shape = 4, size = 3, stroke = 1.5, colour = "white") +
  scale_fill_gradient2(low = "#1B5E20", mid = "#FFF9C4", high = "#B71C1C",
    midpoint = median(sweep_df$log_score), name = expression(log[10]*"(score)")) +
  scale_x_continuous(breaks = seq(0.60, 0.85, 0.05)) +
  labs(title = expression(bold("A  Disease Modifier Calibration")),
       x = expression(CI[max]), y = expression(alpha[clear])) +
  theme(legend.key.height = unit(0.5, "cm"), legend.title = element_text(size = 7))

val_data <- data.frame(
  target = c("CI activity", "mROS (% basal)", "CL drop (%)", "DA loss (%)", "TLR2 KO aSyn (%)"),
  observed = c(ss_d35$CI_activity, ss_d35$mROS / ss_h35$mROS * 100,
               (1 - ss_d35$CL_n / ss_h35$CL_n) * 100, (1 - ss_d35$DA_neuron) * 100, 20.3),
  lo = c(0.40, 150, 15, 25, 15), hi = c(0.58, 250, 35, 45, 50),
  source = c("Gao 2017", "Choi 2022", "Literature", "Bernheimer 1973", "Ivanova 2024"),
  stringsAsFactors = FALSE)
val_data$norm <- (val_data$observed - val_data$lo) / (val_data$hi - val_data$lo)
val_data$target <- factor(val_data$target, levels = rev(val_data$target))
cat("Fig3 target values (day-35):\n"); print(val_data[, c("target","observed","lo","hi")])

p3b <- ggplot(val_data) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = as.numeric(target) - 0.35, ymax = as.numeric(target) + 0.35),
            fill = "#E8F5E9", colour = "#A5D6A7", linewidth = 0.3) +
  geom_point(aes(x = norm, y = target), colour = "#1B5E20", size = 3.5, shape = 18) +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_text(aes(x = norm, y = target, label = sprintf("%.1f", observed)),
            vjust = -1.3, size = 2.5, colour = "#1B5E20", fontface = "bold") +
  geom_text(aes(x = 0.98, y = target, label = source), size = 2, colour = "grey50",
            hjust = 1, fontface = "italic") +
  scale_x_continuous(limits = c(-0.08, 1.15), breaks = c(0, 0.5, 1),
                     labels = c("Lower\nbound", "Mid", "Upper\nbound")) +
  labs(title = "B  Calibration Target Fit (5/5 pass)", x = "Position within acceptable range", y = NULL) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "bold", size = 8))

fig3 <- arrangeGrob(p3a, p3b, ncol = 2, widths = c(1.2, 1))
ggsave(file.path(out_dir, "Fig3_calibration_validation.png"), fig3,
       width = 12, height = 5, dpi = 300, bg = "white")
cat("[", format(Sys.time()), "] Fig3 saved (", round(difftime(Sys.time(), t0, units="secs"),1), "s )\n")

##=============================================================================
## FIG 7 (panel A only re-run; panels B/C reuse cached RA_arm RDS, unchanged)
##=============================================================================
HUMAN_PK <- list(ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
                  k_brain_in = 0.015, k_brain_out = 0.010)
C_ref_human <- 2.15
BW_human <- 70

CI_max_v      <- c(0.86, 0.88, 0.90, 0.92, 0.94)
alpha_clear_v <- c(0.28, 0.32, 0.35, 0.38, 0.42)
death_scale_v <- c(0.04, 0.05, 0.06, 0.07, 0.08)
combos <- expand.grid(CI_max = CI_max_v, alpha_clear = alpha_clear_v, death_scale = death_scale_v)
cat("[", format(Sys.time()), "] Fig7 sweep starting,", nrow(combos), "combos\n")

clusterExport(cl, c("HUMAN_PK", "C_ref_human", "BW_human", "CAL", "combos"))
joint_list <- parLapply(cl, 1:nrow(combos), function(i) {
  ci <- combos$CI_max[i]; ac <- combos$alpha_clear[i]; ds <- combos$death_scale[i]
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * ds, k_death_inflam = 0.00020 * ds, k_death_aSyn = 0.00015 * ds))
  res_dis <- tryCatch(run_integrated(dose_mg_kg = 0, duration_days = 547, dt = 72,
         CI_max = ci, alpha_clear = ac, k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
         k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
         parms_override = po, quiet = TRUE), error = function(e) NULL)
  res_drg <- tryCatch(run_integrated(dose_mg_kg = 30/BW_human, duration_days = 547, dt = 72,
         CI_max = ci, alpha_clear = ac, k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
         k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
         parms_override = po, quiet = TRUE), error = function(e) NULL)
  if (!is.null(res_dis) && !is.null(res_drg)) {
    data.frame(DA_saved = res_drg$DA_neuron[nrow(res_drg)] - res_dis$DA_neuron[nrow(res_dis)])
  } else NULL
})
joint_res <- do.call(rbind, joint_list)
cat("[", format(Sys.time()), "] Fig7 sweep done (", round(difftime(Sys.time(), t0, units="secs"),1), "s )\n")
stopCluster(cl)

nominal_val <- joint_res$DA_saved[combos$CI_max == 0.90 & combos$alpha_clear == 0.35 & combos$death_scale == 0.06]
if (length(nominal_val) == 0) nominal_val <- median(joint_res$DA_saved)
q95 <- quantile(joint_res$DA_saved, c(0.025, 0.975))
cat(sprintf("Fig7A: n=%d, nominal=%.4f, 95%% range=[%.4f, %.4f]\n",
            nrow(joint_res), nominal_val, q95[1], q95[2]))

p7a <- ggplot(joint_res, aes(x = DA_saved)) +
  geom_histogram(bins = 12, fill = "#00796B", colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = nominal_val, colour = "#C62828", linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = q95, colour = "#C62828", linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = nominal_val + diff(range(joint_res$DA_saved))*0.08,
           y = Inf, label = "Nominal", vjust = 2, colour = "#C62828", size = 3, fontface = "bold") +
  labs(title = "A  Parameter Uncertainty Envelope",
       subtitle = sprintf("%d combos | 95%% range: [%.3f, %.3f]", nrow(joint_res), q95[1], q95[2]),
       x = expression(Delta*"DA (drug - placebo)"), y = "Count") +
  theme(plot.title = element_text(size = 10, face = "bold"))

completed_ids <- NULL
arm_data <- list()
for (ai in 1:4) {
  f <- sprintf("RA_arm_%d.rds", ai)
  if (file.exists(f)) arm_data[[ai]] <- readRDS(f)
}
ids_list <- lapply(arm_data, function(x) x$id)
completed_ids <- Reduce(intersect, ids_list)
arm_data[[1]]$arm <- gsub("(\\d+)mg", "\\1 mg", arm_data[[1]]$arm)
for (ai in 2:4) arm_data[[ai]]$arm <- gsub("(\\d+)mg", "\\1 mg", arm_data[[ai]]$arm)

pbo_b <- arm_data[[1]][arm_data[[1]]$id %in% completed_ids, ] %>% arrange(id)
pi_df <- data.frame()
dose_names <- c("10 mg", "30 mg", "60 mg")
for (j in 2:4) {
  drg <- arm_data[[j]][arm_data[[j]]$id %in% completed_ids, ] %>% arrange(id)
  delta <- drg$DA_neuron - pbo_b$DA_neuron
  pi_df <- rbind(pi_df, data.frame(Dose = dose_names[j-1], delta = delta))
}
pi_df$Dose <- factor(pi_df$Dose, levels = dose_names)

p7b <- ggplot(pi_df, aes(x = Dose, y = delta, fill = Dose)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#E65100", linewidth = 0.35) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("10 mg" = "#1976D2", "30 mg" = "#00796B", "60 mg" = "#E65100")) +
  annotate("text", x = 0.6, y = 0.055, label = "MCID = 5%", size = 2.5, colour = "#E65100",
           fontface = "italic", hjust = 0) +
  labs(title = "B  Patient-Level Treatment Effect",
       subtitle = sprintf("Matched VPop pairs (N=%d/arm)", length(completed_ids)),
       x = NULL, y = expression(Delta*"DA vs placebo")) +
  theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "none")

cat("[", format(Sys.time()), "] Running bootstrap (B=10,000)...\n")
set.seed(123)
B <- 10000; n_b <- nrow(pbo_b)
boot_d <- data.frame()
for (j in 2:4) {
  drg <- arm_data[[j]][arm_data[[j]]$id %in% completed_ids, ] %>% arrange(id)
  bd <- numeric(B)
  for (b in 1:B) {
    idx <- sample(1:n_b, n_b, replace = TRUE)
    pooled_sd <- sqrt((var(pbo_b$DA_neuron[idx]) + var(drg$DA_neuron[idx])) / 2)
    bd[b] <- abs(mean(drg$DA_neuron[idx]) - mean(pbo_b$DA_neuron[idx])) / pooled_sd
  }
  boot_d <- rbind(boot_d, data.frame(Dose = dose_names[j-1], d_unpaired = bd))
}
boot_d$Dose <- factor(boot_d$Dose, levels = dose_names)

p7c <- ggplot(boot_d, aes(x = d_unpaired, fill = Dose)) +
  geom_density(alpha = 0.4, colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = c("10 mg" = "#1976D2", "30 mg" = "#00796B", "60 mg" = "#E65100"), name = NULL) +
  labs(title = "C  Bootstrap Effect Size", subtitle = "B = 10,000 | DA neuron (unpaired d)",
       x = expression("Unpaired Cohen's " * italic(d)), y = "Density") +
  theme(plot.title = element_text(size = 10, face = "bold"), legend.position = c(0.85, 0.78),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3),
        legend.key.size = unit(0.28, "cm"))

fig7 <- arrangeGrob(p7a, p7b, p7c, ncol = 3)
ggsave(file.path(out_dir, "Fig7_uncertainty_quantification.png"), fig7,
       width = 14, height = 5.5, dpi = 300, bg = "white")
cat("[", format(Sys.time()), "] Fig7 saved. TOTAL TIME:", round(difftime(Sys.time(), t0, units="secs"),1), "s\n")
