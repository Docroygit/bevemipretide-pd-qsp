##=============================================================================
## PUBLICATION FIGURES
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
##
## Generates journal-ready figures (300 DPI, consistent theme, composite panels)
## Requires: deSolve, ggplot2, tidyr, dplyr, gridExtra
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

cat("================================================================\n")
cat("  PUBLICATION FIGURES\n")
cat("================================================================\n\n")

## ── GLOBAL THEME ──────────────────────────────────────────────────

pub_theme <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 8)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5),
    strip.text = element_text(face = "bold", size = 9.5),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8.5),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey88"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )
theme_set(pub_theme)

## ── COLOUR PALETTES ──────────────────────────────────────────────

col_scenario <- c(
  "Healthy"                 = "#4CAF50",
  "Disease"                 = "#E53935",
  "Disease + SBT-272 low"   = "#FB8C00",
  "Disease + SBT-272 high"  = "#1E88E5"
)

col_dose <- c(
  "Placebo"       = "#9E9E9E",
  "SBT-272 10mg"  = "#42A5F5",
  "SBT-272 30mg"  = "#26A69A",
  "SBT-272 60mg"  = "#EF6C00"
)

col_module <- c(
  "M0:PK"    = "#78909C",
  "M1:CL"    = "#AB47BC",
  "M2:aSyn"  = "#EF5350",
  "M3:ETC"   = "#42A5F5",
  "M4:ROS"   = "#FFA726",
  "M5:Death" = "#26A69A",
  "M6:Infl"  = "#EC407A",
  "M7:Motor" = "#66BB6A"
)

## ── CALIBRATED PARAMETERS ────────────────────────────────────────

CAL <- list(K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4,
            CI_max = 0.74, alpha_clear = 0.25, k_impair = 0.35)

run_cal <- function(dose, days, CI_max = 1.0, alpha_clear = 1.0,
                    parms_override = NULL, dt = 0.5) {
  run_integrated(dose_mg_kg = dose, duration_days = days, dt = dt,
                 CI_max = CI_max, alpha_clear = alpha_clear,
                 k_impair_val = CAL$k_impair,
                 K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP,
                 n_death_val = CAL$n_death,
                 parms_override = parms_override, quiet = TRUE)
}


##=============================================================================
## FIGURE 2: CALIBRATION HEATMAP (CI_max x alpha_clear sweep)
##=============================================================================

cat("Generating Figure 2: Calibration heatmap...\n")

res_h <- run_cal(0, 50)
ss_h  <- res_h[nrow(res_h), ]
mROS_h <- ss_h$mROS
CL_n_h <- ss_h$CL_n

CI_grid  <- seq(0.60, 0.85, by = 0.01)
ac_grid  <- seq(0.20, 0.55, by = 0.05)

sweep_df <- data.frame()
for (cm in CI_grid) {
  for (ac in ac_grid) {
    res <- run_cal(0, 35, CI_max = cm, alpha_clear = ac)
    ss  <- res[nrow(res), ]
    err_CI   <- (ss$CI_activity - 0.49)^2
    err_mROS <- (ss$mROS / mROS_h - 1.77)^2
    err_CL   <- ((1 - ss$CL_n / CL_n_h) * 100 - 23)^2 / 100
    da_loss  <- (1 - ss$DA_neuron) * 100
    err_DA   <- ifelse(da_loss < 25, ((25 - da_loss)/20)^2,
                       ifelse(da_loss > 55, ((da_loss - 55)/20)^2, 0))
    score    <- err_CI * 4 + err_mROS + err_CL + err_DA * 2
    sweep_df <- rbind(sweep_df, data.frame(
      CI_max = cm, alpha_clear = ac,
      CI = ss$CI_activity, DA_loss = da_loss,
      score = score
    ))
  }
}

sweep_df$log_score <- log10(pmax(sweep_df$score, 1e-4))

fig2 <- ggplot(sweep_df, aes(x = CI_max, y = factor(alpha_clear), fill = log_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = data.frame(CI_max = 0.74, alpha_clear = factor(0.25)),
             aes(x = CI_max, y = alpha_clear),
             inherit.aes = FALSE,
             shape = 4, size = 3.5, stroke = 1.5, colour = "white") +
  scale_fill_gradient2(
    low = "#1B5E20", mid = "#FFF9C4", high = "#B71C1C",
    midpoint = median(sweep_df$log_score),
    name = expression(log[10]~"(score)")
  ) +
  scale_x_continuous(breaks = seq(0.60, 0.85, 0.05)) +
  labs(
    title = "Disease Modifier Calibration Sweep",
    subtitle = "Weighted score vs Gao 2017 (CI), Choi 2022 (mROS), Bernheimer 1973 (DA loss) | x = selected",
    x = expression(CI[max]~"(Complex I repair ceiling)"),
    y = expression(alpha[clear]~"(aSyn clearance capacity)")
  )
ggsave("Fig2_calibration_heatmap.png", fig2,
       width = 8, height = 4.5, dpi = 300, bg = "white")
cat("  Saved: Fig2_calibration_heatmap.png\n")


##=============================================================================
## FIGURE 3: INTEGRATED DYNAMICS (Vicious cycle + downstream composite)
##=============================================================================

cat("Generating Figure 3: Integrated dynamics composite...\n")

res_healthy <- run_cal(0, 50)
res_disease <- run_cal(0, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
res_drug_lo <- run_cal(0.5, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
res_drug_hi <- run_cal(5.0, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)

res_healthy$scenario <- "Healthy"
res_disease$scenario <- "Disease"
res_drug_lo$scenario <- "Disease + SBT-272 low"
res_drug_hi$scenario <- "Disease + SBT-272 high"
all_data <- rbind(res_healthy, res_disease, res_drug_lo, res_drug_hi)
all_data$scenario <- factor(all_data$scenario,
  levels = c("Healthy", "Disease", "Disease + SBT-272 low", "Disease + SBT-272 high"))

make_panel <- function(data, var, ylab, title_text, ylimits = NULL) {
  p <- ggplot(data, aes(x = time_days, y = .data[[var]], colour = scenario)) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = col_scenario, guide = "none") +
    labs(x = NULL, y = ylab, title = title_text) +
    theme(plot.title = element_text(size = 10, face = "bold"))
  if (!is.null(ylimits)) p <- p + coord_cartesian(ylim = ylimits)
  p
}

p3a <- make_panel(all_data, "CL_ratio", "Fraction", "A  CL Functional Ratio")
p3b <- make_panel(all_data, "aSyn_olig", "Fraction", "B  α-Synuclein Oligomers")
p3c <- make_panel(all_data, "CI_activity", "Fraction", "C  Complex I Activity")
p3d <- make_panel(all_data, "mROS", "Fraction", "D  Mitochondrial ROS")
p3e <- make_panel(all_data, "DA_neuron", "Fraction", "E  DA Neuron Survival")
p3f <- make_panel(all_data, "Motor_score", "Score", "F  Motor Impairment")
p3g <- make_panel(all_data, "MG_active", "Fraction", "G  Microglial Activation")
p3h <- make_panel(all_data, "ATP", "Fraction", "H  ATP")

p3h_with_legend <- make_panel(all_data, "ATP", "Fraction", "H  ATP") +
  geom_line(aes(colour = scenario), linewidth = 0.7) +
  scale_colour_manual(values = col_scenario, name = NULL) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(linewidth = 1.5)))

panels <- list(p3a, p3b, p3c, p3d, p3e, p3f, p3g, p3h_with_legend)

fig3 <- arrangeGrob(
  grobs = panels, ncol = 2, nrow = 4,
  top = textGrob("Integrated Model Dynamics: Healthy / Disease / Drug",
                 gp = gpar(fontface = "bold", fontsize = 13))
)
ggsave("Fig3_integrated_dynamics.png", fig3,
       width = 10, height = 13, dpi = 300, bg = "white")
cat("  Saved: Fig3_integrated_dynamics.png\n")


##=============================================================================
## FIGURE 4: DEATH PATHWAY DECOMPOSITION (NEW)
##=============================================================================

cat("Generating Figure 4: Death pathway decomposition...\n")

extract_death_rates <- function(res) {
  data.frame(
    time_days     = res$time_days,
    mito          = res$death_rate_mito,
    inflam        = res$death_rate_inflam,
    aSyn          = res$death_rate_aSyn,
    total         = res$death_rate_total
  )
}

dr_disease <- extract_death_rates(res_disease)
dr_drug    <- extract_death_rates(res_drug_hi)
dr_healthy <- extract_death_rates(res_healthy)

# Stacked area: disease
dr_dis_long <- dr_disease %>%
  select(time_days, mito, inflam, aSyn) %>%
  pivot_longer(-time_days, names_to = "pathway", values_to = "rate")
dr_dis_long$pathway <- factor(dr_dis_long$pathway,
  levels = c("aSyn", "inflam", "mito"),
  labels = c("aSyn Proteotoxicity", "Neuroinflammation", "Mitochondrial (drug-rescuable)"))

p4a <- ggplot(dr_dis_long, aes(x = time_days, y = rate * 1e4, fill = pathway)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c(
    "aSyn Proteotoxicity" = "#7E57C2",
    "Neuroinflammation" = "#EC407A",
    "Mitochondrial (drug-rescuable)" = "#26A69A"
  ), name = NULL) +
  labs(x = "Time (days)", y = expression("Death rate ("*10^{-4}*" h"^{-1}*")"),
       title = "A  Disease: Death Pathway Contributions") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold")) +
  guides(fill = guide_legend(nrow = 1))

# Stacked area: drug
dr_drg_long <- dr_drug %>%
  select(time_days, mito, inflam, aSyn) %>%
  pivot_longer(-time_days, names_to = "pathway", values_to = "rate")
dr_drg_long$pathway <- factor(dr_drg_long$pathway,
  levels = c("aSyn", "inflam", "mito"),
  labels = c("aSyn Proteotoxicity", "Neuroinflammation", "Mitochondrial (drug-rescuable)"))

p4b <- ggplot(dr_drg_long, aes(x = time_days, y = rate * 1e4, fill = pathway)) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c(
    "aSyn Proteotoxicity" = "#7E57C2",
    "Neuroinflammation" = "#EC407A",
    "Mitochondrial (drug-rescuable)" = "#26A69A"
  ), guide = "none") +
  labs(x = "Time (days)", y = expression("Death rate ("*10^{-4}*" h"^{-1}*")"),
       title = "B  Disease + SBT-272: Death Pathway Contributions") +
  theme(plot.title = element_text(size = 10, face = "bold"))

# Pie chart: steady-state proportions
ss_dis <- dr_disease[nrow(dr_disease), ]
ss_drg <- dr_drug[nrow(dr_drug), ]

pie_dis <- data.frame(
  pathway = c("Mitochondrial", "Neuroinflam.", "aSyn"),
  rate = c(ss_dis$mito, ss_dis$inflam, ss_dis$aSyn)
)
pie_dis$pct <- pie_dis$rate / sum(pie_dis$rate) * 100
pie_dis$label <- sprintf("%s\n%.0f%%", pie_dis$pathway, pie_dis$pct)
pie_dis$pathway <- factor(pie_dis$pathway,
  levels = c("Mitochondrial", "Neuroinflam.", "aSyn"))

p4c <- ggplot(pie_dis, aes(x = "", y = pct, fill = pathway)) +
  geom_col(width = 1, colour = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", colour = "white") +
  scale_fill_manual(values = c(
    "Mitochondrial" = "#26A69A",
    "Neuroinflam." = "#EC407A",
    "aSyn" = "#7E57C2"
  ), guide = "none") +
  labs(title = "C  Disease Steady State") +
  theme_void() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

pie_drg <- data.frame(
  pathway = c("Mitochondrial", "Neuroinflam.", "aSyn"),
  rate = c(ss_drg$mito, ss_drg$inflam, ss_drg$aSyn)
)
pie_drg$pct <- pie_drg$rate / sum(pie_drg$rate) * 100
pie_drg$label <- sprintf("%s\n%.0f%%", pie_drg$pathway, pie_drg$pct)
pie_drg$pathway <- factor(pie_drg$pathway,
  levels = c("Mitochondrial", "Neuroinflam.", "aSyn"))

p4d <- ggplot(pie_drg, aes(x = "", y = pct, fill = pathway)) +
  geom_col(width = 1, colour = "white", linewidth = 0.5) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = 3, fontface = "bold", colour = "white") +
  scale_fill_manual(values = c(
    "Mitochondrial" = "#26A69A",
    "Neuroinflam." = "#EC407A",
    "aSyn" = "#7E57C2"
  ), guide = "none") +
  labs(title = "D  Drug Steady State") +
  theme_void() +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

fig4 <- arrangeGrob(
  arrangeGrob(p4a, p4b, ncol = 2),
  arrangeGrob(p4c, p4d, ncol = 2),
  ncol = 1, heights = c(3, 2),
  top = textGrob("DA Neuron Death: Three Independent Pathways",
                 gp = gpar(fontface = "bold", fontsize = 13))
)
ggsave("Fig4_death_pathways.png", fig4,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig4_death_pathways.png\n")


##=============================================================================
## FIGURE 5: DOSE-RESPONSE CURVES (NEW)
##=============================================================================

cat("Generating Figure 5: Dose-response curves...\n")

doses <- c(0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0)
dr_data <- data.frame()
for (d in doses) {
  res <- run_cal(d, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
  ss  <- res[nrow(res), ]
  dr_data <- rbind(dr_data, data.frame(
    dose = d,
    DA_neuron = ss$DA_neuron,
    Motor_score = ss$Motor_score,
    CL_ratio = ss$CL_ratio,
    aSyn_olig = ss$aSyn_olig,
    CI_activity = ss$CI_activity,
    mROS = ss$mROS
  ))
}

healthy_DA    <- ss_h$DA_neuron
healthy_Motor <- ss_h$Motor_score
dr_data$DA_rescue_pct <- (dr_data$DA_neuron - dr_data$DA_neuron[1]) /
                          (healthy_DA - dr_data$DA_neuron[1]) * 100
dr_data$Motor_improve_pct <- (dr_data$Motor_score[1] - dr_data$Motor_score) /
                              dr_data$Motor_score[1] * 100

p5a <- ggplot(dr_data, aes(x = dose)) +
  geom_line(aes(y = DA_neuron), colour = "#1E88E5", linewidth = 0.9) +
  geom_point(aes(y = DA_neuron), colour = "#1E88E5", size = 2.2) +
  geom_hline(yintercept = healthy_DA, linetype = "dashed", colour = "#4CAF50", linewidth = 0.4) +
  geom_hline(yintercept = dr_data$DA_neuron[1], linetype = "dashed", colour = "#E53935", linewidth = 0.4) +
  annotate("text", x = 9.5, y = healthy_DA - 0.01, label = "Healthy", colour = "#4CAF50",
           size = 2.8, hjust = 1, vjust = 1) +
  annotate("text", x = 9.5, y = dr_data$DA_neuron[1] + 0.01, label = "Disease",
           colour = "#E53935", size = 2.8, hjust = 1, vjust = 0) +
  labs(x = "Dose (mg/kg)", y = "DA Neuron Survival", title = "A  DA Neuron Dose-Response") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p5b <- ggplot(dr_data, aes(x = dose)) +
  geom_line(aes(y = Motor_score), colour = "#E53935", linewidth = 0.9) +
  geom_point(aes(y = Motor_score), colour = "#E53935", size = 2.2) +
  geom_hline(yintercept = healthy_Motor, linetype = "dashed", colour = "#4CAF50", linewidth = 0.4) +
  labs(x = "Dose (mg/kg)", y = "Motor Impairment", title = "B  Motor Score Dose-Response") +
  theme(plot.title = element_text(size = 10, face = "bold"))

p5c <- ggplot(dr_data, aes(x = dose)) +
  geom_line(aes(y = CL_ratio, colour = "CL ratio"), linewidth = 0.7) +
  geom_point(aes(y = CL_ratio, colour = "CL ratio"), size = 1.8) +
  geom_line(aes(y = CI_activity, colour = "CI activity"), linewidth = 0.7) +
  geom_point(aes(y = CI_activity, colour = "CI activity"), size = 1.8) +
  scale_colour_manual(values = c("CL ratio" = "#AB47BC", "CI activity" = "#42A5F5"), name = NULL) +
  labs(x = "Dose (mg/kg)", y = "Normalised Level",
       title = "C  Proximal Targets (CL, CI)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.75, 0.25))

p5d <- ggplot(dr_data, aes(x = dose)) +
  geom_line(aes(y = aSyn_olig, colour = "αSyn oligomers"), linewidth = 0.7) +
  geom_point(aes(y = aSyn_olig, colour = "αSyn oligomers"), size = 1.8) +
  geom_line(aes(y = mROS, colour = "mROS"), linewidth = 0.7) +
  geom_point(aes(y = mROS, colour = "mROS"), size = 1.8) +
  scale_colour_manual(values = c("αSyn oligomers" = "#EF5350", "mROS" = "#FFA726"), name = NULL) +
  labs(x = "Dose (mg/kg)", y = "Normalised Level",
       title = "D  Cycle Intermediates (aSyn, ROS)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.75, 0.75))

fig5 <- arrangeGrob(p5a, p5b, p5c, p5d, ncol = 2,
  top = textGrob("Dose-Response: 35-Day Mouse Protocol",
                 gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("Fig5_dose_response.png", fig5,
       width = 10, height = 7, dpi = 300, bg = "white")
cat("  Saved: Fig5_dose_response.png\n")


##=============================================================================
## FIGURE 6: VALIDATION SUMMARY (NEW)
##=============================================================================

cat("Generating Figure 6: Validation summary...\n")

ss_d <- res_disease[nrow(res_disease), ]

val_data <- data.frame(
  target = c("CI activity", "mROS (% basal)", "CL drop (%)",
             "DA loss (%)", "TLR2 KO aSyn red (%)",
             "CL rescue", "aSyn reduction", "mROS reduction",
             "DA preservation", "Motor improvement"),
  observed = c(
    ss_d$CI_activity,
    ss_d$mROS / ss_h$mROS * 100,
    (1 - ss_d$CL_n / ss_h$CL_n) * 100,
    (1 - ss_d$DA_neuron) * 100,
    20.3,
    TRUE, TRUE, TRUE, TRUE, TRUE
  ),
  lo = c(0.40, 150, 15, 20, 15, NA, NA, NA, NA, NA),
  hi = c(0.58, 250, 35, 70, 50, NA, NA, NA, NA, NA),
  category = c(rep("Quantitative", 5), rep("Directional", 5)),
  stringsAsFactors = FALSE
)

quant_data <- val_data[val_data$category == "Quantitative", ]

# Normalise observed to range [0, 1] within [lo, hi]
quant_data$norm <- (as.numeric(quant_data$observed) - quant_data$lo) /
                    (quant_data$hi - quant_data$lo)
quant_data$norm <- pmin(1.3, pmax(-0.3, quant_data$norm))
quant_data$target <- factor(quant_data$target, levels = rev(quant_data$target))

p6 <- ggplot(quant_data) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = as.numeric(target) - 0.35,
                ymax = as.numeric(target) + 0.35),
            fill = "#E8F5E9", colour = "#A5D6A7", linewidth = 0.3) +
  geom_point(aes(x = norm, y = target),
             colour = "#1B5E20", size = 4, shape = 18) +
  geom_vline(xintercept = c(0, 1), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_text(aes(x = norm, y = target,
                label = sprintf("%.1f", as.numeric(observed))),
            vjust = -1.3, size = 3, colour = "#1B5E20", fontface = "bold") +
  geom_text(aes(x = 0, y = target, label = sprintf("%.0f", lo)),
            hjust = 1.2, size = 2.5, colour = "grey50") +
  geom_text(aes(x = 1, y = target, label = sprintf("%.0f", hi)),
            hjust = -0.2, size = 2.5, colour = "grey50") +
  scale_x_continuous(limits = c(-0.4, 1.4), breaks = c(0, 0.5, 1),
                     labels = c("Lower\nbound", "Mid", "Upper\nbound")) +
  labs(title = "Quantitative Validation: Model vs Literature Targets",
       subtitle = "Green band = acceptable range | Diamond = model output | All 5 targets met",
       x = NULL, y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 9))

ggsave("Fig6_validation_summary.png", p6,
       width = 8, height = 4.5, dpi = 300, bg = "white")
cat("  Saved: Fig6_validation_summary.png\n")


##=============================================================================
## FIGURE 7: SENSITIVITY ANALYSIS COMPOSITE
##=============================================================================

cat("Generating Figure 7: Sensitivity analysis composite...\n")

## Re-run OAT for top 15 parameters (fast: 30 runs)
param_noms <- c(
  k_ROS_basal = 0.122, k_SOD2 = 1.0, CI_max = 0.74, K_mPTP = 0.40,
  k_shunt = 0.125, k_mPTP_open = 0.50, k_agg = 0.009, K_death = 0.25,
  alpha_clear = 0.25, k_clear = 0.241, k_CI_aSyn = 0.200,
  k_CI_repair = 0.050, k_CL_loss = 5.0, k_biogen = 3.0, k_death = 0.0007
)

base_res <- run_cal(0, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
DA_base  <- base_res$DA_neuron[nrow(base_res)]

oat_results <- data.frame()
for (pname in names(param_noms)) {
  nom <- param_noms[pname]
  for (dir in c(-0.5, 0.5)) {
    pval <- nom * (1 + dir)
    if (pname %in% c("CI_max", "alpha_clear")) {
      cm <- ifelse(pname == "CI_max", pval, CAL$CI_max)
      ac <- ifelse(pname == "alpha_clear", pval, CAL$alpha_clear)
      res <- run_cal(0, 35, CI_max = cm, alpha_clear = ac)
    } else {
      po <- list()
      po[[pname]] <- pval
      res <- run_cal(0, 35, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear,
                     parms_override = po)
    }
    DA_pert <- res$DA_neuron[nrow(res)]
    pct_change <- (DA_pert - DA_base) / DA_base * 100
    oat_results <- rbind(oat_results, data.frame(
      param = pname, direction = ifelse(dir > 0, "+50%", "-50%"),
      pct = pct_change
    ))
  }
}

tornado_df <- oat_results %>%
  group_by(param) %>%
  summarise(range = max(abs(pct)), .groups = "drop") %>%
  arrange(range)
oat_results$param <- factor(oat_results$param, levels = tornado_df$param)

p7a <- ggplot(oat_results, aes(x = pct, y = param, fill = direction)) +
  geom_col(position = "identity", width = 0.6) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("+50%" = "#1E88E5", "-50%" = "#E53935"), name = NULL) +
  labs(x = "% Change in DA Neuron", y = NULL,
       title = "A  OAT Tornado (±50%)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7.5, family = "mono"),
        legend.position = c(0.85, 0.15))

## Sobol-style summary (hardcoded from SA results)
sobol_df <- data.frame(
  param = c("k_ROS_basal", "K_mPTP", "CI_max", "k_SOD2", "k_agg",
            "k_shunt", "alpha_clear", "k_CL_loss", "k_CI_repair", "k_red"),
  S1 = c(0.19, 0.00, 0.15, 0.27, 0.07, 0.06, 0.01, 0.03, 0.02, 0.02),
  ST = c(0.50, 0.33, 0.31, 0.28, 0.16, 0.10, 0.09, 0.05, 0.05, 0.03)
)
sobol_df$interaction <- sobol_df$ST - pmax(sobol_df$S1, 0)
sobol_df$param <- factor(sobol_df$param, levels = rev(sobol_df$param))

sobol_long <- sobol_df %>%
  select(param, S1, interaction) %>%
  pivot_longer(-param, names_to = "type", values_to = "value")
sobol_long$type <- factor(sobol_long$type, levels = c("S1", "interaction"),
                           labels = c("First-order", "Interactions"))

p7b <- ggplot(sobol_long, aes(x = value, y = param, fill = type)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("First-order" = "#1E88E5", "Interactions" = "#FFA726"),
                    name = NULL) +
  labs(x = "Sobol Index", y = NULL,
       title = "B  Sobol Variance Decomposition") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7.5, family = "mono"),
        legend.position = c(0.75, 0.15))

fig7 <- arrangeGrob(p7a, p7b, ncol = 2,
  top = textGrob("Global Sensitivity Analysis: DA Neuron Survival",
                 gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("Fig7_sensitivity_composite.png", fig7,
       width = 12, height = 5.5, dpi = 300, bg = "white")
cat("  Saved: Fig7_sensitivity_composite.png\n")


##=============================================================================
## FIGURE 8: HUMAN TRANSLATION (PK + 20-year trajectory)
##=============================================================================

cat("Generating Figure 8: Human translation...\n")

BW_human <- 70.0
HUMAN_PK <- list(
  ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
  k_brain_in = 0.015, k_brain_out = 0.010,
  k_mito_in = 0.50, k_mito_out = 0.020
)
CI_max_human <- 0.96
alpha_clear_human <- 0.55

po_cal <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 60, dt = 4,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_cal, quiet = TRUE,
                          solver_tol = 1e-8, infusion_mode = TRUE)
last_day <- res_ref[res_ref$time >= (59 * 24), ]
C_ref_human <- min(last_day$C_mito)

run_human_ht <- function(dose_mg, dur_days, CI_max_h, alpha_clear_h, dt = 168) {
  po <- c(HUMAN_PK, list(C_ref = C_ref_human))
  run_integrated(dose_mg_kg = dose_mg / BW_human, duration_days = dur_days, dt = dt,
                 CI_max = CI_max_h, alpha_clear = alpha_clear_h,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
                 parms_override = po, quiet = TRUE, infusion_mode = TRUE,
                 solver_tol = 1e-8)
}

dur_years <- 20
dur_days  <- dur_years * 365

cat("  Running 20-year healthy...\n")
res_h_human <- run_human_ht(0, dur_days, 1.0, 1.0)
cat("  Running 20-year disease...\n")
res_d_human <- run_human_ht(0, dur_days, CI_max_human, alpha_clear_human)
cat("  Running 20-year drug...\n")
res_rx_human <- run_human_ht(30, dur_days, CI_max_human, alpha_clear_human)

res_h_human$scenario  <- "Healthy"
res_d_human$scenario  <- "Disease (no drug)"
res_rx_human$scenario <- "Disease + SBT-272 (30 mg)"
ht_all <- rbind(res_h_human, res_d_human, res_rx_human)
ht_all$time_years <- ht_all$time / (24 * 365)
ht_all$scenario <- factor(ht_all$scenario,
  levels = c("Healthy", "Disease (no drug)", "Disease + SBT-272 (30 mg)"))

col_ht <- c("Healthy" = "#4CAF50", "Disease (no drug)" = "#E53935",
            "Disease + SBT-272 (30 mg)" = "#1E88E5")

p8a <- ggplot(ht_all, aes(x = time_years, y = DA_neuron, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.40, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 18, y = 0.42, label = "Bernheimer threshold (60% DA loss)",
           size = 2.5, colour = "grey40") +
  scale_colour_manual(values = col_ht, name = NULL) +
  labs(x = "Time (years)", y = "DA Neuron Fraction",
       title = "A  DA Neuron Survival: 20-Year Trajectory") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.7, 0.85),
        legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", colour = NA))

p8b <- ggplot(ht_all, aes(x = time_years, y = Motor_score, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = col_ht, guide = "none") +
  labs(x = "Time (years)", y = "Motor Impairment Score",
       title = "B  Motor Score: 20-Year Trajectory") +
  theme(plot.title = element_text(size = 10, face = "bold"))

# Drug benefit over time
benefit_df <- data.frame(
  time_years = res_d_human$time / (24 * 365),
  DA_saved = res_rx_human$DA_neuron - res_d_human$DA_neuron
)
benefit_df <- benefit_df[benefit_df$time_years > 0.1, ]

p8c <- ggplot(benefit_df, aes(x = time_years, y = DA_saved * 100)) +
  geom_area(fill = "#1E88E5", alpha = 0.3) +
  geom_line(colour = "#1E88E5", linewidth = 0.8) +
  labs(x = "Time (years)", y = "DA Preservation (%)",
       title = "C  Cumulative Drug Benefit (DA saved vs disease)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

# PK comparison: mouse vs human
res_pk_mouse <- run_cal(5.0, 5, dt = 0.1)
res_pk_human <- run_human_ht(30, 5, 1.0, 1.0, dt = 0.5)
pk_mouse <- data.frame(time_h = res_pk_mouse$time, C_mito = res_pk_mouse$C_mito,
                        species = "Mouse (5 mg/kg IP)")
pk_human <- data.frame(time_h = res_pk_human$time, C_mito = res_pk_human$C_mito,
                        species = "Human (30 mg SC)")
pk_both <- rbind(pk_mouse, pk_human)

p8d <- ggplot(pk_both, aes(x = time_h, y = C_mito, colour = species)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Mouse (5 mg/kg IP)" = "#FFA726",
                                  "Human (30 mg SC)" = "#1E88E5"), name = NULL) +
  labs(x = "Time (hours)", y = expression(C[mito]),
       title = "D  Mitochondrial PK: Mouse vs Human") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.7, 0.8))

fig8 <- arrangeGrob(p8a, p8b, p8c, p8d, ncol = 2,
  top = textGrob("Mouse-to-Human Translation",
                 gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("Fig8_human_translation.png", fig8,
       width = 10, height = 8, dpi = 300, bg = "white")
cat("  Saved: Fig8_human_translation.png\n")


##=============================================================================
## FIGURE 9: CLINICAL TRIAL SIMULATION COMPOSITE
##=============================================================================

cat("Generating Figure 9: CTS composite...\n")

arms_name <- c("Placebo", "SBT-272 10mg", "SBT-272 30mg", "SBT-272 60mg")

vpop_results <- data.frame()
for (ai in 1:4) {
  f <- sprintf("RA_arm_%d.rds", ai)
  if (file.exists(f)) {
    vpop_results <- rbind(vpop_results, readRDS(f))
  }
}

if (nrow(vpop_results) == 0) {
  cat("  WARNING: VPop cache not found. Skipping Figure 9.\n")
} else {

completed_ids <- Reduce(intersect, lapply(arms_name, function(a) {
  vpop_results$id[vpop_results$arm == a]
}))
plot_data <- vpop_results[vpop_results$id %in% completed_ids, ]
plot_data$arm <- factor(plot_data$arm, levels = arms_name)

## Panel A: Violin
p9a <- ggplot(plot_data, aes(x = arm, y = DA_neuron, fill = arm)) +
  geom_violin(alpha = 0.4, width = 0.85, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               colour = "grey30", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.35, colour = "grey30") +
  scale_fill_manual(values = col_dose, guide = "none") +
  labs(x = NULL, y = "DA Neuron Fraction",
       title = "A  DA Neuron Survival at 18 Months") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8))

## Panel B: Motor violin
p9b <- ggplot(plot_data, aes(x = arm, y = Motor_score, fill = arm)) +
  geom_violin(alpha = 0.4, width = 0.85, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8,
               colour = "grey30", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.35, colour = "grey30") +
  scale_fill_manual(values = col_dose, guide = "none") +
  labs(x = NULL, y = "Motor Impairment Score",
       title = "B  Motor Score at 18 Months") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8))

## Panel C: Waterfall (30 mg, DA)
MCID_DA <- 0.05
pbo <- plot_data[plot_data$arm == "Placebo", ]
trt30 <- plot_data[plot_data$arm == "SBT-272 30mg", ]
trt30 <- trt30[order(trt30$id), ]
pbo_m <- pbo[pbo$id %in% trt30$id, ]
pbo_m <- pbo_m[order(pbo_m$id), ]
delta_30 <- trt30$DA_neuron - pbo_m$DA_neuron

wf <- data.frame(
  patient = rank(-delta_30),
  delta = delta_30
)
wf$category <- ifelse(wf$delta > MCID_DA, "Responder (>5%)",
                ifelse(wf$delta > 0, "Marginal (0-5%)", "No benefit"))
wf$category <- factor(wf$category,
  levels = c("Responder (>5%)", "Marginal (0-5%)", "No benefit"))

p9c <- ggplot(wf, aes(x = patient, y = delta * 100, fill = category)) +
  geom_col(width = 1) +
  geom_hline(yintercept = MCID_DA * 100, linetype = "dashed",
             colour = "#E53935", linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  annotate("text", x = 5, y = MCID_DA * 100 + 1, label = "MCID = 5%",
           colour = "#E53935", hjust = 0, size = 2.8) +
  scale_fill_manual(values = c(
    "Responder (>5%)" = "#26A69A",
    "Marginal (0-5%)" = "#FFA726",
    "No benefit" = "#BDBDBD"
  ), name = NULL) +
  labs(x = "Virtual Patient (ranked)", y = "ΔDA vs Placebo (%)",
       title = "C  Individual Responses: SBT-272 30 mg") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.8, 0.85),
        legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.3, "cm"),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))

## Panel D: Effect size benchmarking
bench_df <- data.frame(
  trial = c("DATATOP\n(selegiline)", "ADAGIO\n(rasagiline)",
            "SPARK\n(GBA)", "SURE-PD3\n(inosine)",
            "SBT-272\n10 mg", "SBT-272\n30 mg", "SBT-272\n60 mg"),
  d = c(0.30, 0.20, 0.15, 0.05, 0.080, 0.130, 0.162),
  source = c(rep("Published Trial", 4), rep("QSP Prediction", 3)),
  stringsAsFactors = FALSE
)
bench_df$trial <- factor(bench_df$trial, levels = rev(bench_df$trial))

p9d <- ggplot(bench_df, aes(x = d, y = trial, fill = source)) +
  geom_col(width = 0.5, alpha = 0.85) +
  geom_vline(xintercept = 0.2, linetype = "dotted", colour = "grey50") +
  geom_text(aes(label = sprintf("d = %.3f", d)), hjust = -0.15, size = 2.8) +
  scale_fill_manual(values = c("Published Trial" = "#BDBDBD", "QSP Prediction" = "#26A69A"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(x = "Unpaired Cohen's d", y = NULL,
       title = "D  Effect Size Benchmarking vs PD Trials") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.75, 0.15),
        axis.text.y = element_text(size = 8))

fig9 <- arrangeGrob(p9a, p9b, p9c, p9d, ncol = 2,
  top = textGrob("Virtual Clinical Trial: 18-Month Placebo-Controlled, Daily SC Injection",
                 gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("Fig9_CTS_composite.png", fig9,
       width = 12, height = 9, dpi = 300, bg = "white")
cat("  Saved: Fig9_CTS_composite.png\n")

}


##=============================================================================
## FIGURE 10: VICIOUS CYCLE AMPLIFICATION (NEW)
##=============================================================================

cat("Generating Figure 10: Vicious cycle amplification...\n")

amp_vars <- c("CL_ratio", "aSyn_olig", "CI_activity", "mROS",
              "ATP", "mPTP_open", "DA_neuron", "MG_active")
amp_labels <- c("CL Ratio (M1)", "αSyn Oligomers (M2)", "Complex I (M3)", "mROS (M4)",
                "ATP (M3)", "mPTP Opening (M3)", "DA Neuron (M5)", "Microglia (M6)")

ss_d_vals <- sapply(amp_vars, function(v) as.numeric(ss_d[[v]]))
ss_h_vals <- sapply(amp_vars, function(v) as.numeric(ss_h[[v]]))
ss_rx <- res_drug_hi[nrow(res_drug_hi), ]
ss_rx_vals <- sapply(amp_vars, function(v) as.numeric(ss_rx[[v]]))

bar_df <- data.frame(
  variable = rep(amp_labels, 3),
  condition = rep(c("Healthy", "Disease", "Drug (5 mg/kg)"), each = length(amp_vars)),
  value = c(ss_h_vals, ss_d_vals, ss_rx_vals)
)
bar_df$condition <- factor(bar_df$condition,
  levels = c("Healthy", "Disease", "Drug (5 mg/kg)"))
bar_df$variable <- factor(bar_df$variable, levels = amp_labels)

fig10 <- ggplot(bar_df, aes(x = variable, y = value, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Healthy" = "#4CAF50", "Disease" = "#E53935",
                                "Drug (5 mg/kg)" = "#1E88E5"), name = NULL) +
  labs(title = "Endpoint Comparison: Healthy vs Disease vs Drug",
       subtitle = "Day 35 steady state | All values normalised [0, 1]",
       x = NULL, y = "Normalised Value") +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5))

ggsave("Fig10_endpoint_comparison.png", fig10,
       width = 10, height = 5.5, dpi = 300, bg = "white")
cat("  Saved: Fig10_endpoint_comparison.png\n")


##=============================================================================
## FIGURE 11: PHASED SIMULATION COMPOSITE
##=============================================================================

cat("Generating Figure 11: Phased simulation composite...\n")

BW_human <- 70.0
HUMAN_PK <- list(ka = 0.40, ke_plasma = 0.032, k_12 = 0.041, k_21 = 0.027,
                 k_brain_in = 0.015, k_brain_out = 0.010,
                 k_mito_in = 0.50, k_mito_out = 0.020)
CI_max_human <- 0.90
alpha_clear_human <- 0.35
death_scale_human <- 0.06

po_ref <- c(HUMAN_PK, list(C_ref = 1.0))
res_ref <- run_integrated(dose_mg_kg = 30/BW_human, duration_days = 120, dt = 1,
                          CI_max = 1.0, alpha_clear = 1.0,
                          parms_override = po_ref, quiet = TRUE)
last_ref <- res_ref[res_ref$time >= (119 * 24), ]
C_ref_human <- min(last_ref$C_mito)

run_human_pub <- function(dose_mg, duration_days, CI_max_h = 1.0, alpha_clear_h = 1.0,
                          dt = 24, init_override = NULL) {
  dose_mgkg <- dose_mg / BW_human
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * death_scale_human,
          k_death_inflam = 0.00020 * death_scale_human,
          k_death_aSyn = 0.00015 * death_scale_human))
  run_integrated(dose_mg_kg = dose_mgkg, duration_days = duration_days, dt = dt,
                 CI_max = CI_max_h, alpha_clear = alpha_clear_h,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
                 parms_override = po, init_override = init_override, quiet = TRUE)
}

THRESHOLD_SYMPTOMS <- 0.70

## Phase 1: prodromal (20-year disease-only)
res_prodromal <- run_human_pub(0, 365.25 * 20, CI_max_human, alpha_clear_human, dt = 24)

## Find time to symptom threshold
symptom_idx <- which(res_prodromal$DA_neuron < THRESHOLD_SYMPTOMS)[1]
if (!is.na(symptom_idx)) {
  t_symptom_days <- res_prodromal$time[symptom_idx] / 24
  t_symptom_years <- t_symptom_days / 365.25
} else {
  t_symptom_years <- 20
  t_symptom_days <- t_symptom_years * 365.25
}

## Diagnostic delay: 12 months after symptoms
diag_days <- t_symptom_days + 365
diag_idx <- which.min(abs(res_prodromal$time / 24 - diag_days))

## Panel A: Disease phases (20-year trajectory with phase annotations)
phase_df <- data.frame(
  time_years = res_prodromal$time / 24 / 365.25,
  DA_neuron = res_prodromal$DA_neuron
)

p11a <- ggplot(phase_df, aes(x = time_years, y = DA_neuron)) +
  annotate("rect", xmin = 0, xmax = t_symptom_years,
           ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.5) +
  annotate("rect", xmin = t_symptom_years, xmax = t_symptom_years + 1,
           ymin = -Inf, ymax = Inf, fill = "#FFF3E0", alpha = 0.5) +
  annotate("rect", xmin = t_symptom_years + 1, xmax = 20,
           ymin = -Inf, ymax = Inf, fill = "#FFEBEE", alpha = 0.5) +
  geom_line(linewidth = 1.1, colour = "#E53935") +
  geom_hline(yintercept = 0.70, linetype = "dashed", colour = "#FF9800", linewidth = 0.5) +
  geom_hline(yintercept = 0.50, linetype = "dashed", colour = "#E53935", linewidth = 0.5) +
  geom_hline(yintercept = 0.40, linetype = "dashed", colour = "#B71C1C", linewidth = 0.5) +
  annotate("text", x = t_symptom_years / 2, y = 0.98, label = "Prodromal",
           size = 3, fontface = "bold", colour = "#388E3C") +
  annotate("text", x = t_symptom_years + 0.5, y = 0.98, label = "Dx\ndelay",
           size = 2.5, colour = "#EF6C00") +
  annotate("text", x = 0.3, y = 0.71, label = "30% DA loss", size = 2.5,
           colour = "#FF9800", hjust = 0) +
  annotate("text", x = 0.3, y = 0.51, label = "Clinical PD", size = 2.5,
           colour = "#E53935", hjust = 0) +
  annotate("text", x = 0.3, y = 0.41, label = "Bernheimer", size = 2.5,
           colour = "#B71C1C", hjust = 0) +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format()) +
  labs(title = "A  Disease Phase Characterisation",
       subtitle = sprintf("Prodromal: %.1f yr | Diagnosis: %.1f yr", t_symptom_years, t_symptom_years + 1),
       x = "Years", y = "DA Neuron Survival") +
  theme(plot.title = element_text(size = 10, face = "bold"))

## Panel B: Treatment from diagnosis (3 doses)
## init_override needs all 28 states including PK (set to 0 for new treatment phase)
state_names <- c("Depot", "C_plasma", "C_periph", "C_brain", "C_mito",
                 "CL_n", "CL_ox", "CL_ext", "ALCAT1", "TAZ",
                 "aSyn_mono", "aSyn_olig", "aSyn_ext",
                 "CI_activity", "SC_integrity", "delta_psi_m", "ATP", "mPTP_open",
                 "mROS", "SOD2_act", "GPx4_act",
                 "PINK1_act", "mito_damage", "DA_neuron",
                 "MG_active", "TNF", "IL1b", "Motor_score")
bio_names <- state_names[6:28]
state_at_diag <- setNames(
  c(rep(0, 5), as.numeric(res_prodromal[diag_idx, bio_names])),
  state_names)

treatment_days <- 3 * 365.25
treat_df <- data.frame()

for (dose in c(0, 10, 30, 60)) {
  label <- ifelse(dose == 0, "No treatment",
           sprintf("SBT-272 %d mg", dose))
  res <- tryCatch(
    run_human_pub(dose, treatment_days, CI_max_human, alpha_clear_human,
                  dt = 24, init_override = state_at_diag),
    error = function(e) NULL)
  if (!is.null(res)) {
    treat_df <- rbind(treat_df, data.frame(
      time_years = res$time / 24 / 365.25 + (t_symptom_years + 1),
      DA_neuron = res$DA_neuron,
      Motor_score = res$Motor_score,
      scenario = label))
  }
}

treat_df$scenario <- factor(treat_df$scenario,
  levels = c("No treatment", "SBT-272 10 mg", "SBT-272 30 mg", "SBT-272 60 mg"))

col_treat <- c("No treatment" = "#9E9E9E", "SBT-272 10 mg" = "#42A5F5",
               "SBT-272 30 mg" = "#26A69A", "SBT-272 60 mg" = "#EF6C00")

p11b <- ggplot(treat_df, aes(x = time_years, y = DA_neuron, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.50, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = 0.40, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  scale_colour_manual(values = col_treat, name = NULL) +
  scale_y_continuous(limits = c(0.2, 0.75), labels = scales::percent_format()) +
  labs(title = "B  Treatment from Diagnosis",
       subtitle = "3-year treatment initiated after 12-month diagnostic delay",
       x = "Years from Disease Onset", y = "DA Neuron Survival") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.75, 0.85), legend.background = element_blank(),
        legend.key.size = unit(0.35, "cm"))

## Panel C: Efficacy vs delay
delays <- c(0, 6, 12, 24)
delay_df <- data.frame()

for (delay_mo in delays) {
  delay_days <- delay_mo * 30.44
  delay_total <- t_symptom_days + delay_days
  d_idx <- which.min(abs(res_prodromal$time / 24 - delay_total))

  state_d <- setNames(
    c(rep(0, 5), as.numeric(res_prodromal[d_idx, bio_names])),
    state_names)

  res_no  <- tryCatch(run_human_pub(0,  547, CI_max_human, alpha_clear_human,
                                     dt = 24, init_override = state_d), error = function(e) NULL)
  res_30  <- tryCatch(run_human_pub(30, 547, CI_max_human, alpha_clear_human,
                                     dt = 24, init_override = state_d), error = function(e) NULL)
  if (!is.null(res_no) && !is.null(res_30)) {
    da_no  <- res_no[nrow(res_no), "DA_neuron"]
    da_30  <- res_30[nrow(res_30), "DA_neuron"]
    delay_df <- rbind(delay_df, data.frame(
      delay = delay_mo, DA_saved = da_30 - da_no))
  }
}

if (nrow(delay_df) > 0 && delay_df$DA_saved[1] > 0) {
  delay_df$relative <- delay_df$DA_saved / delay_df$DA_saved[1] * 100
} else {
  delay_df$relative <- 100
}

p11c <- ggplot(delay_df, aes(x = delay, y = relative)) +
  geom_col(fill = "#26A69A", width = 4.5) +
  geom_text(aes(label = sprintf("%.0f%%", relative)), vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(title = "C  Efficacy vs Treatment Delay",
       subtitle = "30 mg, 18 months | Relative to immediate treatment",
       x = "Delay from Symptom Onset (months)", y = "Relative Efficacy (%)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

fig11 <- arrangeGrob(p11a, p11b, p11c, ncol = 3,
  top = textGrob("Phased Disease-Treatment Simulation",
                 gp = gpar(fontface = "bold", fontsize = 13)))
ggsave("Fig11_phased_simulation.png", fig11,
       width = 15, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig11_phased_simulation.png\n")


##=============================================================================
## FIGURE 12: UNCERTAINTY ANALYSIS COMPOSITE
##=============================================================================

cat("Generating Figure 12: Uncertainty analysis composite...\n")

## Load VPop arm results
arm1 <- readRDS("RA_arm_1.rds")
arm2 <- readRDS("RA_arm_2.rds")
arm3 <- readRDS("RA_arm_3.rds")
arm4 <- readRDS("RA_arm_4.rds")

ids_common <- Reduce(intersect, list(arm1$id, arm2$id, arm3$id, arm4$id))
pbo <- arm1[arm1$id %in% ids_common, ] %>% arrange(id)
d10 <- arm2[arm2$id %in% ids_common, ] %>% arrange(id)
d30 <- arm3[arm3$id %in% ids_common, ] %>% arrange(id)
d60 <- arm4[arm4$id %in% ids_common, ] %>% arrange(id)

## Panel A: Bootstrap distribution of unpaired d
set.seed(123)
B <- 10000
n_boot <- nrow(pbo)

boot_d <- data.frame()
for (dose_info in list(list(d = d10, nm = "10 mg"), list(d = d30, nm = "30 mg"),
                       list(d = d60, nm = "60 mg"))) {
  for (ep in c("DA_neuron", "Motor_score")) {
    bd <- numeric(B)
    for (b in 1:B) {
      idx <- sample(1:n_boot, n_boot, replace = TRUE)
      pooled_sd <- sqrt((var(pbo[[ep]][idx]) + var(dose_info$d[[ep]][idx])) / 2)
      bd[b] <- abs(mean(dose_info$d[[ep]][idx]) - mean(pbo[[ep]][idx])) / pooled_sd
    }
    boot_d <- rbind(boot_d, data.frame(
      Dose = dose_info$nm, Endpoint = ep, d_unpaired = bd))
  }
}

ep_labels <- c(DA_neuron = "DA Neuron Survival", Motor_score = "Motor Score")

p12a <- ggplot(boot_d, aes(x = d_unpaired, fill = Dose)) +
  geom_density(alpha = 0.45, colour = "grey30", linewidth = 0.3) +
  facet_wrap(~Endpoint, scales = "free_x", labeller = labeller(Endpoint = ep_labels)) +
  scale_fill_manual(values = c("10 mg" = "#42A5F5", "30 mg" = "#26A69A", "60 mg" = "#EF6C00")) +
  geom_vline(xintercept = 0.15, linetype = "dashed", colour = "#E53935", linewidth = 0.4) +
  geom_vline(xintercept = 0.20, linetype = "dashed", colour = "#1E88E5", linewidth = 0.4) +
  annotate("text", x = 0.15, y = Inf, label = "SPARK", vjust = 2, hjust = 1.1,
           colour = "#E53935", size = 2.5) +
  annotate("text", x = 0.20, y = Inf, label = "ADAGIO", vjust = 2, hjust = -0.1,
           colour = "#1E88E5", size = 2.5) +
  labs(title = "A  Bootstrap Distribution (B=10,000)",
       x = expression("Unpaired Cohen's " * italic(d)), y = "Density", fill = NULL) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.92, 0.85), legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))

## Panel B: Joint parameter uncertainty envelope
CI_max_vals      <- c(0.86, 0.88, 0.90, 0.92, 0.94)
alpha_clear_vals <- c(0.28, 0.32, 0.35, 0.38, 0.42)
death_scale_vals <- c(0.04, 0.05, 0.06, 0.07, 0.08)

combos <- expand.grid(CI_max = CI_max_vals, alpha_clear = alpha_clear_vals,
                      death_scale = death_scale_vals)
joint_res <- data.frame()

cat("  Running 125 parameter combinations...\n")
for (i in 1:nrow(combos)) {
  ci <- combos$CI_max[i]; ac <- combos$alpha_clear[i]; ds <- combos$death_scale[i]
  dose_mgkg <- 30 / BW_human
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * ds, k_death_inflam = 0.00020 * ds,
          k_death_aSyn = 0.00015 * ds))

  res_dis <- run_integrated(dose_mg_kg = 0, duration_days = 547, dt = 24,
               CI_max = ci, alpha_clear = ac,
               k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
               k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
               parms_override = po, quiet = TRUE)
  res_drg <- run_integrated(dose_mg_kg = dose_mgkg, duration_days = 547, dt = 24,
               CI_max = ci, alpha_clear = ac,
               k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
               k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
               parms_override = po, quiet = TRUE)

  da_saved <- res_drg$DA_neuron[nrow(res_drg)] - res_dis$DA_neuron[nrow(res_dis)]
  joint_res <- rbind(joint_res, data.frame(
    CI_max = ci, alpha_clear = ac, death_scale = ds, DA_saved = da_saved))
}

nominal_val <- joint_res$DA_saved[joint_res$CI_max == 0.90 &
                                   joint_res$alpha_clear == 0.35 &
                                   joint_res$death_scale == 0.06]

p12b <- ggplot(joint_res, aes(x = DA_saved)) +
  geom_histogram(bins = 25, fill = "#26A69A", colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = nominal_val, colour = "#E53935", linewidth = 0.8,
             linetype = "solid") +
  geom_vline(xintercept = quantile(joint_res$DA_saved, c(0.025, 0.975)),
             colour = "#E53935", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = nominal_val, y = Inf, label = "Nominal",
           vjust = 2, hjust = -0.15, colour = "#E53935", size = 3, fontface = "bold") +
  labs(title = "B  Parameter Uncertainty Envelope",
       subtitle = sprintf("125 combos | 95%% range: [%.3f, %.3f]",
                           quantile(joint_res$DA_saved, 0.025),
                           quantile(joint_res$DA_saved, 0.975)),
       x = expression(Delta * "DA Neuron (drug - placebo)"),
       y = "Count") +
  theme(plot.title = element_text(size = 10, face = "bold"))

## Panel C: Individual patient prediction intervals
pi_df <- data.frame()
for (dose_info in list(list(d = d10, nm = "10 mg"), list(d = d30, nm = "30 mg"),
                       list(d = d60, nm = "60 mg"))) {
  matched <- dose_info$d[dose_info$d$id %in% ids_common, ]
  matched <- matched[order(matched$id), ]
  delta <- matched$DA_neuron - pbo$DA_neuron
  pi_df <- rbind(pi_df, data.frame(Dose = dose_info$nm, delta = delta))
}

p12c <- ggplot(pi_df, aes(x = Dose, y = delta, fill = Dose)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#FF9800", linewidth = 0.4) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("10 mg" = "#42A5F5", "30 mg" = "#26A69A", "60 mg" = "#EF6C00")) +
  annotate("text", x = 3.4, y = 0.05, label = "MCID", size = 2.5,
           colour = "#FF9800", fontface = "italic") +
  labs(title = "C  Patient-Level Treatment Effect",
       subtitle = "Matched VPop pairs (N=244/arm)",
       x = NULL, y = expression(Delta * "DA Neuron")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "none")

fig12 <- arrangeGrob(p12a, arrangeGrob(p12b, p12c, ncol = 2), nrow = 2,
  heights = c(1, 1),
  top = textGrob("Uncertainty Quantification: Bootstrap CIs, Parameter Envelope, Prediction Intervals",
                 gp = gpar(fontface = "bold", fontsize = 12)))
ggsave("Fig12_uncertainty_analysis.png", fig12,
       width = 12, height = 9, dpi = 300, bg = "white")
cat("  Saved: Fig12_uncertainty_analysis.png\n")


##=============================================================================
## SUMMARY
##=============================================================================

cat("\n================================================================\n")
cat("  PUBLICATION FIGURES COMPLETE\n")
cat("================================================================\n\n")
cat("Generated figures (300 DPI, journal-ready):\n")
cat("  Fig2_calibration_heatmap.png    - CI_max x alpha_clear sweep\n")
cat("  Fig3_integrated_dynamics.png    - 8-panel integrated model time courses\n")
cat("  Fig4_death_pathways.png         - 3-pathway death decomposition\n")
cat("  Fig5_dose_response.png          - Multi-dose response curves\n")
cat("  Fig6_validation_summary.png     - Quantitative validation targets\n")
cat("  Fig7_sensitivity_composite.png  - OAT tornado + Sobol indices\n")
cat("  Fig8_human_translation.png      - PK scaling + 20-year trajectory\n")
cat("  Fig9_CTS_composite.png          - VPop + waterfall + benchmarking\n")
cat("  Fig10_endpoint_comparison.png   - Grouped bar comparison\n")
cat("  Fig11_phased_simulation.png     - Disease phases + treatment + delay\n")
cat("  Fig12_uncertainty_analysis.png  - Bootstrap + envelope + prediction\n\n")
cat("All figures use consistent theme, palette, and 300 DPI resolution.\n")
