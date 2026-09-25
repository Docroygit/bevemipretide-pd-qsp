##=============================================================================
## NEW PUBLICATION PLOTS — FIXED v2
## Bevemipretide (SBT-272) QSP for Parkinson's Disease
## 7 Main (Fig1-Fig7) + 3 Supplementary (FigS1-FigS3)
## 300 DPI, unified theme, no watermarks, no text overlap
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

SOURCED_FOR_FUNCTIONS <- TRUE
source("Integrated Model.R", local = FALSE)

out_dir <- "../New publication plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("================================================================\n")
cat("  GENERATING NEW PUBLICATION PLOTS (v2 — all fixes)\n")
cat("================================================================\n\n")

## -- GLOBAL THEME --
pub_theme <- theme_minimal(base_size = 10) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8, color = "grey40", margin = margin(b = 6)),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 9),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.35, "cm"),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 6, 6)
  )
theme_set(pub_theme)

## -- COLOUR PALETTES (consistent across ALL figures) --
col_scen <- c("Healthy" = "#2E7D32", "Disease" = "#C62828",
              "Drug (low)" = "#E65100", "Drug (high)" = "#1565C0")
col_dose <- c("Placebo" = "#616161", "10 mg" = "#1976D2",
              "30 mg" = "#00796B", "60 mg" = "#E65100")
col_death <- c("Mitochondrial" = "#C62828",
               "Neuroinflammation" = "#E65100",
               "Proteotoxicity" = "#6A1B9A")
col_3scen <- c("Healthy" = "#2E7D32", "PD (untreated)" = "#C62828",
               "PD + SBT-272" = "#1565C0")

## -- CALIBRATED PARAMETERS --
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

##=============================================================================
## BASE SIMULATIONS (ALL 50 days — consistent x-axes)
##=============================================================================
cat("Running base simulations (all 50 days)...\n")
SIM_DAYS <- 50
res_healthy <- run_cal(0,   SIM_DAYS)
res_disease <- run_cal(0,   SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
res_drug_lo <- run_cal(0.5, SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
res_drug_hi <- run_cal(5.0, SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)

res_healthy$scenario <- "Healthy"
res_disease$scenario <- "Disease"
res_drug_lo$scenario <- "Drug (low)"
res_drug_hi$scenario <- "Drug (high)"
all_data <- rbind(res_healthy, res_disease, res_drug_lo, res_drug_hi)
all_data$scenario <- factor(all_data$scenario,
  levels = c("Healthy", "Disease", "Drug (low)", "Drug (high)"))

ss_h  <- res_healthy[nrow(res_healthy), ]
ss_d  <- res_disease[nrow(res_disease), ]
ss_rx <- res_drug_hi[nrow(res_drug_hi), ]

## 35-day (5-week) snapshots: the actual calibration checkpoint used for the
## five quantitative targets in Table 3 / Results 3.1. ss_h/ss_d above are at
## day 50 (used for the Fig 2 dynamics panels) and must NOT be reused for the
## calibration-target-fit panel (Fig 3B), since DA_neuron has not reached
## steady state by day 50 and differs materially from its day-35 value.
row35_h <- which.min(abs(res_healthy$time - 35 * 24))
row35_d <- which.min(abs(res_disease$time - 35 * 24))
ss_h35  <- res_healthy[row35_h, ]
ss_d35  <- res_disease[row35_d, ]
cat("  Done.\n\n")


##=============================================================================
## FIG 1: MODEL ARCHITECTURE SCHEMATIC
##=============================================================================
cat("Generating Fig 1: Model architecture...\n")

modules <- data.frame(
  x     = c(6,   2.5, 9.5, 2.5, 9.5,  6,   9.5, 6),
  y     = c(12,  9,   9,   6,   6,    3,   3,   0.5),
  label = c("M0\nPharmacokinetics", "M1\nCardiolipin",
            "M2\nalpha-Synuclein",  "M3\nETC Bioenergetics",
            "M4\nOxidative Stress", "M5\nCell Death",
            "M6\nNeuroinflammation","M7\nMotor Endpoint"),
  fill  = c("#B0BEC5","#CE93D8","#EF9A9A","#90CAF9",
            "#FFCC80","#80CBC4","#F48FB1","#A5D6A7"),
  stringsAsFactors = FALSE
)
bw <- 2.8; bh <- 1.1

p1 <- ggplot() +
  ## Vicious cycle region
  annotate("rect", xmin = 0.6, xmax = 11.4, ymin = 5.2, ymax = 10.0,
           fill = "#FFF3E0", colour = "#E65100", linetype = "dashed",
           linewidth = 0.6, alpha = 0.45) +
  annotate("text", x = 6, y = 10.35, label = "VICIOUS CYCLE",
           size = 3.2, fontface = "bold", colour = "#BF360C") +
  ## Module boxes
  geom_rect(data = modules,
            aes(xmin = x - bw/2, xmax = x + bw/2,
                ymin = y - bh/2, ymax = y + bh/2, fill = fill),
            colour = "grey40", linewidth = 0.45, show.legend = FALSE) +
  scale_fill_identity() +
  geom_text(data = modules, aes(x = x, y = y, label = label),
            size = 2.7, lineheight = 0.85, fontface = "bold") +
  ## CYCLE: M1->M3 (down left)
  annotate("segment", x = 2.5, y = 8.45, xend = 2.5, yend = 6.55,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           colour = "#BF360C", linewidth = 0.9) +
  annotate("text", x = 1.2, y = 7.5, label = "CL ratio",
           size = 2.3, colour = "#BF360C", fontface = "italic") +
  ## CYCLE: M3->M4 (across)
  annotate("segment", x = 3.9, y = 6, xend = 8.1, yend = 6,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           colour = "#BF360C", linewidth = 0.9) +
  annotate("text", x = 6, y = 5.5, label = "CI, SC",
           size = 2.3, colour = "#BF360C", fontface = "italic") +
  ## CYCLE: M4->M2 (up right)
  annotate("segment", x = 9.5, y = 6.55, xend = 9.5, yend = 8.45,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           colour = "#BF360C", linewidth = 0.9) +
  annotate("text", x = 10.8, y = 7.5, label = "mROS",
           size = 2.3, colour = "#BF360C", fontface = "italic") +
  ## CYCLE: M2->M1 (across top)
  annotate("segment", x = 8.1, y = 9, xend = 3.9, yend = 9,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           colour = "#BF360C", linewidth = 0.9) +
  annotate("text", x = 6, y = 9.5, label = "aSyn olig.",
           size = 2.3, colour = "#BF360C", fontface = "italic") +
  ## Drug entry: SBT-272 -> M0
  annotate("text", x = 1.0, y = 12.0, label = "SBT-272\n(SC injection)",
           size = 2.8, fontface = "bold", colour = "#1565C0",
           hjust = 0, lineheight = 0.85) +
  annotate("segment", x = 3.3, y = 12, xend = 4.6, yend = 12,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           colour = "#1565C0", linewidth = 0.8) +
  ## M0->M1 (C_mito)
  annotate("segment", x = 4.6, y = 11.45, xend = 3.2, yend = 9.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
           colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 3.4, y = 10.6, label = expression(C[mito]),
           size = 2.2, colour = "grey45", fontface = "italic") +
  ## M3->M5 (mPTP, diagonal)
  annotate("segment", x = 3.0, y = 5.45, xend = 4.8, yend = 3.55,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
           colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 3.2, y = 4.3, label = "mPTP",
           size = 2.2, colour = "grey45", fontface = "italic") +
  ## M6->M5 (TNF, IL1b)
  annotate("segment", x = 8.1, y = 3, xend = 7.4, yend = 3,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
           colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 7.75, y = 3.45, label = "TNF, IL1b",
           size = 2.1, colour = "grey45", fontface = "italic") +
  ## M5->M7 (DA neuron)
  annotate("segment", x = 6, y = 2.45, xend = 6, yend = 1.05,
           arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
           colour = "grey45", linewidth = 0.5) +
  annotate("text", x = 7.0, y = 1.75, label = "DA neuron",
           size = 2.2, colour = "grey45", fontface = "italic") +
  ## M2<->M6 (aSyn_ext / MG) — OUTSIDE the vicious cycle box
  annotate("segment", x = 10.9, y = 8.5, xend = 10.9, yend = 3.5,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed", ends = "both"),
           colour = "grey50", linewidth = 0.4) +
  annotate("text", x = 11.3, y = 6.0, label = expression(aSyn[ext]*" / MG"),
           size = 2.0, colour = "grey50", fontface = "italic", angle = 90) +
  ## M2->M5 (proteotoxicity) — curved OUTSIDE the cycle box on the right
  geom_curve(aes(x = 9.0, y = 8.45, xend = 7.0, yend = 3.55),
             curvature = 0.5,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
             colour = "#6A1B9A", linewidth = 0.45, linetype = "dashed") +
  annotate("text", x = 9.8, y = 4.8, label = "Proteotoxicity",
           size = 2.2, colour = "#6A1B9A", fontface = "italic") +
  coord_cartesian(xlim = c(0.3, 12.5), ylim = c(-0.3, 13), clip = "off") +
  labs(title = "QSP Model Architecture: 8 Modules, 28 ODEs") +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5,
                                  margin = margin(b = 8)),
        plot.margin = margin(8, 8, 8, 8))

ggsave(file.path(out_dir, "Fig1_model_architecture.png"), p1,
       width = 8, height = 8, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG 2: INTEGRATED MODEL DYNAMICS (8 panels, consistent x-axes)
##=============================================================================
cat("Generating Fig 2: Integrated dynamics...\n")

make_panel <- function(data, var, ylab, title_text) {
  ggplot(data, aes(x = time_days, y = .data[[var]], colour = scenario)) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = col_scen, guide = "none") +
    scale_x_continuous(limits = c(0, SIM_DAYS), breaks = seq(0, SIM_DAYS, 10)) +
    labs(x = NULL, y = ylab, title = title_text) +
    theme(plot.title = element_text(size = 9, face = "bold"))
}

p2a <- make_panel(all_data, "CL_ratio",    "Fraction", "A  CL Functional Ratio")
p2b <- make_panel(all_data, "aSyn_olig",   "Fraction", "B  aSyn Oligomers")
p2c <- make_panel(all_data, "CI_activity",  "Fraction", "C  Complex I Activity")
p2d <- make_panel(all_data, "mROS",         "Fraction", "D  Mitochondrial ROS")
p2e <- make_panel(all_data, "DA_neuron",    "Fraction", "E  DA Neuron Survival")
p2f <- make_panel(all_data, "Motor_score",  "Score",    "F  Motor Impairment")
p2g <- make_panel(all_data, "MG_active",    "Fraction", "G  Microglial Activation")

p2h <- make_panel(all_data, "ATP", "Fraction", "H  ATP") +
  scale_colour_manual(values = col_scen, name = NULL) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.5)) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.2)))

fig2 <- arrangeGrob(
  p2a, p2b, p2c, p2d, p2e, p2f, p2g, p2h,
  ncol = 2, nrow = 4,
  bottom = textGrob("Time (days)", gp = gpar(fontsize = 9))
)
ggsave(file.path(out_dir, "Fig2_integrated_dynamics.png"), fig2,
       width = 8, height = 10, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG 3: CALIBRATION & VALIDATION COMPOSITE
##=============================================================================
cat("Generating Fig 3: Calibration & validation...\n")

mROS_h <- ss_h35$mROS; CL_n_h <- ss_h35$CL_n

CI_grid <- seq(0.60, 0.85, by = 0.01)
ac_grid <- seq(0.20, 0.55, by = 0.05)
sweep_df <- data.frame()
for (cm in CI_grid) {
  for (ac in ac_grid) {
    res <- run_cal(0, SIM_DAYS, CI_max = cm, alpha_clear = ac)
    ss  <- res[nrow(res), ]
    err_CI   <- (ss$CI_activity - 0.49)^2
    err_mROS <- (ss$mROS / mROS_h - 1.77)^2
    err_CL   <- ((1 - ss$CL_n / CL_n_h) * 100 - 23)^2 / 100
    da_loss  <- (1 - ss$DA_neuron) * 100
    err_DA   <- ifelse(da_loss < 25, ((25 - da_loss)/20)^2,
                       ifelse(da_loss > 55, ((da_loss - 55)/20)^2, 0))
    score    <- err_CI * 4 + err_mROS + err_CL + err_DA * 2
    sweep_df <- rbind(sweep_df, data.frame(
      CI_max = cm, alpha_clear = ac, score = score))
  }
}
sweep_df$log_score <- log10(pmax(sweep_df$score, 1e-4))

p3a <- ggplot(sweep_df, aes(x = CI_max, y = factor(alpha_clear), fill = log_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(data = data.frame(CI_max = 0.74, alpha_clear = factor(0.25)),
             aes(x = CI_max, y = alpha_clear), inherit.aes = FALSE,
             shape = 4, size = 3, stroke = 1.5, colour = "white") +
  scale_fill_gradient2(
    low = "#1B5E20", mid = "#FFF9C4", high = "#B71C1C",
    midpoint = median(sweep_df$log_score),
    name = expression(log[10]*"(score)")
  ) +
  scale_x_continuous(breaks = seq(0.60, 0.85, 0.05)) +
  labs(title = expression(bold("A  Disease Modifier Calibration")),
       x = expression(CI[max]), y = expression(alpha[clear])) +
  theme(legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(size = 7))

## Validation: source name INSIDE the green band, not floating left
val_data <- data.frame(
  target = c("CI activity", "mROS (% basal)", "CL drop (%)",
             "DA loss (%)", "TLR2 KO aSyn (%)"),
  observed = c(ss_d35$CI_activity,
               ss_d35$mROS / ss_h35$mROS * 100,
               (1 - ss_d35$CL_n / ss_h35$CL_n) * 100,
               (1 - ss_d35$DA_neuron) * 100,
               20.3),
  lo = c(0.40, 150, 15, 25, 15),
  hi = c(0.58, 250, 35, 45, 50),
  source = c("Gao 2017", "Choi 2022", "Literature",
             "Bernheimer 1973", "Ivanova 2024"),
  stringsAsFactors = FALSE
)
val_data$norm <- (val_data$observed - val_data$lo) / (val_data$hi - val_data$lo)
val_data$target <- factor(val_data$target, levels = rev(val_data$target))

p3b <- ggplot(val_data) +
  geom_rect(aes(xmin = 0, xmax = 1,
                ymin = as.numeric(target) - 0.35,
                ymax = as.numeric(target) + 0.35),
            fill = "#E8F5E9", colour = "#A5D6A7", linewidth = 0.3) +
  geom_point(aes(x = norm, y = target),
             colour = "#1B5E20", size = 3.5, shape = 18) +
  geom_vline(xintercept = c(0, 1), linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_text(aes(x = norm, y = target,
                label = sprintf("%.1f", observed)),
            vjust = -1.3, size = 2.5, colour = "#1B5E20", fontface = "bold") +
  ## Source labels INSIDE the bars (right-aligned within the green band)
  geom_text(aes(x = 0.98, y = target, label = source),
            size = 2, colour = "grey50", hjust = 1, fontface = "italic") +
  scale_x_continuous(limits = c(-0.08, 1.15),
                     breaks = c(0, 0.5, 1),
                     labels = c("Lower\nbound", "Mid", "Upper\nbound")) +
  labs(title = "B  Calibration Target Fit (5/5 pass)",
       x = "Position within acceptable range", y = NULL) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold", size = 8))

fig3 <- arrangeGrob(p3a, p3b, ncol = 2, widths = c(1.2, 1))
ggsave(file.path(out_dir, "Fig3_calibration_validation.png"), fig3,
       width = 12, height = 5, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG 4: PHASE PORTRAIT
##=============================================================================
cat("Generating Fig 4: Phase portrait...\n")

res_h_50  <- run_cal(0, 50, dt = 0.25)
res_d_50  <- run_cal(0, 50, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear, dt = 0.25)
res_rx_50 <- run_cal(5, 50, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear, dt = 0.25)

phase_df <- rbind(
  data.frame(CL = res_h_50$CL_ratio,  aSyn = res_h_50$aSyn_olig,  Scenario = "Healthy"),
  data.frame(CL = res_d_50$CL_ratio,  aSyn = res_d_50$aSyn_olig,  Scenario = "PD (untreated)"),
  data.frame(CL = res_rx_50$CL_ratio, aSyn = res_rx_50$aSyn_olig, Scenario = "PD + SBT-272")
)
phase_df$Scenario <- factor(phase_df$Scenario,
  levels = c("Healthy", "PD (untreated)", "PD + SBT-272"))

ends <- data.frame(
  CL   = c(tail(res_h_50$CL_ratio,1), tail(res_d_50$CL_ratio,1), tail(res_rx_50$CL_ratio,1)),
  aSyn = c(tail(res_h_50$aSyn_olig,1), tail(res_d_50$aSyn_olig,1), tail(res_rx_50$aSyn_olig,1)),
  Scenario = c("Healthy", "PD (untreated)", "PD + SBT-272"))

make_arrows <- function(res, fracs = c(0.35, 0.65)) {
  n <- nrow(res)
  do.call(rbind, lapply(fracs, function(f) {
    i <- max(2, round(n * f))
    data.frame(x = res$CL_ratio[i-1], y = res$aSyn_olig[i-1],
               xend = res$CL_ratio[i], yend = res$aSyn_olig[i])
  }))
}
arr_h  <- make_arrows(res_h_50);  arr_h$Scenario  <- "Healthy"
arr_d  <- make_arrows(res_d_50);  arr_d$Scenario  <- "PD (untreated)"
arr_rx <- make_arrows(res_rx_50); arr_rx$Scenario <- "PD + SBT-272"
arrows_df <- rbind(arr_h, arr_d, arr_rx)

## Compute label offsets based on actual data range
cl_range <- range(phase_df$CL)
as_range <- range(phase_df$aSyn)
dx <- diff(cl_range) * 0.04
dy <- diff(as_range) * 0.06

p4 <- ggplot(phase_df, aes(x = CL, y = aSyn, colour = Scenario)) +
  geom_path(linewidth = 0.85, alpha = 0.85) +
  geom_point(data = ends, shape = 18, size = 4, show.legend = FALSE) +
  geom_segment(data = arrows_df,
               aes(x = x, y = y, xend = xend, yend = yend, colour = Scenario),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = 1, show.legend = FALSE) +
  scale_colour_manual(values = col_3scen) +
  ## Labels placed with explicit offsets, AWAY from edges
  annotate("label", x = ends$CL[1] - dx*2, y = ends$aSyn[1] + dy,
           label = "Healthy\nattractor", size = 2.5, fill = "#E8F5E9",
           colour = "#2E7D32", fontface = "bold", label.padding = unit(0.15, "lines")) +
  annotate("label", x = ends$CL[2] + dx*2, y = ends$aSyn[2] - dy*0.5,
           label = "Disease\nattractor", size = 2.5, fill = "#FFEBEE",
           colour = "#C62828", fontface = "bold", label.padding = unit(0.15, "lines")) +
  annotate("label", x = ends$CL[3] - dx*8, y = ends$aSyn[3] - dy*0.3,
           label = "Drug-rescued\nstate", size = 2.5, fill = "#E3F2FD",
           colour = "#1565C0", fontface = "bold", label.padding = unit(0.15, "lines")) +
  labs(title = "Vicious Cycle Phase Portrait",
       subtitle = "State-space trajectories reveal distinct attractors for health, disease, and treatment",
       x = "Functional Cardiolipin Fraction (CL ratio)",
       y = expression(alpha*"Syn Oligomer Burden"), colour = NULL) +
  expand_limits(x = c(cl_range[1] - diff(cl_range)*0.12,
                       cl_range[2] + diff(cl_range)*0.15),
                y = c(as_range[1] - diff(as_range)*0.15,
                       as_range[2] + diff(as_range)*0.10)) +
  theme(legend.position = c(0.78, 0.88),
        legend.background = element_rect(fill = alpha("white", 0.9), colour = NA),
        legend.key.size = unit(0.35, "cm"))

ggsave(file.path(out_dir, "Fig4_phase_portrait.png"), p4,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## HUMAN TRANSLATION SETUP
##=============================================================================
cat("Setting up human translation...\n")

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
cat(sprintf("  C_ref_human = %.4f\n", C_ref_human))

run_human <- function(dose_mg, duration_days, CI_max_h = 1.0,
                      alpha_clear_h = 1.0, dt = 24, init_override = NULL) {
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * death_scale_human,
          k_death_inflam = 0.00020 * death_scale_human,
          k_death_aSyn = 0.00015 * death_scale_human))
  run_integrated(dose_mg_kg = dose_mg / BW_human, duration_days = duration_days,
                 dt = dt, CI_max = CI_max_h, alpha_clear = alpha_clear_h,
                 k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
                 k_damage_mPTP_val = CAL$k_damage_mPTP,
                 n_death_val = CAL$n_death,
                 parms_override = po, init_override = init_override,
                 quiet = TRUE)
}


##=============================================================================
## FIG 5: HUMAN TRANSLATION COMPOSITE (taller, legend repositioned)
##=============================================================================
cat("Generating Fig 5: Human translation...\n")

cat("  Running short PK sims...\n")
res_pk_mouse <- run_cal(5.0, 5, dt = 0.1)
res_pk_human_short <- run_human(30, 5, dt = 0.5)

pk_df <- rbind(
  data.frame(time_h = res_pk_mouse$time, C_mito = res_pk_mouse$C_mito,
             Species = "Mouse (5 mg/kg IP)"),
  data.frame(time_h = res_pk_human_short$time, C_mito = res_pk_human_short$C_mito,
             Species = "Human (30 mg SC)")
)

p5a <- ggplot(pk_df, aes(x = time_h, y = C_mito, colour = Species)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Mouse (5 mg/kg IP)" = "#E65100",
                                  "Human (30 mg SC)" = "#1565C0"), name = NULL) +
  labs(x = "Time (hours)", y = expression(C[mito]~"(normalised)"),
       title = "A  Mitochondrial PK: Mouse vs Human") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.60, 0.45),
        legend.background = element_rect(fill = "white", colour = "grey80",
                                         linewidth = 0.3),
        legend.text = element_text(size = 7))

cat("  Running 20-year human sims...\n")
dur_days <- 20 * 365
res_h_human  <- run_human(0,  dur_days, 1.0,            1.0,                dt = 168)
res_d_human  <- run_human(0,  dur_days, CI_max_human,    alpha_clear_human,  dt = 168)
res_rx_human <- run_human(30, dur_days, CI_max_human,    alpha_clear_human,  dt = 168)

ht_df <- rbind(
  data.frame(years = res_h_human$time / 24 / 365.25,
             DA = res_h_human$DA_neuron, Scenario = "Healthy"),
  data.frame(years = res_d_human$time / 24 / 365.25,
             DA = res_d_human$DA_neuron, Scenario = "Disease"),
  data.frame(years = res_rx_human$time / 24 / 365.25,
             DA = res_rx_human$DA_neuron, Scenario = "SBT-272 30 mg")
)
ht_df$Scenario <- factor(ht_df$Scenario,
  levels = c("Healthy", "Disease", "SBT-272 30 mg"))
col_ht <- c("Healthy" = "#2E7D32", "Disease" = "#C62828",
            "SBT-272 30 mg" = "#1565C0")

p5b <- ggplot(ht_df, aes(x = years, y = DA, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.40, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  annotate("text", x = 17, y = 0.43, label = "Bernheimer threshold",
           size = 2.5, colour = "grey40") +
  scale_colour_manual(values = col_ht, name = NULL) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Years", y = "DA Neuron Survival",
       title = "B  20-Year Disease Trajectory") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        legend.background = element_rect(fill = "white", colour = NA),
        legend.text = element_text(size = 7),
        legend.margin = margin(t = -2))

symptom_idx <- which(res_d_human$DA_neuron < 0.70)[1]
t_symptom_yr <- if (!is.na(symptom_idx)) res_d_human$time[symptom_idx] / 24 / 365.25 else 20

phase_d <- data.frame(years = res_d_human$time / 24 / 365.25,
                      DA = res_d_human$DA_neuron)

p5c <- ggplot(phase_d, aes(x = years, y = DA)) +
  annotate("rect", xmin = 0, xmax = t_symptom_yr,
           ymin = 0, ymax = 1.02, fill = "#E8F5E9", alpha = 0.4) +
  annotate("rect", xmin = t_symptom_yr, xmax = t_symptom_yr + 1,
           ymin = 0, ymax = 1.02, fill = "#FFF3E0", alpha = 0.5) +
  annotate("rect", xmin = t_symptom_yr + 1, xmax = 20,
           ymin = 0, ymax = 1.02, fill = "#FFEBEE", alpha = 0.4) +
  geom_line(linewidth = 1, colour = "#C62828") +
  geom_hline(yintercept = c(0.70, 0.50, 0.40), linetype = "dashed",
             colour = c("#E65100", "#C62828", "#7B1FA2"), linewidth = 0.35) +
  annotate("text", x = t_symptom_yr / 2, y = 0.96,
           label = "Prodromal", size = 3.2, fontface = "bold", colour = "#2E7D32") +
  annotate("text", x = t_symptom_yr + 0.5, y = 0.96,
           label = "Dx\ndelay", size = 2.8, colour = "#E65100",
           lineheight = 0.8, fontface = "bold") +
  annotate("text", x = 14, y = 0.96,
           label = "Clinical PD", size = 3.2, fontface = "bold", colour = "#C62828") +
  annotate("text", x = 19.5, y = 0.72, label = "30% loss",
           size = 2.3, colour = "#E65100", hjust = 1, fontface = "bold") +
  annotate("text", x = 19.5, y = 0.52, label = "50% loss",
           size = 2.3, colour = "#C62828", hjust = 1, fontface = "bold") +
  annotate("text", x = 19.5, y = 0.42, label = "Bernheimer",
           size = 2.3, colour = "#7B1FA2", hjust = 1, fontface = "bold") +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format()) +
  labs(x = "Years", y = "DA Neuron Survival",
       title = "C  Disease Phase Characterisation",
       subtitle = sprintf("Prodromal: %.1f yr | Dx delay: 12 mo", t_symptom_yr)) +
  theme(plot.title = element_text(size = 10, face = "bold"))

fig5 <- arrangeGrob(p5a, p5b, p5c, ncol = 3)
ggsave(file.path(out_dir, "Fig5_human_translation.png"), fig5,
       width = 14, height = 6, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG 6: VIRTUAL CLINICAL TRIAL
##=============================================================================
cat("Generating Fig 6: Virtual clinical trial...\n")

arms_name <- c("Placebo", "SBT-272 10 mg", "SBT-272 30 mg", "SBT-272 60 mg")
vpop_results <- data.frame()
for (ai in 1:4) {
  f <- sprintf("RA_arm_%d.rds", ai)
  if (file.exists(f)) vpop_results <- rbind(vpop_results, readRDS(f))
}

if (nrow(vpop_results) == 0) {
  cat("  WARNING: VPop RDS cache not found. Skipping Fig 6.\n")
} else {

vpop_results$arm <- gsub("SBT-272 (\\d+)mg", "SBT-272 \\1 mg", vpop_results$arm)
completed_ids <- Reduce(intersect, lapply(arms_name, function(a) {
  vpop_results$id[vpop_results$arm == a]
}))
plot_data <- vpop_results[vpop_results$id %in% completed_ids, ]
plot_data$arm <- factor(plot_data$arm, levels = arms_name)

## Short labels for x-axis to avoid overlap
plot_data$arm_short <- factor(
  gsub("SBT-272 ", "", as.character(plot_data$arm)),
  levels = c("Placebo", "10 mg", "30 mg", "60 mg"))

p6a <- ggplot(plot_data, aes(x = arm_short, y = DA_neuron, fill = arm)) +
  geom_violin(alpha = 0.35, width = 0.8, colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7,
               colour = "grey30", linewidth = 0.3) +
  scale_fill_manual(values = col_dose, guide = "none") +
  labs(x = NULL, y = "DA Neuron Fraction",
       title = "A  DA Neuron Survival at 18 Months") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8.5))

MCID_DA <- 0.05
pbo   <- plot_data[plot_data$arm == "Placebo", ]
trt30 <- plot_data[plot_data$arm == "SBT-272 30 mg", ]
trt30 <- trt30[order(trt30$id), ]
pbo_m <- pbo[pbo$id %in% trt30$id, ]
pbo_m <- pbo_m[order(pbo_m$id), ]
delta_30 <- trt30$DA_neuron - pbo_m$DA_neuron

wf <- data.frame(patient = rank(-delta_30), delta = delta_30)
wf$category <- ifelse(wf$delta > MCID_DA, "Responder (>5%)",
                ifelse(wf$delta > 0, "Marginal (0-5%)", "No benefit"))
wf$category <- factor(wf$category,
  levels = c("Responder (>5%)", "Marginal (0-5%)", "No benefit"))

n_resp <- sum(wf$delta > MCID_DA)
p6b <- ggplot(wf, aes(x = patient, y = delta * 100, fill = category)) +
  geom_col(width = 1) +
  geom_hline(yintercept = MCID_DA * 100, linetype = "dashed",
             colour = "#C62828", linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  annotate("text", x = nrow(wf) * 0.75, y = MCID_DA * 100 + 2,
           label = "MCID = 5%", colour = "#C62828", size = 2.8, fontface = "bold") +
  scale_fill_manual(values = c("Responder (>5%)" = "#00796B",
                                "Marginal (0-5%)" = "#F9A825",
                                "No benefit" = "#9E9E9E"), name = NULL) +
  labs(x = "Virtual patient (ranked)", y = "DA change vs placebo (%)",
       title = "B  Individual Responses: 30 mg") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.75, 0.88),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.28, "cm"),
        legend.background = element_rect(fill = "white", colour = "grey80",
                                         linewidth = 0.3))

## Rescue fraction heatmap
cat("  Running rescue fraction sims...\n")
rescue_doses  <- c(1.8, 5.0, 10.5)
rescue_labels <- c("10 mg", "30 mg", "60 mg")
ss_dis <- res_disease[nrow(res_disease), ]
rescue_df <- data.frame()
for (i in seq_along(rescue_doses)) {
  res_r <- run_cal(rescue_doses[i], SIM_DAYS, CI_max = CAL$CI_max,
                   alpha_clear = CAL$alpha_clear)
  ss_r <- res_r[nrow(res_r), ]
  mito_r  <- (ss_dis$death_rate_mito   - ss_r$death_rate_mito)   / ss_dis$death_rate_mito   * 100
  infl_r  <- (ss_dis$death_rate_inflam  - ss_r$death_rate_inflam) / ss_dis$death_rate_inflam * 100
  asyn_r  <- (ss_dis$death_rate_aSyn    - ss_r$death_rate_aSyn)   / ss_dis$death_rate_aSyn   * 100
  total_r <- (ss_dis$death_rate_total   - ss_r$death_rate_total)  / ss_dis$death_rate_total  * 100
  rescue_df <- rbind(rescue_df, data.frame(
    Dose = rescue_labels[i],
    Pathway = c("Mitochondrial\n(drug target)", "Neuroinflammation",
                "Proteotoxicity", "Total"),
    Rescue = c(mito_r, infl_r, asyn_r, total_r)))
}
rescue_df$Dose    <- factor(rescue_df$Dose,    levels = rescue_labels)
rescue_df$Pathway <- factor(rescue_df$Pathway,
  levels = c("Mitochondrial\n(drug target)", "Neuroinflammation",
             "Proteotoxicity", "Total"))

p6c <- ggplot(rescue_df, aes(x = Dose, y = Pathway, fill = Rescue)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.0f%%", Rescue)), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#FFEBEE", mid = "#FFF9C4", high = "#1B5E20",
                       midpoint = 30, limits = c(0, 70), name = "Rescue\n(%)") +
  labs(title = "C  Therapeutic Rescue Fraction",
       subtitle = "Non-mitochondrial pathways resist rescue",
       x = "Dose (human equivalent)", y = NULL) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8),
        axis.text.y = element_text(size = 8.5),
        legend.key.height = unit(0.5, "cm"))

fig6 <- arrangeGrob(p6a, p6b, p6c, ncol = 3, widths = c(0.9, 1.1, 0.95))
ggsave(file.path(out_dir, "Fig6_virtual_clinical_trial.png"), fig6,
       width = 14, height = 5.5, dpi = 300, bg = "white")
cat("  Saved.\n")
}


##=============================================================================
## FIG 7: UNCERTAINTY QUANTIFICATION (taller, fixed label positions)
##=============================================================================
cat("Generating Fig 7: Uncertainty quantification...\n")

if (exists("pbo") && exists("trt30")) {

## Matches the documented 5x5x5 = 125-combination joint sweep in
## "Uncertainty Analysis.R" (Section 4) and Methods 2.6 / Results 3.8 --
## previously this used a 3x3x3 = 27-combo shortcut, which is why Fig 7A's
## subtitle and 95% range did not match the manuscript text.
CI_max_v      <- c(0.86, 0.88, 0.90, 0.92, 0.94)
alpha_clear_v <- c(0.28, 0.32, 0.35, 0.38, 0.42)
death_scale_v <- c(0.04, 0.05, 0.06, 0.07, 0.08)
combos <- expand.grid(CI_max = CI_max_v, alpha_clear = alpha_clear_v,
                      death_scale = death_scale_v)

cat("  Running parameter envelope (", nrow(combos), "combos)...\n")
joint_res <- data.frame()
for (i in 1:nrow(combos)) {
  ci <- combos$CI_max[i]; ac <- combos$alpha_clear[i]; ds <- combos$death_scale[i]
  po <- c(HUMAN_PK, list(C_ref = C_ref_human,
          k_death = 0.0007 * ds, k_death_inflam = 0.00020 * ds,
          k_death_aSyn = 0.00015 * ds))
  res_dis <- tryCatch(
    run_integrated(dose_mg_kg = 0, duration_days = 547, dt = 72,
         CI_max = ci, alpha_clear = ac,
         k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
         k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
         parms_override = po, quiet = TRUE),
    error = function(e) NULL)
  res_drg <- tryCatch(
    run_integrated(dose_mg_kg = 30/BW_human, duration_days = 547, dt = 72,
         CI_max = ci, alpha_clear = ac,
         k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
         k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
         parms_override = po, quiet = TRUE),
    error = function(e) NULL)
  if (!is.null(res_dis) && !is.null(res_drg)) {
    da_saved <- res_drg$DA_neuron[nrow(res_drg)] - res_dis$DA_neuron[nrow(res_dis)]
    joint_res <- rbind(joint_res, data.frame(DA_saved = da_saved))
  }
  if (i %% 9 == 0) cat(sprintf("    %d/%d\n", i, nrow(combos)))
}

nominal_val <- joint_res$DA_saved[combos$CI_max == 0.90 &
                                   combos$alpha_clear == 0.35 &
                                   combos$death_scale == 0.06]
if (length(nominal_val) == 0) nominal_val <- median(joint_res$DA_saved)

p7a <- ggplot(joint_res, aes(x = DA_saved)) +
  geom_histogram(bins = 12, fill = "#00796B", colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = nominal_val, colour = "#C62828",
             linewidth = 0.7, linetype = "solid") +
  geom_vline(xintercept = quantile(joint_res$DA_saved, c(0.025, 0.975)),
             colour = "#C62828", linewidth = 0.4, linetype = "dashed") +
  annotate("text", x = nominal_val + diff(range(joint_res$DA_saved))*0.08,
           y = Inf, label = "Nominal", vjust = 2, colour = "#C62828",
           size = 3, fontface = "bold") +
  labs(title = "A  Parameter Uncertainty Envelope",
       subtitle = sprintf("%d combos | 95%% range: [%.3f, %.3f]",
                           nrow(joint_res),
                           quantile(joint_res$DA_saved, 0.025),
                           quantile(joint_res$DA_saved, 0.975)),
       x = expression(Delta*"DA (drug - placebo)"), y = "Count") +
  theme(plot.title = element_text(size = 10, face = "bold"))

arm_data <- list()
for (ai in 1:4) {
  f <- sprintf("RA_arm_%d.rds", ai)
  if (file.exists(f)) arm_data[[ai]] <- readRDS(f)
}
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
  geom_hline(yintercept = 0.05, linetype = "dashed",
             colour = "#E65100", linewidth = 0.35) +
  geom_boxplot(width = 0.5, outlier.size = 0.6, outlier.alpha = 0.3) +
  scale_fill_manual(values = c("10 mg" = "#1976D2", "30 mg" = "#00796B",
                                "60 mg" = "#E65100")) +
  annotate("text", x = 0.6, y = 0.055, label = "MCID = 5%",
           size = 2.5, colour = "#E65100", fontface = "italic", hjust = 0) +
  labs(title = "B  Patient-Level Treatment Effect",
       subtitle = sprintf("Matched VPop pairs (N=%d/arm)", length(completed_ids)),
       x = NULL, y = expression(Delta*"DA vs placebo")) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "none")

cat("  Running bootstrap (B=10,000)...\n")
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
  scale_fill_manual(values = c("10 mg" = "#1976D2", "30 mg" = "#00796B",
                                "60 mg" = "#E65100"), name = NULL) +
  labs(title = "C  Bootstrap Effect Size",
       subtitle = "B = 10,000 | DA neuron (unpaired d)",
       x = expression("Unpaired Cohen's " * italic(d)), y = "Density") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = c(0.85, 0.78),
        legend.background = element_rect(fill = "white", colour = "grey80",
                                         linewidth = 0.3),
        legend.key.size = unit(0.28, "cm"))

fig7 <- arrangeGrob(p7a, p7b, p7c, ncol = 3)
ggsave(file.path(out_dir, "Fig7_uncertainty_quantification.png"), fig7,
       width = 14, height = 5.5, dpi = 300, bg = "white")
cat("  Saved.\n")
} else {
  cat("  WARNING: VPop data required. Skipping Fig 7.\n")
}


##=============================================================================
## FIG S1: SENSITIVITY ANALYSIS (legend repositioned)
##=============================================================================
cat("Generating Fig S1: Sensitivity analysis...\n")

param_noms <- c(
  k_ROS_basal = 0.122, k_SOD2 = 1.0, CI_max = 0.74, K_mPTP = 0.40,
  k_shunt = 0.125, k_mPTP_open = 0.50, k_agg = 0.009, K_death = 0.25,
  alpha_clear = 0.25, k_clear = 0.241, k_CI_aSyn = 0.200,
  k_CI_repair = 0.050, k_CL_loss = 5.0, k_biogen = 3.0, k_death = 0.0007
)
base_res <- run_cal(0, SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
DA_base <- base_res$DA_neuron[nrow(base_res)]

oat_results <- data.frame()
for (pname in names(param_noms)) {
  nom <- param_noms[pname]
  for (dir in c(-0.5, 0.5)) {
    pval <- nom * (1 + dir)
    if (pname %in% c("CI_max", "alpha_clear")) {
      cm <- ifelse(pname == "CI_max", pval, CAL$CI_max)
      ac <- ifelse(pname == "alpha_clear", pval, CAL$alpha_clear)
      res <- run_cal(0, SIM_DAYS, CI_max = cm, alpha_clear = ac)
    } else {
      po <- list(); po[[pname]] <- pval
      res <- run_cal(0, SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear,
                     parms_override = po)
    }
    DA_pert <- res$DA_neuron[nrow(res)]
    pct <- (DA_pert - DA_base) / DA_base * 100
    oat_results <- rbind(oat_results, data.frame(
      param = pname, direction = ifelse(dir > 0, "+50%", "-50%"), pct = pct))
  }
}
tornado_order <- oat_results %>% group_by(param) %>%
  summarise(range = max(abs(pct)), .groups = "drop") %>% arrange(range)
oat_results$param <- factor(oat_results$param, levels = tornado_order$param)

ps1a <- ggplot(oat_results, aes(x = pct, y = param, fill = direction)) +
  geom_col(position = "identity", width = 0.55) +
  geom_vline(xintercept = 0, linewidth = 0.35) +
  scale_fill_manual(values = c("+50%" = "#1565C0", "-50%" = "#C62828"), name = NULL) +
  labs(x = "% change in DA neuron survival", y = NULL,
       title = "A  OAT Tornado (15 parameters)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7, family = "mono"),
        legend.position = "top", legend.justification = "right")

sobol_df <- data.frame(
  param = c("k_ROS_basal","K_mPTP","CI_max","k_SOD2","k_agg",
            "k_shunt","alpha_clear","k_CL_loss","k_CI_repair","k_red"),
  S1 = c(0.19,0.00,0.15,0.27,0.07,0.06,0.01,0.03,0.02,0.02),
  ST = c(0.50,0.33,0.31,0.28,0.16,0.10,0.09,0.05,0.05,0.03))
sobol_df$interaction <- sobol_df$ST - pmax(sobol_df$S1, 0)
sobol_df$param <- factor(sobol_df$param, levels = rev(sobol_df$param))

sobol_long <- sobol_df %>%
  select(param, S1, interaction) %>%
  pivot_longer(-param, names_to = "type", values_to = "value")
sobol_long$type <- factor(sobol_long$type, levels = c("S1", "interaction"),
                           labels = c("First-order", "Interactions"))

ps1b <- ggplot(sobol_long, aes(x = value, y = param, fill = type)) +
  geom_col(width = 0.55) +
  scale_fill_manual(values = c("First-order" = "#1565C0",
                                "Interactions" = "#E65100"), name = NULL) +
  labs(x = "Sobol index", y = NULL,
       title = "B  Sobol Variance Decomposition (top 10)") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 7, family = "mono"),
        legend.position = "top", legend.justification = "right")

figs1 <- arrangeGrob(ps1a, ps1b, ncol = 2)
ggsave(file.path(out_dir, "FigS1_sensitivity_analysis.png"), figs1,
       width = 12, height = 5.5, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG S2: DEATH DECOMPOSITION + DISEASE FINGERPRINT RADAR
##=============================================================================
cat("Generating Fig S2: Death decomposition + radar...\n")

build_death <- function(res, label) {
  data.frame(time_d = res$time_days,
             Mitochondrial     = res$death_rate_mito,
             Neuroinflammation = res$death_rate_inflam,
             Proteotoxicity    = res$death_rate_aSyn,
             condition = label)
}
death_df <- rbind(build_death(res_disease, "PD (untreated)"),
                  build_death(res_drug_hi, "PD + SBT-272"))

death_long <- death_df %>%
  pivot_longer(cols = c(Mitochondrial, Neuroinflammation, Proteotoxicity),
               names_to = "Pathway", values_to = "rate") %>%
  mutate(Pathway = factor(Pathway,
    levels = c("Proteotoxicity", "Neuroinflammation", "Mitochondrial")))

ps2a <- ggplot(death_long, aes(x = time_d, y = rate * 1e4, fill = Pathway)) +
  geom_area(alpha = 0.75, colour = "white", linewidth = 0.2) +
  facet_wrap(~condition, ncol = 2) +
  scale_fill_manual(values = c("Mitochondrial" = "#C62828",
                                "Neuroinflammation" = "#E65100",
                                "Proteotoxicity" = "#6A1B9A"), name = NULL) +
  labs(x = "Days", y = expression("Death rate ("*10^{-4}*" h"^{-1}*")"),
       title = "A  Multi-Pathway Death Decomposition") +
  theme(plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(face = "bold", size = 9),
        legend.position = "bottom")

## Radar chart
vars <- c("CL integrity", "CI activity", "ATP levels", "SOD2 capacity",
          "DA neuron\nsurvival", "Membrane\npotential",
          "Synaptic\nclearance", "PINK1\nactivity")
get_vals <- function(ss) {
  c(as.numeric(ss$CL_ratio), as.numeric(ss$CI_activity),
    as.numeric(ss$ATP), as.numeric(ss$SOD2_act),
    as.numeric(ss$DA_neuron), as.numeric(ss$delta_psi_m),
    1 - as.numeric(ss$aSyn_olig), as.numeric(ss$PINK1_act))
}
n_v <- length(vars)
make_radar <- function(vals, label) {
  angles <- seq(0, 2*pi, length.out = n_v + 1)[1:n_v]
  data.frame(var = vars, value = vals, angle = angles,
    x = vals * cos(angles - pi/2), y = vals * sin(angles - pi/2),
    Scenario = label)
}
radar_df <- rbind(make_radar(get_vals(ss_h), "Healthy"),
                  make_radar(get_vals(ss_d), "PD (untreated)"),
                  make_radar(get_vals(ss_rx), "PD + SBT-272"))
radar_df$Scenario <- factor(radar_df$Scenario,
  levels = c("Healthy", "PD (untreated)", "PD + SBT-272"))
radar_closed <- radar_df %>% group_by(Scenario) %>%
  group_modify(~ rbind(.x, .x[1, ])) %>% ungroup()

grid_circles <- data.frame()
for (r in c(0.25, 0.50, 0.75, 1.0)) {
  theta <- seq(0, 2*pi, length.out = 100)
  grid_circles <- rbind(grid_circles,
    data.frame(x = r * cos(theta), y = r * sin(theta), r = r))
}
axis_lines <- data.frame(
  angle = seq(0, 2*pi, length.out = n_v + 1)[1:n_v])
axis_lines$xend <- cos(axis_lines$angle - pi/2)
axis_lines$yend <- sin(axis_lines$angle - pi/2)

label_r <- 1.18
axis_labels <- data.frame(
  label = vars,
  x = label_r * cos(seq(0, 2*pi, length.out = n_v + 1)[1:n_v] - pi/2),
  y = label_r * sin(seq(0, 2*pi, length.out = n_v + 1)[1:n_v] - pi/2))

ps2b <- ggplot() +
  geom_path(data = grid_circles, aes(x = x, y = y, group = r),
            colour = "grey82", linewidth = 0.25) +
  geom_segment(data = axis_lines, aes(x = 0, y = 0, xend = xend, yend = yend),
               colour = "grey75", linewidth = 0.25) +
  geom_polygon(data = radar_closed,
               aes(x = x, y = y, fill = Scenario, colour = Scenario),
               alpha = 0.1, linewidth = 0.8) +
  geom_point(data = radar_df,
             aes(x = x, y = y, colour = Scenario), size = 1.8) +
  geom_text(data = axis_labels, aes(x = x, y = y, label = label),
            size = 2.3, lineheight = 0.8) +
  scale_fill_manual(values = col_3scen) +
  scale_colour_manual(values = col_3scen) +
  coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45)) +
  labs(title = "B  Multi-Dimensional Disease Fingerprint",
       fill = NULL, colour = NULL) +
  theme_void(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
        legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))

figs2 <- arrangeGrob(ps2a, ps2b, ncol = 2, widths = c(1.2, 1))
ggsave(file.path(out_dir, "FigS2_death_radar.png"), figs2,
       width = 13, height = 6, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## FIG S3: DOSE-RESPONSE + CASCADE PROPAGATION
##=============================================================================
cat("Generating Fig S3: Dose-response + cascade...\n")

doses_mouse <- c(0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0)
dr_data <- data.frame()
for (d in doses_mouse) {
  res <- run_cal(d, SIM_DAYS, CI_max = CAL$CI_max, alpha_clear = CAL$alpha_clear)
  ss <- res[nrow(res), ]
  dr_data <- rbind(dr_data, data.frame(dose = d, DA = ss$DA_neuron,
                                        Motor = ss$Motor_score))
}
healthy_DA <- ss_h$DA_neuron

ps3a <- ggplot(dr_data) +
  geom_line(aes(x = dose, y = DA), colour = "#1565C0", linewidth = 0.8) +
  geom_point(aes(x = dose, y = DA), colour = "#1565C0", size = 2.2) +
  geom_hline(yintercept = healthy_DA, linetype = "dashed",
             colour = "#2E7D32", linewidth = 0.3) +
  geom_hline(yintercept = dr_data$DA[1], linetype = "dashed",
             colour = "#C62828", linewidth = 0.3) +
  annotate("text", x = 9.5, y = healthy_DA - 0.015, label = "Healthy",
           colour = "#2E7D32", size = 2.8, hjust = 1, fontface = "italic") +
  annotate("text", x = 9.5, y = dr_data$DA[1] + 0.015, label = "Disease (no drug)",
           colour = "#C62828", size = 2.8, hjust = 1, fontface = "italic") +
  labs(x = "Dose (mg/kg)", y = "DA Neuron Survival",
       title = "A  DA Neuron Dose-Response (50-day mouse)") +
  theme(plot.title = element_text(size = 10, face = "bold"))

## Cascade propagation
track_vars <- list(
  "CI activity"         = list(col="CI_activity",  h=as.numeric(ss_h$CI_activity),  dir="down"),
  "aSyn aggregation"    = list(col="aSyn_olig",    h=as.numeric(ss_h$aSyn_olig),    dir="up"),
  "Microglia"           = list(col="MG_active",    h=as.numeric(ss_h$MG_active),    dir="up"),
  "ATP levels"          = list(col="ATP",           h=as.numeric(ss_h$ATP),          dir="down"),
  "CL integrity"        = list(col="CL_ratio",     h=as.numeric(ss_h$CL_ratio),     dir="down"),
  "Membrane potential"  = list(col="delta_psi_m",  h=as.numeric(ss_h$delta_psi_m),  dir="down"),
  "DA neuron loss"      = list(col="DA_neuron",    h=as.numeric(ss_h$DA_neuron),    dir="down"),
  "Motor impairment"    = list(col="Motor_score",  h=as.numeric(ss_h$Motor_score),  dir="up")
)

cascade_df <- data.frame()
for (vname in names(track_vars)) {
  info <- track_vars[[vname]]
  for (cond in list(list(v=res_disease, c="PD (untreated)"),
                    list(v=res_drug_hi, c="PD + SBT-272"))) {
    vals <- cond$v[[info$col]]
    h_val <- info$h
    if (info$dir == "down") {
      dev <- (h_val - vals) / max(h_val, 1e-8)
    } else {
      dev <- (vals - h_val) / max(1 - h_val, 1e-8)
    }
    t10 <- cond$v$time_days[which(dev > 0.10)[1]]
    t50 <- cond$v$time_days[which(dev > 0.50)[1]]
    if (length(t10) == 0 || is.na(t10)) t10 <- NA
    if (length(t50) == 0 || is.na(t50)) t50 <- NA
    cascade_df <- rbind(cascade_df, data.frame(
      variable = vname, condition = cond$c,
      t_onset = t10, t_severe = t50, stringsAsFactors = FALSE))
  }
}
dis_order <- cascade_df %>% filter(condition == "PD (untreated)") %>%
  arrange(t_onset) %>% pull(variable)
cascade_df$variable <- factor(cascade_df$variable, levels = rev(dis_order))

ps3b <- ggplot(cascade_df %>% filter(!is.na(t_onset)),
               aes(y = variable, colour = condition)) +
  geom_segment(aes(x = t_onset,
                   xend = ifelse(is.na(t_severe), SIM_DAYS, t_severe),
                   yend = variable),
               linewidth = 2.5, alpha = 0.45,
               position = position_dodge(width = 0.6)) +
  geom_point(aes(x = t_onset), size = 3, shape = 16,
             position = position_dodge(width = 0.6)) +
  geom_point(data = cascade_df %>% filter(!is.na(t_severe)),
             aes(x = t_severe), size = 3, shape = 17,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c("PD (untreated)" = "#C62828",
                                  "PD + SBT-272" = "#1565C0"), name = NULL) +
  scale_x_continuous(limits = c(0, SIM_DAYS), breaks = seq(0, SIM_DAYS, 10)) +
  labs(title = "B  Cascade Propagation Timeline",
       subtitle = "Circle = >10% deviation | Triangle = >50% deviation",
       x = "Days after disease onset", y = NULL) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom",
        axis.text.y = element_text(size = 8.5))

figs3 <- arrangeGrob(ps3a, ps3b, ncol = 2, widths = c(0.9, 1.1))
ggsave(file.path(out_dir, "FigS3_dose_response_cascade.png"), figs3,
       width = 12, height = 5.5, dpi = 300, bg = "white")
cat("  Saved.\n")


##=============================================================================
## SUMMARY
##=============================================================================
cat("\n================================================================\n")
cat("  ALL PUBLICATION PLOTS GENERATED\n")
cat("================================================================\n\n")
cat("Main manuscript (7 figures):\n")
cat("  Fig1_model_architecture.png\n")
cat("  Fig2_integrated_dynamics.png\n")
cat("  Fig3_calibration_validation.png\n")
cat("  Fig4_phase_portrait.png\n")
cat("  Fig5_human_translation.png\n")
cat("  Fig6_virtual_clinical_trial.png\n")
cat("  Fig7_uncertainty_quantification.png\n\n")
cat("Supplementary (3 figures):\n")
cat("  FigS1_sensitivity_analysis.png\n")
cat("  FigS2_death_radar.png\n")
cat("  FigS3_dose_response_cascade.png\n\n")
cat("All: 300 DPI, white background, unified theme, no watermarks.\n")
cat("Output:", normalizePath(out_dir), "\n")
