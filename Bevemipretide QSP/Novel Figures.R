##=============================================================================
## Novel Publication Figures for Bevemipretide PD QSP Model
## Unique visualisations: phase portrait, death decomposition,
## radar fingerprint, cascade propagation
##=============================================================================

library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(grid)

source("Integrated Model.R", local = FALSE)

pub_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(colour = "grey40", size = 9),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
theme_set(pub_theme)

CAL <- list(k_impair = 0.35, K_death = 0.25, k_damage_mPTP = 4.0, n_death = 4)


##=============================================================================
## Run base simulations (mouse, 50 days for fuller trajectories)
##=============================================================================

cat("Running base simulations...\n")

res_h <- run_integrated(dose_mg_kg = 0, duration_days = 50, dt = 0.25,
            CI_max = 1.0, alpha_clear = 1.0,
            k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
            k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
            quiet = TRUE)

res_d <- run_integrated(dose_mg_kg = 0, duration_days = 50, dt = 0.25,
            CI_max = 0.74, alpha_clear = 0.25,
            k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
            k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
            quiet = TRUE)

res_rx <- run_integrated(dose_mg_kg = 5.0, duration_days = 50, dt = 0.25,
            CI_max = 0.74, alpha_clear = 0.25,
            k_impair_val = CAL$k_impair, K_death_val = CAL$K_death,
            k_damage_mPTP_val = CAL$k_damage_mPTP, n_death_val = CAL$n_death,
            quiet = TRUE)

cat("  Simulations complete.\n")


##=============================================================================
## FIGURE 13: Phase Portrait — Vicious Cycle Attractor Dynamics
##=============================================================================

cat("Generating Fig13: Phase Portrait...\n")

phase_df <- rbind(
  data.frame(CL_ratio = res_h$CL_ratio, aSyn = res_h$aSyn_olig,
             time_d = res_h$time / 24, Scenario = "Healthy"),
  data.frame(CL_ratio = res_d$CL_ratio, aSyn = res_d$aSyn_olig,
             time_d = res_d$time / 24, Scenario = "PD (untreated)"),
  data.frame(CL_ratio = res_rx$CL_ratio, aSyn = res_rx$aSyn_olig,
             time_d = res_rx$time / 24, Scenario = "PD + SBT-272")
)
phase_df$Scenario <- factor(phase_df$Scenario,
  levels = c("Healthy", "PD (untreated)", "PD + SBT-272"))

col_phase <- c("Healthy" = "#4CAF50", "PD (untreated)" = "#E53935",
               "PD + SBT-272" = "#1565C0")

## Arrow positions at ~40% and ~70% of each trajectory
make_arrows <- function(df, fracs = c(0.35, 0.65)) {
  n <- nrow(df)
  do.call(rbind, lapply(fracs, function(f) {
    i <- max(2, round(n * f))
    data.frame(x = df$CL_ratio[i-1], y = df$aSyn_olig[i-1],
               xend = df$CL_ratio[i], yend = df$aSyn_olig[i])
  }))
}

arr_h  <- make_arrows(res_h);  arr_h$Scenario  <- "Healthy"
arr_d  <- make_arrows(res_d);  arr_d$Scenario  <- "PD (untreated)"
arr_rx <- make_arrows(res_rx); arr_rx$Scenario <- "PD + SBT-272"
arrows <- rbind(arr_h, arr_d, arr_rx)

## Start and end markers
ends <- data.frame(
  CL_ratio = c(tail(res_h$CL_ratio,1), tail(res_d$CL_ratio,1), tail(res_rx$CL_ratio,1)),
  aSyn     = c(tail(res_h$aSyn_olig,1), tail(res_d$aSyn_olig,1), tail(res_rx$aSyn_olig,1)),
  Scenario = c("Healthy", "PD (untreated)", "PD + SBT-272")
)

p13 <- ggplot(phase_df, aes(x = CL_ratio, y = aSyn, colour = Scenario)) +
  geom_path(linewidth = 0.9, alpha = 0.85) +
  geom_point(data = ends, aes(x = CL_ratio, y = aSyn, colour = Scenario),
             shape = 18, size = 4.5, show.legend = FALSE) +
  geom_segment(data = arrows,
               aes(x = x, y = y, xend = xend, yend = yend, colour = Scenario),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               linewidth = 1.1, show.legend = FALSE) +
  scale_colour_manual(values = col_phase) +
  annotate("label", x = tail(res_h$CL_ratio,1), y = tail(res_h$aSyn_olig,1) + 0.015,
           label = "Healthy\nsteady state", size = 2.8, fill = "#E8F5E9",
           colour = "#2E7D32", label.size = 0.3, fontface = "bold") +
  annotate("label", x = tail(res_d$CL_ratio,1) + 0.015, y = tail(res_d$aSyn_olig,1),
           label = "Disease\nattractor", size = 2.8, fill = "#FFEBEE",
           colour = "#C62828", label.size = 0.3, fontface = "bold") +
  annotate("label", x = tail(res_rx$CL_ratio,1), y = tail(res_rx$aSyn_olig,1) - 0.018,
           label = "Drug-rescued\nstate", size = 2.8, fill = "#E3F2FD",
           colour = "#0D47A1", label.size = 0.3, fontface = "bold") +
  annotate("segment", x = 0.70, xend = 0.63, y = 0.08, yend = 0.15,
           arrow = arrow(length = unit(0.12, "cm")),
           colour = "grey50", linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = 0.72, y = 0.07, label = "Vicious cycle\npulls system here",
           size = 2.5, colour = "grey45", fontface = "italic") +
  labs(title = "Vicious Cycle Phase Portrait",
       subtitle = "State-space trajectories reveal distinct attractors for health, disease, and treatment",
       x = "Functional Cardiolipin Fraction (CL ratio)",
       y = expression(alpha*"-Synuclein Oligomer Burden"),
       colour = NULL) +
  theme(legend.position = c(0.82, 0.92),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_line(colour = "grey90"))

ggsave("Fig13_phase_portrait.png", p13, width = 7, height = 6, dpi = 300, bg = "white")
cat("  Saved: Fig13_phase_portrait.png\n")


##=============================================================================
## FIGURE 14: Death Pathway Decomposition Over Disease Course
##=============================================================================

cat("Generating Fig14: Death Pathway Decomposition...\n")

## Build decomposition from both disease and drug simulations
build_death_df <- function(res, label) {
  data.frame(
    time_d = res$time / 24,
    Mitochondrial = res$death_rate_mito,
    Neuroinflammation = res$death_rate_inflam,
    Proteotoxicity = res$death_rate_aSyn,
    condition = label
  )
}

death_df <- rbind(build_death_df(res_d, "PD (untreated)"),
                  build_death_df(res_rx, "PD + SBT-272"))

death_long <- death_df %>%
  pivot_longer(cols = c(Mitochondrial, Neuroinflammation, Proteotoxicity),
               names_to = "Pathway", values_to = "rate") %>%
  mutate(Pathway = factor(Pathway,
    levels = c("Mitochondrial", "Neuroinflammation", "Proteotoxicity")))

col_death <- c("Mitochondrial" = "#E53935", "Neuroinflammation" = "#FF9800",
               "Proteotoxicity" = "#7B1FA2")

p14 <- ggplot(death_long, aes(x = time_d, y = rate * 1000, fill = Pathway)) +
  geom_area(alpha = 0.75, colour = "white", linewidth = 0.3) +
  facet_wrap(~condition, ncol = 2) +
  scale_fill_manual(values = col_death) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(title = "Multi-Pathway Death Decomposition",
       subtitle = "Relative contributions shift over disease course; SBT-272 selectively attenuates the dominant mitochondrial pathway",
       x = "Days", y = expression("Death Rate ("*10^{-3}*" h"^{-1}*")"),
       fill = NULL) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.major = element_line(colour = "grey90"))

ggsave("Fig14_death_decomposition.png", p14, width = 10, height = 5, dpi = 300, bg = "white")
cat("  Saved: Fig14_death_decomposition.png\n")


##=============================================================================
## FIGURE 15: Multi-Dimensional Disease Fingerprint (Radar Chart)
##=============================================================================

cat("Generating Fig15: Disease Fingerprint Radar...\n")

## Extract steady-state values (last row)
ss_h  <- tail(res_h, 1)
ss_d  <- tail(res_d, 1)
ss_rx <- tail(res_rx, 1)

## Variables: orient so outer = healthy (normalise each)
vars <- c("CL\nintegrity", "CI\nactivity", "ATP\nlevels",
          "Antioxidant\n(SOD2)", "DA neuron\nsurvival",
          "Membrane\npotential", "Synaptic\nclearance", "Mitophagy\n(PINK1)")

## Raw values: for "good" variables use as-is; these are all in [0,1]
get_vals <- function(ss) {
  c(ss$CL_ratio, ss$CI_activity, ss$ATP, ss$SOD2_act,
    ss$DA_neuron, ss$delta_psi_m,
    1 - ss$aSyn_olig,   # invert: low oligomers = good clearance
    ss$PINK1_act)
}

vals_h  <- get_vals(ss_h)
vals_d  <- get_vals(ss_d)
vals_rx <- get_vals(ss_rx)

n_vars <- length(vars)

## Build radar data: close the polygon by repeating the first point
make_radar <- function(vals, label) {
  angles <- seq(0, 2 * pi, length.out = n_vars + 1)[1:n_vars]
  data.frame(
    var = factor(vars, levels = vars),
    value = vals,
    angle = angles,
    x = vals * cos(angles - pi/2),
    y = vals * sin(angles - pi/2),
    Scenario = label
  )
}

radar_df <- rbind(
  make_radar(vals_h, "Healthy"),
  make_radar(vals_d, "PD (untreated)"),
  make_radar(vals_rx, "PD + SBT-272")
)
radar_df$Scenario <- factor(radar_df$Scenario,
  levels = c("Healthy", "PD (untreated)", "PD + SBT-272"))

## Close the polygons
close_polygon <- function(df) {
  rbind(df, df[1, ])
}
radar_closed <- radar_df %>% group_by(Scenario) %>%
  group_modify(~ close_polygon(.x)) %>% ungroup()

col_radar <- c("Healthy" = "#4CAF50", "PD (untreated)" = "#E53935",
               "PD + SBT-272" = "#1565C0")

## Grid circles
grid_circles <- data.frame()
for (r in c(0.25, 0.50, 0.75, 1.0)) {
  theta <- seq(0, 2*pi, length.out = 100)
  grid_circles <- rbind(grid_circles,
    data.frame(x = r * cos(theta), y = r * sin(theta), r = r))
}

## Axis lines
axis_lines <- data.frame(
  angle = seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars]
)
axis_lines$xend <- cos(axis_lines$angle - pi/2)
axis_lines$yend <- sin(axis_lines$angle - pi/2)

## Axis labels
label_r <- 1.15
axis_labels <- data.frame(
  label = vars,
  x = label_r * cos(seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars] - pi/2),
  y = label_r * sin(seq(0, 2*pi, length.out = n_vars + 1)[1:n_vars] - pi/2)
)

p15 <- ggplot() +
  geom_path(data = grid_circles, aes(x = x, y = y, group = r),
            colour = "grey82", linewidth = 0.3) +
  geom_segment(data = axis_lines, aes(x = 0, y = 0, xend = xend, yend = yend),
               colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = radar_closed, aes(x = x, y = y, fill = Scenario, colour = Scenario),
               alpha = 0.12, linewidth = 0.9) +
  geom_point(data = radar_df, aes(x = x, y = y, colour = Scenario), size = 2.2) +
  geom_text(data = axis_labels, aes(x = x, y = y, label = label),
            size = 2.7, lineheight = 0.85) +
  scale_fill_manual(values = col_radar) +
  scale_colour_manual(values = col_radar) +
  coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
  labs(title = "Multi-Dimensional Disease Fingerprint",
       subtitle = "Outer ring = optimal function | Disease collapses multiple axes; SBT-272 partially restores",
       fill = NULL, colour = NULL) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.subtitle = element_text(colour = "grey40", size = 9, hjust = 0.5),
        legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 10))

ggsave("Fig15_disease_fingerprint.png", p15, width = 7.5, height = 7.5, dpi = 300, bg = "white")
cat("  Saved: Fig15_disease_fingerprint.png\n")


##=============================================================================
## FIGURE 16: Cascade Propagation — Temporal Sequence of Dysfunction
##=============================================================================

cat("Generating Fig16: Cascade Propagation...\n")

## Track when each variable deviates >10% from its healthy steady-state
## and when it hits 50% deviation (severe)

track_vars <- list(
  "CL integrity"      = list(col = "CL_ratio",    healthy = as.numeric(ss_h$CL_ratio),    dir = "down"),
  "CI activity"       = list(col = "CI_activity",  healthy = as.numeric(ss_h$CI_activity),  dir = "down"),
  "SC integrity"      = list(col = "SC_integrity", healthy = as.numeric(ss_h$SC_integrity), dir = "down"),
  "ATP levels"        = list(col = "ATP",           healthy = as.numeric(ss_h$ATP),          dir = "down"),
  "Membrane potential" = list(col = "delta_psi_m",  healthy = as.numeric(ss_h$delta_psi_m),  dir = "down"),
  "mROS elevation"    = list(col = "mROS",          healthy = as.numeric(ss_h$mROS),         dir = "up"),
  "aSyn aggregation"  = list(col = "aSyn_olig",     healthy = as.numeric(ss_h$aSyn_olig),    dir = "up"),
  "Microglia activation" = list(col = "MG_active",  healthy = as.numeric(ss_h$MG_active),    dir = "up"),
  "mPTP opening"      = list(col = "mPTP_open",     healthy = as.numeric(ss_h$mPTP_open),    dir = "up"),
  "DA neuron loss"    = list(col = "DA_neuron",      healthy = as.numeric(ss_h$DA_neuron),    dir = "down"),
  "Motor impairment"  = list(col = "Motor_score",    healthy = as.numeric(ss_h$Motor_score),  dir = "up")
)

cascade_df <- data.frame()
for (vname in names(track_vars)) {
  info <- track_vars[[vname]]
  vals_dis <- res_d[[info$col]]
  vals_drg <- res_rx[[info$col]]
  h_val <- info$healthy

  for (cond_info in list(list(v = vals_dis, cond = "PD (untreated)"),
                         list(v = vals_drg, cond = "PD + SBT-272"))) {
    vals <- cond_info$v

    if (info$dir == "down") {
      deviation <- (h_val - vals) / max(h_val, 1e-8)
    } else {
      deviation <- (vals - h_val) / max(1 - h_val, 1e-8)
    }

    t10 <- res_d$time[which(deviation > 0.10)[1]] / 24
    t50 <- res_d$time[which(deviation > 0.50)[1]] / 24
    if (is.na(t10)) t10 <- NA
    if (is.na(t50)) t50 <- NA

    cascade_df <- rbind(cascade_df, data.frame(
      variable = vname, condition = cond_info$cond,
      t_onset = t10, t_severe = t50,
      stringsAsFactors = FALSE
    ))
  }
}

## Order by onset time (disease trajectory)
dis_order <- cascade_df %>% filter(condition == "PD (untreated)") %>%
  arrange(t_onset) %>% pull(variable)
cascade_df$variable <- factor(cascade_df$variable, levels = rev(dis_order))

col_cascade <- c("PD (untreated)" = "#E53935", "PD + SBT-272" = "#1565C0")

p16 <- ggplot(cascade_df %>% filter(!is.na(t_onset)),
              aes(y = variable, colour = condition)) +
  geom_segment(aes(x = t_onset, xend = ifelse(is.na(t_severe), max(cascade_df$t_onset, na.rm=T)*1.1, t_severe),
                   yend = variable),
               linewidth = 3, alpha = 0.5,
               position = position_dodge(width = 0.6)) +
  geom_point(aes(x = t_onset), size = 3.5, shape = 16,
             position = position_dodge(width = 0.6)) +
  geom_point(data = cascade_df %>% filter(!is.na(t_severe)),
             aes(x = t_severe), size = 3.5, shape = 17,
             position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = col_cascade) +
  scale_x_continuous(breaks = seq(0, 50, 5)) +
  labs(title = "Cascade Propagation Timeline",
       subtitle = expression("Sequential dysfunction onset ("*circle*" >10% deviation) to severe ("*triangle*" >50% deviation) | SBT-272 delays propagation"),
       x = "Days after disease onset", y = NULL, colour = NULL) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.major.x = element_line(colour = "grey90"),
        axis.text.y = element_text(size = 9))

ggsave("Fig16_cascade_propagation.png", p16, width = 10, height = 6, dpi = 300, bg = "white")
cat("  Saved: Fig16_cascade_propagation.png\n")


##=============================================================================
## SUMMARY
##=============================================================================

cat("\n================================================================\n")
cat("  NOVEL FIGURES COMPLETE\n")
cat("================================================================\n\n")
cat("Generated (300 DPI, journal-ready):\n")
cat("  Fig13_phase_portrait.png        - Vicious cycle attractor dynamics\n")
cat("  Fig14_death_decomposition.png   - Multi-pathway death over disease course\n")
cat("  Fig15_disease_fingerprint.png   - Radar chart: multi-dimensional state\n")
cat("  Fig16_cascade_propagation.png   - Sequential dysfunction timeline\n\n")
