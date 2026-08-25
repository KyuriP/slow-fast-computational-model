# ============================================================
# 05_figures_make_all.R
# ============================================================
# Purpose
# -------
# Create the main manuscript figures from the final saved results:
#   1) Faceted representative trajectories (A/B × shocks)
#   2) Faceted bifurcation overlays (A/B × shocks)
#   3) Faceted ensemble trajectories
#   4) Heatmap figure from saved no-shock regime scan
#   5) Hysteresis figure from saved ramp simulation
#
# Required inputs
# ---------------
#   res/main/res_A_clean.rds
#   res/main/res_AS_clean.rds
#   res/main/res_B0_clean.rds
#   res/main/res_BS_clean.rds
#   res/heatmaps/heat_zoom_noshock_hyst2.rds
#   res/hysteresis/sim_B2.rds
#
# Dependencies
# ------------
#   R/utils/utils_fastlayer.R
#   R/utils/utils_slowfast.R
#   R/utils/utils_diagnostics.R
#   R/utils/utils_plotting.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
  library(scales)
  library(grid)
  library(cowplot)
})

# ------------------------------------------------------------
# Source shared utilities
# ------------------------------------------------------------
source("R/utils/utils_fastlayer.R")
source("R/utils/utils_slowfast.R")
source("R/utils/utils_diagnostics.R")
source("R/utils/utils_plotting.R")

# ------------------------------------------------------------
# Load final scenario results
# ------------------------------------------------------------
res_A  <- readRDS("res/main/res_A_clean.rds")
res_AS <- readRDS("res/main/res_AS_clean.rds")
res_B0 <- readRDS("res/main/res_B0_clean.rds")
res_BS <- readRDS("res/main/res_BS_clean.rds")

# ------------------------------------------------------------
# Shared settings
# ------------------------------------------------------------
burn_final <- 500

pal <- c(
  "Low baseline context"  = "#7ABD7D",
  "High baseline context" = "#FF5733"
)

scenario_labels <- c(
  A = "M (monostable)",
  B = "B (bistable)"
)

shock_labels <- c(
  No  = "No shocks",
  Yes = "Shocks"
)

# ------------------------------------------------------------
# Representative replicate extraction
# ------------------------------------------------------------
pick_representative_reps <- function(metrics_df) {
  stopifnot(all(c("group", "rep", "mean_m") %in% names(metrics_df)))
  
  if ("scenario" %in% names(metrics_df)) {
    metrics_df %>%
      group_by(scenario, group) %>%
      mutate(target = mean(mean_m, na.rm = TRUE)) %>%
      slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(scenario, group, rep)
  } else {
    metrics_df %>%
      group_by(group) %>%
      mutate(target = mean(mean_m, na.rm = TRUE)) %>%
      slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(group, rep)
  }
}

get_rep_sim <- function(res_obj) {
  rep_pick <- pick_representative_reps(res_obj$metrics) %>%
    select(group, rep)
  
  res_obj$sim_all %>%
    inner_join(rep_pick, by = c("group", "rep"))
}

sim_A_rep  <- get_rep_sim(res_A)
sim_AS_rep <- get_rep_sim(res_AS)
sim_B0_rep <- get_rep_sim(res_B0)
sim_BS_rep <- get_rep_sim(res_BS)

# ============================================================
# 1) FACETED REPRESENTATIVE TRAJECTORIES
# ============================================================

scenario_levels_traj <- c("M", "B")
scenario_labels_traj <- c("M (monostable)", "B (bistable)")
shock_levels_traj    <- c("No", "Yes")
shock_labels_traj    <- c("No shocks", "Shocks")

sim_rep_all <- bind_rows(
  sim_A_rep  %>% mutate(scenario = "M", shocks = "No"),
  sim_AS_rep %>% mutate(scenario = "M", shocks = "Yes"),
  sim_B0_rep %>% mutate(scenario = "B", shocks = "No"),
  sim_BS_rep %>% mutate(scenario = "B", shocks = "Yes")
) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_levels_traj, labels = scenario_labels_traj),
    shocks   = factor(shocks,   levels = shock_levels_traj,    labels = shock_labels_traj),
    group    = factor(group, levels = names(pal))
  )

long_ts <- sim_rep_all %>%
  select(t, group, scenario, shocks, P, m) %>%
  pivot_longer(c(P, m), names_to = "series", values_to = "value") %>%
  mutate(
    series = factor(series, levels = c("P", "m"), labels = c("P[t]", "m[t]"))
  )

shock_lines <- NULL
has_shock_cols <- any(c("shock_any", "shock", "shock_exo", "shock_endo") %in% names(sim_rep_all))

if (has_shock_cols) {
  shock_lines <- sim_rep_all %>%
    group_by(scenario, shocks) %>%
    group_modify(~ shock_times(.x, burn_in = burn_final, by_group = TRUE)) %>%
    ungroup() %>%
    tidyr::crossing(series = levels(long_ts$series)) %>%
    mutate(
      series = factor(series, levels = levels(long_ts$series)),
      group  = factor(group, levels = names(pal))
    )
}

burn_label_df <- tibble(
  scenario = factor(scenario_labels_traj[1], levels = levels(long_ts$scenario)),
  shocks   = factor(shock_labels_traj[1],    levels = levels(long_ts$shocks)),
  series   = factor("P[t]", levels = levels(long_ts$series)),
  x        = burn_final * 0.35,
  y        = Inf,
  lab      = "Burn-in"
)

fig_traj_faceted <- ggplot(long_ts, aes(t, value, colour = group)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = burn_final,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.15
  ) +
  geom_text(
    data = burn_label_df,
    aes(x = x, y = y, label = lab),
    colour = "grey35",
    fontface = 2,
    size = 3.5,
    vjust = 1.2,
    hjust = 0
  ) +
  geom_line(alpha = 0.90, linewidth = 0.6, key_glyph = "path") +
  {
    if (!is.null(shock_lines) && nrow(shock_lines) > 0) {
      geom_vline(
        data = shock_lines,
        aes(xintercept = t, colour = group),
        linetype = "dashed",
        linewidth = 0.35,
        alpha = 0.9,
        show.legend = FALSE
      )
    }
  } +
  ggh4x::facet_nested(
    rows = vars(scenario, series),
    cols = vars(shocks),
    scales = "free_y",
    labeller = labeller(series = label_parsed)
  ) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(x = "Time step", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = 0.4),
    strip.text.y = element_text(size = 11, face = "bold"),
    strip.text.x = element_text(size = 11, face = "bold")
  )

fig_traj_faceted

# ggsave("img/example_ts.png", fig_traj_faceted, width = 10, height = 9)

# ============================================================
# 2) FACETED BIFURCATION OVERLAYS
# ============================================================

strip_titles <- function(p) p + labs(title = NULL, subtitle = NULL)
kill_axis_titles <- function(p) p + labs(x = NULL, y = NULL)

legend_dot_fix <- function(p, alpha = 0.9, size = 2.5) {
  p + guides(colour = guide_legend(override.aes = list(alpha = alpha, size = size)))
}

plot_margin_fix <- theme(plot.margin = margin(8, 8, 8, 8))

axis_base  <- theme(axis.title.x = element_blank(), axis.title.y = element_blank())
axis_top   <- theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank())
axis_right <- theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank())

strip_cell <- function(label = "", rotate = 0,
                       fill = "grey95", line_col = "grey30",
                       text_size = 11, fontface = "bold",
                       border = TRUE) {
  grid::grobTree(
    grid::rectGrob(gp = gpar(
      fill = if (border) fill else NA,
      col  = if (border) line_col else NA,
      lwd  = if (border) 0.8 else 0
    )),
    grid::textGrob(label, rot = rotate,
             gp = grid::gpar(fontsize = text_size, fontface = fontface))
  )
}

make_overlay_panel <- function(bf, bif_info, sim_df, ribbon_alpha) {
  p <- plot_branch_overlay(
    bf = bf,
    bif_info = bif_info,
    sim_df = sim_df,
    title = NULL,
    burn_in = burn_final,
    cols = pal,
    ribbon_alpha = ribbon_alpha,
    max_points_per_group = 8000,
    legend_pos = "none",
    base_size = 14
  )
  
  p %>%
    strip_titles() %>%
    kill_axis_titles() %>%
    legend_dot_fix() +
    plot_margin_fix
}

p_ov_A  <- make_overlay_panel(res_A$bf,  res_A$bif_info,  sim_A_rep,  ribbon_alpha = 0.00)
p_ov_AS <- make_overlay_panel(res_AS$bf, res_AS$bif_info, sim_AS_rep, ribbon_alpha = 0.00)
p_ov_B0 <- make_overlay_panel(res_B0$bf, res_B0$bif_info, sim_B0_rep, ribbon_alpha = 0.12)
p_ov_BS <- make_overlay_panel(res_BS$bf, res_BS$bif_info, sim_BS_rep, ribbon_alpha = 0.12)

pA  <- p_ov_A  + axis_base + axis_top
pAS <- p_ov_AS + axis_base + axis_top + axis_right
pB0 <- p_ov_B0 + axis_base
pBS <- p_ov_BS + axis_base + axis_right

xlab_grob <- wrap_elements(
  full = textGrob("Context load  P", gp = gpar(fontsize = 11))
)

row_M  <- wrap_elements(full = strip_cell("M (monostable)", rotate = 270))
row_B  <- wrap_elements(full = strip_cell("B (bistable)",   rotate = 270))

col_no <- wrap_elements(full = strip_cell("No shocks"))
col_sh <- wrap_elements(full = strip_cell("Shocks"))

ga <- guide_area()

strip_w  <- 0.05
header_h <- 0.075
xlab_h   <- 0.07
leg_h    <- 0.16

facet_overlay_right_strips <-
  wrap_plots(
    a = col_no, b = col_sh, c = plot_spacer(),
    d = pA,     e = pAS,    f = row_M,
    g = pB0,    h = pBS,    i = row_B,
    j = xlab_grob, k = plot_spacer(),
    l = ga, m = plot_spacer(),
    design = "
      abc
      def
      ghi
      jjk
      llm
    ",
    widths  = c(1, 1, strip_w),
    heights = c(header_h, 1, 1, xlab_h, leg_h)
  ) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Bifurcation / equilibrium mapping with simulated trajectories",
    theme = theme(plot.title = element_text(face = "bold"))
  ) &
  theme(legend.position = "bottom")

y_lab <- grid::textGrob(
  "Mean symptom activation  m",
  rot = 90,
  gp = grid::gpar(fontsize = 12)
)

fig_overlay_faceted <- cowplot::ggdraw() +
  cowplot::draw_grob(y_lab, x = 0.03, y = 0.5, width = 0.04, height = 0.1) +
  cowplot::draw_plot(facet_overlay_right_strips, x = 0.06, y = 0, width = 0.94, height = 1)

fig_overlay_faceted

# ggsave("img/example_overlay_M.pdf", fig_overlay_faceted, width = 10, height = 7)

# ============================================================
# 3) FACETED ENSEMBLE FIGURE
# ============================================================

meta_ens <- tibble::tibble(
  key      = c("M0", "MS", "B0", "BS"),
  res      = list(res_A, res_AS, res_B0, res_BS),
  scenario = c("M (monostable)", "M (monostable)", "B (bistable)", "B (bistable)"),
  shocks   = c("No shocks", "Shocks", "No shocks", "Shocks")
) %>%
  mutate(
    scenario = factor(scenario, levels = c("M (monostable)", "B (bistable)")),
    shocks   = factor(shocks,   levels = c("No shocks", "Shocks"))
  )

extract_sim_ensemble <- function(res_obj, scenario, shocks) {
  res_obj$sim_all %>%
    select(any_of(c("t", "group", "rep", "P", "m"))) %>%
    mutate(
      scenario = scenario,
      shocks   = shocks,
      group    = factor(group, levels = names(pal))
    )
}

sim_ens <- pmap_dfr(meta_ens, \(key, res, scenario, shocks) {
  extract_sim_ensemble(res, scenario, shocks)
}) %>%
  filter(t > burn_final)

long_ens <- sim_ens %>%
  pivot_longer(c(P, m), names_to = "var_code", values_to = "value") %>%
  mutate(
    var_code = factor(var_code, levels = c("P", "m"), labels = c("P[t]", "m[t]"))
  )

summ_ens <- long_ens %>%
  group_by(scenario, shocks, group, var_code, t) %>%
  summarise(
    mu  = mean(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

ylims <- long_ens %>%
  group_by(var_code) %>%
  summarise(
    lo = quantile(value, 0.01, na.rm = TRUE),
    hi = quantile(value, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

P_lo <- -1
P_hi <- 3
m_lo <- ylims$lo[ylims$var_code == "m[t]"]
m_hi <- ylims$hi[ylims$var_code == "m[t]"]

p_ens <- ggplot(long_ens, aes(t, value, colour = group)) +
  geom_ribbon(
    data = summ_ens,
    aes(x = t, ymin = q10, ymax = q90, fill = group),
    inherit.aes = FALSE,
    alpha = 0.12
  ) +
  geom_line(
    aes(group = interaction(rep, group)),
    alpha = 0.06,
    linewidth = 0.25
  ) +
  geom_line(
    data = summ_ens,
    aes(y = mu, group = group),
    linewidth = 0.9
  ) +
  ggh4x::facet_nested(
    rows = vars(scenario, var_code),
    cols = vars(shocks),
    scales = "free_y",
    labeller = labeller(var_code = label_parsed)
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      "P[t]" = scale_y_continuous(
        limits = c(P_lo, P_hi),
        oob = scales::squish,
        minor_breaks = NULL
      ),
      "m[t]" = scale_y_continuous(
        limits = c(m_lo, m_hi),
        oob = scales::squish,
        minor_breaks = NULL
      )
    )
  ) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_fill_manual(values = pal, name = NULL) +
  labs(x = "Time step", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = 0.4),
    strip.text.y = element_text(size = 11, face = "bold"),
    strip.text.x = element_text(size = 11, face = "bold"),
    panel.spacing.y = unit(0.6, "lines"),
    panel.spacing.x = unit(0.6, "lines")
  ) +
  guides(
    colour = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1)),
    fill   = "none"
  )

p_ens

# ggsave("img/ensemble_ts2_M.pdf", p_ens, width = 10, height = 9)

# ============================================================
# 4) HEATMAP FIGURE 
# ============================================================

heat_zoom <- readRDS("res/heatmaps/heat_zoom_noshock_hyst2.rds") |>
  mutate(
    betaJ  = as.numeric(betaJ),
    P_base = as.numeric(P_base)
  )

base_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    plot.title = element_text(face = "bold", hjust = 0, size = 14),
    plot.subtitle = element_text(hjust = 0, size = 10),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    plot.margin = margin(6, 6, 6, 6)
  )

cb_pr <- guide_colorbar(
  barwidth = unit(9, "cm"),
  barheight = unit(0.35, "cm"),
  title.position = "left"
)

cb_sw <- guide_colorbar(
  barwidth = unit(9, "cm"),
  barheight = unit(0.35, "cm"),
  title.position = "left"
)

p_frac <- ggplot(heat_zoom, aes(betaJ, P_base)) +
  geom_tile(aes(fill = Pr_high)) +
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0, 1),
    name = expression(f[high]),
    guide = cb_pr
  ) +
  labs(
    title = "High-state occupancy",
    subtitle = "Fraction of time labeled high",
    x = expression(beta * J),
    y = expression(P[base])
  ) +
  base_theme +
  theme(
    legend.title = element_text(margin = margin(r = 12))
  )

p_sw <- ggplot(heat_zoom, aes(betaJ, P_base)) +
  geom_tile(aes(fill = Switches)) +
  scale_fill_viridis_c(
    option = "viridis",
    trans = "sqrt",
    name = expression(N[switch]),
    guide = cb_sw
  ) +
  labs(
    title = "Switching frequency",
    subtitle = expression("Mean number of switches (" * N[switch] * ")"),
    x = expression(beta * J),
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(margin = margin(r = 12))
  )

fig_heat2 <- (p_frac | p_sw) +
  plot_layout(
    widths = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.margin = margin(6, 6, 6, 6),
    panel.spacing = unit(2, "mm")
  )

fig_heat2

# ggsave("img/Figure_heatmaps.png", fig_heat2, width = 10, height = 6)
# ggsave("img/Figure_heatmaps_2panel_plasma_clean.pdf", fig_heat2, width = 10, height = 6)


# ============================================================
# 5) HYSTERESIS FIGURE
# ============================================================

sim_B2 <- readRDS("res/hysteresis/sim_B2.rds")

dfB <- sim_B2 %>%
  mutate(
    t = as.numeric(t),
    P = as.numeric(P),
    m = as.numeric(m),
    P_base = as.numeric(P_base)
  )

Tup <- dfB$t[which.max(dfB$P_base)]

dfB <- dfB %>%
  mutate(
    phase = if_else(t <= Tup, "Ramp up", "Ramp down"),
    phase = factor(phase, levels = c("Ramp up", "Ramp down"))
  )

detect_jumps <- function(x, t, thr = 0.40, min_gap = 30) {
  idx_raw <- which(abs(diff(x)) > thr) + 1L
  if (length(idx_raw) == 0) return(integer(0))
  
  keep <- c(TRUE, diff(t[idx_raw]) >= min_gap)
  idx_raw[keep]
}

jump_idx <- detect_jumps(dfB$m, dfB$t, thr = 0.40, min_gap = 30)
jump_t   <- dfB$t[jump_idx]

jump_pts <- dfB %>%
  slice(jump_idx) %>%
  select(t, P, m, phase)

col_up    <- "#0072B2"
col_down  <- "#D55E00"
col_phase <- c("Ramp up" = col_up, "Ramp down" = col_down)

col_base <- "grey40"
col_P    <- "grey10"
col_jump <- "grey45"

phase_rect <- tibble::tibble(
  xmin = c(min(dfB$t, na.rm = TRUE), Tup),
  xmax = c(Tup, max(dfB$t, na.rm = TRUE)),
  ymin = -Inf,
  ymax = Inf,
  phase = factor(c("Ramp up", "Ramp down"), levels = c("Ramp up", "Ramp down"))
)

base_theme_hyst <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 11)
  )

pP_hyst <- ggplot(dfB, aes(x = t)) +
  geom_rect(
    data = phase_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phase),
    inherit.aes = FALSE,
    alpha = 0.06,
    colour = NA
  ) +
  geom_line(aes(y = P, linetype = "realized"), colour = col_P, linewidth = 0.9) +
  geom_line(aes(y = P_base, linetype = "base"), colour = col_base, linewidth = 0.9) +
  geom_vline(xintercept = Tup, linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  annotate(
    "text",
    x = Tup,
    y = max(c(dfB$P, dfB$P_base), na.rm = TRUE) + 0.04,
    label = "Ramp reversal",
    colour = "grey35",
    size = 3.1,
    vjust = 0
  ) +
  scale_fill_manual(values = col_phase, guide = "none") +
  scale_linetype_manual(
    name = NULL,
    values = c(realized = "solid", base = "dashed"),
    breaks = c("realized", "base"),
    labels = c(
      realized = expression("Realized " * P[t]),
      base     = expression("Imposed " * P[base](t))
    )
  ) +
  guides(
    linetype = guide_legend(
      keywidth = unit(1.0, "cm"),
      keyheight = unit(0.4, "cm"),
      override.aes = list(
        colour = c(col_P, col_base),
        linewidth = c(0.9, 0.9)
      ),
      order = 1
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Context ramp protocol and realized context",
    subtitle = expression("Dashed line: imposed " * P[base](t) * "; solid line: realized " * P[t] * "."),
    x = NULL,
    y = expression(P[t])
  ) +
  base_theme_hyst +
  theme(plot.margin = margin(12, 8, 6, 6))

pm_hyst <- ggplot(dfB, aes(x = t, y = m)) +
  geom_rect(
    data = phase_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phase),
    inherit.aes = FALSE,
    alpha = 0.06,
    colour = NA
  ) +
  geom_line(
    data = dfB %>% filter(phase == "Ramp up"),
    aes(colour = phase),
    linewidth = 0.8
  ) +
  geom_line(
    data = dfB %>% filter(phase == "Ramp down"),
    aes(colour = phase),
    linewidth = 0.8
  ) +
  geom_vline(xintercept = jump_t, linetype = "dashed", colour = col_jump, linewidth = 0.45, alpha = 0.8) +
  geom_vline(xintercept = Tup, linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  scale_fill_manual(values = col_phase, guide = "none") +
  scale_colour_manual(
    values = col_phase,
    breaks = c("Ramp up", "Ramp down"),
    name = NULL
  ) +
  labs(
    title = "Symptom trajectory under slow ramp",
    subtitle = "Vertical dashed lines mark abrupt jumps between branches.",
    x = "Time step",
    y = expression(m[t]),
    colour = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  coord_cartesian(clip = "off") +
  base_theme_hyst +
  theme(plot.margin = margin(6, 8, 6, 6))

ploop_hyst <- ggplot(dfB, aes(x = P, y = m, colour = phase)) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  geom_point(data = jump_pts, aes(P, m), size = 2.2, stroke = 0.5, fill = "white", shape = 21) +
  scale_colour_manual(
    values = c("Ramp up" = col_up, "Ramp down" = col_down),
    breaks = c("Ramp up", "Ramp down"),
    guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Hysteresis loop in state space",
    subtitle = "Same context level can map to different symptom levels depending on history.",
    x = expression(P[t]),
    y = expression(m[t])
  ) +
  base_theme_hyst

left_col_hyst <- pP_hyst / pm_hyst + plot_layout(heights = c(1, 1))

hyst_fig <- (left_col_hyst | ploop_hyst) +
  plot_layout(widths = c(1.1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12)
  )

hyst_fig

# ggsave("img/Figure_hysteresis.pdf", hyst_fig, width = 13, height = 7)
