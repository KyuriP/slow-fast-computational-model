
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggh4x)
  library(scales)   # for squish(), pretty breaks, etc.
  library(grid)     # for grobs in the overlay-strip layout
  library(cowplot)
})


# ============================================================
# Faceted representative trajectories figure (A/B × shocks)
# - Columns: No shocks | Shocks
# - Rows: Scenario A then B, each with sub-rows P[t] and m[t]
# - Burn-in: shaded region + single "Burn-in" label (top-left only)
# - Shocks: vertical dashed lines colored by group in panels
#           + neutral grey dashed "Shock" legend key
# - Y scales: fixed per series
#     P[t] limits [-1.5, 1.5], breaks {1.5, 0, -0.5}
#     m[t] limits [0, 1.0],   breaks {1.0, 0.5, 0.0}
# - Clean legends: Shock (linetype) + Groups (colour) with line-only keys
# ============================================================


# ---- facet labels ----
scenario_levels <- c("A","B")
scenario_labels <- c("A (monostable)","B (bistable)")
shock_levels    <- c("No","Yes")
shock_labels    <- c("No shocks","Shocks")

# ---- combine representative sims ----
sim_all <- bind_rows(
  sim_A_rep  %>% mutate(scenario="A", shocks="No"),
  sim_AS_rep %>% mutate(scenario="A", shocks="Yes"),
  sim_B0_rep %>% mutate(scenario="B", shocks="No"),
  sim_BS_rep %>% mutate(scenario="B", shocks="Yes")
) %>%
  mutate(
    scenario = factor(scenario, levels=scenario_levels, labels=scenario_labels),
    shocks   = factor(shocks,   levels=shock_levels,    labels=shock_labels),
    group    = factor(group, levels = names(pal))
  )

# ---- long format: P and m ----
# Use parseable math strings for strip labels: P[t], m[t]
long_ts <- sim_all %>%
  select(t, group, scenario, shocks, P, m) %>%
  pivot_longer(c(P, m), names_to="series", values_to="value") %>%
  mutate(
    series = factor(series,
                    levels = c("P","m"),
                    labels = c("P[t]", "m[t]"))
  )

# ---- fold lines: only for P[t] ----
folds_all <- bind_rows(
  folds_to_df(res_A$bif_info)   %>% mutate(scenario=scenario_labels[1], shocks=shock_labels[1]),
  folds_to_df(res_A_S$bif_info) %>% mutate(scenario=scenario_labels[1], shocks=shock_labels[2]),
  folds_to_df(res_B0$bif_info)  %>% mutate(scenario=scenario_labels[2], shocks=shock_labels[1]),
  folds_to_df(res_B_S$bif_info) %>% mutate(scenario=scenario_labels[2], shocks=shock_labels[2])
) %>%
  mutate(
    scenario = factor(scenario, levels = scenario_labels),
    shocks   = factor(shocks,   levels = shock_labels),
    series   = factor("P[t]", levels = levels(long_ts$series))
  )

# ---- shock lines: compute per (scenario, shocks), colored by group ----
has_shock_cols <- any(c("shock_any","shock","shock_exo","shock_endo") %in% names(sim_all))
shock_lines <- NULL
if (has_shock_cols) {
  shock_lines <- sim_all %>%
    group_by(scenario, shocks) %>%
    group_modify(~ shock_times(.x, burn_in = burnin, by_group = TRUE)) %>%  # must return columns incl. t and group
    ungroup() %>%
    tidyr::crossing(series = levels(long_ts$series)) %>%                    # show on both P[t] and m[t]
    mutate(
      series = factor(series, levels = levels(long_ts$series)),
      group  = factor(group, levels = names(pal))
    )
}

# ---- Burn-in label ONLY once (top-left facet: A, P[t], No shocks) ----
burn_label_df <- tibble(
  scenario = factor(scenario_labels[1], levels = levels(long_ts$scenario)),
  shocks   = factor(shock_labels[1],    levels = levels(long_ts$shocks)),
  series   = factor("P[t]",             levels = levels(long_ts$series)),
  x        = burnin * 0.35,
  y        = Inf,
  lab      = "Burn-in"
)

# ---- Plot ----
p <- ggplot(long_ts, aes(t, value, colour = group)) +
  
  # Burn-in shading across facets (no legend)
  annotate(
    "rect",
    xmin = -Inf, xmax = burnin,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.15
  ) +
  
  # Burn-in text only once
  geom_text(
    data = burn_label_df,
    aes(x = x, y = y, label = lab),
    colour = "grey35",
    fontface =2,
    size = 3.5,
    vjust = 1.2,
    hjust = 0
  ) +
  
  # Trajectories (force legend glyph to be a simple line)
  geom_line(alpha = 0.90, linewidth = 0.6, key_glyph = "path") +
  
  # # Fold lines (P[t] only)
  # geom_hline(
  #   data = folds_all,
  #   aes(yintercept = P_fold),
  #   colour = "grey40",
  #   linetype = "dotted",
  #   linewidth = 0.6,
  #   show.legend = FALSE
  # ) +
  
  # Shocks: coloured by group in panels (no legend)
  { if (!is.null(shock_lines) && nrow(shock_lines) > 0)
    geom_vline(
      data = shock_lines,
      aes(xintercept = t, colour = group),
      linetype = "dashed",
      linewidth = 0.35,
      alpha = 0.9,
      show.legend = FALSE
    )
  } +
  
  # # Dummy layer: neutral grey dashed "Shock" legend key only
  # geom_vline(
  #   data = tibble(x = Inf),
  #   aes(xintercept = x, linetype = "Shock"),
  #   colour = "grey35",
  #   linewidth = 0.35,
  #   show.legend = TRUE
  # ) +
  # 
  # Facets (nested row strips with big A/B label, series sublabel)
  ggh4x::facet_nested(
    rows = vars(scenario, series),
    cols = vars(shocks),
    scales = "free_y",
    labeller = labeller(series = label_parsed)
  ) +
  
  # # Fixed y limits/breaks per series (reduces clutter)
  # ggh4x::facetted_pos_scales(
  #   y = list(
  #     "P[t]" = scale_y_continuous(
  #       limits = c(-1.5, 1.5),
  #       breaks = c(1.5, 0, -0.5),
  #       minor_breaks = NULL,
  #       labels = function(x) format(x, trim = TRUE)  
  #     ),
  #     "m[t]" = scale_y_continuous(
  #       limits = c(0, 1.0),
  #       breaks = c(1.0, 0.5, 0.0),
  #       minor_breaks = NULL,
  #       labels = function(x) format(x, trim = TRUE)
  #     )
  #   )
  # ) +
  
  # Scales + labels
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(name = NULL, values = c("Shock" = "dashed"), breaks = "Shock") +
  labs(x = "Time step", y = NULL) +
  
  # Theme + boxed strip labels
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", colour = "grey30", linewidth = 0.4),
    #ggh4x.strip.nest.text.y = element_text(face = "bold", size = 13),
    strip.text.y = element_text(size = 11, face = "bold"),
    strip.text.x = element_text(size = 11, face = "bold")
  ) 

p


# ggsave("img/example_ts.png",p, width=10,height=9)











# ============================================================
# FINAL (FIXED): Overlay facet-style 2×2 (A/B × shocks)
# - Facet-like strips with fixed size + aligned
# - Legend in its own row (prevents B-strip "overhang")
# - Outer axes only (ticks only on bottom row / left col)
# - No per-panel titles/subtitles
# - Global title/subtitle + global axis titles as insets
# ============================================================

# --------------------------
# Helpers
# --------------------------

strip_titles <- function(p) p + labs(title = NULL, subtitle = NULL)
kill_axis_titles <- function(p) p + labs(x = NULL, y = NULL)

legend_dot_fix <- function(p, alpha = 0.9, size = 2.5) {
  p + guides(colour = guide_legend(override.aes = list(alpha = alpha, size = size)))
}

plot_margin_fix <- theme(plot.margin = margin(8, 8, 8, 8))

# outer-axes-only rules
axis_base  <- theme(axis.title.x = element_blank(), axis.title.y = element_blank())
axis_top   <- theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank())
axis_right <- theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank())

# strip cell (corner is borderless)
strip_cell <- function(label = "", rotate = 0,
                       fill = "grey95", line_col = "grey30",
                       text_size = 11, fontface = "bold",
                       border = TRUE) {
  grobTree(
    rectGrob(gp = gpar(
      fill = if (border) fill else NA,
      col  = if (border) line_col else NA,
      lwd  = if (border) 0.8 else 0
    )),
    textGrob(label, rot = rotate,
             gp = gpar(fontsize = text_size, fontface = fontface))
  )
}

# --------------------------
# 1) Build four overlay panels
# --------------------------

make_overlay_panel <- function(bf, bif_info, sim_df, ribbon_alpha) {
  p <- plot_branch_overlay(
    bf = bf, bif_info = bif_info, sim_df = sim_df,
    title = NULL,
    burn_in = burnin,
    cols = pal,
    ribbon_alpha = ribbon_alpha,
    max_points_per_group = 8000,
    legend_pos = "none",
    base_size = 14
  )
  
  # enforce: no titles/subtitles, no axis titles, legend readable, consistent margins
  p <- strip_titles(p)
  p <- kill_axis_titles(p)
  p <- legend_dot_fix(p)
  p <- p + plot_margin_fix
  p
}

p_ov_A  <- make_overlay_panel(res_A$bf,   res_A$bif_info,   sim_A_rep,  ribbon_alpha = 0.00)
p_ov_AS <- make_overlay_panel(res_A_S$bf, res_A_S$bif_info, sim_AS_rep, ribbon_alpha = 0.00)
p_ov_B0 <- make_overlay_panel(res_B0$bf,  res_B0$bif_info,  sim_B0_rep, ribbon_alpha = 0.12)
p_ov_BS <- make_overlay_panel(res_B_S$bf, res_B_S$bif_info, sim_BS_rep, ribbon_alpha = 0.12)

# outer-axes-only versions (facet-like)
pA  <- p_ov_A  + axis_base + axis_top
pAS <- p_ov_AS + axis_base + axis_top + axis_right
pB0 <- p_ov_B0 + axis_base
pBS <- p_ov_BS + axis_base + axis_right

# --- axis label grobs ---
# y-label ONCE (we'll place it in the A row and make it visually centered)
ylab_once <- wrap_elements(
  full = textGrob(
    "Mean symptom activation  m",
    rot = 90,
    gp = gpar(fontsize = 11)
  )
)

xlab_grob <- wrap_elements(
  full = textGrob("Context load  P", gp = gpar(fontsize = 11))
)

# --- right-side strips: flip reading direction bottom-to-top (rotate 270) ---
row_A  <- wrap_elements(full = strip_cell("A (monostable)", rotate = 270))
row_B  <- wrap_elements(full = strip_cell("B (bistable)",   rotate = 270))

# --- column headers ---
col_no <- wrap_elements(full = strip_cell("No shocks"))
col_sh <- wrap_elements(full = strip_cell("Shocks"))

# --- spacers + guide area ---
ga <- guide_area()

# ---- remove y-label column: 3 columns only (No-shocks | Shocks | right strip) ----
# widths now: plot, plot, strip
strip_w  <- 0.05
header_h <- 0.075
xlab_h   <- 0.07
leg_h    <- 0.16

facet_overlay_right_strips <-
  wrap_plots(
    # Row 1: column headers
    a = col_no, b = col_sh, c = plot_spacer(),
    
    # Row 2: A row
    d = pA,     e = pAS,    f = row_A,
    
    # Row 3: B row
    g = pB0,    h = pBS,    i = row_B,
    
    # Row 4: x-label row (spans plot columns)
    j = xlab_grob, k = plot_spacer(),
    
    # Row 5: legend row
    l = ga,       m = plot_spacer(),
    
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
    # subtitle = expression(
    #   "Grey curves show equilibrium " * m^"*" * "(P) (solid = stable, dashed = unstable); " *
    #     "points show post-burn-in trajectory samples " * (P[t] * "," * m[t]) * "."
    # ),
    theme = theme(plot.title = element_text(face = "bold"))
  ) &
  theme(legend.position = "bottom")

facet_overlay_right_strips




# your existing plot
p_main <- facet_overlay_right_strips

# y label as a grob (outside)
y_lab <- grid::textGrob(
  "Mean symptom activation  m",
  rot = 90,
  gp = grid::gpar(fontsize = 12)
)

# draw label in the left margin area + plot on the right
p_with_ylab <- cowplot::ggdraw() +
  cowplot::draw_grob(y_lab, x = 0.03, y = 0.5, width = 0.04, height = 0.1) +
  cowplot::draw_plot(p_main, x = 0.06, y = 0, width = 0.94, height = 1)

p_with_ylab

# --------------------------
# 4) Save
# --------------------------
# ggsave("img/example_overlay_only.png", p_with_ylab, width = 10, height = 7)

ggsave("img/example_overlay.pdf", p_with_ylab, width = 10, height = 7)





















#| label: results-A-B-ensemble
#| include: true
#| fig-width: 12
#| fig-height: 10

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# ---- 1) Build ensemble long + summary ----

extract_sim <- function(res_obj, scenario_label){
  res_obj$sim_all %>%
    dplyr::select(dplyr::any_of(c("t","group","rep","P","m"))) %>%
    dplyr::mutate(
      scenario = scenario_label,
      group = dplyr::recode(group, !!!group_recode),
      group = factor(group, levels = names(pal))
    )
}

# Explicit scenario order (edit labels here if you want)
scenario_levels <- c(
  "A (βJ<4) no shocks",
  "A (βJ<4) both shocks",
  "B0 (βJ>4) no shocks",
  "B (βJ>4) both shocks"
)

res_list <- list(
  "A (βJ<4) no shocks"   = res_A,
  "A (βJ<4) both shocks" = res_A_S,
  "B0 (βJ>4) no shocks"  = res_B0,
  "B (βJ>4) both shocks" = res_B_S
)

sim_ens <- imap_dfr(res_list, extract_sim) %>%
  filter(t > burnin) %>%
  mutate(scenario = factor(scenario, levels = scenario_levels))

long_ens <- sim_ens %>%
  pivot_longer(c(m, P), names_to = "var_code", values_to = "value") %>%
  mutate(var_code = factor(var_code, levels = c("m","P")))

summ_ens <- long_ens %>%
  group_by(scenario, group, var_code, t) %>%
  summarise(
    mu  = mean(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

# ---- 2) Plot builder (one column: either m(t) or P(t)) ----

build_col <- function(var_keep = c("m","P"), ylim = NULL, title = NULL) {
  var_keep <- match.arg(var_keep, choices = c("m","P"))
  
  dat  <- long_ens %>% filter(var_code == var_keep)
  summ <- summ_ens %>% filter(var_code == var_keep)
  
  ggplot(dat, aes(t, value, color = group)) +
    geom_ribbon(
      data = summ,
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
      data = summ,
      aes(y = mu, group = group),
      linewidth = 0.9
    ) +
    facet_grid(rows = vars(scenario), scales = "fixed") +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    labs(x = NULL, y = NULL, title = title, color = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      # tighter facet spacing
      panel.spacing.y = unit(0.6, "lines")
    ) +
    { if (!is.null(ylim)) coord_cartesian(ylim = ylim) }
}
# Give both plots the same x label so axis_titles="collect" can pick it up
p_m <- p_m + labs(x = "Time step")
p_P <- p_P + labs(x = "Time step")

ens_fig <-
  (p_m | p_P) +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(
    theme = theme(
      plot.margin = margin(6, 6, 6, 6)
    )
  )

# Apply global theme LAST (using &)
ens_fig <- ens_fig & theme(
  legend.position = "bottom",
  axis.title.x = element_text(margin = margin(t = 8))
)

ens_fig









#| label: results-A-B-ensemble-faceted
#| include: true
#| fig-width: 12
#| fig-height: 9


# ----------------------------
# 1) Map your result objects to (scenario, shocks)
# ----------------------------
meta <- tibble::tibble(
  key      = c("A0", "AS", "B0", "BS"),
  res      = list(res_A, res_A_S, res_B0, res_B_S),
  scenario = c("A (monostable)", "A (monostable)", "B (bistable)", "B (bistable)"),
  shocks   = c("No shocks",      "Shocks",        "No shocks",    "Shocks")
) %>%
  mutate(
    scenario = factor(scenario, levels = c("A (monostable)", "B (bistable)")),
    shocks   = factor(shocks,   levels = c("No shocks", "Shocks"))
  )

extract_sim <- function(res_obj, scenario, shocks){
  res_obj$sim_all %>%
    select(any_of(c("t","group","rep","P","m"))) %>%
    mutate(
      scenario = scenario,
      shocks   = shocks,
      group = recode(group, !!!group_recode),
      group = factor(group, levels = names(pal))
    )
}

sim_ens <- pmap_dfr(meta, \(key, res, scenario, shocks) {
  extract_sim(res, scenario, shocks)
}) %>%
  filter(t > burnin)

# ----------------------------
# 2) Long + summary
# ----------------------------
long_ens <- sim_ens %>%
  pivot_longer(c(P, m), names_to = "var_code", values_to = "value") %>%
  mutate(
    var_code = factor(var_code, levels = c("P","m"), labels = c("P[t]","m[t]"))
  )

summ_ens <- long_ens %>%
  group_by(scenario, shocks, group, var_code, t) %>%
  summarise(
    mu  = mean(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    .groups = "drop"
  )



# 1) Robust global limits for comparability (same across A/B and shocks)
ylims <- long_ens %>%
  group_by(var_code) %>%
  summarise(
    lo = quantile(value, 0.01, na.rm = TRUE),
    hi = quantile(value, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

P_lo <- -1 #ylims$lo[ylims$var_code == "P[t]"]
P_hi <- 3 #ylims$hi[ylims$var_code == "P[t]"]
m_lo <- ylims$lo[ylims$var_code == "m[t]"]
m_hi <- ylims$hi[ylims$var_code == "m[t]"]

# ----------------------------
# 3) Faceted ensemble plot (matches your desired layout)
# ----------------------------
p_ens <- ggplot(long_ens, aes(t, value, colour = group)) +
  # 10–90% ribbon
  geom_ribbon(
    data = summ_ens,
    aes(x = t, ymin = q10, ymax = q90, fill = group),
    inherit.aes = FALSE,
    alpha = 0.12
  ) +
  # faint individual reps
  geom_line(
    aes(group = interaction(rep, group)),
    alpha = 0.06,
    linewidth = 0.25
  ) +
  
  # mean line
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
    # make legend dots visible despite low alpha in the plot
    colour = guide_legend(override.aes = list(alpha = 0.9, linewidth = 1)),
    fill   = "none"
  )

p_ens

# ggsave("img/ensemble_ts2.pdf",p_ens, width=10,height=9)






#| label: heatmaps-hyst2-noshock-refined
#| include: true
#| fig-width: 12
#| fig-height: 8.5

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# ---- load saved grid results ----
heat_zoom <- readRDS("res_clean/heat_zoom_noshock_hyst2.rds")

# ---- sanity checks ----
req <- c("betaJ","P_base","Pr_high","Switches","MeanHighDur")
stopifnot(all(req %in% names(heat_zoom)))

# ---- clean + order ----
heat_zoom <- heat_zoom %>%
  mutate(
    betaJ  = as.numeric(betaJ),
    P_base = as.numeric(P_base)
  ) %>%
  arrange(betaJ, P_base)

# ---- common theme ----
heat_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 11),
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 9),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold")
  )

# ---- common axes / styling helper ----
heat_axes <- list(
  scale_x_continuous(
    expand = c(0, 0),
    breaks = scales::pretty_breaks(5)
  ),
  scale_y_continuous(
    expand = c(0, 0),
    breaks = scales::pretty_breaks(5)
  ),
  coord_cartesian(expand = FALSE)
)

# ---- optional visual reference line at betaJ = 4 ----
betaJ_threshold_layer <- list(
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.45, colour = "grey20")
)

# ============================================================
# Panel A: Pr(high)  (diverging, midpoint = 0.5)
# ============================================================
p_frac <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = Pr_high)) +
  geom_tile() +
  geom_contour(
    aes(z = Pr_high),
    breaks = 0.5,
    colour = "black",
    linewidth = 0.45
  ) +
  betaJ_threshold_layer +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426",
    midpoint = 0.5, limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = scales::label_number(accuracy = 0.01),
    name = "Pr(high)"
  ) +
  labs(
    title = "A. Occupancy of high-symptom mode",
    x = expression(beta*J),
    y = expression(P[base])
  ) +
  heat_axes +
  heat_theme +
  guides(fill = guide_colorbar(
    barwidth = unit(8, "cm"),
    barheight = unit(0.35, "cm"),
    title.position = "top"
  ))

# ============================================================
# Panel B: Switches (sqrt-scaled sequential)
# ============================================================
p_sw <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = Switches)) +
  geom_tile() +
  betaJ_threshold_layer +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",   # helps if highly skewed
    name = "Switches\n(sqrt scale)"
  ) +
  labs(
    title = "B. Switching frequency",
    x = expression(beta*J),
    y = NULL
  ) +
  heat_axes +
  heat_theme +
  guides(fill = guide_colorbar(
    barwidth = unit(4.2, "cm"),
    barheight = unit(0.35, "cm"),
    title.position = "top"
  ))

# ============================================================
# Panel C: Mean high duration (sequential)
# ============================================================
p_dur <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = MeanHighDur)) +
  geom_tile() +
  betaJ_threshold_layer +
  scale_fill_viridis_c(
    option = "plasma",
    name = "Mean high\nduration"
  ) +
  labs(
    title = "C. Mean duration of high episodes",
    x = expression(beta*J),
    y = NULL
  ) +
  heat_axes +
  heat_theme +
  guides(fill = guide_colorbar(
    barwidth = unit(4.2, "cm"),
    barheight = unit(0.35, "cm"),
    title.position = "top"
  ))

# ============================================================
# Layout: large left panel + stacked right panels
# ============================================================
right_col <- p_sw / p_dur + plot_layout(heights = c(1, 1))

heat_fig <- (p_frac | right_col) +
  plot_layout(widths = c(1.25, 1), guides = "collect") +
  plot_annotation(
    title = "Hysteresis heatmaps (no shocks): occupancy, switching, and dwell-time structure over " %+%
      "parameter space",
    subtitle = "Contour in panel A marks Pr(high) = 0.5. Dashed vertical line indicates the fast-layer threshold βJ = 4."
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10)
  )

heat_fig












library(dplyr)
library(ggplot2)
library(patchwork)

heat_zoom <- readRDS("res_clean/heat_zoom_noshock_hyst2.rds") |>
  mutate(betaJ = as.numeric(betaJ),
         P_base = as.numeric(P_base))

base_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    plot.title = element_text(face = "bold", hjust = 0, size = 14),
    plot.subtitle = element_text(hjust = 0, size = 10),
    axis.title = element_text(size = 13),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    legend.position = "bottom",
    plot.margin = margin(6, 6, 6, 6)
  )

# Longer colorbars
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

# ---- Panel A: high-state occupancy ----
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

# ---- Panel B: switching frequency ----
p_sw <- ggplot(heat_zoom, aes(betaJ, P_base)) +
  geom_tile(aes(fill = Switches)) +
  scale_fill_viridis_c(
    option = "plasma",
    trans = "sqrt",
    name = expression(N[switch]),
    guide = cb_sw
  ) +
  labs(
    title = "Switching frequency",
    subtitle = expression("Mean number of switches ("*N[switch]*")"),
    x = expression(beta * J),
    y = NULL
  ) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(margin = margin(r = 12))
  )


# ---- Combine panels ----
fig_heat2 <- (p_frac | p_sw) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.margin = margin(6, 6, 6, 6)
  ) +
  plot_layout(
    widths = c(1, 1),
    guides = "collect"
  )

# Reduce visual separation between panels
fig_heat2 <- fig_heat2 & theme(panel.spacing = unit(2, "mm"))

# Print
fig_heat2

# ggsave("img/Figure_heatmaps.png", fig_heat2, width = 10, height = 6)
# ggsave("img/Figure_heatmaps_2panel_plasma_clean.pdf", fig_heat2, width = 10, height = 6)



















#| label: scenario-B2-hysteresis-paper-figure
#| include: true
#| fig-width: 12
#| fig-height: 7.5

library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(grid)

# -----------------------------
# 1) Data prep (single trajectory)
# -----------------------------
# sim_B2 assumed already loaded:
sim_B2 <- readRDS("res/sim_B2.rds")

dfB <- sim_B2 %>%
  mutate(
    t = as.numeric(t),
    P = as.numeric(P),
    m = as.numeric(m),
    P_base = as.numeric(P_base)
  )

# Turning point of ramp (more robust than hard-coding 0.5*T)
Tup <- dfB$t[which.max(dfB$P_base)]

dfB <- dfB %>%
  mutate(
    phase = if_else(t <= Tup, "Ramp up", "Ramp down"),
    phase = factor(phase, levels = c("Ramp up", "Ramp down"))
  )

# -----------------------------
# 2) Jump detection (robust helper)
# -----------------------------
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

# -----------------------------
# 3) Styling
# -----------------------------
# Colorblind-safe-ish phase colors (Okabe-Ito family)
col_up   <- "#0072B2"   # blue
col_down <- "#D55E00"   # vermillion
col_phase <- c("Ramp up" = col_up, "Ramp down" = col_down)

col_base <- "grey40"    # imposed P_base(t)
col_P    <- "grey10"    # realized P(t)
col_jump <- "grey45"

phase_rect <- tibble::tibble(
  xmin = c(min(dfB$t, na.rm = TRUE), Tup),
  xmax = c(Tup, max(dfB$t, na.rm = TRUE)),
  ymin = -Inf,
  ymax = Inf,
  phase = factor(c("Ramp up", "Ramp down"), levels = c("Ramp up", "Ramp down"))
)

base_theme <- theme_minimal(base_size = 12) +
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

# -----------------------------
# 4) Panel A: P(t) and imposed ramp
# -----------------------------
pP <- ggplot(dfB, aes(x = t)) +
  geom_rect(
    data = phase_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phase),
    inherit.aes = FALSE, alpha = 0.06, colour = NA
  ) +
  geom_line(aes(y = P, linetype = "realized"), colour = col_P, linewidth = 0.9) +
  geom_line(aes(y = P_base, linetype = "base"), colour = col_base, linewidth = 0.9) +
  geom_vline(xintercept = Tup, linetype = "dotted", colour = "grey55", linewidth = 0.5) +
  annotate(
    "text",
    x = Tup, y = max(c(dfB$P, dfB$P_base), na.rm = TRUE) + 0.04,
    label = "Ramp reversal",
    colour = "grey35", size = 3.1, vjust = 0
  ) +
  scale_fill_manual(values = col_phase, guide = "none") +
  scale_linetype_manual(
    name = NULL,
    values = c(realized = "solid", base = "dashed"),
    breaks = c("realized", "base"),
    labels = c(
      realized = expression("Realized " * P[t]),
      base = expression("Imposed " * P[base](t))
    )
  ) +
  guides(
    linetype = guide_legend(
      keywidth = unit(1.0, "cm"),   # <- longer key so dash is obvious
      keyheight = unit(0.4, "cm"),
      override.aes = list(
        colour = c(col_P, col_base),
        linewidth = c(0.9, 0.9)
      ),
      order = 1
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +  # top headroom
  coord_cartesian(clip = "off") +
  labs(
    title = "Context ramp protocol and realized context",
    subtitle = expression("Dashed line: imposed " * P[base](t) * "; solid line: realized " * P[t] * "."),
    x = NULL,
    y = expression(P[t])
  ) +
  base_theme +
  theme(plot.margin = margin(12, 8, 6, 6))

# -----------------------------
# 5) Panel B: m(t) with phase and jump markers
# -----------------------------
pm <- ggplot(dfB, aes(x = t, y = m)) +
  geom_rect(
    data = phase_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phase),
    inherit.aes = FALSE,
    alpha = 0.06,
    colour = NA
  ) +
  # plot separate phase segments for a clean phase legend
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
  guides(
    linetype = guide_legend(
      keywidth = unit(1.0, "cm"),   # <- longer key so dash is obvious
      keyheight = unit(0.4, "cm"),
      override.aes = list(
        colour = c(col_P, col_base),
        linewidth = c(0.9, 0.9)
      ),
      order = 1
    )
  )+
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
  theme(plot.margin = margin(6, 8, 6, 6))+
  base_theme

# -----------------------------
# 6) Panel C: hysteresis loop (P, m)
# -----------------------------
# Split path by phase so line colors map cleanly to legend
loop_up   <- dfB %>% filter(phase == "Ramp up")
loop_down <- dfB %>% filter(phase == "Ramp down")

# A few directional arrows to emphasize traversal (optional)
arrow_idx_up <- round(seq(200, nrow(loop_up) - 200, length.out = 3))
arrow_idx_dn <- round(seq(200, nrow(loop_down) - 200, length.out = 3))

arrows_up <- if (nrow(loop_up) > 250) {
  tibble::tibble(
    x = loop_up$P[pmax(1, arrow_idx_up - 15)],
    y = loop_up$m[pmax(1, arrow_idx_up - 15)],
    xend = loop_up$P[arrow_idx_up],
    yend = loop_up$m[arrow_idx_up],
    phase = "Ramp up"
  )
} else tibble::tibble(x=numeric(),y=numeric(),xend=numeric(),yend=numeric(),phase=character())

arrows_dn <- if (nrow(loop_down) > 250) {
  tibble::tibble(
    x = loop_down$P[pmax(1, arrow_idx_dn - 15)],
    y = loop_down$m[pmax(1, arrow_idx_dn - 15)],
    xend = loop_down$P[arrow_idx_dn],
    yend = loop_down$m[arrow_idx_dn],
    phase = "Ramp down"
  )
} else tibble::tibble(x=numeric(),y=numeric(),xend=numeric(),yend=numeric(),phase=character())

arrow_df <- bind_rows(arrows_up, arrows_dn) %>%
  mutate(phase = factor(phase, levels = c("Ramp up", "Ramp down")))

ploop <- ggplot(dfB, aes(x = P, y = m, colour = phase)) +
  geom_path(linewidth = 0.8, alpha = 0.8) +
  geom_point(data = jump_pts, aes(P, m), size = 2.2, stroke = 0.5, fill = "white", shape = 21) +
  scale_colour_manual(
    values = c("Ramp up" = col_up, "Ramp down" = col_down),
    breaks = c("Ramp up", "Ramp down"),
    guide = "none"   # <- HARD OFF (won't be revived by patchwork theme)
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
  base_theme

# -----------------------------
# 7) Assemble: left stack (A,B), right loop (C)
# -----------------------------
left_col <- pP / pm + plot_layout(heights = c(1, 1))

hyst_fig <- (left_col | ploop) +
  plot_layout(widths = c(1.1, 1), guides = "collect") +
  plot_annotation(
    # title = "Hysteresis under a slow exogenous context ramp (no feedback)",
    # subtitle = "Fast layer bistable-capable (βJ > 4), with feedback off (b = 0). Path dependence appears as branch-dependent jumps and a loop in (P, m).",
    tag_levels = "A"
  ) &
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12)
  )

hyst_fig


# ggsave("img/Figure_hysteresis.pdf", hyst_fig, width = 13, height = 7)


