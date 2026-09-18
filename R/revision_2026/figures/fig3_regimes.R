# ============================================================
# R/revision_2026/figures/fig3_regimes.R
# ============================================================
# Figure 3 (reframed): feedback, recovery, and history-dependent regimes.
#
# Replaces fig3_recovery_feedback.R as the main-text Figure 3, and
# promotes the history-dependence regime check (05_supp_regime_history_
# dependence.R / figS_history_dependence.R) from supplementary evidence to
# a main theoretical figure. That's a real change in what the paper
# claims, not just layout -- decided over the earlier, more cautious
# "keep it supplementary" call, because otherwise the figure sequence
# reads as four separate demonstrations with no central model insight
# (the role the old mean-field bistability figure used to play, before it
# got dropped).
#
# The new claim, made only from the revised (0/1) heterogeneous-threshold
# model, never from the old mean-field approximation:
#   - under the locked moderate-feedback regime (b=0.50), the system
#     recovers -- low/high context perturbations decay back out
#   - under a stronger feedback regime (b=1.0, still non-runaway -- max|P|
#     ~2 there vs. a saturation threshold of 6), the same architecture
#     becomes history-dependent: identical dynamics started from a low-
#     vs. high-burden state stay separated instead of reconverging
#   - that separation shows up cleanly across a b-grid, not just at one
#     cherry-picked value
#
# 2026-08-28: merged from three panels to two. The original design had a
# separate panel A (feedback off vs. the locked b=0.50 only) and panel B
# (the same shock-and-recovery design swept across a wider b-grid,
# deliberately excluding b=0/0.50 to avoid duplicating panel A). Once
# panel A's window matched panel B's, the two panels were plotting the
# same experiment on the same timescale, so panel A wasn't earning its own
# space anymore. Folded its two lines (b=0, b=0.50) directly into the grid
# panel below -- one shock-and-recovery panel now covers the full range
# from "no feedback" to "runaway/plateau," plus the regime-index panel.
# b=0 and b=0.50 keep the individual-chain spaghetti layer (the two
# "headline" lines); the other grid points show mean + ribbon only, since
# spaghetti across more than ~4 lines gets illegible.
#
# Two panels:
#   (A) Shock-and-recovery trajectories across the feedback-strength grid
#       (b = 0, 0.5, 0.75, 0.9, 1.0, 1.3; from
#       03c_sim_feedback_shock_grid_extended.R), one line per b: recovers
#       (b=0, b=0.50) to settles at an elevated plateau as b increases.
#   (B) History-dependence index across feedback strength (full b-grid
#       from the regime-check script), reconvergent vs.
#       initial-state-dependent, with the locked main-text b marked.
#
# Panel B's first version, dropped: raw low-/high-initial-state
# trajectories at b=1.0 from the regime-check script's no-shock design
# (chains started directly at a low- vs. high-burden state, not perturbed
# mid-run). That tests a real but different question -- initial-condition
# dependence rather than "does a shock recover" -- and never shared panel
# A's shock-at-t=0 structure, so it read as unmotivated alongside it.
# Individual per-chain noise was also comparable in size to the
# between-group separation, so it looked like two overlapping noisy
# clouds even though the underlying means were cleanly separated.
#
# Both panels plot symptom burden (m) / a burden-derived index, sharing
# one visual language. Legend: colour (+ linetype for the off/on pair)
# identifies the compared conditions. The spaghetti layer for b=0/0.50 is
# included in the legend's aesthetic mapping (not show.legend=FALSE), so
# the legend key shows both the thin/faint individual-trajectory style and
# the bold mean-line style together, making explicit that the figure is
# built from many raw trajectories, not one computed curve.
#
# Outputs
# -------
#   figs/revision_2026/Figure3_regimes.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# RColorBrewer used (namespaced, not attached) for panel A's Set3 palette --
# install.packages("RColorBrewer") if not already installed.
source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_*, pal_feedback

roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
base_sz <- base_size_panel_title  # 10.5, from theme_publication.R's global tokens
interval_mult <- 1.96

# Local (this figure only) darker gray in place of theme_publication.R's
# col_off ("#7F7F7F") wherever this script encodes "feedback off"/
# "reconvergent" -- the shared mid-gray read as too washed-out against
# white under the ribbon/spaghetti alpha layers. Not changed in
# theme_publication.R itself since col_off is reused elsewhere (e.g.
# Figure 4's "Fixed P" arm) where this hasn't been an issue.
gray_off_dark <- "#4D4D4D"

# Text sizes bumped throughout this figure specifically -- theme_pub()'s
# shared sizes are fixed regardless of the base_size argument, so bumping
# them needs local theme() overrides per panel (pGrid/pB below). Only this
# figure is affected; sizes settled after a couple of rounds of feedback
# on legibility.
axis_text_sz  <- 12
axis_title_sz <- 13
panel_title_sz <- 15
legend_text_sz <- 12.5

# ------------------------------------------------------------------------
# Panel A: shock-and-recovery trajectories across the feedback-strength
# grid (03c_sim_feedback_shock_grid_extended.R), b = 0, 0.5, 0.75, 0.9,
# 1.0, 1.3 -- one line per b, showing the shift from "recovers" (b=0,
# b=0.50) to "settles at an elevated plateau" as b increases, all in one
# panel.
#
# b=0 ("feedback off") and b=0.50 (the locked main-text value) are folded
# back into this grid (previously shown separately, see the merge note
# above). Re-sourced from 03c_sim_feedback_shock_grid_extended.R,
# windowed to 0-2500 steps so the off/on lines show actual reconvergence
# rather than cutting off mid-recovery.
#
# Caveat: the peak/end-of-window NUMBERS quoted in the Results paragraph
# (e.g. "4.34 active symptoms," "2.57 active symptoms") still come from
# the original locked Simulation 2/3 run, not from this 03c companion run --
# same parameters, different random draws, so the b=0/b=0.50 curves here
# are qualitatively identical but not pixel-identical to those exact numbers.
#
# b=1.10 was checked against the rerun (2500-step window) and dropped:
# unlike b=1.30 (visibly flattens, rising only ~0.03 over the last 1000
# steps) and b=1.00 (~0.02 residual drift, essentially flat), b=1.10 was
# still rising by ~0.08 over the same span with no sign of leveling off --
# plausibly critical slowing near the transition itself, but not something
# to show as "resolved" next to genuinely plateaued curves. Left out for
# the same reason 1.25/1.50 were left out of the original 1500-step version.
display_b_grid <- c(0, 0.50, 0.75, 0.90, 1.00, 1.30)

sim3c <- read_csv("res/revision_2026/sim3c/sim3c_summary.csv", show_col_types = FALSE)
sim3c_raw <- readRDS("res/revision_2026/sim3c/sim3c_raw.rds")$traj

plot_window_grid <- c(-100, 2500)
N <- 9

grid_w <- sim3c |> filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2]) |>
  arrange(b, time_since_shock) |>
  group_by(b) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N)) |>
  ungroup() |>
  mutate(mean_M_smooth = mean_m_smooth * N, se_M_smooth = se_m_smooth * N)  # count-scale, see header note

b_levels <- sort(unique(grid_w$b))
# Colour ramps (RGB, then Lab-space) and the Set3 qualitative palette were
# all tried here and still left some adjacent b's too close in hue to
# tell apart at a glance. Switched to a hand-picked set of 6 colours for
# maximum pairwise separation (gray, green, blue, violet, amber, red)
# instead of any interpolated/automatic palette -- a deliberate break from
# the "no rainbow" restraint used elsewhere, because with 6 ordered lines
# sharing one panel, legibility matters more than palette minimalism here.
grid_pal <- c(
  "0"    = gray_off_dark,  # feedback off
  "0.5"  = col_low,        # locked main-text value
  "0.75" = col_on,
  "0.9"  = "#8E6FAE",
  "1"    = "#D4A017",
  "1.3"  = col_high
)
grid_pal <- grid_pal[as.character(b_levels)]  # order/subset to match whatever b's actually ended up in the data
# In-panel end-of-line labels just say "b = 0" for the zero line, matching
# the plain "b = X.XX" style of the other five -- "Feedback off (b = 0)"
# cluttered the plot area. The fuller phrasing is kept for the legend
# only, where there's room and it reads better to someone scanning cold.
grid_labels_plot <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_plot[as.character(0)] <- "b = 0"

grid_labels_legend <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_legend[as.character(0)] <- "Feedback off (b = 0)"

# Individual-chain spaghetti kept ONLY for b=0 and b=0.50 (the two
# "headline" lines) -- with 6 lines now in this panel, spaghetti across
# all of them would repeat the illegibility problem from earlier drafts.
# The other four grid points show mean + ribbon only.
set.seed(3)
n_sample_chains_headline <- 30L
set.seed(11)
n_sample_chains_grid <- 10L

sim3c_raw_grid <- sim3c_raw |>
  filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2])

headline_chains <- sim3c_raw_grid |>
  filter(b %in% c(0, 0.5)) |>
  group_by(b) |>
  filter(chain %in% sample(unique(chain), n_sample_chains_headline)) |>
  ungroup()

other_grid_chains <- sim3c_raw_grid |>
  filter(!b %in% c(0, 0.5)) |>
  group_by(b) |>
  filter(chain %in% sample(unique(chain), n_sample_chains_grid)) |>
  ungroup()

grid_chains <- bind_rows(headline_chains, other_grid_chains) |>
  arrange(b, chain, time_since_shock) |>
  group_by(b, chain) |>
  mutate(m = roll_mean(m)) |>
  ungroup() |>
  filter(!is.na(m)) |>
  mutate(M = m * N)  # count-scale version for plotting, see grid_w note above

# grid_w includes b=0 again, so the pre-shock baseline is computed
# directly from those rows -- feedback doesn't change pre-shock dynamics,
# so b=0 is as valid a source for this reference line as any other b.
m_ref_grid <- mean(grid_w$mean_m_smooth[grid_w$b == 0 & grid_w$time_since_shock < 0], na.rm = TRUE)
M_ref_grid <- m_ref_grid * N
x_max_grid <- plot_window_grid[2]

grid_end <- grid_w |> filter(!is.na(mean_M_smooth)) |> group_by(b) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup() |>
  arrange(mean_M_smooth)
# Spread out end-labels that would otherwise collide (several b's plateau
# close together). label_gap is the original m-scale value (0.045) x N=9,
# now that grid_end is on the count scale.
label_gap <- 0.045 * N
if (nrow(grid_end) > 1) {
  y <- grid_end$mean_M_smooth
  for (i in 2:length(y)) {
    if (y[i] - y[i - 1] < label_gap) y[i] <- y[i - 1] + label_gap
  }
  grid_end$label_y <- y
} else {
  grid_end$label_y <- grid_end$mean_M_smooth
}

# Local tweaks, this panel only: mean lines thinned slightly from the
# shared main_line_width (6 lines sharing the panel reads better a touch
# thinner); end-of-line "b = ..." labels enlarged for legibility.
grid_line_width <- main_line_width * 0.82
grid_label_size <- 4.4

# y-axis zoomed to the data's actual range rather than ggplot's default
# expansion -- b=0 and b=0.50 both recover to nearly the same baseline,
# and with 6 lines sharing one axis (including b=1.30's much higher
# plateau) that small gap was getting visually compressed. The range
# below includes both the mean+CI band (grid_w) and the individual-chain
# spaghetti (grid_chains): computing it from grid_w alone let raw
# per-chain noise swing outside the mean's CI and get silently clipped by
# coord_cartesian(ylim=...) at the top/bottom of the panel.
grid_y_lo <- min(
  grid_w$mean_M_smooth - interval_mult * grid_w$se_M_smooth,
  grid_chains$M,
  na.rm = TRUE
)
grid_y_hi <- max(
  grid_w$mean_M_smooth + interval_mult * grid_w$se_M_smooth,
  grid_chains$M,
  na.rm = TRUE
)
grid_y_pad <- 0.05 * (grid_y_hi - grid_y_lo)
grid_ylim <- c(grid_y_lo - grid_y_pad, grid_y_hi + grid_y_pad)

pGrid <- ggplot(grid_w, aes(x = time_since_shock, y = mean_M_smooth, colour = factor(b), group = factor(b))) +
  geom_hline(yintercept = M_ref_grid, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = grid_chains,
    aes(x = time_since_shock, y = M, group = interaction(chain, b), colour = factor(b)),
    inherit.aes = FALSE, alpha = spaghetti_alpha, linewidth = spaghetti_width
  ) +
  geom_ribbon(
    aes(ymin = mean_M_smooth - interval_mult * se_M_smooth, ymax = mean_M_smooth + interval_mult * se_M_smooth, fill = factor(b)),
    alpha = ribbon_alpha, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = grid_line_width) +
  geom_text(
    data = grid_end, aes(x = x_max_grid + 30, y = label_y, label = grid_labels_plot[as.character(b)], colour = factor(b)),
    hjust = 0, size = grid_label_size, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(xlim = c(plot_window_grid[1], x_max_grid), ylim = grid_ylim, clip = "off") +
  scale_colour_manual(values = grid_pal, labels = grid_labels_legend, name = NULL) +
  scale_fill_manual(values = grid_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.11))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +  # 6 entries now: wrap to 2 rows so labels don't crowd
  labs(
    title = panel_title("A", "Feedback strength changes whether the system recovers"),
    x = "Steps since perturbation", y = "Mean number of active symptoms"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    # Right margin bumped up: end-of-line labels are considerably larger
    # now and need headroom to avoid clipping, especially the long
    # "Feedback off (b = 0)" label -- check the rendered PNG for clipping.
    plot.margin = margin(5.5, 88, 5.5, 5.5),
    legend.position = "bottom",
    legend.key.width = unit(0.9, "cm"),
    legend.text = element_text(size = legend_text_sz),
    plot.title = element_text(size = panel_title_sz),
    axis.title = element_text(size = axis_title_sz),
    axis.text = element_text(size = axis_text_sz)
  )

sep <- read_csv("res/revision_2026/supp_history/history_separation_summary.csv", show_col_types = FALSE)

# gray_off_dark in place of col_off, same fix as panel A -- "reconvergent"
# points were rendering in the same washed-out gray that was hard to
# distinguish from the panel background.
regime_pal <- c(convergent = gray_off_dark, `history-dependent` = col_high, `runaway/saturation` = "black")

# Locate the locked main-text b (0.50) by nearest value, not exact
# equality -- the grid is generated via seq(0, 1.5, by = 0.1), and exact
# `== 0.5` against a seq()-built float risks silently matching nothing due
# to floating-point drift.
locked_b_row <- sep[which.min(abs(sep$b - 0.5)), ]

pB <- ggplot(sep, aes(x = b, y = separation_m)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  # Connecting line, drawn in a neutral colour underneath the per-point
  # regime colouring, so the point colours still carry the
  # convergent/history-dependent classification -- turns the panel from a
  # scatter of dots into an actual continuous regime-transition curve now
  # that the grid is dense (16 points).
  geom_line(colour = "grey45", linewidth = 0.6) +
  geom_pointrange(
    aes(ymin = separation_m - interval_mult * se_separation_m, ymax = separation_m + interval_mult * se_separation_m, colour = regime_flag),
    linewidth = 0.75, size = 0.5
  ) +
  annotate("point", x = locked_b_row$b, y = locked_b_row$separation_m, shape = 21, size = 5.5, colour = col_on, stroke = 1.0) +
  annotate("text", x = locked_b_row$b, y = locked_b_row$separation_m + 0.07, label = "main simulation\nvalue", size = 3.8, colour = col_on, hjust = 0.5, vjust = 0, lineheight = 0.85) +
  # Relabeled at the display level only -- regime_flag's underlying values
  # ("convergent"/"history-dependent", from
  # 05_supp_regime_history_dependence.R) are unchanged, no need to rerun
  # that script. "Reconvergent"/"initial-state-dependent" reads clearer
  # for a psych audience and matches the manuscript prose.
  scale_colour_manual(
    values = regime_pal, name = NULL,
    labels = c(
      convergent = "reconvergent",
      `history-dependent` = "initial-state-dependent",
      `runaway/saturation` = "runaway/saturation"
    )
  ) +
  labs(
    title = panel_title("B", "Initial-state dependence emerges at higher feedback"),
    x = "Feedback strength (b)", y = "Difference in late mean symptom activation"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = legend_text_sz),
    plot.title = element_text(size = panel_title_sz),
    axis.title = element_text(size = axis_title_sz),
    axis.text = element_text(size = axis_text_sz)
  )

# ------------------------------------------------------------------------
# Combine + save -- A (shock-and-recovery across the full b-grid,
# including b=0/b=0.50), B (the regime-index summary tying the grid
# together). Each panel keeps its own local legend rather than one shared
# legend -- A/B's colour scales are two different variables, not levels of one.
# ------------------------------------------------------------------------
fig3 <- (pGrid / pB) +
  plot_layout(heights = c(1.7, 1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_regimes.pdf", fig3, width =12, height = 10)
ggsave("figs/revision_2026/Figure3_regimes.png", fig3, width = 7.8, height = 8.4, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_regimes.pdf (+ .png)\n")
cat("\nThis REPLACES fig3_recovery_feedback.R / FigureS_history_dependence.R\n")
cat("as main-text Figure 3. Those two scripts are kept (not deleted) but are\n")
cat("now superseded for main-text purposes.\n")
