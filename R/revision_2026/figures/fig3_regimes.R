# ============================================================
# R/revision_2026/figures/fig3_regimes.R
# ============================================================
# Figure 3 (REFRAMED): Feedback, recovery, and history-dependent regimes.
#
# This REPLACES fig3_recovery_feedback.R as the main-text Figure 3, and
# promotes the history-dependence regime check (05_supp_regime_history_
# dependence.R / figS_history_dependence.R) from supplementary evidence
# to a main theoretical figure. That is a real change in what the paper
# claims, not just a layout change -- explicitly decided over the
# earlier, more cautious "keep it supplementary" guidance, because the
# figure sequence otherwise reads as four separate demonstrations with
# no central model insight (the role the old mean-field bistability
# figure used to play, before it was dropped per Denny's feedback).
#
# The new claim, made only from the revised (0/1) heterogeneous-
# threshold model (never from the old mean-field approximation):
#   - under the LOCKED moderate-feedback regime (b=0.50), the system
#     recovers -- low/high context perturbations decay back out
#   - under a STRONGER feedback regime (b=1.0, still non-runaway --
#     max|P| ~2 there, vs. a saturation threshold of 6), the same
#     architecture becomes history-dependent: identical dynamics started
#     from a low- vs. high-burden state stay separated instead of
#     reconverging
#   - that separation appears cleanly across a b-grid, not just at one
#     cherry-picked value
#
# 2026-08-28: MERGED FROM THREE PANELS TO TWO. The original design had a
# separate panel A (feedback off vs. the locked b=0.50 case only) and
# panel B (the same shock-and-recovery design swept across a wider
# b-grid, deliberately excluding b=0/0.50 to avoid duplicating panel A).
# Once panel A's window was extended to match panel B's (see git history
# for that intermediate version), the two panels were plotting the
# identical experiment on the identical timescale with heavy overlap in
# purpose -- per review, panel A wasn't earning its own space anymore.
# Folded panel A's two lines (b=0 "feedback off", b=0.50 "the locked
# main-text value") directly INTO the grid panel below, so there is now
# one shock-and-recovery panel covering the full range from "no feedback"
# to "runaway/plateau," plus the regime-index panel. b=0 and b=0.50 keep
# panel A's individual-chain spaghetti layer (the two "headline" lines);
# the other grid points show mean + ribbon only, to avoid the
# already-flagged problem of individual-chain spaghetti becoming
# illegible past ~4 simultaneous lines.
#
# Two panels:
#   (A) Shock-and-recovery trajectories across the feedback-strength grid
#       (b = 0, 0.5, 0.75, 0.9, 1.0, 1.3; from
#       03c_sim_feedback_shock_grid_extended.R), one line per b, showing
#       the shift from "recovers" (b=0, b=0.50) to "settles at an
#       elevated plateau" as b increases.
#   (B) History-dependence index across feedback strength (the full
#       b-grid from the regime-check script), reconvergent vs.
#       initial-state-dependent, with the locked main-text b explicitly
#       marked.
#
# FIRST VERSION OF PANEL B (pre-merge), DROPPED: raw low-/high-initial-
# state trajectories at b=1.0, from the regime-check script's NO-SHOCK
# design (chains started directly at a low- vs. high-burden state, not
# perturbed mid-run). That tests initial-condition dependence, a real
# question, but a genuinely different one from the shock-recovery panel's
# "does a shock recover" -- it never shared that panel's shock-at-t=0
# structure, which made it read as unclear/unmotivated alongside it. It
# also had a second, independent problem: individual per-chain noise was
# comparable in size to the actual between-group separation, so the
# raw-trajectory view read as two overlapping noisy clouds rather than
# two separated regimes, even though the underlying means were cleanly
# separated.
#
# Both panels plot symptom burden (m) / a burden-derived index, sharing
# one visual language.
#
# Legend: colour (+ linetype for the off/on pair) identify the compared
# conditions. Per review, the spaghetti (individual-chain) layer for
# b=0/0.50 is included in the legend-generating aesthetic mapping (not
# show.legend=FALSE), so the legend key shows both the thin/faint
# individual-trajectory style and the bold mean-line style together --
# making explicit, via the legend itself, that the figure is built from
# many raw trajectories, not just a single computed curve.
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

# 2026-08-27: local (THIS FIGURE ONLY) darker gray, used in place of the
# shared theme_publication.R col_off ("#7F7F7F") wherever this script
# encodes "feedback off"/"reconvergent" -- the shared mid-gray read as too
# washed-out against white, especially under the ribbon/spaghetti alpha
# layers. Not changed in theme_publication.R itself, since col_off is
# reused elsewhere (e.g. Figure 4's "Fixed P" arm) where this hasn't been
# flagged as an issue -- no reason to risk changing figures that work.
gray_off_dark <- "#4D4D4D"

# 2026-08-27: text sizes bumped throughout this figure specifically (axis
# text/titles, panel titles, legend text) -- theme_pub()'s shared sizes
# (axis text 9pt, axis title 10pt, panel title 10.5pt) are FIXED
# regardless of the base_size argument (see theme_publication.R's own
# note), so bumping them requires local theme() overrides per panel,
# applied below in pGrid/pB. Only this figure is affected.
# bumped 2026-08-28 (axis_text 11->13, axis_title 12->14, legend_text
# 10.5->12.5, panel_title 13->15) per feedback that text was too small,
# then axis text/title dialed back slightly (13->12, 14->13) per
# follow-up feedback that it had overshot.
axis_text_sz  <- 12
axis_title_sz <- 13
panel_title_sz <- 15
legend_text_sz <- 12.5

# ------------------------------------------------------------------------
# Panel A: shock-and-recovery trajectories across the feedback-strength
# grid (03c_sim_feedback_shock_grid_extended.R), b = 0, 0.5, 0.75, 0.9,
# 1.0, 1.3 -- one line per b, so the reader sees the shift from "recovers"
# (b=0, b=0.50) to "settles at an elevated plateau" as b increases, all
# in one panel.
#
# 2026-08-28: b=0 ("feedback off") and b=0.50 (the locked main-text
# value) folded back INTO this grid (previously shown separately in a
# now-removed panel A, then deliberately excluded from this panel to
# avoid duplicating it -- see the top-of-file merge note). Re-sourced
# from 03c_sim_feedback_shock_grid_extended.R throughout, windowed to
# 0-2500 steps so the off/on lines are shown actually reconverging
# rather than cut off mid-recovery.
#
# Caveat, stated plainly: the peak/end-of-window NUMBERS quoted in the
# Results paragraph (e.g. "4.34 active symptoms," "2.57 active symptoms")
# still come from the ORIGINAL locked Simulation 2/3 run, not from this
# 03c companion run -- 03c uses identical parameters but is a separate
# simulation (different random draws), so the b=0/b=0.50 curves here will
# be qualitatively identical but not pixel-identical to those exact
# quoted numbers.
#
# b=1.10 was checked against the rerun (2500-step window) and DROPPED:
# unlike b=1.30 (which visibly flattens, rising only ~0.03 over the last
# 1000 steps) and b=1.00 (~0.02 residual drift, essentially flat),
# b=1.10 was still rising by ~0.08 over the same span with no sign of
# leveling off -- plausibly critical slowing down near the transition
# itself, but not something we can currently show as "resolved" next to
# genuinely plateaued curves without overclaiming. Left out for the same
# reason 1.25/1.50 were left out of the original 1500-step version.
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
# 2026-08-28: colour ramps (RGB, then Lab-space) and the Set3 qualitative
# palette were all tried here and all still left some adjacent b's too
# close in hue to tell apart at a glance. Switched to a HAND-PICKED set of
# 6 colours chosen for maximum pairwise separation around the colour
# wheel (gray, green, blue, violet, amber, red) rather than any
# interpolated/automatic palette -- deliberately breaks from the "no
# rainbow" restraint used elsewhere in this figure set, because with 6
# ordered lines sharing one panel, legibility here matters more than
# palette minimalism.
grid_pal <- c(
  "0"    = gray_off_dark,  # feedback off
  "0.5"  = col_low,        # locked main-text value
  "0.75" = col_on,
  "0.9"  = "#8E6FAE",
  "1"    = "#D4A017",
  "1.3"  = col_high
)
grid_pal <- grid_pal[as.character(b_levels)]  # order/subset to match whatever b's actually ended up in the data
# Two label sets (2026-08-28): the in-panel end-of-line labels now just
# say "b = 0" for the zero line, matching the plain "b = X.XX" style of
# the other five lines -- "Feedback off (b = 0)" was cluttering the plot
# area. The fuller "Feedback off (b = 0)" phrasing is kept for the LEGEND
# only, where there's room for it and it's more helpful to a reader
# scanning the legend cold.
grid_labels_plot <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_plot[as.character(0)] <- "b = 0"

grid_labels_legend <- setNames(sprintf("b = %.2f", b_levels), b_levels)
grid_labels_legend[as.character(0)] <- "Feedback off (b = 0)"

# Individual-chain spaghetti: kept ONLY for b=0 and b=0.50 (the two
# "headline" lines, formerly panel A's dedicated off/on comparison) --
# with 6 lines now in this panel, spaghetti across ALL of them would
# repeat the earlier-established problem (illegible past ~4 simultaneous
# lines). The other four grid points show mean + ribbon only.
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

# grid_w now DOES contain b=0 again -- pre-shock baseline computed
# directly from those rows (feedback doesn't change pre-shock dynamics,
# so b=0 is as valid a source for this reference line as any other b).
m_ref_grid <- mean(grid_w$mean_m_smooth[grid_w$b == 0 & grid_w$time_since_shock < 0], na.rm = TRUE)
M_ref_grid <- m_ref_grid * N
x_max_grid <- plot_window_grid[2]

grid_end <- grid_w |> filter(!is.na(mean_M_smooth)) |> group_by(b) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup() |>
  arrange(mean_M_smooth)
# Spread out end-labels that would otherwise collide (several b's plateau
# close together) -- same fix-class as panel A's off/on label collision.
# label_gap is the original m-scale value (0.045) x N=9, now that
# grid_end is on the count scale.
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

# Local tweaks (2026-08-28, this panel only): mean lines thinned slightly
# from the shared main_line_width (now that 6 lines share the panel, a
# touch thinner reads as less heavy/overlapping); end-of-line "b = ..."
# labels enlarged from 2.8 to 3.4 for legibility.
grid_line_width <- main_line_width * 0.82
grid_label_size <- 4.4  # bumped again 2026-08-28 (was 3.4, then 2.8 originally)

# y-axis zoomed to the data's actual range (2026-08-28), rather than
# ggplot's default expansion -- b=0 and b=0.50 both recover to nearly the
# same baseline, and with 6 lines sharing one axis (including b=1.30's
# much higher plateau) that small gap was getting visually compressed.
# Tightening the range + padding, combined with giving this panel more
# vertical room in the final layout below, makes better use of the
# panel's height for exactly the part of the range where the lines are
# close together.
#
# FIX (2026-08-28): the range above was computed only from the mean+CI
# band (grid_w), not the individual-chain spaghetti (grid_chains) drawn
# underneath it -- raw per-chain noise swings well outside the mean's CI,
# so coord_cartesian(ylim=...) was silently clipping the spaghetti lines
# flat at the top/bottom of the panel instead of showing their actual
# peaks/troughs. Range now also includes grid_chains$M so the axis covers
# everything actually drawn, not just the mean lines.
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
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +  # 6 entries now (was 4) -- wrap to 2 rows so labels don't crowd/shrink in one row
  labs(
    title = panel_title("A", "Feedback strength changes whether the system recovers"),
    x = "Steps since perturbation", y = "Mean number of active symptoms"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    # right margin bumped 62 -> 72 -> 88 (2026-08-28): end-of-line labels
    # are now considerably larger (grid_label_size 4.4, was 2.8 originally)
    # and need more headroom to avoid clipping, especially the long
    # "Feedback off (b = 0)" label -- check the rendered PNG for clipping,
    # this is an estimate, not a measured value.
    plot.margin = margin(5.5, 88, 5.5, 5.5),
    legend.position = "bottom",
    legend.key.width = unit(0.9, "cm"),
    legend.text = element_text(size = legend_text_sz),
    plot.title = element_text(size = panel_title_sz),
    axis.title = element_text(size = axis_title_sz),
    axis.text = element_text(size = axis_text_sz)
  )

sep <- read_csv("res/revision_2026/supp_history/history_separation_summary.csv", show_col_types = FALSE)

# gray_off_dark in place of col_off (2026-08-27, same fix as panel A) --
# "reconvergent" points were rendering in the same washed-out gray that
# was hard to distinguish from the panel background.
regime_pal <- c(convergent = gray_off_dark, `history-dependent` = col_high, `runaway/saturation` = "black")

# Locate the locked main-text b (0.50) by nearest value, not exact
# equality -- the grid is now generated via seq(0, 1.5, by = 0.1), and
# exact `== 0.5` comparisons against a seq()-built float are a real risk
# of silently matching nothing due to floating-point drift.
locked_b_row <- sep[which.min(abs(sep$b - 0.5)), ]

pB <- ggplot(sep, aes(x = b, y = separation_m)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "grey55", linewidth = 0.4) +
  # Connecting line -- now that the grid is dense (16 points), this turns
  # the panel from a scatter of dots into an actual continuous regime-
  # transition curve. Drawn in a single neutral colour underneath the
  # per-point regime colouring, so the point colours still carry the
  # convergent/history-dependent classification.
  geom_line(colour = "grey45", linewidth = 0.6) +
  geom_pointrange(
    aes(ymin = separation_m - interval_mult * se_separation_m, ymax = separation_m + interval_mult * se_separation_m, colour = regime_flag),
    linewidth = 0.75, size = 0.5
  ) +
  # Reverted to the circled-point + text callout per review -- the plain
  # vertical-line version was tried but the circle+text combo was
  # preferred after all.
  annotate("point", x = locked_b_row$b, y = locked_b_row$separation_m, shape = 21, size = 5.5, colour = col_on, stroke = 1.0) +
  # text size bumped 2.7 -> 3.8 (2026-08-28) per feedback -- check the
  # rendered PNG for collision with nearby points/gridlines now that the
  # label is noticeably bigger.
  annotate("text", x = locked_b_row$b, y = locked_b_row$separation_m + 0.07, label = "main simulation\nvalue", size = 3.8, colour = col_on, hjust = 0.5, vjust = 0, lineheight = 0.85) +
  # 2026-08-27: relabeled at the DISPLAY level only -- regime_flag's
  # underlying values ("convergent"/"history-dependent", from
  # 05_supp_regime_history_dependence.R) are unchanged, so this doesn't
  # require rerunning that script. "Reconvergent"/"initial-state-
  # dependent" reads clearer for a psych audience than "convergent"/
  # "history-dependent", and matches the manuscript prose now avoiding
  # "history-dependence" language.
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
# together). Each panel keeps its OWN local legend rather than collecting
# into one shared legend -- A/B's colour scales are two different
# variables, not levels of one.
# ------------------------------------------------------------------------
# height ratio bumped 1.15 -> 1.5 (2026-08-28), and canvas height 7.6 ->
# 8.4, so panel A's now-tightened y-range (grid_ylim above) gets more
# vertical pixels to resolve the close b=0/b=0.50 gap -- check the
# rendered PNG; if panel A now looks disproportionately tall next to B,
# scale this back down a little.
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
