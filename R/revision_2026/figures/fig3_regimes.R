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
# Three panels (a first attempt at panel B -- raw low-/high-initial-state
# trajectories at b=1.0, no shock -- was tried and dropped; see note
# below; this is the SECOND version of panel B, built differently):
#   (A) Moderate feedback (b=0.50): shock response and recovery
#       -- reuses Sim 2 (feedback off) vs. Sim 3 (feedback on) symptom
#       burden, same data/smoothing as the old Figure 3 panel D.
#   (B) The SAME shock-and-recovery design as panel A, swept across the
#       full feedback-strength grid (b = 0, 0.5, 0.75, 1.0, 1.25, 1.5;
#       03c_sim_feedback_shock_grid_extended.R) instead of just off/on --
#       one line per b, all sharing panel A's shock marker and visual
#       language, showing the same experiment shift from "recovers" to
#       "settles at an elevated plateau" as b increases.
#   (C) History-dependence index across feedback strength (the full
#       b-grid from the regime-check script), convergent vs.
#       history-dependent, with the locked main-text b explicitly marked.
#
# FIRST VERSION OF PANEL B, DROPPED: raw low-/high-initial-state
# trajectories at b=1.0, from the regime-check script's NO-SHOCK design
# (chains started directly at a low- vs. high-burden state, not
# perturbed mid-run). That tests initial-condition dependence, a real
# question, but a genuinely different one from panel A's "does a shock
# recover" -- it never shared panel A's shock-at-t=0 structure, which
# made it read as unclear/unmotivated next to panel A. It also had a
# second, independent problem: individual per-chain noise was comparable
# in size to the actual between-group separation, so the raw-trajectory
# view read as two overlapping noisy clouds rather than two separated
# regimes, even though the underlying means were cleanly separated.
# Replaced with the current panel B, which reuses panel A's exact
# experimental design (shock, not initial-condition) at multiple b, so
# it is now a direct visual extension of panel A rather than a
# differently-designed companion.
#
# All three panels plot symptom burden (m) / a burden-derived index,
# sharing one visual language rather than mixing P_t and m_t panels as
# the old 2x2 layout did.
#
# Legend: colour (+ linetype in A) identify the compared conditions in
# each panel. Per review, the spaghetti (individual-chain) layer in
# panel A is included in the legend-generating aesthetic mapping (not
# show.legend=FALSE), so the legend key shows both the thin/faint
# individual-trajectory style and the bold mean-line style together --
# making explicit, via the legend itself, that the figure is built from
# many raw trajectories, not just a single computed curve. Panel B now
# also carries this treatment (added once the panel was trimmed to 4
# lines instead of 6 -- with 6 conditions individual-chain spaghetti
# stopped being legible, the lesson from the first panel-B attempt; 4 is
# fine).
#
# Outputs
# -------
#   figs/revision_2026/Figure3_regimes.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_*, pal_feedback

roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
base_sz <- base_size_panel_title  # 10.5, from theme_publication.R's global tokens
interval_mult <- 1.96

# ------------------------------------------------------------------------
# Panel A data: moderate feedback (b=0.50), symptom burden, off vs on
# ------------------------------------------------------------------------
sim2 <- read_csv("res/revision_2026/sim2/sim2_summary.csv", show_col_types = FALSE)
sim3 <- read_csv("res/revision_2026/sim3/sim3_summary.csv", show_col_types = FALSE)

plot_window_a <- c(-100, 750)
N <- 9

sim2_w <- sim2 |> filter(time_since_shock >= plot_window_a[1], time_since_shock <= plot_window_a[2]) |>
  mutate(feedback = "off")
sim3_on_w <- sim3 |> filter(feedback == "on", time_since_shock >= plot_window_a[1], time_since_shock <= plot_window_a[2])

a_summary <- bind_rows(
  sim2_w |> select(time_since_shock, feedback, mean_m, se_M),
  sim3_on_w |> select(time_since_shock, feedback, mean_m, se_M)
) |>
  mutate(feedback = factor(feedback, levels = c("off", "on"))) |>
  arrange(feedback, time_since_shock) |>
  group_by(feedback) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N)) |>
  ungroup()

m_ref_a <- mean(a_summary$mean_m_smooth[a_summary$feedback == "off" & a_summary$time_since_shock < 0], na.rm = TRUE)

set.seed(3)
n_sample_chains <- 30L
sim2_raw <- readRDS("res/revision_2026/sim2/sim2_raw.rds")$traj |>
  filter(time_since_shock >= plot_window_a[1], time_since_shock <= plot_window_a[2]) |>
  mutate(feedback = "off")
sim3_raw <- readRDS("res/revision_2026/sim3/sim3_raw.rds")$traj |>
  filter(feedback == "on", time_since_shock >= plot_window_a[1], time_since_shock <= plot_window_a[2])

a_chains <- bind_rows(
  sim2_raw |> filter(chain %in% sample(unique(sim2_raw$chain), n_sample_chains)) |>
    select(time_since_shock, chain, feedback, m),
  sim3_raw |> filter(chain %in% sample(unique(sim3_raw$chain), n_sample_chains)) |>
    select(time_since_shock, chain, feedback, m)
) |>
  mutate(feedback = factor(feedback, levels = c("off", "on"))) |>
  arrange(feedback, chain, time_since_shock) |>
  group_by(feedback, chain) |>
  mutate(m = roll_mean(m)) |>
  ungroup() |>
  filter(!is.na(m))
# ^ raw m is a per-sweep discrete symptom COUNT (0-9 active symptoms,
# resampled every step), not a smooth continuous process like raw P --
# at the single-sweep level it is high-frequency noise, so 30 overlaid
# raw chains rendered a dense vertical hash instead of legible
# trajectories. Applying the same rolling-window smoothing used for the
# aggregate mean to each INDIVIDUAL chain first fixes this: each line
# still shows genuine chain-to-chain variability (it is not the same
# curve repeated), just without the per-sweep sampling noise dominating
# the visual.

a_end <- a_summary |> filter(!is.na(mean_m_smooth)) |> group_by(feedback) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup() |>
  arrange(mean_m_smooth)

# The off/on end labels collide when the two curves end up close together
# (here: off ~0.29, on ~0.32) -- nudge the label Y positions apart a
# little without moving the actual lines, same fix-class as Figure 2's
# peak-label collision earlier.
if (nrow(a_end) == 2 && diff(range(a_end$mean_m_smooth)) < 0.05) {
  a_end$label_y <- a_end$mean_m_smooth + c(-0.028, 0.028)
} else {
  a_end$label_y <- a_end$mean_m_smooth
}

feedback_labels <- c(off = "Feedback off", on = "Feedback on")
x_max_a <- plot_window_a[2]

# REVERTED (2026-08-25): x-axis-matching to panel B's 0-1500 window was
# tried, but with panel B now genuinely covering the b=0/0.50 contrast
# itself, panel A's line visibly stopping around step 750 while the axis
# ran to 1500 read as broken/incomplete rather than "intentionally
# focused." Panel A is the zoomed, precise main-simulation comparison --
# it keeps its own native window (0-750) rather than matching panel B.
x_axis_max_a <- x_max_a

pA <- ggplot(a_summary, aes(x = time_since_shock, y = mean_m_smooth, colour = feedback, linetype = feedback)) +
  geom_hline(yintercept = m_ref_a, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = a_chains,
    aes(x = time_since_shock, y = m, group = interaction(chain, feedback), colour = feedback),
    inherit.aes = FALSE, alpha = spaghetti_alpha, linewidth = spaghetti_width
  ) +
  geom_ribbon(
    aes(ymin = mean_m_smooth - interval_mult * se_m_smooth, ymax = mean_m_smooth + interval_mult * se_m_smooth, fill = feedback),
    alpha = ribbon_alpha, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = main_line_width) +
  geom_text(
    data = a_end, aes(x = x_max_a + 25, y = label_y, label = feedback_labels[as.character(feedback)]),
    hjust = 0, size = 3.0, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(xlim = c(plot_window_a[1], x_axis_max_a), clip = "off") +
  scale_colour_manual(values = pal_feedback, labels = feedback_labels) +
  scale_fill_manual(values = pal_feedback, guide = "none") +
  scale_linetype_manual(values = c(off = "dashed", on = "solid"), labels = feedback_labels) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.16))) +
  labs(
    title = panel_title("A", "Moderate feedback slows recovery"),
    x = "Steps since shock", y = "Mean symptom activation (m)", colour = NULL, linetype = NULL
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    plot.margin = margin(5.5, 60, 5.5, 5.5),
    legend.position = "bottom",
    legend.key.width = unit(1.1, "cm"),
    legend.text = element_text(size = base_sz - 2),
    legend.margin = margin(t = -4)
  )

# ------------------------------------------------------------------------
# Panel B: shock-and-recovery trajectories across the feedback-strength
# grid (03c_sim_feedback_shock_grid_extended.R) -- same shock design and
# visual language as panel A, just with one line per b instead of two
# (off/on), so the reader can see the SAME experiment shift from
# "recovers" to "settles at an elevated plateau" as b increases.
#
# DISPLAY SET revised 2026-08-25: b=0/0.50 dropped -- panel A now covers
# that exact contrast directly, so repeating it here was redundant.
# Instead shows b = 0.75, 0.90, 1.00, 1.10, 1.30: points strictly between
# the established "recovers" and "clear elevated plateau" regimes, plus
# one point further into the history-dependent side. b=1.25/1.50 stay out
# of the DISPLAY set -- b=1.50 in particular is still the least-resolved
# point in the grid. 03c's window was extended (1500 -> 2500 steps) and
# a b=1.30 point was added specifically so this set could include
# something past b=1.00 that is ACTUALLY plateaued, not just labeled as
# if it were -- check the rendered PNG (does the b=1.30 line visibly
# flatten before the window ends?) before trusting that it belongs here.
# If b=1.30 is still rising at step 2500, drop it back out, same as
# b=1.25/1.50 were dropped from the original 1500-step version.
# ------------------------------------------------------------------------
sim3c <- read_csv("res/revision_2026/sim3c/sim3c_summary.csv", show_col_types = FALSE)

# b=1.10 was checked against the rerun (2500-step window) and DROPPED:
# unlike b=1.30 (which visibly flattens, rising only ~0.03 over the last
# 1000 steps) and b=1.00 (~0.02 residual drift, essentially flat),
# b=1.10 was still rising by ~0.08 over the same span with no sign of
# leveling off -- plausibly critical slowing down near the transition
# itself, but not something we can currently show as "resolved" next to
# genuinely plateaued curves without overclaiming. Left out for the same
# reason 1.25/1.50 were left out of the original 1500-step version.
display_b_grid <- c(0.75, 0.90, 1.00, 1.30)

plot_window_grid <- c(-100, 2500)
grid_w <- sim3c |> filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2]) |>
  arrange(b, time_since_shock) |>
  group_by(b) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N)) |>
  ungroup()

b_levels <- sort(unique(grid_w$b))
grid_pal <- setNames(colorRampPalette(c(col_off, col_on, col_high))(length(b_levels)), b_levels)
grid_labels <- setNames(sprintf("b = %.2f", b_levels), b_levels)

# Individual-chain spaghetti, same per-chain-smoothing fix used for panel
# A (raw m is a per-sweep discrete count -- smooth EACH chain first, not
# just the aggregate mean, or it renders as noise-hash rather than
# legible trajectories). display_b_grid has 4 lines -- enough visual room
# for thinned-down individual chains without them turning into a block.
set.seed(11)
# Thinned from 25 to 10 chains/b per review -- with panel A now carrying
# the precise, focused off/on comparison, panel B's job is to show the
# broader shape of the transition across b, not to re-carry the same
# individual-trajectory detail. Fewer, fainter background chains keep
# the multi-line panel from getting visually noisy.
n_sample_chains_grid <- 10L
sim3c_raw <- readRDS("res/revision_2026/sim3c/sim3c_raw.rds")$traj |>
  filter(b %in% display_b_grid, time_since_shock >= plot_window_grid[1], time_since_shock <= plot_window_grid[2])

grid_chains <- sim3c_raw |>
  group_by(b) |>
  filter(chain %in% sample(unique(chain), n_sample_chains_grid)) |>
  ungroup() |>
  arrange(b, chain, time_since_shock) |>
  group_by(b, chain) |>
  mutate(m = roll_mean(m)) |>
  ungroup() |>
  filter(!is.na(m))

# grid_w no longer contains b=0 (dropped from the display set, see above),
# so the old "pre-shock baseline from the b=0 rows" trick would silently
# return NaN here. Reuse panel A's m_ref_a instead -- same pre-shock
# baseline concept (feedback doesn't change pre-shock dynamics, so it
# doesn't matter which b the reference line is computed from).
m_ref_grid <- m_ref_a
x_max_grid <- plot_window_grid[2]

grid_end <- grid_w |> filter(!is.na(mean_m_smooth)) |> group_by(b) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup() |>
  arrange(mean_m_smooth)
# Spread out end-labels that would otherwise collide (several b's plateau
# close together) -- same fix-class as panel A's off/on label collision.
label_gap <- 0.045
if (nrow(grid_end) > 1) {
  y <- grid_end$mean_m_smooth
  for (i in 2:length(y)) {
    if (y[i] - y[i - 1] < label_gap) y[i] <- y[i - 1] + label_gap
  }
  grid_end$label_y <- y
} else {
  grid_end$label_y <- grid_end$mean_m_smooth
}

pGrid <- ggplot(grid_w, aes(x = time_since_shock, y = mean_m_smooth, colour = factor(b), group = factor(b))) +
  geom_hline(yintercept = m_ref_grid, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = grid_chains,
    aes(x = time_since_shock, y = m, group = interaction(chain, b), colour = factor(b)),
    inherit.aes = FALSE, alpha = spaghetti_alpha, linewidth = spaghetti_width
  ) +
  geom_ribbon(
    aes(ymin = mean_m_smooth - interval_mult * se_m_smooth, ymax = mean_m_smooth + interval_mult * se_m_smooth, fill = factor(b)),
    alpha = ribbon_alpha, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = main_line_width) +
  geom_text(
    data = grid_end, aes(x = x_max_grid + 30, y = label_y, label = grid_labels[as.character(b)], colour = factor(b)),
    hjust = 0, size = 2.8, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(xlim = c(plot_window_grid[1], x_max_grid), clip = "off") +
  scale_colour_manual(values = grid_pal, labels = grid_labels, name = NULL) +
  scale_fill_manual(values = grid_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.11))) +
  labs(
    title = panel_title("B", "Recovery changes across feedback strengths"),
    x = "Steps since shock", y = "Mean symptom activation (m)"
  ) +
  theme_pub(base_size = base_sz) +
  theme(
    plot.margin = margin(5.5, 62, 5.5, 5.5),
    legend.position = "bottom",
    legend.key.width = unit(0.9, "cm"),
    legend.text = element_text(size = base_sz - 3)
  )

# ------------------------------------------------------------------------
# Panel C: history-dependence index across the full feedback-strength grid
# ------------------------------------------------------------------------
sep <- read_csv("res/revision_2026/supp_history/history_separation_summary.csv", show_col_types = FALSE)

regime_pal <- c(convergent = col_off, `history-dependent` = col_high, `runaway/saturation` = "black")

# Locate the locked main-text b (0.50) by nearest value, not exact
# equality -- the grid is now generated via seq(0, 1.5, by = 0.1), and
# exact `== 0.5` comparisons against a seq()-built float are a real risk
# of silently matching nothing due to floating-point drift.
locked_b_row <- sep[which.min(abs(sep$b - 0.5)), ]

pC <- ggplot(sep, aes(x = b, y = separation_m)) +
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
  annotate("text", x = locked_b_row$b, y = locked_b_row$separation_m + 0.06, label = "main simulation\nvalue", size = 2.7, colour = col_on, hjust = 0.5, vjust = 0, lineheight = 0.85) +
  scale_colour_manual(values = regime_pal, name = NULL) +
  labs(
    title = panel_title("C", "History-dependence emerges at higher feedback"),
    x = "Feedback strength (b)", y = expression(Delta ~ "mean symptom activation (late)")
  ) +
  theme_pub(base_size = base_sz) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------
# Combine + save -- A (locked b=0.50 recovery, precise/high-chain-count),
# B (same shock design swept across the full b-grid, shared visual
# language with A), C (the regime-index summary tying the grid together).
# Each panel keeps its OWN local legend rather than collecting into one
# shared legend -- A/B/C's colour scales are three different variables,
# not levels of one; collecting them earlier produced one crammed,
# clipped row that also read as though unrelated categories (e.g.
# "Feedback off" and "convergent") were comparable.
# ------------------------------------------------------------------------
fig3 <- (pA / pGrid / pC) +
  plot_layout(heights = c(1, 1, 1.05))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_regimes.pdf", fig3, width = 7.8, height = 10.6)
ggsave("figs/revision_2026/Figure3_regimes.png", fig3, width = 7.8, height = 10.6, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_regimes.pdf (+ .png)\n")
cat("\nThis REPLACES fig3_recovery_feedback.R / FigureS_history_dependence.R\n")
cat("as main-text Figure 3. Those two scripts are kept (not deleted) but are\n")
cat("now superseded for main-text purposes.\n")
