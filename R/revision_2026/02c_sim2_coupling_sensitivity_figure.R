# ============================================================
# R/revision_2026/02c_sim2_coupling_sensitivity_figure.R
# ============================================================
# Re-plots the coupling-sensitivity recovery figure from the raw
# per-chain output already saved by 02b_sim_stress_recovery_coupling_
# sensitivity.R -- no new simulation, just a cleaner final figure for
# the appendix. Requires that 02b has been run (with c=0.3 in its
# c_grid -- rerun 02b first if sim2_coupling_sensitivity_raw.rds
# predates that).
#
# Changes from 02b's own figure:
#   1. Three conditions only: c in {0.3, 1.0 (reference), 3.0}. The
#      full six/seven-point grid showed recovery is flat from c=0.3 to
#      c=2.0 and only pulls away at c=2.5-3.0 -- three points (a low
#      end, the reference, and the clearly-diverged high end) tell
#      that story without the clutter of the intermediate, redundant
#      points.
#   2. Individual chains as light spaghetti under a bold mean line,
#      matching the Simulation 1 burn-in trace figure. Each chain is
#      smoothed with a short rolling mean (spag_smooth_window steps)
#      before computing its own R(t) -- raw per-step chain values are
#      as noisy as the whole condition's peak-minus-baseline window,
#      so without smoothing the normalization blows that noise up into
#      a solid haze instead of a light texture.
#   3. Colour: fixed, fully-saturated categorical blue/grey/red rather
#      than a distance-scaled diverging ramp. With only 3 categories,
#      scaling saturation by how far c is from 1 made asymmetric grids
#      (like 0.3 vs 3, which are not equidistant from 1) look washed
#      out on the closer side -- fixed full-strength colours at each
#      rank position (lowest / reference / highest) read clearly
#      regardless of the actual numeric spacing.
#
# Requires: res/revision_2026/sim2/sim2_coupling_sensitivity_raw.rds
# Output:   figs/revision_2026/fig_sim2_coupling_sensitivity.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub()

raw <- readRDS("res/revision_2026/sim2/sim2_coupling_sensitivity_raw.rds")
traj_all <- raw$traj
post_shock_steps <- raw$params$post_shock_steps

c_subset <- c(0.3, 1.0, 3.0)
missing_c <- setdiff(c_subset, unique(traj_all$c_coupling))
if (length(missing_c) > 0) {
  stop(sprintf(
    "c = %s not found in the saved raw data. Rerun 02b_sim_stress_recovery_coupling_sensitivity.R first (its c_grid now includes 0.3) before sourcing this script.",
    paste(missing_c, collapse = ", ")
  ))
}

traj <- traj_all |> filter(c_coupling %in% c_subset)

summary_tbl <- traj |>
  group_by(c_coupling, step, time_since_shock) |>
  summarise(mean_M = mean(M), .groups = "drop")

baseline_tbl <- summary_tbl |>
  filter(time_since_shock < 0, time_since_shock >= -50) |>
  group_by(c_coupling) |>
  summarise(mean_M_pre = mean(mean_M), .groups = "drop")

peak_tbl <- summary_tbl |>
  filter(time_since_shock >= 0, time_since_shock < 20) |>
  group_by(c_coupling) |>
  summarise(mean_M_peak = max(mean_M), .groups = "drop")

recovery_mean <- summary_tbl |>
  left_join(baseline_tbl, by = "c_coupling") |>
  left_join(peak_tbl, by = "c_coupling") |>
  mutate(R = (mean_M - mean_M_pre) / (mean_M_peak - mean_M_pre))

# A light, reproducible sample of individual chains per condition --
# chains are iid by construction, so the first n of them is as good a
# sample as a random draw, and this keeps the figure deterministic
# without introducing a new seed.
n_spag_chains <- 15L
spag_smooth_window <- 21L  # centered rolling mean, in sweeps

recovery_spag <- traj |>
  filter(chain <= n_spag_chains) |>
  arrange(c_coupling, chain, step) |>
  group_by(c_coupling, chain) |>
  mutate(M_smooth = as.numeric(stats::filter(M, rep(1 / spag_smooth_window, spag_smooth_window), sides = 2))) |>
  ungroup() |>
  filter(!is.na(M_smooth)) |>
  left_join(baseline_tbl, by = "c_coupling") |>
  left_join(peak_tbl, by = "c_coupling") |>
  mutate(R = (M_smooth - mean_M_pre) / (mean_M_peak - mean_M_pre))

c_lab_str <- sprintf("c = %.2g%s", c_subset, ifelse(c_subset == 1, " (reference)", ""))
c_labels <- setNames(c_lab_str, as.character(c_subset))

# --- fixed, fully-saturated categorical colour by rank, not magnitude ---
below <- sort(c_subset[c_subset < 1], decreasing = TRUE)  # nearest-to-1 first
above <- sort(c_subset[c_subset > 1])                      # nearest-to-1 first
pos_of <- c(setNames(-seq_along(below) / max(1, length(below)), as.character(below)),
            setNames(seq_along(above) / max(1, length(above)), as.character(above)))
if (1 %in% c_subset) pos_of <- c(pos_of, setNames(0, "1"))

diverging_ramp <- grDevices::colorRampPalette(c("#6166AC", "grey50", "#CF746C"))(201)
get_col <- function(c_val) {
  pos <- pos_of[[as.character(c_val)]]
  idx <- round((pos + 1) / 2 * 200) + 1
  diverging_ramp[idx]
}
pal_coupling <- setNames(sapply(c_subset, get_col), c_lab_str)

recovery_mean <- recovery_mean |> mutate(c_label = factor(c_labels[as.character(c_coupling)], levels = c_lab_str))
recovery_spag <- recovery_spag |> mutate(c_label = factor(c_labels[as.character(c_coupling)], levels = c_lab_str))

spag_width <- 0.14
spag_alpha <- 0.10
mean_width <- 0.25  # narrower than the usual main_line_width, so close conditions stay distinguishable
mean_alpha <- 1  # slightly transparent so overlapping mean lines blend instead of fully occluding each other

p_recovery <- ggplot(
  recovery_mean |> filter(time_since_shock >= -20, time_since_shock <= post_shock_steps),
  aes(time_since_shock, R, colour = c_label)
) +
  geom_hline(yintercept = c(0, 1), linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_line(data = recovery_spag |> filter(time_since_shock >= -20, time_since_shock <= post_shock_steps),
            aes(group = interaction(c_label, chain)),
            linewidth = spag_width, alpha = spag_alpha) +
  geom_line(linewidth = mean_width, alpha = mean_alpha) +
  scale_colour_manual(values = pal_coupling, name = NULL) +
  coord_cartesian(ylim = c(-0.3, 1.3)) +
  labs(x = "Steps since shock", y = "R(t)  (normalized excess activation)",
       title = "Sensitivity of recovery to symptom-coupling strength",
       subtitle = expression(omega^{(c)} == c %.% omega * ",  topology fixed -- thin lines are individual chains (smoothed), bold = mean of 1000")) +
  theme_pub()

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/fig_sim2_coupling_sensitivity.pdf", p_recovery, width = 8, height = 5.5)

cat("Done: figs/revision_2026/fig_sim2_coupling_sensitivity.pdf\n")
cat(sprintf("Conditions plotted: %s\n", paste(c_subset, collapse = ", ")))
cat(sprintf("Spaghetti: %d chains per condition, smoothed over %d steps, alpha=%.2f, linewidth=%.2f\n",
            n_spag_chains, spag_smooth_window, spag_alpha, spag_width))
cat(sprintf("Mean line width: %.2f\n", mean_width))
