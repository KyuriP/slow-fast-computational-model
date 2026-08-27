# ============================================================
# R/revision_2026/figures/fig3_perturbation_recovery.R
# ============================================================
# Figure for Simulation 2 (2026-08-27, replaces fig3_stress_recovery.R):
# a perturbation to the slow field is described entirely in prose in the
# manuscript, with no accompanying figure -- the numbers (P_t jumps to
# ~1.00, symptom activation peaks at ~4.73, both recover) read as
# "floating" without a time-course panel. This script builds that figure:
# mean P_t and mean number of active symptoms over time, feedback off
# (b=0, matching Simulation 2's design), each as a bold mean trajectory
# with a light 95% ribbon and a SMALL sample of individual chains.
#
# Renamed from fig3_stress_recovery.R -- Simulation 2 is not modeling
# stressor occurrence, it's a one-time perturbation to P_t, so
# "perturbation recovery" matches the manuscript's own framing
# ("perturbation", "perturbation onset") rather than "stress".
#
# 2026-08-27, three rounds of revision on panel B specifically (panel A's
# mean + light ribbon + small chain sample worked from the start and is
# unchanged below):
#   v1: ALL 200 chains shown as spaghetti, matching
#       fig3_recovery_feedback.R's original panel B. Individual symptom
#       counts are an integer 0-9 series -- 200 overlapping lines
#       rendered as a solid orange block that obscured the mean recovery
#       curve entirely.
#   v2 (overcorrection): spaghetti removed, mean + ribbon only. Read as
#       too sterile/flat.
#   v3: a small sample of 12 chains, faint and thin. Still didn't work --
#       discrete integer-valued step trajectories don't read as "faint
#       texture" the way continuous P_t trajectories do; even 12 of them
#       looked like a barcode rather than noise around a trend.
#   v4 (this version): panel B is now a heatmap over (time, symptom
#       count) with a bold mean line overlaid, not a spaghetti/ribbon
#       plot at all. This fits the data better than any line-based
#       approach: symptom count is discrete (0-9) and the quantity of
#       real interest is how the FULL distribution across chains shifts
#       and reshapes after the perturbation, not individual chain paths.
#       Color intensity = proportion of the 200 chains at each symptom
#       count, within a time bin. Panel A keeps the small-sample-
#       spaghetti approach, since P_t is continuous and that version
#       already read cleanly.
#
# Uses raw symptom counts (mean_M), not the mean_m fraction, so the
# figure's y-axis matches the units actually cited in the Results prose
# ("an average of 2.48 active symptoms", "peak of about 4.73 active
# symptoms") -- no unit conversion for the reader to do in their head.
#
# Outputs
# -------
#   figs/revision_2026/Figure3_perturbation_recovery.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_P, col_M

# ------------------------------------------------------------------------
# 0. Load data (Simulation 2 only -- feedback off, b=0)
# ------------------------------------------------------------------------
sim2 <- read_csv("res/revision_2026/sim2/sim2_summary.csv", show_col_types = FALSE)

plot_window <- c(-100, 750)
sim2_w <- sim2 |> filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])

base_sz <- 11

# ------------------------------------------------------------------------
# 0a. Rolling-window smoothing for the symptom-count panel (display only,
#     same as fig3_recovery_feedback.R). Two-sided moving average; edges
#     become NA and are silently dropped by geom_line.
# ------------------------------------------------------------------------
roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))

sim2_w <- sim2_w |>
  arrange(time_since_shock) |>
  mutate(mean_M_smooth = roll_mean(mean_M), se_M_smooth = roll_mean(se_M))

# ------------------------------------------------------------------------
# 0b. A SMALL sample of individual-chain trajectories (see header note --
#     12, not all 200) for a faint spaghetti layer behind the mean in
#     both panels. Individual-chain M is NOT rolling-smoothed -- the raw
#     per-chain step noise is part of the point of showing a few of them
#     (visible texture, not another smoothed summary curve).
# ------------------------------------------------------------------------
set.seed(3)
n_sample_chains <- 25L

sim2_raw <- readRDS("res/revision_2026/sim2/sim2_raw.rds")$traj |>
  filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])
sim2_chains_w <- sim2_raw |> filter(chain %in% sample(unique(sim2_raw$chain), n_sample_chains))

# ------------------------------------------------------------------------
# 0c. Panel B heatmap data: bin time into the same effective window used
#     for the rolling-mean smoothing (k=15, see 0a) so the heatmap's time
#     resolution matches the mean line's, then compute the proportion of
#     all 200 chains at each symptom count (0-9) within each time bin.
#     Uses the FULL sim2_raw (all chains), not the 12-chain sample --
#     the heatmap needs the full distribution to estimate proportions
#     sensibly, unlike the spaghetti layer which deliberately only shows
#     a few individual paths.
# ------------------------------------------------------------------------
time_bin_width <- 15

heat_df <- sim2_raw |>
  mutate(time_bin = round(time_since_shock / time_bin_width) * time_bin_width) |>
  count(time_bin, M, name = "n_chains") |>
  group_by(time_bin) |>
  mutate(prop = n_chains / sum(n_chains)) |>
  ungroup()

# ------------------------------------------------------------------------
# 1. y-limits and interval settings
# ------------------------------------------------------------------------
interval_mult <- 1.96  # 95% intervals

range_pad <- function(x, pad = 0.06, lower = -Inf, upper = Inf) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (d == 0) d <- 1
  c(max(lower, r[1] - pad * d), min(upper, r[2] + pad * d))
}

p_ylim <- range_pad(
  c(sim2_w$mean_P - interval_mult * sim2_w$se_P, sim2_w$mean_P + interval_mult * sim2_w$se_P),
  pad = 0.06
)

p_ref <- mean(sim2_w$mean_P[sim2_w$time_since_shock < 0], na.rm = TRUE)
# Note: M_ylim/M_ref (ribbon range + pre-perturbation reference line) are
# no longer needed -- panel B replaced its ribbon+hline design with the
# heatmap below, which uses a fixed 0-9 axis (the full possible range of
# symptom counts) instead.

# ------------------------------------------------------------------------
# 2. Panel A: slow field P_t -- bold mean + light ribbon + a small
#    faint sample of individual chains
# ------------------------------------------------------------------------
pA <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_P)) +
  geom_hline(yintercept = p_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = sim2_chains_w,
    aes(x = time_since_shock, y = P, group = chain),
    inherit.aes = FALSE,
    colour = col_P, alpha = 0.22, linewidth = 0.18
  ) +
  geom_ribbon(
    aes(ymin = mean_P - interval_mult * se_P, ymax = mean_P + interval_mult * se_P),
    fill = col_P, alpha = 0.12
  ) +
  geom_line(colour = col_P, linewidth = 0.55) +
  coord_cartesian(ylim = p_ylim, clip = "off") +
  labs(title = panel_title("A", "Slow context rises and recovers"), x = NULL, y = "Slow field (P)") +
  theme_pub(base_size = base_sz) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ------------------------------------------------------------------------
# 3. Panel B: heatmap of the symptom-count DISTRIBUTION over time, with
#    the bold mean line overlaid. geom_tile is drawn first (its own data/
#    aes, inherit.aes=FALSE) so it sits behind the mean line. The mean
#    line is drawn twice -- a thicker white "halo" first, then the actual
#    dark line on top -- so it stays legible against both the light
#    (low-proportion) and saturated (high-proportion) ends of the fill
#    scale, rather than picking one line colour that only contrasts with
#    part of the heatmap.
# ------------------------------------------------------------------------
pB <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_M_smooth)) +
  geom_tile(
    data = heat_df,
    aes(x = time_bin, y = M, fill = prop),
    inherit.aes = FALSE,
    width = time_bin_width, height = 1
  ) +
  scale_fill_gradient(
    low = "white", high = col_M, name = "Proportion of chains",
    labels = scales::percent_format(accuracy = 1),
    breaks = scales::breaks_width(0.1),
    guide = guide_colorbar(
      barwidth = unit(5.5, "cm"), barheight = unit(0.42, "cm"),
      title.position = "top", title.hjust = 0.5
    )
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
  geom_line(colour = "white", linewidth = 1.5, lineend = "round") +
  geom_line(colour = "grey15", linewidth = 0.8) +
  scale_y_continuous(breaks = 0:9, limits = c(-0.5, 9.5), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = panel_title("B", "Symptom activation shifts upward and recovers"),
       x = "Steps since perturbation", y = "Number of active symptoms") +
  theme_pub(base_size = base_sz)

# ------------------------------------------------------------------------
# 4. Combine + save -- vertical stack, shared x-axis (perturbation-onset
#    line lines up visually between panels). Heights slightly unequal
#    (0.9 : 1.1) -- panel A alone looked a touch too tall relative to B
#    once B carries the heatmap + its own legend row underneath.
# ------------------------------------------------------------------------
fig_perturbation_recovery <- pA / pB + plot_layout(heights = c(0.9, 1.1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_perturbation_recovery.pdf", fig_perturbation_recovery, width = 6.4, height = 6.4)
ggsave("figs/revision_2026/Figure3_perturbation_recovery.png", fig_perturbation_recovery, width = 6.4, height = 6.4, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_perturbation_recovery.pdf (+ .png)\n")
cat("\nNOTE: this is a NEW figure, inserted before the existing Figure 3\n")
cat("(regimes, fig3_regimes.R / Figure3_regimes.pdf) in document order --\n")
cat("that figure and Figure 4 (network trio) will renumber automatically\n")
cat("via LaTeX's \\ref, no manuscript label changes needed beyond adding\n")
cat("this new figure environment. See chat for the exact LaTeX block and\n")
cat("Results paragraph to insert in the Simulation 2 section.\n")
