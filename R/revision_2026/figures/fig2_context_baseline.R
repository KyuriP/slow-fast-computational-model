# ============================================================
# R/revision_2026/figures/fig2_context_baseline.R
# ============================================================
# Figure 2: same symptom-symptom coupling, different slow context,
# different symptom burden.
#
# REBUILT as a two-panel figure. The 3-panel version (schematic + PMF +
# summary) was cut for the same reason as Figure 4's original summary
# panel: it tried to do too many jobs (illustrate the fixed network,
# show the distribution, show the summary) and read as an "analysis
# report panel" rather than a single clean visual claim. The fixed
# network is now stated in the caption/Methods instead of drawn as a
# schematic -- a placeholder 7-node cartoon risked implying a wrong N
# anyway (flagged in an earlier round), so dropping it removes both the
# clutter and that accuracy risk at once.
#
#   (A) Burden distributions, wider, no legend -- lines labeled directly
#       at their right end (Low / Middle / High context).
#   (B) Compact summary by context, narrower and tighter than the
#       previous version (labels were drifting too far from the points).
#
# Column-name notes (this repo's actual structure):
#   - sim1_raw.rds is a LIST (burden / symptom_activation / params /
#     design), not a flat data frame -- the burden table is at $burden.
#   - burden columns are: condition, P, chain, sweep, M, m (m = M/N).
#   - There is no valid precomputed SE in sim1_summary.csv (sd_M there is
#     sweep-level, and sequential Gibbs sweeps within a chain are
#     autocorrelated -- treating all 200 chains x 200 sweeps as
#     independent understates the true SE). Panel B recomputes SE
#     properly from per-CHAIN means (200 independent chains).
#
# Outputs
# -------
#   figs/revision_2026/Figure2_context_baseline.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_*, pal_context

# ------------------------------------------------------------------------
# 0. Load data
# ------------------------------------------------------------------------
raw_list <- readRDS("res/revision_2026/sim1/sim1_raw.rds")
burden   <- raw_list$burden                              # condition, P, chain, sweep, M, m

# ------------------------------------------------------------------------
# 1. Three simulated conditions, labeled for display
# ------------------------------------------------------------------------
condition_labels <- c(low = "Low context", middle = "Middle context", high = "High context")

burden_show <- burden |>
  mutate(P_group = factor(condition_labels[as.character(condition)],
                           levels = unname(condition_labels[c("low", "middle", "high")])))

# Named palette keyed to the DISPLAY labels (pal_context from
# theme_publication.R is keyed to low/middle/high; ggplot's scale needs
# the labels actually used in the data).
pal_context_display <- setNames(pal_context, condition_labels[names(pal_context)])

# ------------------------------------------------------------------------
# 2. Panel A: burden distributions (probability-mass line/point plot --
#    M is a discrete 0-9 count, not continuous, so this stays a PMF plot
#    rather than a density). No legend -- each line is labeled directly,
#    so the reader doesn't have to cross-reference a legend key.
#
#    Labels are placed at each curve's own PEAK, not at a shared x
#    position (e.g. the right end, M=9) -- all three distributions
#    converge to ~0% by M=9, so a right-edge label collided all three
#    text strings on top of each other there. The peaks sit at different
#    M values (low~1, middle~2, high~3-4) and are naturally well
#    separated vertically, so labeling there is both readable and a more
#    informative anchor than the flat right tail.
# ------------------------------------------------------------------------
prop_tbl <- burden_show |>
  count(P_group, M, name = "n") |>
  group_by(P_group) |>
  mutate(prop = n / sum(n)) |>
  ungroup() |>
  complete(P_group, M = 0:9, fill = list(n = 0, prop = 0))

peak_labels <- prop_tbl |>
  group_by(P_group) |>
  slice_max(prop, n = 1, with_ties = FALSE) |>
  ungroup()

pA <- ggplot(prop_tbl, aes(x = M, y = prop, colour = P_group, group = P_group)) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 2.0) +
  # geom_label (not geom_text): the three curves cross each other near
  # each peak, so plain text sat directly on top of a differently-colored
  # line behind it and read as "crossed with the graph." A borderless,
  # semi-opaque white label punches a clean halo behind the text so it
  # stays legible regardless of what's plotted underneath.
  geom_label(
    data = peak_labels,
    aes(x = M, y = prop, label = P_group),
    vjust = -0.65, size = 3.1, fontface = "plain", show.legend = FALSE,
    label.size = 0, label.padding = unit(0.12, "lines"),
    fill = scales::alpha("white", 0.82)
  ) +
  scale_colour_manual(values = pal_context_display, guide = "none") +
  scale_x_continuous(breaks = 0:9, expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0.02, 0.11))) +
  labs(title = panel_title("A", "Context level shifts symptom activation"),
       x = "Number of active symptoms (M)", y = "Probability") +
  theme_pub(base_size = 10.5)

# ------------------------------------------------------------------------
# 3. Panel B: compact summary contrast. Mean burden as point + 95%
#    interval (per-chain SE, not sweep-level), high-burden probability as
#    a short text label. Tightened relative to the earlier version: less
#    horizontal travel between the point, its mean label, and the Pr(M>=5)
#    label, so the panel reads as one compact strip instead of a row of
#    scattered text.
# ------------------------------------------------------------------------
chain_summary <- burden_show |>
  group_by(condition, P, P_group, chain) |>
  summarise(
    chain_mean_M = mean(M),
    chain_pr_high = mean(M >= 5),
    .groups = "drop"
  )

summ_corrected <- chain_summary |>
  group_by(condition, P, P_group) |>
  summarise(
    mean_M = mean(chain_mean_M),
    se_M = sd(chain_mean_M) / sqrt(n()),
    pr_high = mean(chain_pr_high),
    se_pr_high = sd(chain_pr_high) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    P_group = factor(
      P_group,
      levels = rev(unname(condition_labels[c("low", "middle", "high")]))
    ),
    mean_label = sprintf("%.2f", mean_M),
    high_label = sprintf("5+ active: %.0f%%", 100 * pr_high)
  )

# CORRECTED after checking the rendered PNG (2026-08-25): the fixed
# label_x=3.95 / x_axis_max=4.55 tried above clipped the label text mid-
# word ("5+ activ") on every row, and on High context (mean_M=3.56) the
# label collided directly with the "3.56" mean-value text -- fixed limits
# don't adapt to the actual data, and scale_x_continuous(limits=...)
# hard-clips anything (including text) past the edge, so a too-narrow
# guess doesn't just look cramped, it silently truncates the label.
# Back to a DATA-DRIVEN position/limit, like the earlier working version,
# but keeping the new shorter "5+ active: %" text:
#   - label_x sits a fixed gap past the widest mean_label, so it can never
#     collide with a mean_M value no matter which condition is largest
#   - x_axis_max leaves generous room (not just enough for the shortest
#     label) for the longest label text ("5+ active: 100%" case)
# Check the rendered PNG again after rerunning -- if it still clips,
# widen the multiplier below rather than shrinking the text.
label_x <- max(summ_corrected$mean_M) + 0.95
x_axis_max <- label_x + 2.3
sep_x   <- max(summ_corrected$mean_M) + 0.62   # light divider between the
                                                # mean-value column and the
                                                # Pr(5+ active) column

# The x-SCALE genuinely runs to x_axis_max (kept, proven not to clip the
# label text -- see the 2026-08-25 fix note above), but only breaks/ticks
# 0:4 are drawn, so the space past 4 was unlabeled and read as "why does
# the axis just stop at 4." Rather than change the scale mechanics again,
# the fix here is to stop leaving that region unexplained: a light
# separator plus explicit "Pr(5+ active)" column header now marks it as a
# deliberate second field, not a truncated axis.
n_grp <- nlevels(summ_corrected$P_group)

pB <- ggplot(summ_corrected, aes(y = P_group, x = mean_M, colour = P_group)) +
  geom_segment(
    aes(x = 0, xend = mean_M, y = P_group, yend = P_group),
    linewidth = 0.9,
    alpha = 0.28
  ) +
  geom_segment(
    aes(
      x = mean_M - 1.96 * se_M,
      xend = mean_M + 1.96 * se_M,
      y = P_group,
      yend = P_group
    ),
    linewidth = 0.75
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(x = mean_M + 0.10, label = mean_label),
    hjust = 0,
    size = 3.0,
    colour = "black"
  ) +
  geom_vline(xintercept = sep_x, colour = "grey85", linewidth = 0.4) +
  geom_text(
    aes(x = label_x, label = high_label),
    hjust = 0,
    size = 3.0,
    colour = "grey35"
  ) +
  # Column headers, drawn once (not per-row) just above the top data row --
  # these are what actually answer "what does panel B mean": the left
  # column is the mean +/- 95% CI, the right column is a separate quantity
  # (probability of high activation), not a continuation of the same axis.
  annotate("text", x = 0.15, y = n_grp + 0.62,
           label = "Mean (95% CI)", hjust = 0, size = 2.7, colour = "grey35",
           fontface = "italic") +
  annotate("text", x = label_x, y = n_grp + 0.62,
           label = "Pr(5+ active)", hjust = 0, size = 2.7, colour = "grey35",
           fontface = "italic") +
  scale_colour_manual(values = pal_context_display, guide = "none") +
  scale_x_continuous(
    limits = c(0, x_axis_max),
    breaks = 0:4,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.0))) +
  labs(
    title = panel_title("B", "Mean activation and probability of high burden"),
    x = "Mean number of active symptoms",
    y = NULL
  ) +
  theme_pub(base_size = 10.5) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10, colour = "black"),
    plot.margin = margin(14, 10, 8, 8)
  )

# ------------------------------------------------------------------------
# 4. Combine + save -- A wider (it's the main visual claim), B narrower
#    and treated as compact numeric support, not a co-equal panel.
# ------------------------------------------------------------------------
fig2 <- pA + pB + plot_layout(widths = c(1.5, 1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/Figure2_context_baseline.pdf", fig2, width = 9.6, height = 4.1)
ggsave("figs/revision_2026/Figure2_context_baseline.png", fig2, width = 9.6, height = 4.1, dpi = 300)

cat("Note: SEs (mean_M and pr_high) are computed from per-chain summaries\n")
cat("(n=200 independent chains per condition), not from sweep-level samples.\n\n")
print(summ_corrected)

cat("\nDone. Files:\n")
cat("  figs/revision_2026/Figure2_context_baseline.pdf (+ .png)\n")
