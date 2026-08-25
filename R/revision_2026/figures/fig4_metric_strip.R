# ============================================================
# R/revision_2026/figures/fig4_metric_strip.R
# ============================================================
# Compact quantitative strip for Figure 4.
# Use with the network trio:
#   A true network | B symptom-only estimate | C context-adjusted estimate
#   D excess global strength | E error on truly absent edges
#
# Replaces fig4_summary_panels.R as the quantitative half of Figure 4.
# That earlier version tried to do three jobs in one figure (design
# schematic + summary stats + replicate-level robustness check), which is
# why it read as sparse/unfocused no matter how the individual panels
# were polished. This version keeps only the two panels that belong in
# the main text:
#   - the design-schematic panel (old panel A) moves to the caption/
#     Methods instead of taking up figure space
#   - the paired-replicate jitter plot (old panel D) moves to supplement
#     as a robustness check, not main-figure content
#   - the remaining "excess global strength over true" panel replaces the
#     old absolute-scale panel B, which was forced to a huge y-range by
#     including the "true" reference as a separate value on the same
#     axis rather than differencing against it directly
#
# Meant to sit alongside fig4_network_trio.R's output as Figure 4's full
# panel set (network trio on top, this strip on the bottom), assembled in
# LaTeX/Overleaf rather than composited in R -- qgraph's base-R graphics
# and this file's ggplot objects aren't directly combinable via patchwork
# without extra wiring (ggplotify::as.ggplot() / wrap_elements()), which
# remains deferred, lower-priority work.
#
# Outputs
# -------
#   figs/revision_2026/Figure4_metric_strip.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")   # theme_pub(), panel_title(), pal_arm

by_arm <- read.csv("res/revision_2026/sim4/sim4_replicated_by_arm.csv") |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")))

se <- function(x) sd(x) / sqrt(length(x))

# ------------------------------------------------------------------------
# 2026-08-25 readability pass, revised: this strip is now placed at the
# SAME LaTeX width as the network trio (\textwidth, not 0.82\textwidth as
# before) so panels D/E print noticeably larger and their error bars are
# easier to read at a glance -- per review, page-width match was preferred
# over a narrower, taller-relative-to-width strip. With both figures now
# sharing the same width-based scale factor (approx 6.5/9.2 ~= 0.71, using
# the trio's current 9.2in native width as the reference), the earlier
# 13.3pt title override (calibrated for the narrower 0.67x placement) would
# now overshoot and print LARGER than the trio's titles. Sizes below are
# recalibrated down for the new, more favorable scale factor while still
# landing at a comparable final size to the trio (~9pt titles).
# Canvas native size bumped slightly (8.0x2.6 -> 8.4x2.75) so that, once
# both figures share the same print width, this strip's height comes out
# close to but slightly less than the trio's (native aspect ratio
# 2.75/8.4 ~= 0.327 vs. the trio's current 3.1/9.2 ~= 0.337).
title_size_strip <- 12.0
axis_title_strip <- 11.2
axis_text_strip  <- 10.2
label_size_strip <- 3.4   # geom_text size is in mm, not pt

# Matches the locked Sim 4 results text (10 edges, weights 0.20-0.45);
# not recomputed here to keep this script's only input the replicated CSV.
true_gs <- 3.05

arm_labels <- c(
  baseline = "Fixed P",
  naive = "P omitted",
  adjusted = "P included"
)

summary_by_arm <- by_arm |>
  group_by(arm) |>
  summarise(
    gs_excess_mean = mean(global_strength_est - true_gs),
    gs_excess_se   = se(global_strength_est - true_gs),
    mtz_mean       = mean(mae_true_zero),
    mtz_se         = se(mae_true_zero),
    .groups = "drop"
  )

# ------------------------------------------------------------------------
# Panel D: excess global strength over true value
#
# NOTE: y-axis limits/breaks below (both panels) are FIXED per the
# design-pass spec rather than computed from summary_by_arm, to remove
# dead space above/below the points. If a future data refresh shifts
# gs_excess_mean/mtz_mean meaningfully, check the rendered PNG for
# clipped points/error bars before trusting these fixed ranges again.
# ------------------------------------------------------------------------
pD <- ggplot(summary_by_arm, aes(x = arm, y = gs_excess_mean, colour = arm)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.45) +
  annotate("text", x = 0.58, y = 0, label = "No inflation", size = 2.9,
           colour = "grey45", hjust = 0, vjust = -0.6) +
  geom_pointrange(
    aes(ymin = gs_excess_mean - gs_excess_se,
        ymax = gs_excess_mean + gs_excess_se),
    linewidth = 0.85,
    size = 0.75,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.2f", gs_excess_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  scale_y_continuous(
    limits = c(-0.1, 3.75),
    breaks = 0:3,
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    title = panel_title("D", "Excess estimated connectivity"),
    x = NULL,
    y = "Estimated − true global strength"
  ) +
  theme_pub(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_size_strip),
    axis.title = element_text(size = axis_title_strip),
    axis.text = element_text(size = axis_text_strip),
    axis.text.x = element_text(size = axis_text_strip)
  )

# ------------------------------------------------------------------------
# Panel E: deviation on truly absent edges
# ------------------------------------------------------------------------
pE <- ggplot(summary_by_arm, aes(x = arm, y = mtz_mean, colour = arm)) +
  geom_pointrange(
    aes(ymin = mtz_mean - mtz_se,
        ymax = mtz_mean + mtz_se),
    linewidth = 0.65,
    size = 0.55,
    alpha = 0.85
  ) +
  geom_text(
    aes(label = sprintf("%.3f", mtz_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  scale_y_continuous(
    limits = c(0.083, 0.112),
    breaks = c(0.09, 0.10, 0.11),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    title = panel_title("E", "Deviation on absent edges"),
    x = NULL,
    y = "MAE on absent edges"
  ) +
  theme_pub(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_size_strip),
    axis.title = element_text(size = axis_title_strip),
    axis.text = element_text(size = axis_text_strip),
    axis.text.x = element_text(size = axis_text_strip)
  )

fig4_metric_strip <- pD | pE

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.pdf",
  fig4_metric_strip,
  width = 8.4,
  height = 2.75
)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.png",
  fig4_metric_strip,
  width = 8.4,
  height = 2.75,
  dpi = 300
)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure4_metric_strip.pdf (+ .png)\n")
cat("\nIMPORTANT: this is now meant to be placed at \\textwidth in LaTeX, same\n")
cat("as the network trio -- update the \\includegraphics line for this file\n")
cat("from width=0.82\\textwidth to width=\\textwidth (or drop the width arg\n")
cat("if the trio's own includegraphics doesn't specify one explicitly).\n")
cat("\nCombine with fig4_network_trio.R's output (Figure4_network_trio.pdf,\n")
cat("panels A-C) as Figure 4's full panel set -- assembled in LaTeX/Overleaf,\n")
cat("not composited in R (qgraph base-R graphics + ggplot don't combine\n")
cat("directly via patchwork without ggplotify/wrap_elements wiring).\n")
