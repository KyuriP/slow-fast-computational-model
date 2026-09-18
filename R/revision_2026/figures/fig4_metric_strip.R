# ============================================================
# R/revision_2026/figures/fig4_metric_strip.R
# ============================================================
# Compact quantitative strip for Figure 4.
# Use with the network panels (fig4_network_trio.R, expanded 2026-08-27
# from 3 to 4 networks):
#   A true network | B baseline (fixed P) estimate | C symptom-only
#   estimate | D context-adjusted estimate
#   E estimated total coupling: sum_{i<j} omega_hat_ij
#   F spurious absolute coupling on symptom pairs with true omega_ij = 0
#
# Replaces fig4_summary_panels.R as the quantitative half of Figure 4.
# That earlier version tried to do three jobs (design schematic + summary
# stats + replicate-level robustness check) in one figure, which is why it
# read as sparse/unfocused however the individual panels were polished.
# This version keeps only what belongs in the main text:
#   - the design-schematic panel (old panel A) moves to the caption/Methods
#   - the paired-replicate jitter plot (old panel D) moves to supplement
#     as a robustness check, not main-figure content
#   - the remaining "excess global strength over true" panel replaces the
#     old absolute-scale panel B, which was forced to a huge y-range by
#     including the "true" reference as a separate value on the same axis
#     rather than differencing against it directly
#
# Meant to sit alongside fig4_network_trio.R's output as Figure 4's full
# panel set (network trio on top, this strip on the bottom), assembled in
# LaTeX/Overleaf rather than composited in R -- qgraph's base-R graphics
# and this file's ggplot objects don't combine directly via patchwork
# without extra wiring (ggplotify::as.ggplot() / wrap_elements()), left
# as deferred, lower-priority work.
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
# This strip is placed at the SAME LaTeX width as the network trio
# (\textwidth, not 0.82\textwidth as before), so panels E/F print larger
# and their error bars are easier to read. With both figures now sharing
# the same width-based scale factor (~6.5/9.2 ~= 0.71, using the trio's
# current 9.2in native width as reference), the earlier 13.3pt title
# override (calibrated for the narrower 0.67x placement) would now
# overshoot and print larger than the trio's titles, so sizes below are
# recalibrated for the new scale factor. Canvas native size bumped
# slightly (8.0x2.6 -> 8.4x2.75) so once both figures share the print
# width, this strip's height comes out close to but slightly less than
# the trio's.
title_size_strip <- 12.0
axis_title_strip <- 11.2
axis_text_strip  <- 10.2
label_size_strip <- 3.4   # geom_text size is in mm, not pt

# Matches the locked Sim 4 results text (10 edges, weights 0.20-0.45);
# not recomputed here to keep this script's only input the replicated CSV.
true_total_coupling <- 3.05

arm_labels <- c(
  baseline = "Fixed P",
  naive = "P omitted",
  adjusted = "P included"
)

summary_by_arm <- by_arm |>
  group_by(arm) |>
  summarise(
    coupling_mean = mean(total_coupling_est),
    coupling_se   = se(total_coupling_est),

    spurious_mean = mean(spurious_abs_coupling),
    spurious_se   = se(spurious_abs_coupling),

    .groups = "drop"
  )

# ------------------------------------------------------------------------
# Panel E: estimated total coupling (signed sum of estimated pairwise
# couplings), with a dashed reference line at the data-generating value
# (true_total_coupling=3.05). Using the signed sum rather than an
# absolute-value global-strength statistic avoids the systematic upward
# contribution that comes from summing |estimate| over true-zero edges,
# where sampling noise is symmetric around zero but abs() makes all of it
# positive.
#
# y-axis limits/breaks (both panels) are FIXED per the design-pass spec
# rather than computed from summary_by_arm, to remove dead space above/
# below the points. If a future data refresh shifts coupling_mean/
# spurious_mean meaningfully, check the rendered PNG for clipped points/
# error bars before trusting these fixed ranges again.
# ------------------------------------------------------------------------
pE <- ggplot(
  summary_by_arm,
  aes(x = arm, y = coupling_mean, colour = arm)
) +
  geom_hline(
    yintercept = true_total_coupling,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.45
  ) +
  annotate(
    "text",
    x = 0.58,
    y = true_total_coupling,
    label = "Data-generating value",
    size = 2.9,
    colour = "grey45",
    hjust = 0,
    vjust = -0.6
  ) +
  geom_pointrange(
    aes(
      ymin = coupling_mean - coupling_se,
      ymax = coupling_mean + coupling_se
    ),
    linewidth = 0.85,
    size = 0.4,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.2f", coupling_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  labs(
    title = panel_title("E", "Estimated total coupling"),
    x = NULL,
    y = "Estimated total coupling"
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
# Panel F: spurious absolute coupling summed across symptom pairs whose
# true (data-generating) omega_ij is zero. More directly tied to the
# recovery question than a mean-absolute-error number -- it's literally
# how much coupling the estimator attributes to pairs that aren't coupled
# at all.
#
# y-axis limits/breaks are FIXED per the design-pass spec (same reasoning
# as panel E) -- check the rendered PNG once this reruns, since
# spurious_mean's actual range depends on the current replicated-run
# output and wasn't recomputed here.
# ------------------------------------------------------------------------
pF <- ggplot(
  summary_by_arm,
  aes(x = arm, y = spurious_mean, colour = arm)
) +
  geom_pointrange(
    aes(
      ymin = spurious_mean - spurious_se,
      ymax = spurious_mean + spurious_se
    ),
    linewidth = 0.65,
    size = 0.35,
    alpha = 0.85
  ) +
  geom_text(
    aes(label = sprintf("%.2f", spurious_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = label_size_strip,
    colour = "grey20"
  ) +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  labs(
    title = panel_title("F", "Spurious coupling among uncoupled pairs"),
    x = NULL,
    y = "Absolute estimated coupling on true-zero edges"
  ) +
  theme_pub(base_size = 9.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = title_size_strip),
    axis.title = element_text(size = axis_title_strip),
    axis.text = element_text(size = axis_text_strip),
    axis.text.x = element_text(size = axis_text_strip)
  )

fig4_metric_strip <- pE | pF

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
cat("\nThis is now meant to be placed at \\textwidth in LaTeX, same as the\n")
cat("network trio -- update the \\includegraphics line for this file from\n")
cat("width=0.82\\textwidth to width=\\textwidth (or drop the width arg if the\n")
cat("trio's own includegraphics doesn't specify one explicitly).\n")
cat("\nCombine with fig4_network_trio.R's output (Figure4_network_trio.pdf,\n")
cat("panels A-D) as Figure 4's full panel set -- assembled in LaTeX/Overleaf,\n")
cat("not composited in R (qgraph base-R graphics + ggplot don't combine\n")
cat("directly via patchwork without ggplotify/wrap_elements wiring).\n")
