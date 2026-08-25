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
  annotate("text", x = 0.58, y = 0, label = "No inflation", size = 2.4,
           colour = "grey45", hjust = 0, vjust = -0.6) +
  geom_pointrange(
    aes(ymin = gs_excess_mean - gs_excess_se,
        ymax = gs_excess_mean + gs_excess_se),
    linewidth = 0.65,
    size = 0.55,
    alpha = 0.85
  ) +
  geom_text(
    aes(label = sprintf("%.2f", gs_excess_mean)),
    nudge_x = 0.10,
    hjust = 0,
    size = 3.0,
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
  theme(legend.position = "none", axis.text.x = element_text(size = 8.3))

# ------------------------------------------------------------------------
# Panel E: error on truly absent edges
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
    size = 3.0,
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
    title = panel_title("E", "Error on absent edges"),
    x = NULL,
    y = "MAE on absent edges"
  ) +
  theme_pub(base_size = 9.5) +
  theme(legend.position = "none", axis.text.x = element_text(size = 8.3))

fig4_metric_strip <- pD | pE

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.pdf",
  fig4_metric_strip,
  width = 8.0,
  height = 2.6
)

ggsave(
  "figs/revision_2026/Figure4_metric_strip.png",
  fig4_metric_strip,
  width = 8.0,
  height = 2.6,
  dpi = 300
)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure4_metric_strip.pdf (+ .png)\n")
cat("\nCombine with fig4_network_trio.R's output (Figure4_network_trio.pdf,\n")
cat("panels A-C) as Figure 4's full panel set -- assembled in LaTeX/Overleaf,\n")
cat("not composited in R (qgraph base-R graphics + ggplot don't combine\n")
cat("directly via patchwork without ggplotify/wrap_elements wiring).\n")
