# ============================================================
# R/revision_2026/figures/fig4_summary_panels.R
# ============================================================
# SUPERSEDED as main-figure content by fig4_metric_strip.R -- this
# 4-panel version tries to do three jobs at once (design schematic +
# summary stats + replicate-level robustness check), which is why it
# reads as sparse/unfocused regardless of per-panel polish. Kept only as
# a development/robustness-check figure (panel D here = the paired
# naive-adjusted jitter plot, useful as a supplement item) -- not to be
# used as main-text Figure 4. See fig4_metric_strip.R + fig4_network_trio.R
# for the current main-figure panels.
#
# Figure 4's quantitative summary panels -- the companion to
# fig4_network_trio.R's qualitative True/Symptom-only/Context-adjusted
# diagrams. Built from the REPLICATED Simulation 4 run (04b, 30 reps),
# the locked quantitative result, not the single-run pilot.
#
# Four panels:
#   A. Design schematic -- what each arm actually is (fixed vs pooled
#      context, P omitted vs included in the estimator). Not a data plot;
#      orients the reader before B-D.
#   B. Estimated global strength by arm, mean +/- SE across 30 replicates,
#      with the true value marked as a reference line.
#   C. Recovery error on true-zero edges by arm (the "phantom edge"
#      metric) -- this is the sharpest single number for "did omitting
#      context fabricate apparent coupling."
#   D. Paired naive-minus-adjusted difference, per replicate -- same data
#      as fig_sim4_replicated_diff.pdf (built during calibration), but
#      restyled with theme_pub() and with independent y-scales per facet
#      so the MAE panel (small values) doesn't get crushed by the global-
#      strength panel (large values) sharing one axis, which was a real
#      bug in the calibration-stage version.
#
# Outputs
# -------
#   figs/revision_2026/Figure4_summary_panels.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")   # theme_pub(), pal_arm

by_arm <- read.csv("res/revision_2026/sim4/sim4_replicated_by_arm.csv") |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")))
diff_all <- read.csv("res/revision_2026/sim4/sim4_replicated_naive_minus_adjusted.csv")

true_gs <- 3.05   # true global strength (10 edges, weights 0.20-0.45), matches
                   # locked Sim 4 results text; not recomputed here to keep
                   # this script's only inputs the two replicated-run CSVs.

se <- function(x) sd(x) / sqrt(length(x))

summary_by_arm <- by_arm |>
  group_by(arm) |>
  summarise(
    gs_mean = mean(global_strength_est), gs_se = se(global_strength_est),
    mtz_mean = mean(mae_true_zero), mtz_se = se(mae_true_zero),
    .groups = "drop"
  )

arm_labels <- c(baseline = "Fixed-context\nbaseline", naive = "Symptom-only\n(naive)",
                 adjusted = "Context-\nadjusted")

# ------------------------------------------------------------------------
# Panel A: design schematic (not a data plot)
# ------------------------------------------------------------------------
# Built with geom_label() rather than hand-computed geom_rect() widths --
# geom_label() auto-sizes its box to its own text, so there is no manual
# "does this text fit in that rectangle" arithmetic to get wrong. The
# previous version of this panel hand-computed rect widths, got them
# wrong, and text overlapped -- geom_label() makes that whole class of
# bug structurally impossible instead of something to get right by hand.
#
# What was still wrong after that fix: the panel's fixed x/y limits were
# sized for a worst-case label, but the actual boxes+arrows only fill the
# middle of that space, leaving a large, unbalanced blank margin (most
# visible on the right) that the data rows below (B|C, full-bleed to the
# panel edges) don't have -- so panel A reads as sparse/unfinished next
# to them. Fix: a full-bleed, per-row background tint (xmin/xmax pinned
# to the SAME xlim_A variable used by scale_x_continuous, so they cannot
# drift out of sync) turns the panel into three deliberate colored rows
# instead of text floating in a void, without hand-measuring label widths.
xlim_A <- c(-1.7, 2.7)
ylim_A <- c(0.3, 3.7)

schematic <- tibble(
  row = c(3, 2, 1),
  arm = c("baseline", "naive", "adjusted"),
  row_label = c("Fixed-context\nbaseline", "Symptom-only\n(naive)", "Context-\nadjusted"),
  context_label = c("everyone at\nfixed P = 0", "P[i] varies,\nper person", "P[i] varies,\nper person"),
  estimator_label = c("network\nestimated", "P omitted\nfrom estimator", "P included\nin estimator")
) |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")),
         ymin = row - 0.42, ymax = row + 0.42)

pA <- ggplot(schematic) +
  geom_rect(aes(xmin = xlim_A[1], xmax = xlim_A[2], ymin = ymin, ymax = ymax, fill = arm),
            alpha = 0.07) +
  geom_segment(aes(x = 1, xend = 2, y = row, yend = row), colour = "grey60", linewidth = 0.3,
               arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_label(aes(x = 0, y = row, label = row_label, colour = arm), fill = "white",
             fontface = "bold", size = 3.1, label.size = 0, hjust = 1, lineheight = 0.9) +
  geom_label(aes(x = 1, y = row, label = context_label), fill = "grey97", colour = "grey20",
             size = 3, label.padding = unit(0.35, "lines"), lineheight = 0.9) +
  geom_label(aes(x = 2, y = row, label = estimator_label, colour = arm, fill = arm), alpha = 0.15,
             fontface = "bold", size = 3, label.padding = unit(0.35, "lines"), lineheight = 0.9) +
  scale_fill_manual(values = pal_arm, guide = "none") +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_continuous(limits = xlim_A, expand = c(0, 0)) +
  scale_y_continuous(limits = ylim_A, expand = c(0, 0)) +
  labs(title = panel_title("A", "What each arm is")) +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "plain", hjust = 0, margin = margin(b = 10)))

# ------------------------------------------------------------------------
# Panel B: global strength by arm
# ------------------------------------------------------------------------
# Points alone read as small and lonely against the wide y-range the
# dashed "true" reference forces (3.05 to ~6.5) -- adding the numeric
# mean next to each point (same move that fixed Figure 2's panel C) gives
# each point real content instead of just a floating dot, and a bigger
# point/line weight makes the three-way comparison read as the main
# content of the panel rather than incidental markers on empty axes.
pB <- ggplot(summary_by_arm, aes(x = arm, y = gs_mean, colour = arm)) +
  geom_hline(yintercept = true_gs, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.55, y = true_gs, label = "true", size = 3.2, colour = "grey40", hjust = 0, vjust = -0.6) +
  geom_pointrange(aes(ymin = gs_mean - gs_se, ymax = gs_mean + gs_se), linewidth = 1.1, size = 0.85) +
  geom_text(aes(label = sprintf("%.2f", gs_mean)), nudge_x = 0.20, hjust = 0, size = 3.4, colour = "grey20") +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(title = panel_title("B", "Estimated global strength"), x = NULL, y = "Global strength (mean ± SE)") +
  theme_pub() +
  theme(legend.position = "none")

# ------------------------------------------------------------------------
# Panel C: recovery error on true-zero edges by arm
# ------------------------------------------------------------------------
pC <- ggplot(summary_by_arm, aes(x = arm, y = mtz_mean, colour = arm)) +
  geom_pointrange(aes(ymin = mtz_mean - mtz_se, ymax = mtz_mean + mtz_se), linewidth = 1.1, size = 0.85) +
  geom_text(aes(label = sprintf("%.3f", mtz_mean)), nudge_x = 0.20, hjust = 0, size = 3.4, colour = "grey20") +
  scale_colour_manual(values = pal_arm, guide = "none") +
  scale_x_discrete(labels = arm_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(title = panel_title("C", "Error on truly-absent edges"), x = NULL, y = "MAE, true-zero edges (mean ± SE)") +
  theme_pub() +
  theme(legend.position = "none")

# ------------------------------------------------------------------------
# Panel D: paired naive-adjusted difference across replicates
# ------------------------------------------------------------------------
diff_long <- diff_all |>
  select(rep, diff_global_strength, diff_mae_true_zero) |>
  pivot_longer(-rep, names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                          diff_global_strength = "Global strength",
                          diff_mae_true_zero = "MAE, true-zero edges"))

pD <- ggplot(diff_long, aes(x = 1, y = value)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.4, colour = pal_arm[["naive"]]) +
  stat_summary(fun = mean, geom = "point", size = 2.6, colour = "grey10") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, colour = "grey10") +
  facet_wrap(~metric, scales = "free_y") +
  scale_x_continuous(breaks = NULL) +
  labs(title = panel_title("D", "Naive − adjusted, per replicate (n=30)"), x = NULL, y = "Naive − adjusted") +
  theme_pub() +
  theme(legend.position = "none")

# ------------------------------------------------------------------------
# Combine -- full-width rows (schematic, then two side-by-side point-range
# panels, then the replicate-level panel), rather than a left/right split.
# Gives panel A room to breathe and reads like a standard stacked-row
# journal multi-panel figure instead of a cramped 2x2 grid.
# ------------------------------------------------------------------------
fig4_summary <- pA / (pB | pC) / pD +
  patchwork::plot_layout(heights = c(0.75, 1.1, 1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/Figure4_summary_panels.pdf", fig4_summary, width = 9, height = 9.5)
ggsave("figs/revision_2026/Figure4_summary_panels.png", fig4_summary, width = 9, height = 9.5, dpi = 200)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure4_summary_panels.pdf (+ .png)\n")
cat("\nCombine with fig4_network_trio.R's output (Figure4_network_trio.pdf) as\n")
cat("Figure 4's full panel set -- network trio + these four summary panels --\n")
cat("or split across main text / appendix as decided.\n")
