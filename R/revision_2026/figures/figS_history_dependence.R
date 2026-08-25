# ============================================================
# R/revision_2026/figures/figS_history_dependence.R
# ============================================================
# Supplementary figure: history-dependence / tipping-like regime check.
# This figure should be included only if the result is clean and not merely
# runaway saturation.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
})

source("R/revision_2026/figures/theme_publication.R")

sim_hist <- readRDS("res/revision_2026/supp_history/history_dependence_raw.rds")
sep <- read_csv("res/revision_2026/supp_history/history_separation_summary.csv", show_col_types = FALSE)

# Choose representative b values:
# - current main feedback level
# - first history-dependent non-runaway regime, if available
candidate_b <- sep |>
  filter(regime_flag == "history-dependent") |>
  arrange(b) |>
  pull(b) |>
  first()

if (is.na(candidate_b)) {
  candidate_b <- max(sep$b)
}

b_show <- unique(c(0.50, candidate_b))

plot_dat <- sim_hist |>
  filter(b %in% b_show) |>
  group_by(b, init, step) |>
  summarise(
    mean_m = mean(m),
    se_m = sd(m) / sqrt(n()),
    mean_P = mean(P),
    se_P = sd(P) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(
    b_label = paste0("b = ", b),
    init = factor(init, levels = c("low", "high"),
                  labels = c("Low initial state", "High initial state"))
  )

init_cols <- c(
  "Low initial state" = "#4C9A6A",
  "High initial state" = "#C85C5C"
)

# ------------------------------------------------------------------------
# Panel A: symptom burden trajectories
# ------------------------------------------------------------------------
pA <- ggplot(plot_dat, aes(x = step, y = mean_m, colour = init, fill = init)) +
  geom_ribbon(aes(ymin = mean_m - 1.96 * se_m, ymax = mean_m + 1.96 * se_m),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ b_label, nrow = 1) +
  scale_colour_manual(values = init_cols) +
  scale_fill_manual(values = init_cols, guide = "none") +
  labs(
    title = panel_title("A", "Late symptom activation can depend on initial state"),
    x = "Simulation step",
    y = expression("Mean symptom activation" ~ (m[t])),
    colour = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

# ------------------------------------------------------------------------
# Panel B: separation index across feedback strengths
# ------------------------------------------------------------------------
pB <- ggplot(sep, aes(x = b, y = separation_m)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = 0.15, linetype = "dotted", colour = "grey55") +
  geom_pointrange(
    aes(ymin = separation_m - 1.96 * se_separation_m,
        ymax = separation_m + 1.96 * se_separation_m,
        colour = regime_flag),
    linewidth = 0.9,
    size = 0.75
  ) +
  scale_colour_manual(
    values = c(
      "convergent" = "grey45",
      "history-dependent" = "#C85C5C",
      "runaway/saturation" = "black"
    )
  ) +
  labs(
    title = panel_title("B", "History-dependence index"),
    x = "Feedback strength (b)",
    y = expression(Delta ~ "mean symptom activation (late)"),
    colour = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

figS <- pA / pB + plot_layout(heights = c(1.2, 1))

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave(
  "figs/revision_2026/FigureS_history_dependence.pdf",
  figS,
  width = 9,
  height = 7
)
ggsave(
  "figs/revision_2026/FigureS_history_dependence.png",
  figS,
  width = 9,
  height = 7,
  dpi = 300
)

cat("Done. Files:\n")
cat("  figs/revision_2026/FigureS_history_dependence.pdf (+ .png)\n")
