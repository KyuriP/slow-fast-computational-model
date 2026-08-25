# ============================================================
# R/revision_2026/figures/fig3_recovery_feedback.R
# ============================================================
# Figure 3: Simulations 2 + 3, as a 2x2 grid of panels
# (A: Slow field, Sim 2 (feedback off)  | C: Slow field, Sim 3 (off vs on))
# (B: Symptom burden, Sim 2             | D: Symptom burden, Sim 3 (off vs on))
#
# STYLING PASS (closest of the four to publication level; this pass is
# about visual hierarchy, not content):
#   - main mean line thinner (0.95 -> 0.70 -> 0.55, two rounds) so it
#     doesn't visually dominate the spaghetti/ribbon underneath it
#   - spaghetti pushed further towards actually-visible (alpha up to
#     ~0.20-0.22, more sample chains shown: 15 -> 35) after the first
#     fade pass was judged still too faint -- the goal is that a reader
#     can see this is built from many individual noisy trajectories, not
#     just infer it from a thin haze behind one dominant mean line
#   - ribbons lightened further, especially around the symptom-burden
#     curves, now that the ribbon has to sit quietly under both the
#     spaghetti and the mean line
#   - direct end-labels ("Feedback on/off") given more horizontal room
#     and de-bolded slightly so they read as integrated labels, not
#     shouted annotations
#   - shock-onset dashed line lightened so it doesn't pull attention
#     away from the trajectories
#   - slightly smaller panel titles/axis text (theme_pub(base_size=11)
#     instead of the default 12) -- small, but this is one of the things
#     that makes a figure read as journal-typeset rather than R-default
#   - left/right panels within each row already share one y-scale
#     (p_ylim / m_ylim computed jointly across Sim 2 + Sim 3) and are
#     combined via patchwork's | / operators, which enforces identical
#     panel geometry -- no change needed there, already a true pair
#
# Outputs
# -------
#   figs/revision_2026/Figure3_recovery_feedback.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("R/revision_2026/figures/theme_publication.R")  # theme_pub(), panel_title(), col_*, pal_feedback

# ------------------------------------------------------------------------
# 0. Load data
# ------------------------------------------------------------------------
sim2 <- read_csv("res/revision_2026/sim2/sim2_summary.csv", show_col_types = FALSE)
sim3 <- read_csv("res/revision_2026/sim3/sim3_summary.csv", show_col_types = FALSE)

plot_window <- c(-100, 750)
sim2_w <- sim2 |> filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])
sim3_w <- sim3 |> filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])

N <- 9   # symptom count, for converting se_M -> se on the m scale
base_sz <- 11

# ------------------------------------------------------------------------
# 0a. Rolling-window smoothing for the m_t panels (display only). Two-
#     sided moving average; edges become NA and are silently dropped by
#     geom_line, trimming only a sliver of the flat baseline/far tail.
# ------------------------------------------------------------------------
roll_mean <- function(x, k = 15) as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))

sim2_w <- sim2_w |>
  arrange(time_since_shock) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N))

sim3_w <- sim3_w |>
  arrange(feedback, time_since_shock) |>
  group_by(feedback) |>
  mutate(mean_m_smooth = roll_mean(mean_m), se_m_smooth = roll_mean(se_M / N)) |>
  ungroup()

# ------------------------------------------------------------------------
# 0b. Raw single-chain trajectories for panels A/C -- thinner and more
#     opaque than the previous pass (see header note).
# ------------------------------------------------------------------------
set.seed(3)
n_sample_chains <- 35L

sim2_raw <- readRDS("res/revision_2026/sim2/sim2_raw.rds")$traj |>
  filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])
sim2_chains_w <- sim2_raw |> filter(chain %in% sample(unique(sim2_raw$chain), n_sample_chains))

sim3_raw <- readRDS("res/revision_2026/sim3/sim3_raw.rds")$traj |>
  filter(time_since_shock >= plot_window[1], time_since_shock <= plot_window[2])
sim3_chains_w <- sim3_raw |> filter(chain %in% sample(unique(sim3_raw$chain), n_sample_chains))

# ------------------------------------------------------------------------
# 1. Shared y-limits and interval settings
# ------------------------------------------------------------------------
interval_mult <- 1.96  # 95% intervals

add_se_m <- function(df) {
  if (!"se_m" %in% names(df)) df <- df |> mutate(se_m = se_M / N)
  df
}
sim2_w <- add_se_m(sim2_w)
sim3_w <- add_se_m(sim3_w)

range_pad <- function(x, pad = 0.06, lower = -Inf, upper = Inf) {
  r <- range(x, na.rm = TRUE)
  d <- diff(r)
  if (d == 0) d <- 1
  c(max(lower, r[1] - pad * d), min(upper, r[2] + pad * d))
}

p_ylim <- range_pad(
  c(
    sim2_w$mean_P - interval_mult * sim2_w$se_P,
    sim2_w$mean_P + interval_mult * sim2_w$se_P,
    sim3_w$mean_P - interval_mult * sim3_w$se_P,
    sim3_w$mean_P + interval_mult * sim3_w$se_P
  ),
  pad = 0.06
)

m_ylim <- range_pad(
  c(
    sim2_w$mean_m_smooth - interval_mult * sim2_w$se_m_smooth,
    sim2_w$mean_m_smooth + interval_mult * sim2_w$se_m_smooth,
    sim3_w$mean_m_smooth - interval_mult * sim3_w$se_m_smooth,
    sim3_w$mean_m_smooth + interval_mult * sim3_w$se_m_smooth
  ),
  pad = 0.08,
  lower = 0,
  upper = 1
)

p_ref <- mean(sim2_w$mean_P[sim2_w$time_since_shock < 0], na.rm = TRUE)
m_ref <- mean(sim2_w$mean_m_smooth[sim2_w$time_since_shock < 0], na.rm = TRUE)

# C/D use colour+linetype to encode FEEDBACK CONDITION (grey dashed = off,
# blue solid = on) -- a different encoding than A/B, where colour encodes
# the QUANTITY (blue = slow field, orange = symptom burden) because there
# is only one condition (feedback off) to show. That mismatch was
# confusing without an explicit key, so C/D now also carry a shared,
# bottom-of-figure legend (collected via patchwork) spelling out colour +
# linetype together, in addition to the direct end-labels.
feedback_labels <- c(off = "Feedback off", on = "Feedback on")

# x-range padded on the right so the direct off/on labels in C/D have
# more breathing room than the previous pass
x_max_plot <- plot_window[2]
x_label_pos <- x_max_plot + 30
x_expand_right <- c(0.01, 0.17)

# ------------------------------------------------------------------------
# 2. Panel A: slow field, Simulation 2
# ------------------------------------------------------------------------
pA <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_P)) +
  geom_hline(yintercept = p_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = sim2_chains_w,
    aes(x = time_since_shock, y = P, group = chain),
    inherit.aes = FALSE,
    colour = col_P, alpha = 0.22, linewidth = 0.16
  ) +
  geom_ribbon(
    aes(ymin = mean_P - interval_mult * se_P, ymax = mean_P + interval_mult * se_P),
    fill = col_P, alpha = 0.08
  ) +
  geom_line(colour = col_P, linewidth = 0.55) +
  coord_cartesian(ylim = p_ylim, clip = "off") +
  labs(title = panel_title("A", "Slow field after shock"), x = NULL, y = "Slow field (P)") +
  theme_pub(base_size = base_sz) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# ------------------------------------------------------------------------
# 3. Panel B: symptom burden, Simulation 2 (rolling mean)
# ------------------------------------------------------------------------
pB <- ggplot(sim2_w, aes(x = time_since_shock, y = mean_m_smooth)) +
  geom_hline(yintercept = m_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = mean_m_smooth - interval_mult * se_m_smooth, ymax = mean_m_smooth + interval_mult * se_m_smooth),
    fill = col_M, alpha = 0.08
  ) +
  geom_line(colour = col_M, linewidth = 0.55) +
  coord_cartesian(ylim = m_ylim) +
  labs(title = panel_title("B", "Mean symptom activation after shock"), x = "Steps since shock", y = "Mean symptom activation (m)") +
  theme_pub(base_size = base_sz)

# ------------------------------------------------------------------------
# 4. Panel C: slow field, Simulation 3 (off vs on), direct end-labels
# ------------------------------------------------------------------------
sim3_end <- sim3_w |> group_by(feedback) |> filter(time_since_shock == max(time_since_shock)) |> ungroup()

pC <- ggplot(
  sim3_w,
  aes(x = time_since_shock, y = mean_P, colour = feedback, fill = feedback, linetype = feedback)
) +
  geom_hline(yintercept = p_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_line(
    data = sim3_chains_w,
    aes(x = time_since_shock, y = P, group = interaction(chain, feedback), colour = feedback),
    inherit.aes = FALSE, alpha = 0.20, linewidth = 0.15, show.legend = FALSE
  ) +
  geom_ribbon(
    aes(ymin = mean_P - interval_mult * se_P, ymax = mean_P + interval_mult * se_P),
    alpha = 0.06, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.55) +
  geom_text(
    data = sim3_end, aes(x = x_label_pos, y = mean_P, label = ifelse(feedback == "off", "Feedback off", "Feedback on")),
    hjust = 0, size = 3.0, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(ylim = p_ylim, xlim = c(plot_window[1], x_max_plot), clip = "off") +
  scale_colour_manual(values = pal_feedback, labels = feedback_labels) +
  scale_fill_manual(values = pal_feedback, guide = "none") +
  scale_linetype_manual(values = c(off = "dashed", on = "solid"), labels = feedback_labels) +
  scale_x_continuous(expand = expansion(mult = x_expand_right)) +
  labs(title = panel_title("C", "Feedback slows field recovery"), x = NULL, y = NULL, colour = NULL, linetype = NULL) +
  theme_pub(base_size = base_sz) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.margin = margin(5.5, 68, 5.5, 5.5)
  )

# ------------------------------------------------------------------------
# 5. Panel D: symptom burden, Simulation 3 (off vs on, rolling mean),
#    direct end-labels
# ------------------------------------------------------------------------
sim3_end_m <- sim3_w |> filter(!is.na(mean_m_smooth)) |> group_by(feedback) |>
  filter(time_since_shock == max(time_since_shock)) |> ungroup()

pD <- ggplot(
  sim3_w,
  aes(x = time_since_shock, y = mean_m_smooth, colour = feedback, fill = feedback, linetype = feedback)
) +
  geom_hline(yintercept = m_ref, colour = "grey78", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = mean_m_smooth - interval_mult * se_m_smooth, ymax = mean_m_smooth + interval_mult * se_m_smooth),
    alpha = 0.06, colour = NA, show.legend = FALSE
  ) +
  geom_line(linewidth = 0.55) +
  geom_text(
    data = sim3_end_m, aes(x = x_label_pos, y = mean_m_smooth, label = ifelse(feedback == "off", "Feedback off", "Feedback on")),
    hjust = 0, size = 3.0, fontface = "plain", show.legend = FALSE
  ) +
  coord_cartesian(ylim = m_ylim, xlim = c(plot_window[1], x_max_plot), clip = "off") +
  scale_colour_manual(values = pal_feedback, labels = feedback_labels) +
  scale_fill_manual(values = pal_feedback, guide = "none") +
  scale_linetype_manual(values = c(off = "dashed", on = "solid"), labels = feedback_labels) +
  scale_x_continuous(expand = expansion(mult = x_expand_right)) +
  labs(title = panel_title("D", "Feedback prolongs symptom activation"), x = "Steps since shock", y = NULL, colour = NULL, linetype = NULL) +
  theme_pub(base_size = base_sz) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 68, 5.5, 5.5)
  )

# ------------------------------------------------------------------------
# 6. Combine + save -- shared legend (colour + linetype merge into one
#    key automatically since both scales use the same "feedback" labels)
#    collected once at the bottom, on top of the direct end-labels
#    already on C/D.
# ------------------------------------------------------------------------
fig3 <- (pA | pC) / (pB | pD) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

ggsave("figs/revision_2026/Figure3_recovery_feedback.pdf", fig3, width = 9.8, height = 6.6)
ggsave("figs/revision_2026/Figure3_recovery_feedback.png", fig3, width = 9.8, height = 6.6, dpi = 300)

cat("Done. Files:\n")
cat("  figs/revision_2026/Figure3_recovery_feedback.pdf (+ .png)\n")
