# ============================================================
# Hysteresis phase GIF
# - animated trajectory in (P, m)
# - color indicates ramp-up vs ramp-down phase
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gganimate)
})

# ------------------------------------------------------------
# 1) Load trajectory
# ------------------------------------------------------------
sim_B2 <- readRDS("res/hysteresis/sim_B2.rds")

# ------------------------------------------------------------
# 2) Prepare animation data
# ------------------------------------------------------------
# If P_base is available, use its peak to define the ramp reversal.
# Otherwise fall back to halfway through the time axis.
if ("P_base" %in% names(sim_B2)) {
  t_turn <- sim_B2$t[which.max(sim_B2$P_base)]
} else {
  t_turn <- max(sim_B2$t, na.rm = TRUE) / 2
}

thin_step <- 5

dfB <- sim_B2 %>%
  mutate(
    phase = if_else(t <= t_turn, "ramp up", "ramp down"),
    phase = factor(phase, levels = c("ramp up", "ramp down"))
  ) %>%
  slice(seq(1, n(), by = thin_step))

# ------------------------------------------------------------
# 3) Animated phase plot
# ------------------------------------------------------------
p_gif <- ggplot(dfB, aes(P, m)) +
  geom_path(aes(colour = phase), alpha = 0.85, linewidth = 1.0) +
  geom_point(aes(colour = phase), size = 2.5) +
  scale_colour_manual(
    values = c(
      "ramp up"   = "#2C7BB6",
      "ramp down" = "#D7191C"
    )
  ) +
  labs(
    title = "Scenario B: hysteresis as a trajectory in (P, m)",
    subtitle = "Frame: t = {frame_time}",
    x = "Context P",
    y = "Mean symptom activation m",
    colour = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  transition_reveal(along = t)

# ------------------------------------------------------------
# 4) Render + save
# ------------------------------------------------------------
anim <- animate(
  p_gif,
  nframes = nrow(dfB),
  fps = 25,
  width = 1200,
  height = 500,
  renderer = gifski_renderer()
)

# anim_save("img/B2_hysteresis_phase.gif", anim)