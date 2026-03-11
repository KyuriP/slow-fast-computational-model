library(gganimate)


# 1) Subsample to keep file size reasonable
dfB <- sim_B2 |>
  mutate(phase = if_else(t <= max(t)/2, "ramp up", "ramp down")) |>
  slice(seq(1, n(), by = 5))  # every 5th step; adjust to 10 if still big

# 2) Animated phase plot: moving point + trailing path
p_gif <- ggplot(dfB, aes(P, m)) +
  geom_path(aes(colour = phase), alpha = 0.85, linewidth = 1.0) +
  geom_point(aes(colour = phase), size = 2.5) +
  scale_colour_manual(values = c("ramp up" = "#2C7BB6", "ramp down" = "#D7191C")) +
  labs(
    title = "Scenario B: hysteresis as a trajectory in (P, m)",
    subtitle = "Frame: t = {frame_time}",
    x = "Context  P",
    y = "Mean symptom activation  m",
    colour = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(legend.position = "top") +
  transition_reveal(along = t)  # reveals path as time progresses

# 3) Render
anim <- animate(
  p_gif,
  nframes = nrow(dfB),   # one frame per sampled row
  fps = 25,
  width = 1200,
  height = 500,
  renderer = gifski_renderer()
)

anim_save("img/B2_hysteresis_phase.gif", anim)




