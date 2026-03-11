# ============================================================
# 01_figure_bistability_geometry.R
# ============================================================
# Figure 1: Geometric intuition for bistability (βJ > 4)
# 0/1 Curie–Weiss mean-field map
#
# Panels
#   A: Fixed-point geometry for βJ < 4 vs βJ > 4 at h_eff = 0
#   B: Bifurcation diagram m*(P) for a bistable-capable regime
#
# Output
#   fig_bistability_geometry
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

theme_set(theme_minimal(base_size = 13))

# ------------------------------------------------------------
# 1) Mean-field helpers
# ------------------------------------------------------------
inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

# Mean-field update map:
#   m_{t+1} = F(m_t; h_eff)
cw01_map <- function(m, beta, J, h_eff) {
  inv_logit(beta * (J * (m - 0.5) + h_eff))
}

# Slope dF/dm evaluated at m
cw01_map_slope <- function(m, beta, J, h_eff) {
  p <- cw01_map(m, beta, J, h_eff)
  beta * J * p * (1 - p)
}

# ------------------------------------------------------------
# 2) Fixed points in [0, 1]
# ------------------------------------------------------------
# Solve F(m) - m = 0 numerically on a dense grid.
find_fixed_points_01 <- function(beta, J, h_eff, grid_n = 4001) {
  f <- function(m) cw01_map(m, beta, J, h_eff) - m
  
  grid <- seq(0, 1, length.out = grid_n)
  vals <- f(grid)
  
  idx <- which(diff(sign(vals)) != 0)
  
  roots <- c()
  for (k in idx) {
    r <- try(stats::uniroot(f, c(grid[k], grid[k + 1]))$root, silent = TRUE)
    if (!inherits(r, "try-error")) {
      roots <- c(roots, r)
    }
  }
  
  sort(unique(roots))
}

# ------------------------------------------------------------
# 3) Fixed-point table across P
# ------------------------------------------------------------
# For each P, compute roots and label them low / mid / high.
fixed_points_for_P <- function(P, beta, J, h0, gammaP, grid_n = 2001) {
  h_eff <- h0 + gammaP * P
  roots <- find_fixed_points_01(beta, J, h_eff, grid_n = grid_n)
  
  if (length(roots) == 0) {
    return(tibble())
  }
  
  tibble(
    P = P,
    m = roots,
    stable = cw01_map_slope(roots, beta, J, h_eff) < 1
  ) %>%
    arrange(m) %>%
    group_by(P) %>%
    mutate(
      k = row_number(),
      K = n(),
      branch = case_when(
        K == 1 ~ ifelse(m < 0.5, "low", "high"),
        K == 3 & k == 1 ~ "low",
        K == 3 & k == 2 ~ "mid",
        K == 3 & k == 3 ~ "high",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
}

build_bifurcation_table <- function(P_grid, beta, J, h0, gammaP) {
  purrr::map_dfr(
    P_grid,
    ~ fixed_points_for_P(.x, beta = beta, J = J, h0 = h0, gammaP = gammaP)
  )
}

# ------------------------------------------------------------
# 4) Parameters
# ------------------------------------------------------------
# Panel A: compare βJ < 4 vs βJ > 4 at the same effective bias h_eff = 0
h_eff_A <- 0

par_low <- list(
  beta = 1.0,
  J    = 3.0    # βJ = 3 < 4
)

par_high <- list(
  beta = 2.0,
  J    = 3.0    # βJ = 6 > 4
)

# Panel B: bistable-capable regime across P
par_bif <- list(
  beta   = 2.0,
  J      = 3.0,
  h0     = 0.0,
  gammaP = 1.0
)

P_grid <- seq(-2, 2, length.out = 400)

# ------------------------------------------------------------
# 5) Panel A: fixed-point geometry
# ------------------------------------------------------------
m_grid <- seq(0, 1, length.out = 900)

df_panel_A <- bind_rows(
  tibble(
    m = m_grid,
    y = cw01_map(m_grid, par_low$beta, par_low$J, h_eff_A),
    case = "BJ < 4"
  ),
  tibble(
    m = m_grid,
    y = cw01_map(m_grid, par_high$beta, par_high$J, h_eff_A),
    case = "BJ > 4"
  ),
  tibble(
    m = m_grid,
    y = m_grid,
    case = "y = m"
  )
)

# At h_eff = 0, m = 0.5 is exactly a fixed point for any beta, J.
m_mid <- 0.5

ann_low <- tibble(
  m = m_mid,
  y = m_mid,
  case = "BJ < 4",
  stability = "Stable fixed point"
)

fp_high <- find_fixed_points_01(par_high$beta, par_high$J, h_eff_A)

# Keep only the outer roots; treat the middle one analytically as exactly 0.5
fp_high_outer <- fp_high[abs(fp_high - 0.5) > 5e-3]

ann_high_outer <- tibble(
  m = fp_high_outer,
  y = fp_high_outer,
  case = "BJ > 4",
  stability = "Stable fixed point"
)

ann_high_mid <- tibble(
  m = m_mid,
  y = m_mid,
  case = "BJ > 4",
  stability = "Unstable fixed point"
)

ann_panel_A <- bind_rows(
  ann_low,
  ann_high_outer,
  ann_high_mid
)

panel_A <- ggplot() +
  geom_line(
    data = df_panel_A %>% filter(case != "y = m"),
    aes(m, y, colour = case),
    linewidth = 1.1
  ) +
  geom_line(
    data = df_panel_A %>% filter(case == "y = m"),
    aes(m, y),
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.9
  ) +
  geom_point(
    data = ann_panel_A %>% filter(stability == "Stable fixed point"),
    aes(m, y, colour = case, shape = stability),
    size = 3.2,
    stroke = 1.0
  ) +
  geom_point(
    data = ann_panel_A %>% filter(stability == "Unstable fixed point"),
    aes(m, y, colour = case, shape = stability),
    size = 3.8,
    stroke = 1.2,
    fill = "white"
  ) +
  scale_colour_manual(
    values = c("BJ < 4" = "#2C7BB6", "BJ > 4" = "#D7191C"),
    breaks = c("BJ < 4", "BJ > 4"),
    labels = c(
      expression(beta * J < 4),
      expression(beta * J > 4)
    )
  ) +
  scale_shape_manual(
    values = c(
      "Stable fixed point"   = 16,
      "Unstable fixed point" = 0
    ),
    breaks = c("Stable fixed point", "Unstable fixed point"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = expression(paste("(A) Fixed-point geometry: ", beta * J < 4, " vs ", beta * J > 4)),
    subtitle = expression(paste("Fixed points satisfy ", m == F(m * ";" * P), ".")),
    x = "Mean burden m",
    y = "Update map F(m;P)",
    colour = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top")

# ------------------------------------------------------------
# 6) Panel B: bifurcation diagram
# ------------------------------------------------------------
bf <- build_bifurcation_table(
  P_grid  = P_grid,
  beta    = par_bif$beta,
  J       = par_bif$J,
  h0      = par_bif$h0,
  gammaP  = par_bif$gammaP
)

# Stable low/high branches only, one row per P
band <- bf %>%
  filter(stable, branch %in% c("low", "high")) %>%
  select(P, branch, m) %>%
  distinct(P, branch, .keep_all = TRUE) %>%
  pivot_wider(names_from = branch, values_from = m) %>%
  drop_na(low, high) %>%
  arrange(P)

panel_B <- ggplot() +
  geom_ribbon(
    data = band,
    aes(x = P, ymin = low, ymax = high),
    fill = "grey80",
    alpha = 0.45
  ) +
  geom_line(
    data = bf %>% filter(!stable, branch == "mid"),
    aes(P, m, linetype = "Unstable equilibrium", group = 1),
    colour = "grey45",
    linewidth = 0.9
  ) +
  geom_line(
    data = bf %>% filter(stable, branch %in% c("low", "high")),
    aes(P, m, linetype = "Stable equilibrium", group = branch),
    colour = "#D7191C",
    linewidth = 1.1
  ) +
  scale_linetype_manual(
    values = c(
      "Stable equilibrium"   = "solid",
      "Unstable equilibrium" = "dashed"
    ),
    name = NULL
  ) +
  labs(
    title = expression("(B) Bifurcation diagram " * m^"*" * "(P)  (" * beta * J > 4 * ")"),
    subtitle = "Grey band: bistable window. Dashed: unstable equilibrium.",
    x = "Slow state P (threshold shift)",
    y = expression("Equilibrium mean burden " * m^"*")
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top")

# ------------------------------------------------------------
# 7) Combine
# ------------------------------------------------------------
fig_bistability_geometry <- (panel_A | panel_B) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold"))
  )

fig_bistability_geometry

# ------------------------------------------------------------
# 8) Save
# ------------------------------------------------------------
# ggsave(
#   filename = "img/Figure_01_bistability_geometry.pdf",
#   plot     = fig_bistability_geometry,
#   width    = 15,
#   height   = 6.5,
#   dpi      = 600
# )