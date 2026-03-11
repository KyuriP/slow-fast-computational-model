## ============================================================
## Figure 1: Geometric intuition for bistability (βJ > 4)
## (0/1 Curie–Weiss mean-field map)
## ============================================================

library(dplyr)
library(ggplot2)
library(purrr)
library(tibble)
library(patchwork)

theme_set(theme_minimal(base_size = 13))

inv_logit <- function(x) 1/(1 + exp(-x))

# Mean-field update map
F_map <- function(m, beta, J, h_eff) {
  inv_logit(beta * (J * (m - 0.5) + h_eff))
}

# Slope of the map at m
Fprime <- function(m, beta, J, h_eff) {
  p <- F_map(m, beta, J, h_eff)
  beta * J * p * (1 - p)
}

# Roots of F(m)-m = 0 in [0,1]
fixed_points_01 <- function(beta, J, h_eff, grid_n = 4001) {
  f <- function(m) F_map(m, beta, J, h_eff) - m
  
  grid <- seq(0, 1, length.out = grid_n)
  vals <- f(grid)
  idx  <- which(diff(sign(vals)) != 0)
  
  roots <- c()
  for (k in idx) {
    r <- try(uniroot(f, c(grid[k], grid[k + 1]))$root, silent = TRUE)
    if (!inherits(r, "try-error")) roots <- c(roots, r)
  }
  sort(unique(roots))
}

# Fixed points for each P + branch labels low/mid/high
fixed_points_for_P <- function(P, beta, J, h0, gammaP, grid_n = 2001) {
  h_eff <- h0 + gammaP * P
  roots <- fixed_points_01(beta, J, h_eff, grid_n = grid_n)
  if (length(roots) == 0) return(tibble())
  
  out <- tibble(
    P = P,
    m = roots,
    stable = Fprime(roots, beta, J, h_eff) < 1
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
  
  out
}

bifurcation_table <- function(P_grid, beta, J, h0, gammaP) {
  purrr::map_dfr(P_grid, ~ fixed_points_for_P(.x, beta, J, h0, gammaP))
}

## ============================================================
## Parameters
## ============================================================

# Panel A: compare βJ<4 vs βJ>4 at same effective bias h_eff=0
h_eff_A <- 0
par_low  <- list(beta = 1.0, J = 3.0)  # βJ = 3 (<4)
par_high <- list(beta = 2.0, J = 3.0)  # βJ = 6 (>4)

# Panel B: bifurcation in P for βJ>4 (match your paper’s B-family)
beta_B  <- 2.0
J_B     <- 3.0
h0_B    <- 0.0
gammaP_B<- 1.0
P_grid  <- seq(-2, 2, length.out = 400)

## ============================================================
## Panel A: intersection geometry
## ============================================================

m_grid <- seq(0, 1, length.out = 900)

dfA <- bind_rows(
  tibble(m = m_grid, y = F_map(m_grid, par_low$beta,  par_low$J,  h_eff_A),
         case = "BJ < 4"),
  tibble(m = m_grid, y = F_map(m_grid, par_high$beta, par_high$J, h_eff_A),
         case = "BJ > 4"),
  tibble(m = m_grid, y = m_grid, case = "y = m")
)

# --- Fixed points for plotting ---
# Key fact: when h_eff = 0, m = 0.5 is EXACTLY a fixed point for any β,J.
m_mid <- 0.5

# Low case: only show the middle fixed point (it is stable when βJ<4)
ann_low <- tibble(
  m = m_mid,
  y = m_mid,
  case = "BJ < 4",
  stability = "Stable fixed point"
)

# High case: compute the outer fixed points numerically + force the exact middle point
fp_high <- fixed_points_01(par_high$beta, par_high$J, h_eff_A)

# remove the numerically-found middle root (near 0.5) and keep outer ones
fp_high_outer <- fp_high[abs(fp_high - 0.5) > 5e-3]

ann_high_outer <- tibble(
  m = fp_high_outer,
  y = fp_high_outer,  # fixed point => y=m
  case = "BJ > 4",
  stability = "Stable fixed point"
)

ann_high_mid <- tibble(
  m = m_mid,
  y = m_mid,
  case = "BJ > 4",
  stability = "Unstable fixed point"
)

# Combine markers. Order matters: draw blue circle, then red square on top.
annA <- bind_rows(
  ann_low,
  ann_high_outer,
  ann_high_mid
)

panel_A <- ggplot() +
  geom_line(
    data = dfA %>% filter(case != "y = m"),
    aes(m, y, colour = case),
    linewidth = 1.1
  ) +
  geom_line(
    data = dfA %>% filter(case == "y = m"),
    aes(m, y),
    linetype = "dashed", colour = "grey40", linewidth = 0.9
  ) +
  # stable points (filled circles)
  geom_point(
    data = annA %>% filter(stability == "Stable fixed point"),
    aes(m, y, colour = case, shape = stability),
    size = 3.2, stroke = 1.0
  ) +
  # unstable point (open square) — only for βJ>4, drawn on top
  geom_point(
    data = annA %>% filter(stability == "Unstable fixed point"),
    aes(m, y, colour = case, shape = stability),
    size = 3.8, stroke = 1.2, fill = "white"
  ) +
  scale_colour_manual(
    values = c("BJ < 4" = "#2C7BB6", "BJ > 4" = "#D7191C"),
    breaks = c("BJ < 4", "BJ > 4"),
    labels = c(
      expression(beta*J < 4),
      expression(beta*J > 4)
    )
  ) +
  scale_shape_manual(
    values = c("Stable fixed point" = 16, "Unstable fixed point" = 0),
    breaks = c("Stable fixed point", "Unstable fixed point"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = expression(paste("(A) Fixed-point geometry: ", beta*J < 4, " vs ", beta*J > 4)),
    subtitle = expression(
      paste("Fixed points satisfy  ", m == F(m~";"~P),".")),
    x = "Mean burden m",
    y = "Update map F(m;P)",
    colour = NULL
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top")

## ============================================================
## Panel B: bifurcation diagram m*(P) for beta*J > 4  (fixed)
## ============================================================

bf <- bifurcation_table(P_grid, beta_B, J_B, h0_B, gammaP_B)


# 1) Make sure stable branches are only low/high and are unique per P
band <- bf %>%
  filter(stable, branch %in% c("low", "high")) %>%
  select(P, branch, m) %>%
  # (important) force exact matching by P; if you see numerical jitter, uncomment:
  # mutate(P = round(P, 6)) %>%
  distinct(P, branch, .keep_all = TRUE) %>%
  pivot_wider(names_from = branch, values_from = m) %>%
  drop_na(low, high) %>%
  arrange(P)

panel_B <- ggplot() +
  geom_ribbon(
    data = band,
    aes(x = P, ymin = low, ymax = high),
    fill = "grey80", alpha = 0.45
  ) +
  geom_line(
    data = bf %>% filter(!stable, branch == "mid"),
    aes(P, m, linetype = "Unstable equilibrium", group = 1),
    colour = "grey45", linewidth = 0.9
  ) +
  geom_line(
    data = bf %>% filter(stable, branch %in% c("low","high")),
    aes(P, m, linetype = "Stable equilibrium", group = branch),
    colour = "#D7191C", linewidth = 1.1
  ) +
  scale_linetype_manual(
    values = c("Stable equilibrium" = "solid",
               "Unstable equilibrium" = "dashed"),
    name = NULL
  ) +
  labs(
    title = expression(paste("(B) Bifurcation diagram ", m^"*", "(P)  (", beta*J, " > 4)")),
    subtitle = "Grey band: bistable window. Dashed: unstable equilibrium.",
    x = "Slow state P (threshold shift)",
    y = expression(paste("Equilibrium mean burden ", m^"*"))
  ) +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top")

panel_B


## ============================================================
## Combine
## ============================================================

figAB <- (panel_A | panel_B) +
  plot_annotation(
    #title = expression(paste("Geometric intuition for the bistability condition ", beta*J, "> 4")),
    theme = theme(plot.title = element_text(face = "bold"))
  )


figAB

# ggsave("Figure_AB_fixedpoints_bifurcation.pdf", figAB, width = 15, height = 6.5, dpi = 600)

