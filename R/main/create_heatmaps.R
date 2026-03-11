# ============================================================
# create_heatmaps.R
# ============================================================
# Purpose
# -------
# Generate no-shock regime heatmaps over (betaJ, P_base) using
# hysteresis-style state labeling with memory.
#
# Output
# ------
#   res_clean/heat_zoom_noshock_hyst2.rds
#
# Notes
# -----
# - This script uses the final no-shock slow parameters from
#   res/final_params.rds.
# - Heatmap summaries are based on a separate hysteresis-style
#   binary state labeling rule:
#       state = 0 if m <= low
#       state = 1 if m >= high
#       otherwise retain previous state
# - This labeling is used only for these heatmaps, not for the
#   branch-based switching summaries reported in the main scenario
#   diagnostics.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(future)
  library(furrr)
  library(ggplot2)
})

source("R/utils_fastlayer.R")
source("R/utils_slowfast.R")
source("R/utils_diagnostics.R")

# ------------------------------------------------------------
# Reproducible parallel RNG
# ------------------------------------------------------------
set.seed(123)
RNGkind("L'Ecuyer-CMRG")
plan(multisession, workers = max(1, parallel::detectCores() - 2))

dir.create("res_clean", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) Grids
# ------------------------------------------------------------
P_vals  <- seq(-0.30, 0.30, by = 0.01)
JB_vals <- seq(0.50, 6.50, by = 0.10)

labels_P <- sprintf("P=%+.2f", P_vals)
group_df <- tibble(
  group  = labels_P,
  P_base = P_vals
)

# ------------------------------------------------------------
# 2) Run settings
# ------------------------------------------------------------
T_steps   <- 5000
burn_in   <- 50
n_rep     <- 3
seed_base <- 20000
P0_init   <- NULL
m0_init   <- 0.05

# ------------------------------------------------------------
# 3) Baseline no-shock parameters
# ------------------------------------------------------------
final_params  <- readRDS("res/final_params.rds")

par_slow_base <- final_params$best_ps_B0
par_slow_base$shock_mode <- "none"

fast_common <- list(
  n_nodes = 100,
  sweeps  = 40,
  h0      = 0.0,
  gammaP  = 1.0
)

# ------------------------------------------------------------
# 4) Hysteresis-style episode metrics
# ------------------------------------------------------------
# State labeling:
#   - state = 0 if m <= low
#   - state = 1 if m >= high
#   - otherwise keep previous state
episode_metrics_hyst2 <- function(m, low = 0.30, high = 0.70) {
  state <- rep(NA_integer_, length(m))
  state[m <= low]  <- 0L
  state[m >= high] <- 1L
  
  init_state <- ifelse(m[1] >= (low + high) / 2, 1L, 0L)
  
  for (i in seq_along(state)) {
    if (is.na(state[i])) {
      state[i] <- if (i == 1) init_state else state[i - 1]
    }
  }
  
  r <- rle(state)
  switches <- max(0L, length(r$values) - 1L)
  
  high_lens <- r$lengths[r$values == 1L]
  entered_high <- length(high_lens) > 0
  
  tibble(
    mean_m            = mean(m, na.rm = TRUE),
    frac_high         = mean(state == 1L, na.rm = TRUE),
    switches          = switches,
    entered_high      = as.integer(entered_high),
    n_high_episodes   = length(high_lens),
    mean_dur_high     = if (entered_high) mean(high_lens) else 0,
    mean_dur_high_cond = if (entered_high) mean(high_lens) else NA_real_
  )
}

# ------------------------------------------------------------
# 5) Run one betaJ slice
# ------------------------------------------------------------
# beta is fixed at 1.0, so betaJ = J
run_one_JB_noshock <- function(betaJ) {
  par_fast <- modifyList(fast_common, list(
    beta = 1.0,
    J    = betaJ
  ))
  
  res <- run_scenario_diagnostics_v3(
    scenario_name = paste0("HM_noshock_betaJ_", sprintf("%.2f", betaJ)),
    par_fast      = par_fast,
    par_slow      = par_slow_base,
    P_bases       = P_vals,
    labels        = labels_P,
    T_steps       = T_steps,
    n_rep         = n_rep,
    seed_base     = seed_base,
    burn_in       = burn_in,
    parallel      = TRUE,
    P0_init       = P0_init,
    m0_init       = m0_init
  )
  
  per_rep <- res$sim_all %>%
    filter(t > burn_in) %>%
    group_by(group, rep) %>%
    group_modify(~ episode_metrics_hyst2(.x$m, low = 0.30, high = 0.70)) %>%
    ungroup() %>%
    left_join(group_df, by = "group") %>%
    mutate(betaJ = betaJ)
  
  per_rep %>%
    group_by(betaJ, P_base) %>%
    summarise(
      Mean_m            = mean(mean_m, na.rm = TRUE),
      Pr_high           = mean(frac_high, na.rm = TRUE),
      Switches          = mean(switches, na.rm = TRUE),
      EnteredHigh       = mean(entered_high, na.rm = TRUE),
      N_high_episodes   = mean(n_high_episodes, na.rm = TRUE),
      MeanHighDur       = mean(mean_dur_high, na.rm = TRUE),
      MeanHighDur_cond  = mean(mean_dur_high_cond, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# 6) Run full grid
# ------------------------------------------------------------
# Outer loop stays sequential; parallelism happens inside
# run_scenario_diagnostics_v3().
heat_zoom <- purrr::map_dfr(
  seq_along(JB_vals),
  function(i) {
    cat(sprintf("\n[%d/%d] Running betaJ = %.2f\n", i, length(JB_vals), JB_vals[i]))
    run_one_JB_noshock(JB_vals[i])
  },
  .progress = TRUE
)

saveRDS(heat_zoom, "res_clean/heat_zoom_noshock_hyst2.rds")

# ------------------------------------------------------------
# 7) Sanity checks
# ------------------------------------------------------------
stopifnot(all(c("betaJ", "P_base", "Pr_high", "Switches", "MeanHighDur") %in% names(heat_zoom)))

print(summary(heat_zoom$Pr_high))
print(summary(heat_zoom$Switches))

# ------------------------------------------------------------
# 8) Quick plots
# ------------------------------------------------------------
p_frac <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = Pr_high)) +
  geom_tile() +
  geom_contour(
    aes(z = Pr_high),
    breaks = 0.5,
    colour = "black",
    linewidth = 0.4
  ) +
  labs(
    x = expression(beta * J),
    y = expression(P[base]),
    fill = "Pr(high)\n(hysteresis)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

p_sw <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = Switches)) +
  geom_tile() +
  labs(
    x = expression(beta * J),
    y = expression(P[base]),
    fill = "Switches\n(hysteresis)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

p_dur <- ggplot(heat_zoom, aes(x = betaJ, y = P_base, fill = MeanHighDur)) +
  geom_tile() +
  labs(
    x = expression(beta * J),
    y = expression(P[base]),
    fill = "Mean high\nduration"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p_frac)
print(p_sw)
print(p_dur)