# ============================================================
# create_finalparams.R
# ============================================================
# Purpose
# -------
# Reproducible parameter selection for the final B-family scenarios:
#   1) Select a no-shock B0 setting
#   2) Select a both-shocks B+S setting conditional on B0
#   3) Save the final tuning object to res/final_params.rds
#
# Output
# ------
#   res/final_params.rds
#   res/b0_table.rds
#   res/bb_table.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(future)
  library(furrr)
})

source("R/utils/utils_fastlayer.R")
source("R/utils/utils_slowfast.R")
source("R/utils/utils_diagnostics.R")

# ------------------------------------------------------------
# Reproducible parallel RNG
# ------------------------------------------------------------
set.seed(123)
RNGkind("L'Ecuyer-CMRG")
plan(multisession, workers = max(1, parallel::detectCores() - 3))

# dir.create("res", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 1) Fast-layer parameterizations
# ------------------------------------------------------------
par_fast_JBlt4 <- list(
  n_nodes = 100,
  sweeps  = 40,
  beta    = 1.0,
  J       = 2.5,   # beta*J = 2.5 < 4
  h0      = 0.0,
  gammaP  = 1.0
)

par_fast_JBgt4 <- list(
  n_nodes = 100,
  sweeps  = 40,
  beta    = 2.0,
  J       = 3.0,   # beta*J = 6.0 > 4
  h0      = 0.0,
  gammaP  = 1.0
)

# ------------------------------------------------------------
# 2) Shared slow parameters (baseline, no shocks)
# ------------------------------------------------------------
par_slow_shared <- list(
  dt       = 0.02,
  kappa    = 0.25,
  b        = 0.25,
  m_star   = 0.25,
  sigmaP   = 0.05,
  lambda_m = 0.001,
  shock_mode = "none"
)

# ------------------------------------------------------------
# 3) Compute folds for the B-family fast layer
# ------------------------------------------------------------
bf_B <- cw01_bifurcation_curve(
  P_grid  = seq(-2, 2, length.out = 400),
  beta    = par_fast_JBgt4$beta,
  J       = par_fast_JBgt4$J,
  h0      = par_fast_JBgt4$h0,
  gammaP  = par_fast_JBgt4$gammaP
)

bif_info_B <- extract_bifurcation_info(bf_B)

P_low_fold  <- bif_info_B$P_low_fold
P_high_fold <- bif_info_B$P_high_fold

stopifnot(is.finite(P_low_fold), is.finite(P_high_fold))

# ------------------------------------------------------------
# 4) Group placement helper
# ------------------------------------------------------------
# Low baseline group is placed left of the lower fold.
# High baseline group is placed left of the upper fold.
make_groups_BC <- function(red_offset = 0.05, green_offset = 0.30) {
  tibble(
    group  = c("Low baseline context", "High baseline context"),
    P_base = c(P_low_fold - green_offset, P_high_fold - red_offset)
  )
}

# ------------------------------------------------------------
# 5) Search settings
# ------------------------------------------------------------
T_steps_scan <- 6000
n_rep_scan   <- 3
burn_in_scan <- 0
m0_init_scan <- 0.05
P0_init_scan <- NULL

# ------------------------------------------------------------
# 6) B0 grid: placement x slow parameters
# ------------------------------------------------------------
grid_place <- tidyr::crossing(
  green_offset = c(0.20, 0.30, 0.40),
  red_offset   = c(0.03, 0.06, 0.10)
) %>%
  mutate(tag = row_number())

grid_B0 <- tidyr::crossing(
  sigmaP = c(0.02, 0.03, 0.04),
  kappa  = c(0.20, 0.30, 0.40),
  b      = c(0.15, 0.25, 0.35),
  m_star = c(0.25)
) %>%
  mutate(tag = row_number())

score_B0 <- function(ms) {
  hi <- ms %>% filter(group == "High baseline context")
  lo <- ms %>% filter(group == "Low baseline context")
  
  if (nrow(hi) == 0 || nrow(lo) == 0) return(-Inf)
  
  red_hi   <- hi$frac_high_abs
  green_hi <- lo$frac_high_abs
  
  pen_red_extreme <- abs(red_hi - 0.5) * 2
  pen_green       <- 10 * green_hi
  
  1 - pen_red_extreme - pen_green
}

tasks_B0 <- tidyr::crossing(
  place_i = seq_len(nrow(grid_place)),
  b0_i    = seq_len(nrow(grid_B0))
)

run_B0_task <- function(place_i, b0_i) {
  gp <- grid_place[place_i, ]
  g0 <- grid_B0[b0_i, ]
  
  groups <- make_groups_BC(
    red_offset   = gp$red_offset,
    green_offset = gp$green_offset
  )
  
  ps <- modifyList(par_slow_shared, list(
    sigmaP = g0$sigmaP,
    kappa  = g0$kappa,
    b      = g0$b,
    m_star = g0$m_star,
    shock_mode = "none"
  ))
  
  scen_id   <- (place_i - 1L) * nrow(grid_B0) + b0_i
  seed_base <- 50000 + 50L * scen_id
  
  res <- run_scenario_diagnostics_v3(
    scenario_name = paste0("B0_place", gp$tag, "_b0", g0$tag),
    par_fast      = par_fast_JBgt4,
    par_slow      = ps,
    P_bases       = groups$P_base,
    labels        = groups$group,
    T_steps       = T_steps_scan,
    n_rep         = n_rep_scan,
    seed_base     = seed_base,
    burn_in       = burn_in_scan,
    parallel      = FALSE,
    P0_init       = P0_init_scan,
    m0_init       = m0_init_scan
  )
  
  ms <- res$metrics_summary
  sc <- score_B0(ms)
  
  tibble(
    place_tag      = gp$tag,
    b0_tag         = g0$tag,
    red_offset     = gp$red_offset,
    green_offset   = gp$green_offset,
    sigmaP         = g0$sigmaP,
    kappa          = g0$kappa,
    b              = g0$b,
    m_star         = g0$m_star,
    score          = sc,
    red_PrHigh     = ms %>% filter(group == "High baseline context") %>% pull(frac_high_abs),
    green_PrHigh   = ms %>% filter(group == "Low baseline context") %>% pull(frac_high_abs)
  )
}

b0_table <- furrr::future_pmap_dfr(
  list(tasks_B0$place_i, tasks_B0$b0_i),
  run_B0_task,
  .options  = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
)

# saveRDS(b0_table, "res/b0_table.rds")

best_row_B0 <- b0_table %>%
  arrange(desc(score)) %>%
  slice(1)

best_groups <- make_groups_BC(
  red_offset   = best_row_B0$red_offset,
  green_offset = best_row_B0$green_offset
)

best_ps_B0 <- modifyList(par_slow_shared, list(
  sigmaP = best_row_B0$sigmaP,
  kappa  = best_row_B0$kappa,
  b      = best_row_B0$b,
  m_star = best_row_B0$m_star,
  shock_mode = "none"
))

# ------------------------------------------------------------
# 7) B+S grid: both shocks conditional on best B0
# ------------------------------------------------------------
grid_Bboth <- tidyr::crossing(
  lambda_exo_lo = c(0.01, 0.02),
  lambda_exo_hi = c(0.03, 0.05),
  lambda1       = c(0.10, 0.15, 0.20),
  m_crit        = c(0.60, 0.70),
  mu            = c(0.35, 0.45, 0.60)
) %>%
  mutate(tag = row_number())

score_Bboth <- function(ms) {
  hi <- ms %>% filter(group == "High baseline context")
  lo <- ms %>% filter(group == "Low baseline context")
  
  if (nrow(hi) == 0 || nrow(lo) == 0) return(-Inf)
  
  green_hi <- lo$frac_high_abs
  red_hi   <- hi$frac_high_abs
  red_sh   <- hi$n_shock
  
  pen_redlock <- 5 * pmax(0, red_hi - 0.95)
  pen_shocks  <- ifelse(is.na(red_sh), 0, 0.10 * pmax(0, red_sh - 30))
  
  10 * green_hi - pen_redlock - pen_shocks
}

run_Bboth_task <- function(i) {
  r <- grid_Bboth[i, ]
  
  ps_lo <- modifyList(best_ps_B0, list(
    shock_mode    = "both",
    lambda_exo    = r$lambda_exo_lo,
    lambda1       = r$lambda1,
    m_crit        = r$m_crit,
    shock_mu_exo  = r$mu,
    shock_sd_exo  = 0.15 * r$mu,
    shock_mu_endo = r$mu,
    shock_sd_endo = 0.15 * r$mu
  ))
  
  ps_hi <- modifyList(best_ps_B0, list(
    shock_mode    = "both",
    lambda_exo    = r$lambda_exo_hi,
    lambda1       = r$lambda1,
    m_crit        = r$m_crit,
    shock_mu_exo  = r$mu,
    shock_sd_exo  = 0.15 * r$mu,
    shock_mu_endo = r$mu,
    shock_sd_endo = 0.15 * r$mu
  ))
  
  par_slow_by_group <- list(
    "Low baseline context"  = ps_lo,
    "High baseline context" = ps_hi
  )
  
  seed_base <- 70000 + 50L * i
  
  res <- run_scenario_diagnostics_v3(
    scenario_name = paste0("Bboth_", r$tag),
    par_fast      = par_fast_JBgt4,
    par_slow      = par_slow_by_group,
    P_bases       = best_groups$P_base,
    labels        = best_groups$group,
    T_steps       = T_steps_scan,
    n_rep         = n_rep_scan,
    seed_base     = seed_base,
    burn_in       = burn_in_scan,
    parallel      = FALSE,
    P0_init       = P0_init_scan,
    m0_init       = m0_init_scan
  )
  
  ms <- res$metrics_summary
  sc <- score_Bboth(ms)
  
  tibble(
    tag           = r$tag,
    score         = sc,
    lambda_exo_lo = r$lambda_exo_lo,
    lambda_exo_hi = r$lambda_exo_hi,
    lambda1       = r$lambda1,
    m_crit        = r$m_crit,
    mu            = r$mu,
    green_PrHigh  = ms %>% filter(group == "Low baseline context") %>% pull(frac_high_abs),
    red_PrHigh    = ms %>% filter(group == "High baseline context") %>% pull(frac_high_abs),
    red_Shocks    = ms %>% filter(group == "High baseline context") %>% pull(n_shock)
  )
}

bb_table <- furrr::future_map_dfr(
  seq_len(nrow(grid_Bboth)),
  run_Bboth_task,
  .options  = furrr::furrr_options(seed = TRUE),
  .progress = TRUE
)

# saveRDS(bb_table, "res/bb_table.rds")

best_bb <- bb_table %>%
  arrange(desc(score)) %>%
  slice(1)

# ------------------------------------------------------------
# 8) Save final tuning object
# ------------------------------------------------------------
final_params <- list(
  timestamp   = as.character(Sys.time()),
  par_fast    = par_fast_JBgt4,
  best_row_B0 = best_row_B0,
  best_groups = best_groups,
  best_ps_B0  = best_ps_B0,
  best_bb     = best_bb
)

# saveRDS(final_params, "res/final_params.rds")

cat("\nBest B0:\n")
print(best_row_B0)

cat("\nBest B+S:\n")
print(best_bb)