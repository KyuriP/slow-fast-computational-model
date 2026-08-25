# ============================================================
# 02_main_sim.R
# ============================================================
# Purpose
# -------
# Generate the four "final story" scenario objects used in the main Results:
#   A    : βJ < 4, no shocks
#   A+S  : βJ < 4, both shocks (exo + endo)
#   B0   : βJ > 4, no shocks
#   B+S  : βJ > 4, both shocks (exo + endo)
#
# Outputs (RDS)
# -------------
#   res/main/res_A_clean.rds
#   res/main/res_AS_clean.rds
#   res/main/res_B0_clean.rds
#   res/main/res_BS_clean.rds
#
# Design guarantees
# -----------------
# 1) Fast-layer parameters are fixed within regime family (A-family vs B-family).
# 2) Group baselines P_base are fixed across all four scenarios.
# 3) Slow parameters are fixed across all scenarios except for shock toggles.
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
source("R/utils/utils_plotting.R")
source("R/utils/utils_diagnostics.R")
source("R/utils/utils_gridsearch.R")

# ------------------------------------------------------------
# Reproducible parallel RNG
# ------------------------------------------------------------
set.seed(123)
RNGkind("L'Ecuyer-CMRG")
plan(multisession, workers = max(1, parallel::detectCores() - 2))

# ------------------------------------------------------------
# 0) Load final tuned parameters
# ------------------------------------------------------------
final_params <- readRDS("res/tuned/final_params.rds")

# ------------------------------------------------------------
# 1) Canonical parameters held fixed across all scenarios
# ------------------------------------------------------------

# Slow parameters (baseline; shocks off by default)
par_slow_base <- final_params$best_ps_B0
par_slow_base$shock_mode <- "none"

# Group baselines fixed across all scenarios
P_bases_all <- final_params$best_groups$P_base
labels_all  <- final_params$best_groups$group

# Fast-layer parameters: only (beta, J) differ between A and B regimes
fast_common <- list(
  n_nodes = 100,
  sweeps  = 40,
  h0      = 0.0,
  gammaP  = 1.0
)

par_fast_A <- modifyList(fast_common, list(beta = 1.0, J = 2.5))  # βJ = 2.5 < 4
par_fast_B <- modifyList(fast_common, list(beta = 2.0, J = 3.0))  # βJ = 6.0 > 4

# ------------------------------------------------------------
# 2) Shock parameterization (used in A+S and B+S)
# ------------------------------------------------------------
# In the manuscript, the group-specific baseline shock hazard is denoted
# lambda_0. In the code, that baseline hazard is implemented as lambda_exo.

shock_common <- list(
  shock_mode    = "both",
  lambda1       = final_params$best_bb$lambda1,
  m_crit        = final_params$best_bb$m_crit,
  shock_mu_exo  = final_params$best_bb$mu,
  shock_sd_exo  = 0.15 * final_params$best_bb$mu,
  shock_mu_endo = final_params$best_bb$mu,
  shock_sd_endo = 0.15 * final_params$best_bb$mu
)

# Group-specific baseline exogenous hazard
par_slow_shock_by_group <- list(
  "Low baseline context" = modifyList(
    modifyList(par_slow_base, shock_common),
    list(lambda_exo = final_params$best_bb$lambda_exo_lo)
  ),
  "High baseline context" = modifyList(
    modifyList(par_slow_base, shock_common),
    list(lambda_exo = final_params$best_bb$lambda_exo_hi)
  )
)

# ------------------------------------------------------------
# 3) Simulation protocol (main results)
# ------------------------------------------------------------
T_final    <- 8000
burn_final <- 500
nrep_final <- 50

# Common initial conditions used in the final paper runs
P0_init <- 0
m0_init <- 0.05

# ------------------------------------------------------------
# 4) Run the 2 x 2 scenario set
# ------------------------------------------------------------
res_A <- run_scenario_diagnostics_v3(
  scenario_name = "A (βJ<4) no shocks",
  par_fast      = par_fast_A,
  par_slow      = par_slow_base,
  P_bases       = P_bases_all,
  labels        = labels_all,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 20000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

res_AS <- run_scenario_diagnostics_v3(
  scenario_name = "A (βJ<4) both shocks",
  par_fast      = par_fast_A,
  par_slow      = par_slow_shock_by_group,
  P_bases       = P_bases_all,
  labels        = labels_all,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 21000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

res_B0 <- run_scenario_diagnostics_v3(
  scenario_name = "B0 (βJ>4) no shocks",
  par_fast      = par_fast_B,
  par_slow      = par_slow_base,
  P_bases       = P_bases_all,
  labels        = labels_all,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 22000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

res_BS <- run_scenario_diagnostics_v3(
  scenario_name = "B (βJ>4) both shocks",
  par_fast      = par_fast_B,
  par_slow      = par_slow_shock_by_group,
  P_bases       = P_bases_all,
  labels        = labels_all,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 23000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

# ------------------------------------------------------------
# 5) Save artifacts for the manuscript plotting code
# ------------------------------------------------------------
# saveRDS(res_A,  "res/main/res_A_clean.rds")
# saveRDS(res_AS, "res/main/res_AS_clean.rds")
# saveRDS(res_B0, "res/main/res_B0_clean.rds")
# saveRDS(res_BS, "res/main/res_BS_clean.rds")
