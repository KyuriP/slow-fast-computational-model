# ============================================================
# inspect_finalparams.R
# ============================================================
# Purpose
# -------
# Re-run and inspect the best B0 and B+S settings stored in
# res/final_params.rds.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(future)
  library(furrr)
  library(ggplot2)
  library(patchwork)
})

source("R/utils_fastlayer.R")
source("R/utils_slowfast.R")
source("R/utils_plotting.R")
source("R/utils_diagnostics.R")

set.seed(123)
RNGkind("L'Ecuyer-CMRG")
plan(multisession, workers = max(1, parallel::detectCores() - 2))

final_params <- readRDS("res/final_params.rds")

par_fast    <- final_params$par_fast
best_groups <- final_params$best_groups
best_ps_B0  <- final_params$best_ps_B0
best_bb     <- final_params$best_bb

pal <- c(
  "Low baseline context"  = "#7ABD7D",
  "High baseline context" = "#FF5733"
)

make_ps_both <- function(best_ps, best_bb_row, which_group = c("lo", "hi")) {
  which_group <- match.arg(which_group)
  
  lambda_exo <- if (which_group == "lo") {
    best_bb_row$lambda_exo_lo
  } else {
    best_bb_row$lambda_exo_hi
  }
  
  modifyList(best_ps, list(
    shock_mode    = "both",
    lambda_exo    = lambda_exo,
    lambda1       = best_bb_row$lambda1,
    m_crit        = best_bb_row$m_crit,
    shock_mu_exo  = best_bb_row$mu,
    shock_sd_exo  = 0.15 * best_bb_row$mu,
    shock_mu_endo = best_bb_row$mu,
    shock_sd_endo = 0.15 * best_bb_row$mu
  ))
}

T_final    <- 8000
nrep_final <- 50
burn_final <- 50
P0_init    <- 0
m0_init    <- 0.05

res_B0 <- run_scenario_diagnostics_v3(
  scenario_name = "B0_final",
  par_fast      = par_fast,
  par_slow      = best_ps_B0,
  P_bases       = best_groups$P_base,
  labels        = best_groups$group,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 90000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

ps_lo <- make_ps_both(best_ps_B0, best_bb, "lo")
ps_hi <- make_ps_both(best_ps_B0, best_bb, "hi")

par_slow_Bboth <- list(
  "Low baseline context"  = ps_lo,
  "High baseline context" = ps_hi
)

res_Bboth <- run_scenario_diagnostics_v3(
  scenario_name = "Bboth_final",
  par_fast      = par_fast,
  par_slow      = par_slow_Bboth,
  P_bases       = best_groups$P_base,
  labels        = best_groups$group,
  T_steps       = T_final,
  n_rep         = nrep_final,
  seed_base     = 91000,
  burn_in       = burn_final,
  parallel      = TRUE,
  P0_init       = P0_init,
  m0_init       = m0_init
)

rep_pick_B0 <- pick_representative_reps(res_B0$metrics)
sim_B0_rep  <- res_B0$sim_all %>%
  inner_join(rep_pick_B0, by = c("scenario", "group", "rep"))

rep_pick_Bb <- pick_representative_reps(res_Bboth$metrics)
sim_Bb_rep  <- res_Bboth$sim_all %>%
  inner_join(rep_pick_Bb, by = c("scenario", "group", "rep"))

p_ts_B0 <- plot_timeseries_with_folds(
  sim_df   = sim_B0_rep,
  folds    = res_B0$bif_info,
  title    = "B0: representative trajectories",
  burn_in  = burn_final,
  cols     = pal,
  show_shocks = TRUE
)

p_ov_B0 <- plot_branch_overlay(
  bf       = res_B0$bf,
  bif_info = res_B0$bif_info,
  sim_df   = sim_B0_rep,
  title    = "B0: trajectory on bifurcation diagram",
  burn_in  = burn_final,
  cols     = pal
)

p_ts_Bb <- plot_timeseries_with_folds(
  sim_df   = sim_Bb_rep,
  folds    = res_Bboth$bif_info,
  title    = "B+S: representative trajectories",
  burn_in  = burn_final,
  cols     = pal,
  show_shocks = TRUE
)

p_ov_Bb <- plot_branch_overlay(
  bf       = res_Bboth$bf,
  bif_info = res_Bboth$bif_info,
  sim_df   = sim_Bb_rep,
  title    = "B+S: trajectory on bifurcation diagram",
  burn_in  = burn_final,
  cols     = pal
)

print(p_ts_B0)
print(p_ov_B0)
print(p_ts_Bb)
print(p_ov_Bb)

cat("\nB0 summary:\n")
print(res_B0$metrics_summary)

cat("\nB+S summary:\n")
print(res_Bboth$metrics_summary)