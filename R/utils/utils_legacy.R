# ============================================================
# Legacy tuning utilities
# ============================================================
# Purpose:
#   Retain the older simulation/diagnostics helpers used during
#   parameter tuning and grid search, including the workflow that
#   produced res/final_params.rds.
#
# Note:
#   These functions are kept for reproducibility of the parameter-
#   selection stage. They are not part of the final manuscript
#   simulation pipeline, which uses:
#     - simulate_slowfast_cw01_v3()
#     - run_scenario_diagnostics_v3()
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(furrr)
})

# ------------------------------------------------------------
# Extract bistability info from a bifurcation table
# ------------------------------------------------------------
# Input:
#   bf with columns P, m, stable
#
# Output:
#   compact summary of:
#     - stable lower / upper branches
#     - whether a bistable window exists
#     - fold locations (if present)
extract_bifurcation_info_legacy <- function(bf) {
  
  stab_by_P <- bf |>
    dplyr::filter(stable) |>
    dplyr::group_by(P) |>
    dplyr::summarise(
      n_stable = dplyr::n(),
      m_low    = min(m),
      m_high   = max(m),
      .groups  = "drop"
    )
  
  bistable <- stab_by_P |>
    dplyr::filter(n_stable == 2) |>
    dplyr::arrange(P)
  
  if (nrow(bistable) == 0) {
    list(
      bistable        = FALSE,
      P_low_fold      = NA_real_,
      P_high_fold     = NA_real_,
      stab_by_P       = stab_by_P,
      bistable_window = bistable
    )
  } else {
    list(
      bistable        = TRUE,
      P_low_fold      = min(bistable$P),
      P_high_fold     = max(bistable$P),
      stab_by_P       = stab_by_P,
      bistable_window = bistable
    )
  }
}

# ------------------------------------------------------------
# Interpolate low / high stable branches m_low(P), m_high(P)
# ------------------------------------------------------------
make_branch_interpolators_legacy <- function(stab_by_P) {
  df <- dplyr::arrange(stab_by_P, P)
  
  f_low  <- function(Px) stats::approx(df$P, df$m_low,  xout = Px, rule = 2)$y
  f_high <- function(Px) stats::approx(df$P, df$m_high, xout = Px, rule = 2)$y
  
  list(m_low = f_low, m_high = f_high)
}

# ------------------------------------------------------------
# Assign each simulated point (P_t, m_t) to the nearest stable branch
# ------------------------------------------------------------
assign_branch_state_legacy <- function(sim_df, branch_funs, burn_in = 500) {
  sim_df |>
    dplyr::mutate(
      m_low_hat  = branch_funs$m_low(P),
      m_high_hat = branch_funs$m_high(P),
      state = dplyr::case_when(
        abs(m_high_hat - m_low_hat) < 1e-6 ~ "single",
        abs(m - m_low_hat) <= abs(m - m_high_hat) ~ "low",
        TRUE ~ "high"
      ),
      after_burn = t > burn_in
    )
}

# ------------------------------------------------------------
# Compute switching / dwell metrics from a labeled trajectory
# ------------------------------------------------------------
# Metrics:
#   mean_m
#   mean_P
#   frac_high_abs
#   frac_high_branch
#   n_switch
#   med_high_dwell
switch_metrics_legacy <- function(df_one, burn_in = 500) {
  
  d <- dplyr::filter(df_one, t > burn_in)
  
  frac_high_abs <- mean(d$m > 0.5, na.rm = TRUE)
  
  d2 <- dplyr::filter(d, state %in% c("low", "high"))
  if (nrow(d2) < 2) {
    return(tibble::tibble(
      mean_m           = mean(d$m, na.rm = TRUE),
      mean_P           = mean(d$P, na.rm = TRUE),
      frac_high_abs    = frac_high_abs,
      frac_high_branch = NA_real_,
      n_switch         = 0,
      med_high_dwell   = NA_real_
    ))
  }
  
  st <- d2$state
  n_switch <- sum(st[-1] != st[-length(st)])
  frac_high_branch <- mean(st == "high")
  
  r <- rle(st)
  high_runs <- r$lengths[r$values == "high"]
  med_high_dwell <- ifelse(length(high_runs) == 0, 0, stats::median(high_runs))
  
  tibble::tibble(
    mean_m           = mean(d$m, na.rm = TRUE),
    mean_P           = mean(d$P, na.rm = TRUE),
    frac_high_abs    = frac_high_abs,
    frac_high_branch = frac_high_branch,
    n_switch         = n_switch,
    med_high_dwell   = med_high_dwell
  )
}

# ------------------------------------------------------------
# Parallel helper: simulate multiple groups (optionally in parallel)
# ------------------------------------------------------------
# This helper uses the older no-shock simulator simulate_slowfast_cw01().
run_scenario_parallel_legacy <- function(
    par_fast,
    par_slow,
    P_bases,
    labels,
    T_steps,
    seeds,
    m0s = NULL,
    parallel_groups = TRUE
){
  stopifnot(
    length(P_bases) == length(labels),
    length(labels)  == length(seeds)
  )
  
  if (is.null(m0s)) m0s <- rep(0.1, length(labels))
  stopifnot(length(m0s) == length(labels))
  
  param_list <- purrr::pmap(
    list(P_base = P_bases, label = labels, seed = seeds, m0 = m0s),
    function(P_base, label, seed, m0){
      list(
        args = c(
          list(
            T_steps = T_steps,
            P_base  = P_base,
            P0      = P_base,
            m0      = m0,
            seed    = seed
          ),
          par_fast,
          par_slow
        ),
        group = label
      )
    }
  )
  
  if (isTRUE(parallel_groups)) {
    sim_list <- furrr::future_map(
      param_list,
      function(cfg){
        out <- do.call(simulate_slowfast_cw01, cfg$args)
        out$group <- cfg$group
        out
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    sim_list <- purrr::map(
      param_list,
      function(cfg){
        out <- do.call(simulate_slowfast_cw01, cfg$args)
        out$group <- cfg$group
        out
      }
    )
  }
  
  dplyr::bind_rows(sim_list)
}

# ------------------------------------------------------------
# Legacy wrapper: run scenario diagnostics
# ------------------------------------------------------------
# This is the older pipeline used during tuning:
#   1) compute fast-layer bifurcation structure
#   2) simulate replicates using simulate_slowfast_cw01()
#   3) assign branch labels
#   4) compute switching metrics
#   5) optionally produce overlay / timeseries plots
run_scenario_diagnostics_legacy <- function(
    scenario_name,
    par_fast,
    par_slow,
    P_bases,
    labels,
    T_steps,
    n_rep = 30,
    seed_base = 1000,
    burn_in = 500,
    P_grid = seq(-2, 2, length.out = 400),
    m0s = NULL,
    parallel_reps = TRUE,
    parallel_groups = TRUE
){
  # 1) bifurcation curve for the fast-layer parameterization
  bf <- cw01_bifurcation_curve(
    P_grid = P_grid,
    beta   = par_fast$beta,
    J      = par_fast$J,
    h0     = par_fast$h0,
    gammaP = par_fast$gammaP
  )
  
  bif_info    <- extract_bifurcation_info_legacy(bf)
  branch_funs <- make_branch_interpolators_legacy(bif_info$stab_by_P)
  
  rep_ids <- seq_len(n_rep)
  
  # 2) simulate replicates
  if (isTRUE(parallel_reps)) {
    sim_all <- furrr::future_map_dfr(
      rep_ids,
      function(r){
        seeds <- seed_base + c(2 * r, 2 * r + 1)
        out <- run_scenario_parallel_legacy(
          par_fast = par_fast,
          par_slow = par_slow,
          P_bases  = P_bases,
          labels   = labels,
          T_steps  = T_steps,
          seeds    = seeds,
          m0s      = m0s,
          parallel_groups = parallel_groups
        )
        out$rep <- r
        out
      },
      .progress = TRUE,
      .options  = furrr::furrr_options(seed = TRUE)
    )
  } else {
    sim_all <- purrr::map_dfr(
      rep_ids,
      function(r){
        seeds <- seed_base + c(2 * r, 2 * r + 1)
        out <- run_scenario_parallel_legacy(
          par_fast = par_fast,
          par_slow = par_slow,
          P_bases  = P_bases,
          labels   = labels,
          T_steps  = T_steps,
          seeds    = seeds,
          m0s      = m0s,
          parallel_groups = parallel_groups
        )
        out$rep <- r
        out
      }
    )
  }
  
  sim_all <- dplyr::mutate(sim_all, scenario = scenario_name)
  
  # 3) assign branch labels
  sim_lab <- sim_all |>
    dplyr::group_by(rep, group) |>
    dplyr::group_modify(~ assign_branch_state_legacy(.x, branch_funs, burn_in = burn_in)) |>
    dplyr::ungroup()
  
  # 4) metrics per replicate
  metrics <- sim_lab |>
    dplyr::group_by(scenario, rep, group) |>
    dplyr::group_modify(~ switch_metrics_legacy(.x, burn_in = burn_in)) |>
    dplyr::ungroup()
  
  metrics_summary <- metrics |>
    dplyr::group_by(scenario, group) |>
    dplyr::summarise(
      mean_m           = mean(mean_m),
      mean_P           = mean(mean_P),
      frac_high_abs    = mean(frac_high_abs, na.rm = TRUE),
      frac_high_branch = mean(frac_high_branch, na.rm = TRUE),
      n_switch         = mean(n_switch),
      med_high_dwell   = stats::median(med_high_dwell, na.rm = TRUE),
      .groups = "drop"
    )
  
  # representative replicate (closest mean_m to group-average mean_m)
  rep_pick <- metrics |>
    dplyr::group_by(scenario, group) |>
    dplyr::mutate(target = mean(mean_m)) |>
    dplyr::slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::select(rep, group)
  
  sim_rep <- dplyr::inner_join(sim_all, rep_pick, by = c("rep", "group"))
  
  # 5) optional plots (requires plotting helpers to be sourced)
  p_ts <- plot_timeseries_with_folds(
    sim_rep, bif_info,
    title   = paste0(scenario_name, ": representative time series (thresholds marked)"),
    burn_in = burn_in
  )
  
  p_overlay <- plot_branch_overlay(
    bf, bif_info,
    dplyr::filter(sim_rep, t > burn_in),
    title = paste0(scenario_name, ": trajectory on bifurcation diagram")
  )
  
  list(
    bf              = bf,
    bif_info        = bif_info,
    sim_all         = sim_all,
    sim_labeled     = sim_lab,
    metrics         = metrics,
    metrics_summary = metrics_summary,
    fig_timeseries  = p_ts,
    fig_overlay     = p_overlay
  )
}

# ------------------------------------------------------------
# Backward-compatible aliases
# ------------------------------------------------------------
# These aliases let older scripts keep working after moving the
# tuning stack into utils_legacy.R.

extract_bifurcation_info_old   <- extract_bifurcation_info_legacy
make_branch_interpolators_old  <- make_branch_interpolators_legacy
assign_branch_state_old        <- assign_branch_state_legacy
switch_metrics                 <- switch_metrics_legacy
run_scenario_parallel          <- run_scenario_parallel_legacy
run_scenario_diagnostics_old   <- run_scenario_diagnostics_legacy
