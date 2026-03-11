# ============================================================
# Parameter-tuning utilities
# ============================================================
# These functions were used to tune parameters and generate
# res/final_params.rds.
#
# They rely on the legacy tuning stack defined in:
#   - R/utils_legacy.R
#
# They are retained for parameter-search reproducibility, but are
# not part of the final manuscript simulation pipeline.
# ============================================================

source("R/utils_legacy.R")

# ------------------------------------------------------------
# Scalar coalesce: replace NA/NaN with fallback
# ------------------------------------------------------------
`%||%` <- function(x, y) {
  if (length(x) == 0L || is.na(x) || is.nan(x)) y else x
}

# ------------------------------------------------------------
# Scoring rules (operate on one-row summaries)
# ------------------------------------------------------------
score_C1 <- function(ms){
  # prefer: high fraction in high mode, few switches, long high dwell
  (ms$frac_high_abs %||% 0) * 2 -
    (ms$n_switch %||% 0) * 0.05 +
    (ms$med_high_dwell %||% 0) * 0.002
}

score_C2 <- function(ms){
  # prefer: substantial time high AND frequent switching (but not ultra-long dwell)
  (ms$frac_high_abs %||% 0) * 1 +
    (ms$n_switch %||% 0) * 0.08 -
    (ms$med_high_dwell %||% 0) * 0.0005
}

# ------------------------------------------------------------
# Cache bifurcation / branch structure once per fast-layer setup
# ------------------------------------------------------------
make_diag_cache <- function(par_fast, P_grid = seq(-2, 2, length.out = 250)){
  bf <- cw01_bifurcation_curve(
    P_grid = P_grid,
    beta   = par_fast$beta,
    J      = par_fast$J,
    h0     = par_fast$h0,
    gammaP = par_fast$gammaP
  )
  
  bif_info    <- extract_bifurcation_info_legacy(bf)
  branch_funs <- make_branch_interpolators_legacy(bif_info$stab_by_P)
  
  list(
    bf = bf,
    bif_info = bif_info,
    branch_funs = branch_funs
  )
}

# ------------------------------------------------------------
# Lightweight evaluator used inside grid search
# ------------------------------------------------------------
# Returns a summary table (one row per group) suitable for scoring.
run_scenario_diagnostics_cached <- function(
    scenario_name,
    par_fast,
    par_slow,
    P_bases,
    labels,
    T_steps,
    cache,
    n_rep = 4,
    seed_base = 6000,
    burn_in = 200,
    m0s = NULL
){
  stopifnot(length(P_bases) == length(labels))
  
  # Simulate a few replicates quickly (avoid nested parallelism)
  sim_all <- purrr::map_dfr(seq_len(n_rep), function(r){
    seeds <- seed_base + c(2 * r, 2 * r + 1)
    
    out <- run_scenario_parallel_legacy(
      par_fast = par_fast,
      par_slow = par_slow,
      P_bases  = P_bases,
      labels   = labels,
      T_steps  = T_steps,
      seeds    = seeds,
      m0s      = m0s,
      parallel_groups = FALSE
    )
    
    out$rep <- r
    out
  }) |>
    dplyr::mutate(scenario = scenario_name)
  
  # Assign branch state (low / high / single) using cached interpolators
  sim_lab <- sim_all |>
    dplyr::group_by(rep, group) |>
    dplyr::group_modify(~ assign_branch_state_legacy(.x, cache$branch_funs, burn_in = burn_in)) |>
    dplyr::ungroup()
  
  # Compute switching metrics per (rep, group)
  metrics <- sim_lab |>
    dplyr::group_by(scenario, rep, group) |>
    dplyr::group_modify(~ switch_metrics_legacy(.x, burn_in = burn_in)) |>
    dplyr::ungroup()
  
  # Summarise across replicates (one row per group)
  metrics |>
    dplyr::group_by(scenario, group) |>
    dplyr::summarise(
      mean_m         = mean(mean_m, na.rm = TRUE),
      mean_P         = mean(mean_P, na.rm = TRUE),
      frac_high_abs  = mean(frac_high_abs, na.rm = TRUE),
      n_switch       = mean(n_switch, na.rm = TRUE),
      med_high_dwell = median(med_high_dwell, na.rm = TRUE),
      .groups = "drop"
    )
}

# ------------------------------------------------------------
# Grid-search driver
# ------------------------------------------------------------
# Returns one row per grid combination (metrics for the "high" group).
grid_search_C <- function(
    par_fast,
    P_bases,
    labels,
    grid,
    cache,
    scenario_prefix,
    dt = 0.02,
    lambda_m = 0.001,
    T_steps = 2000,
    n_rep   = 4,
    burn_in = 200,
    seed_base = 6000
){
  stopifnot(nrow(grid) >= 1)
  
  grid_rows <- split(grid, seq_len(nrow(grid)))
  
  furrr::future_map_dfr(
    seq_along(grid_rows),
    function(i){
      row <- grid_rows[[i]]
      
      par_slow <- list(
        dt       = dt,
        kappa    = row$kappa,
        b        = row$b,
        m_star   = row$m_star,
        sigmaP   = row$sigmaP,
        lambda_m = lambda_m
      )
      
      ms <- run_scenario_diagnostics_cached(
        scenario_name = paste0(
          scenario_prefix,
          " k=", row$kappa,
          " b=", row$b,
          " s=", row$sigmaP,
          " m*=", row$m_star
        ),
        par_fast  = par_fast,
        par_slow  = par_slow,
        P_bases   = P_bases,
        labels    = labels,
        T_steps   = T_steps,
        cache     = cache,
        n_rep     = n_rep,
        seed_base = seed_base,
        burn_in   = burn_in
      )
      
      # Score using the higher-baseline context group
      ms_hi <- ms |>
        dplyr::filter(group == labels[2]) |>
        dplyr::slice(1)
      
      tibble::tibble(
        kappa = row$kappa,
        b = row$b,
        sigmaP = row$sigmaP,
        m_star = row$m_star,
        mean_m = ms_hi$mean_m,
        mean_P = ms_hi$mean_P,
        frac_high_abs  = ms_hi$frac_high_abs,
        n_switch       = ms_hi$n_switch,
        med_high_dwell = ms_hi$med_high_dwell
      )
    },
    .options  = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )
}

