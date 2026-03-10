# ============================================================
# Diagnostics utilities (branch labeling + switching metrics)
# ============================================================
# Pipeline overview:
#   1) From the fast-layer bifurcation table bf(P,m): find folds and stable branches
#   2) Turn stable branches into interpolators m_low(P), m_high(P)
#   3) Use those interpolators to label each simulated point (P_t, m_t) as:
#        "single" (monostable) or "low"/"high" (bistable branches)
#   4) From labeled trajectories: compute occupancy, switching, and dwell time metrics
#   5) Wrap everything into a scenario runner that simulates replicates and summarizes
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(furrr)
})

# ------------------------------------------------------------
# 1) Bifurcation info: folds + stable branch summaries per P
# ------------------------------------------------------------
# bf is the output of cw01_bifurcation_curve(), with columns:
#   P: context grid value
#   m: equilibrium burden at that P
#   stable: TRUE/FALSE for stability of that equilibrium
#
# Goal: reduce bf to a compact description:
#   - for each P: how many stable equilibria exist (1 vs 2)
#   - stable lower/upper branches m_low(P), m_high(P)
#   - fold locations: where the number of stable equilibria changes
extract_bifurcation_info <- function(bf) {
  stopifnot(all(c("P", "m", "stable") %in% names(bf)))
  
  # For each P, summarize the stable equilibria:
  # - monostable: n_stable = 1 and m_low == m_high
  # - bistable:  n_stable = 2 with distinct low/high equilibria
  stab_by_P <- bf %>%
    filter(stable) %>%
    group_by(P) %>%
    summarise(
      n_stable = n(),
      m_low    = min(m),
      m_high   = max(m),
      .groups  = "drop"
    )
  
  # Bistable window = P values where two stable equilibria coexist
  bistable_window <- stab_by_P %>%
    filter(n_stable == 2) %>%
    arrange(P)
  
  # Fold locations are the endpoints of the bistable window
  if (nrow(bistable_window) == 0) {
    list(
      bistable        = FALSE,
      P_low_fold      = NA_real_,
      P_high_fold     = NA_real_,
      stab_by_P       = stab_by_P,
      bistable_window = bistable_window
    )
  } else {
    list(
      bistable        = TRUE,
      P_low_fold      = min(bistable_window$P),
      P_high_fold     = max(bistable_window$P),
      stab_by_P       = stab_by_P,
      bistable_window = bistable_window
    )
  }
}


# ------------------------------------------------------------
# 1b) Interpolators m_low(P), m_high(P)
# ------------------------------------------------------------
# These allow us to evaluate the *stable* branch values at arbitrary P_t
# in a simulation (P_t is continuous, not restricted to the P_grid used
# to compute bf).
#
# Note: outside the computed P-grid, rule=2 keeps the endpoint values.
make_branch_interpolators <- function(stab_by_P) {
  stopifnot(all(c("P", "m_low", "m_high") %in% names(stab_by_P)))
  
  df <- arrange(stab_by_P, P)
  
  f_low  <- function(Px) approx(df$P, df$m_low,  xout = Px, rule = 2)$y
  f_high <- function(Px) approx(df$P, df$m_high, xout = Px, rule = 2)$y
  
  list(m_low = f_low, m_high = f_high)
}


# ------------------------------------------------------------
# 2) Label each simulated time point: single/low/high
# ------------------------------------------------------------
# Input: a trajectory with columns t, P, m
# Output: same trajectory with:
#   m_low_hat(P_t), m_high_hat(P_t): predicted stable branch values at that P_t
#   state:
#     - "single" if the two branches coincide (monostable at that P)
#     - otherwise "low" or "high" depending on which stable branch is closer to m_t
#
# This gives a branch-based state label that is more robust than a fixed
# threshold (e.g., m > 0.5), because the location of branches shifts with P.
assign_branch_state <- function(sim_df, branch_funs, burn_in = 500) {
  stopifnot(all(c("t", "P", "m") %in% names(sim_df)))
  
  sim_df %>%
    mutate(
      m_low_hat  = branch_funs$m_low(P),
      m_high_hat = branch_funs$m_high(P),
      
      # If branches coincide, we are effectively monostable at that P_t
      state = case_when(
        abs(m_high_hat - m_low_hat) < 1e-6 ~ "single",
        
        # Otherwise assign to the nearest stable branch
        abs(m - m_low_hat) <= abs(m - m_high_hat) ~ "low",
        TRUE ~ "high"
      ),
      
      after_burn = t > burn_in
    )
}


# ------------------------------------------------------------
# 3) Per-replicate metrics from a labeled trajectory
# ------------------------------------------------------------
# We compute three families of summaries:
#   (A) mean_m, mean_P: simple post-burn-in averages
#   (B) frac_high_abs: fraction of post-burn-in time with m > 0.5
#       (a simple absolute-burden summary, not the heatmap hysteresis label)
#   (C) branch-based metrics when state ∈ {low, high}:
#         - frac_high_branch: fraction of branch-labeled time on the high branch
#         - n_switch: number of low <-> high transitions after removing
#                     monostable ("single") periods
#         - med_high_dwell: median run length of consecutive high-branch time
#                           in simulation steps
#
# Edge cases:
#   - If trajectory never visits {low, high} (only "single"), branch metrics
#     are returned as NA/0 where appropriate.

switch_metrics_one <- function(df_one, burn_in = 500) {
  d <- df_one %>% filter(t > burn_in)
  
  if (nrow(d) == 0) {
    return(tibble(
      mean_m = NA_real_,
      mean_P = NA_real_,
      frac_high_abs = NA_real_,
      frac_high_branch = NA_real_,
      n_switch = NA_real_,
      med_high_dwell = NA_real_
    ))
  }
  
  # Absolute (threshold) high occupancy — always defined
  frac_high_abs <- mean(d$m > 0.5, na.rm = TRUE)
  
  # Branch-based metrics require a two-branch region
  d2 <- d %>% filter(state %in% c("low", "high"))
  if (nrow(d2) < 2) {
    return(tibble(
      mean_m = mean(d$m, na.rm = TRUE),
      mean_P = mean(d$P, na.rm = TRUE),
      frac_high_abs = frac_high_abs,
      frac_high_branch = NA_real_,
      n_switch = 0,
      med_high_dwell = NA_real_
    ))
  }
  
  st <- d2$state
  n_switch <- sum(st[-1] != st[-length(st)])
  frac_high_branch <- mean(st == "high")
  
  r <- rle(st)
  high_runs <- r$lengths[r$values == "high"]
  med_high_dwell <- if (length(high_runs) == 0) 0 else stats::median(high_runs)
  
  tibble(
    mean_m = mean(d$m, na.rm = TRUE),
    mean_P = mean(d$P, na.rm = TRUE),
    frac_high_abs = frac_high_abs,
    frac_high_branch = frac_high_branch,
    n_switch = n_switch,
    med_high_dwell = med_high_dwell
  )
}


# ------------------------------------------------------------
# 3b) Apply labeling + metrics to a full simulation table
# ------------------------------------------------------------
# sim_all must contain: t, P, m, scenario, group, rep
# Returns:
#   - sim_labeled: sim_all with state labels added
#   - metrics: per (scenario, group, rep) summary metrics
compute_switch_metrics <- function(sim_all, bif_info, burn_in = 500) {
  req <- c("t", "P", "m", "scenario", "group", "rep")
  stopifnot(all(req %in% names(sim_all)))
  
  branch_funs <- make_branch_interpolators(bif_info$stab_by_P)
  
  sim_labeled <- sim_all %>%
    group_by(scenario, group, rep) %>%
    group_modify(~ assign_branch_state(.x, branch_funs, burn_in = burn_in)) %>%
    ungroup()
  
  metrics <- sim_labeled %>%
    group_by(scenario, group, rep) %>%
    group_modify(~ switch_metrics_one(.x, burn_in = burn_in)) %>%
    ungroup()
  
  list(sim_labeled = sim_labeled, metrics = metrics)
}


# ------------------------------------------------------------
# 4) Scenario runner (v3): run group×rep tasks and compute metrics
# ------------------------------------------------------------
# Purpose:
#   - simulate n_rep replicates for each group baseline P_base
#   - compute fast-layer bifurcation structure for the given par_fast
#   - label trajectories by stable branches and compute metrics
#   - add shock counts (v3 outputs shock_any/exo/endo columns)
#
# Inputs:
#   par_fast: list (beta, J, h0, gammaP, n_nodes, sweeps, ...)
#   par_slow: either:
#       (i) shared list used for all groups, or
#       (ii) named list of lists keyed by group label (per-group slow params)
#   P0_init:
#       NULL -> P0 = P_base for each group
#       numeric -> common initial P0 for all groups (helps comparability)
#   m0_init:
#       common initial symptom burden for all groups/replicates
run_scenario_diagnostics_v3 <- function(
    scenario_name,
    par_fast,
    par_slow,
    P_bases,
    labels,
    T_steps    = 6000,
    n_rep      = 30,
    seed_base  = 12000,
    burn_in    = 500,
    parallel   = TRUE,
    P0_init    = NULL,
    m0_init    = 0.05,
    P_grid     = seq(-2, 2, length.out = 400),
    sim_fn     = simulate_slowfast_cw01_v3
){
  stopifnot(length(P_bases) == length(labels))
  
  # ---- allow per-group slow params via a named list keyed by label ----
  get_ps <- function(label) {
    if (is.list(par_slow) && !is.null(names(par_slow)) && label %in% names(par_slow)) {
      par_slow[[label]]
    } else {
      par_slow
    }
  }
  
  # ---- task table: each (group, rep) is an independent job ----
  group_df <- tibble(group = labels, P_base = P_bases)
  tasks <- crossing(group_df, rep = seq_len(n_rep)) %>%
    mutate(seed = seed_base + row_number())
  
  # ---- simulate one job (one group, one replicate) ----
  run_one <- function(P_base, label, rep_i, seed_i) {
    ps <- get_ps(label)
    
    # Prevent slow param list from overwriting per-task fixed values
    reserved <- c("T_steps", "P_base", "P0", "m0", "seed")
    if (!is.null(ps) && !is.null(names(ps))) ps <- ps[setdiff(names(ps), reserved)]
    
    # Only pass slow params that sim_fn actually accepts
    valid <- names(formals(sim_fn))
    ps <- ps[names(ps) %in% valid]
    
    # Avoid collisions with fixed args or par_fast
    fixed_names <- c("T_steps", "P_base", "P0", "m0", "seed", names(par_fast))
    ps <- ps[setdiff(names(ps), fixed_names)]
    if (any(duplicated(names(ps)))) ps <- ps[!duplicated(names(ps), fromLast = TRUE)]
    
    P0_use <- if (is.null(P0_init)) P_base else P0_init
    
    # Common initial conditions for this task
    sim <- do.call(sim_fn, c(
      list(T_steps = T_steps, P_base = P_base, P0 = P0_use, m0 = m0_init, seed = seed_i),
      par_fast,
      ps
    ))
    
    sim$scenario <- scenario_name
    sim$group    <- label
    sim$rep      <- rep_i
    sim$seed     <- seed_i
    sim
  }
  
  # ---- run all jobs (parallel over group×rep) ----
  if (isTRUE(parallel)) {
    sim_all <- future_pmap_dfr(
      list(tasks$P_base, tasks$group, tasks$rep, tasks$seed),
      run_one,
      .options  = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    sim_all <- pmap_dfr(
      list(tasks$P_base, tasks$group, tasks$rep, tasks$seed),
      run_one
    )
  }
  
  # ---- compute fast-layer equilibrium structure for this par_fast ----
  bf <- cw01_bifurcation_curve(
    P_grid  = P_grid,
    beta    = par_fast$beta,
    J       = par_fast$J,
    h0      = par_fast$h0,
    gammaP  = par_fast$gammaP
  )
  bif_info <- extract_bifurcation_info(bf)
  
  # ---- label trajectories + compute switching metrics ----
  lab_and_metrics <- compute_switch_metrics(sim_all, bif_info = bif_info, burn_in = burn_in)
  sim_labeled <- lab_and_metrics$sim_labeled
  metrics     <- lab_and_metrics$metrics
  
  # ---- shock counts (v3 always outputs shock_any/exo/endo columns) ----
  shock_rep <- sim_all %>%
    filter(t > burn_in) %>%
    group_by(scenario, group, rep) %>%
    summarise(
      n_shock      = sum(shock_any == 1L, na.rm = TRUE),
      n_shock_exo  = sum(abs(shock_exo)  > 0, na.rm = TRUE),
      n_shock_endo = sum(abs(shock_endo) > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  metrics <- metrics %>%
    left_join(shock_rep, by = c("scenario", "group", "rep"))
  
  # ---- summarise over replicates (scenario × group) ----
  metrics_summary <- metrics %>%
    group_by(scenario, group) %>%
    summarise(
      mean_m           = mean(mean_m, na.rm = TRUE),
      mean_P           = mean(mean_P, na.rm = TRUE),
      frac_high_abs    = mean(frac_high_abs, na.rm = TRUE),
      frac_high_branch = mean(frac_high_branch, na.rm = TRUE),
      n_switch         = mean(n_switch, na.rm = TRUE),
      med_high_dwell   = stats::median(med_high_dwell, na.rm = TRUE),
      n_shock          = mean(n_shock, na.rm = TRUE),
      n_shock_exo      = mean(n_shock_exo, na.rm = TRUE),
      n_shock_endo     = mean(n_shock_endo, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- pick representative replicate per group (closest mean_m to group mean) ----
  rep_pick <- metrics %>%
    group_by(scenario, group) %>%
    mutate(target = mean(mean_m, na.rm = TRUE)) %>%
    slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(scenario, group, rep)
  
  sim_rep <- sim_all %>%
    inner_join(rep_pick, by = c("scenario", "group", "rep"))
  
  list(
    sim_all         = sim_all,
    sim_labeled     = sim_labeled,
    sim_rep         = sim_rep,
    bf              = bf,
    bif_info        = bif_info,
    metrics         = metrics,
    metrics_summary = metrics_summary,
    rep_pick        = rep_pick
  )
}














# ------------------------------------------------------------
# Previous version (v1, v2): retained only for record
# and should not be used for final manuscript results.
# ------------------------------------------------------------
run_scenario_diagnostics <- function(
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
  
  bif_info    <- extract_bifurcation_info(bf)
  branch_funs <- make_branch_interpolators(bif_info$stab_by_P)
  
  rep_ids <- seq_len(n_rep)
  
  # 2) simulate replicates
  if (isTRUE(parallel_reps)) {
    sim_all <- furrr::future_map_dfr(
      rep_ids,
      function(r){
        seeds <- seed_base + c(2*r, 2*r + 1)
        out <- run_scenario_parallel(
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
        seeds <- seed_base + c(2*r, 2*r + 1)
        out <- run_scenario_parallel(
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
    dplyr::group_by(scenario, rep, group) |>
    dplyr::group_modify(~ assign_branch_state(.x, branch_funs, burn_in = burn_in)) |>
    dplyr::ungroup()

  
  # 4) metrics per replicate
  metrics <- sim_lab |>
    dplyr::group_by(scenario, rep, group) |>
    dplyr::group_modify(~ switch_metrics(.x, burn_in = burn_in)) |>
    dplyr::ungroup()
  
  metrics_summary <- metrics %>%
    dplyr::group_by(scenario, group) %>%
    dplyr::summarise(
      mean_m           = mean(mean_m, na.rm = TRUE),
      mean_P           = mean(mean_P, na.rm = TRUE),
      frac_high_abs    = mean(frac_high_abs, na.rm = TRUE),
      frac_high_branch = mean(frac_high_branch, na.rm = TRUE),
      n_switch         = mean(n_switch, na.rm = TRUE),
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
  
  sim_rep <- dplyr::inner_join(sim_all, rep_pick, by = c("rep","group"))
  
  # 5) plots
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
    bf            = bf,
    bif_info      = bif_info,
    sim_all       = sim_all,
    sim_labeled   = sim_lab,
    metrics       = metrics,
    metrics_summary = metrics_summary,
    fig_timeseries  = p_ts,
    fig_overlay     = p_overlay
  )
}








# ------------------------------------------------------------
# Scenario runner v2:
# - supports per-group slow params (named list keyed by group label)
# - runs (group x rep) as independent tasks (good for shocks)
# - robust to duplicated/overlapping args in do.call()
# - computes switching metrics + optional shock summaries
# ------------------------------------------------------------
run_scenario_diagnostics_v2 <- function(
    scenario_name,
    par_fast,
    par_slow,
    P_bases,
    labels,
    T_steps    = 6000,
    n_rep      = 30,
    seed_base  = 12000,
    burn_in    = 500,
    parallel   = TRUE,
    P_grid     = seq(-2, 2, length.out = 400)
){
  stopifnot(length(P_bases) == length(labels))
  
  # ---- allow per-group slow params (named list keyed by label) ----
  # par_slow can be:
  #   (1) a single list of params (same for all groups)
  #   (2) a named list of lists keyed by group label
  get_par_slow_for_group <- function(label){
    if (is.list(par_slow) && !is.null(names(par_slow)) && label %in% names(par_slow)) {
      return(par_slow[[label]])
    }
    par_slow
  }
  
  # ---- build task grid: every (group, rep) is a separate job ----
  group_df <- tibble::tibble(
    group  = labels,
    P_base = P_bases
  )
  
  tasks <- tidyr::crossing(group_df, rep = seq_len(n_rep)) |>
    dplyr::mutate(seed = seed_base + dplyr::row_number())
  
  run_one <- function(P_base, label, rep_i, seed_i){
    
    ps <- get_par_slow_for_group(label)
    
    # Prevent slow-params lists from overriding per-task setup
    reserved <- c("T_steps","P_base","P0","m0","seed")
    if (!is.null(ps) && !is.null(names(ps))) {
      ps <- ps[setdiff(names(ps), reserved)]
    }
    
    # ---- CRITICAL: sanitize slow params to avoid do.call collisions ----
    # keep only valid formals
    valid <- names(formals(simulate_slowfast_cw01_v2))
    ps <- ps[names(ps) %in% valid]
    
    # drop any names that would collide with fixed args or par_fast
    fixed_names <- c("T_steps","P_base","P0","m0","seed", names(par_fast))
    ps <- ps[setdiff(names(ps), fixed_names)]
    
    # de-duplicate names defensively
    if (any(duplicated(names(ps)))) ps <- ps[!duplicated(names(ps), fromLast = TRUE)]
    # -------------------------------------------------------------------
    
    sim <- do.call(simulate_slowfast_cw01_v2, c(
      list(
        T_steps = T_steps,
        P_base  = P_base,
        P0      = P_base,
        m0      = 0.05,
        seed    = seed_i
      ),
      par_fast,
      ps
    ))
    
    sim$scenario <- scenario_name
    sim$group    <- label
    sim$rep      <- rep_i
    sim$seed     <- seed_i
    sim
  }
  
  # ---- run sims (parallel over tasks: group x rep) ----
  if (isTRUE(parallel)) {
    sim_all <- furrr::future_pmap_dfr(
      list(
        P_base = tasks$P_base,
        label  = tasks$group,
        rep_i  = tasks$rep,
        seed_i = tasks$seed
      ),
      run_one,
      .options  = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  } else {
    sim_all <- purrr::pmap_dfr(
      list(
        P_base = tasks$P_base,
        label  = tasks$group,
        rep_i  = tasks$rep,
        seed_i = tasks$seed
      ),
      run_one
    )
  }
  
  # ---- fast-layer bifurcation structure ----
  bf <- cw01_bifurcation_curve(
    P_grid  = P_grid,
    beta    = par_fast$beta,
    J       = par_fast$J,
    h0      = par_fast$h0,
    gammaP  = par_fast$gammaP
  )
  bif_info <- extract_bifurcation_info(bf)
  
  # ---- switching/dwell metrics ----
  metrics <- compute_switch_metrics(sim_all, bif_info = bif_info, burn_in = burn_in)
  
  # ---- OPTIONAL: shock diagnostics (works with sim v2 outputs shock_*) ----
  has_shock_any  <- ("shock_any"  %in% names(sim_all))
  has_shock_exo  <- ("shock_exo"  %in% names(sim_all))
  has_shock_endo <- ("shock_endo" %in% names(sim_all))
  
  if (has_shock_any || has_shock_exo || has_shock_endo) {
    
    shock_rep <- sim_all |>
      dplyr::filter(t > burn_in) |>
      dplyr::group_by(scenario, rep, group) |>
      dplyr::summarise(
        n_shock = dplyr::if_else(
          has_shock_any,
          sum(shock_any == 1L, na.rm = TRUE),
          sum(
            (abs(dplyr::coalesce(shock_exo,  0)) > 0) |
              (abs(dplyr::coalesce(shock_endo, 0)) > 0),
            na.rm = TRUE
          )
        ),
        n_shock_exo  = if (has_shock_exo)  sum(abs(shock_exo)  > 0, na.rm = TRUE) else NA_real_,
        n_shock_endo = if (has_shock_endo) sum(abs(shock_endo) > 0, na.rm = TRUE) else NA_real_,
        mean_abs_shock = mean(
          abs(dplyr::coalesce(shock_exo, 0) + dplyr::coalesce(shock_endo, 0)),
          na.rm = TRUE
        ),
        .groups = "drop"
      )
    
    metrics <- metrics |>
      dplyr::left_join(shock_rep, by = c("scenario","rep","group"))
  }
  
  # ---- summary over replicates ----
  metrics_summary <- metrics |>
    dplyr::group_by(scenario, group) |>
    dplyr::summarise(
      mean_m           = mean(mean_m, na.rm = TRUE),
      mean_P           = mean(mean_P, na.rm = TRUE),
      frac_high_abs    = mean(frac_high_abs, na.rm = TRUE),
      frac_high_branch = mean(frac_high_branch, na.rm = TRUE),
      n_switch         = mean(n_switch, na.rm = TRUE),
      med_high_dwell   = stats::median(med_high_dwell, na.rm = TRUE),
      
      n_shock        = if ("n_shock" %in% names(metrics)) mean(n_shock, na.rm = TRUE) else NA_real_,
      n_shock_exo    = if ("n_shock_exo" %in% names(metrics)) mean(n_shock_exo, na.rm = TRUE) else NA_real_,
      n_shock_endo   = if ("n_shock_endo" %in% names(metrics)) mean(n_shock_endo, na.rm = TRUE) else NA_real_,
      mean_abs_shock = if ("mean_abs_shock" %in% names(metrics)) mean(mean_abs_shock, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  list(
    sim_all         = sim_all,
    bf              = bf,
    bif_info        = bif_info,
    metrics         = metrics,
    metrics_summary = metrics_summary
  )
}






