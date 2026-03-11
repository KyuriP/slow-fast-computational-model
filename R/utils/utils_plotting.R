# ------------------------------------------------------------
# Plotting helpers for slow–fast scenarios (uniform)
# ------------------------------------------------------------

# ---- Small utility: safely extract fold info ----
folds_to_df <- function(folds) {
  if (is.null(folds)) return(tibble::tibble(P_fold = numeric(0), lab = character(0)))
  
  P_low  <- folds$P_low_fold
  P_high <- folds$P_high_fold
  
  tibble::tibble(
    P_fold = c(P_low, P_high),
    lab    = c("Lower fold", "Upper fold")
  ) %>%
    dplyr::filter(is.finite(P_fold))
}

# ---- Representative replicate selection: one rep per group (closest to group-mean mean_m) ----
pick_representative_reps <- function(metrics_df) {
  # expects columns: group, rep, mean_m (scenario column optional)
  stopifnot(all(c("group", "rep", "mean_m") %in% names(metrics_df)))
  
  # if scenario exists, pick per scenario x group; else per group
  if ("scenario" %in% names(metrics_df)) {
    metrics_df %>%
      dplyr::group_by(scenario, group) %>%
      dplyr::mutate(target = mean(mean_m, na.rm = TRUE)) %>%
      dplyr::slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(scenario, group, rep)
  } else {
    metrics_df %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(target = mean(mean_m, na.rm = TRUE)) %>%
      dplyr::slice_min(abs(mean_m - target), n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(group, rep)
  }
}

# ---- Prep sim df for plotting (optional recode, burn-in, t-window) ----
prep_sim_for_plot <- function(sim_df, burn_in = 0, t_window = NULL, group_map = NULL) {
  df <- sim_df
  
  if (!is.null(group_map)) {
    df <- df %>% dplyr::mutate(group = dplyr::recode(group, !!!group_map))
  }
  
  if (burn_in > 0 && "t" %in% names(df)) {
    df <- df %>% dplyr::filter(t > burn_in)
  }
  
  if (!is.null(t_window) && "t" %in% names(df)) {
    t_max <- max(df$t, na.rm = TRUE)
    df <- df %>% dplyr::filter(t >= (t_max - t_window))
  }
  
  df
}

# ------------------------------------------------------------
# Time-series panel (P(t) and m(t)) with optional fold lines in P(t)
# ------------------------------------------------------------
plot_timeseries_with_folds <- function(
    sim_df,
    folds = NULL,
    title = NULL,
    burn_in = 0,
    t_window = NULL,
    show_burn = FALSE,
    legend_pos = "bottom",
    cols = NULL,
    base_size = 12,
    # ---- NEW: shock lines ----
    show_shocks = TRUE,
    shocks_by_group = TRUE,
    shock_linetype = "dashed",
    shock_linewidth = 0.35,
    shock_alpha = 0.9,
    shock_subtitle = "Dashed vertical lines mark jump shocks in P(t)."
) {
  
  sim_df <- prep_sim_for_plot(sim_df, burn_in = 0, t_window = t_window, group_map = NULL)
  folds_df <- folds_to_df(folds)
  
  base_theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_pos,
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # ---- NEW: compute shock times once ----
  st <- NULL
  has_shock_cols <- any(c("shock_any","shock","shock_exo","shock_endo") %in% names(sim_df))
  if (isTRUE(show_shocks) && has_shock_cols) {
    st <- shock_times(sim_df, burn_in = burn_in, by_group = shocks_by_group)
  }
  
  pP <- ggplot2::ggplot(sim_df, ggplot2::aes(t, P, colour = group)) +
    ggplot2::geom_line(alpha = 0.90, linewidth = 0.6) +
    { if (nrow(folds_df) > 0)
      ggplot2::geom_hline(
        data = folds_df,
        ggplot2::aes(yintercept = P_fold),
        linetype = "dotted",
        colour   = "grey40",
        linewidth = 0.6
        )
    } +
    { if (show_burn && burn_in > 0)
      ggplot2::geom_vline(
        xintercept = burn_in,
        linetype = "dashed",
        colour = "grey50",
        linewidth = 0.6
      )
    } +
    # ---- NEW: shock lines on P(t) ----
  { if (!is.null(st) && nrow(st) > 0)
    ggplot2::geom_vline(
      data = st,
      ggplot2::aes(xintercept = t, colour = if (isTRUE(shocks_by_group)) group else NULL),
      linetype = shock_linetype,
      linewidth = shock_linewidth,
      alpha = shock_alpha,
      show.legend = FALSE
    )
  } +
    ggplot2::labs(title = title, x = NULL, y = "Context  P(t)", colour = NULL) +
    { if (!is.null(cols)) ggplot2::scale_colour_manual(values = cols) } +
    base_theme
  
  pm <- ggplot2::ggplot(sim_df, ggplot2::aes(t, m, colour = group)) +
    ggplot2::geom_line(alpha = 0.90, linewidth = 0.6) +
    { if (show_burn && burn_in > 0)
      ggplot2::geom_vline(
        xintercept = burn_in,
        linetype = "dashed",
        colour = "grey50",
        linewidth = 0.6
      )
    } +
    # ---- NEW: shock lines on m(t) ----
  { if (!is.null(st) && nrow(st) > 0)
    ggplot2::geom_vline(
      data = st,
      ggplot2::aes(xintercept = t, colour = if (isTRUE(shocks_by_group)) group else NULL),
      linetype = shock_linetype,
      linewidth = shock_linewidth,
      alpha = shock_alpha,
      show.legend = FALSE
    )
  } +
    ggplot2::labs(
      x = "Time step",
      y = "Mean symptoms  m(t)",
      colour = NULL,
      subtitle = if (!is.null(st) && nrow(st) > 0) shock_subtitle else NULL
    ) +
    { if (!is.null(cols)) ggplot2::scale_colour_manual(values = cols) } +
    base_theme
  
  patchwork::wrap_plots(pP, pm, ncol = 1, heights = c(1, 1)) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}




# ------------------------------------------------------------
# Bifurcation overlay: equilibrium m*(P) + trajectory points
# ------------------------------------------------------------
plot_branch_overlay <- function(
    bf, bif_info, sim_df, title = NULL,
    burn_in = 0,
    point_alpha = 0.10,
    point_size  = 0.6,
    max_points_per_group = 3000,
    legend_pos = "none",
    cols = NULL,
    ribbon_alpha = 0.12,
    base_size = 12
) {
  
  if (burn_in > 0 && "t" %in% names(sim_df)) {
    sim_df <- dplyr::filter(sim_df, t > burn_in)
  }
  
  sim_plot <- sim_df |>
    dplyr::group_by(group) |>
    dplyr::group_modify(\(.x, .g) {
      n_keep <- min(nrow(.x), max_points_per_group)
      .x[sample.int(nrow(.x), n_keep), , drop = FALSE]
    }) |>
    dplyr::ungroup()
  
  bw <- tibble::tibble(P = numeric(0), m_low = numeric(0), m_high = numeric(0))
  
  if (isTRUE(bif_info$bistable) &&
      !is.null(bif_info$bistable_window) &&
      nrow(bif_info$bistable_window) > 0) {
    
    bw <- bif_info$bistable_window |>
      dplyr::filter(is.finite(P), is.finite(m_low), is.finite(m_high)) |>
      dplyr::arrange(P) |>
      dplyr::group_by(P) |>
      dplyr::summarise(
        m_low  = min(m_low),
        m_high = max(m_high),
        .groups = "drop"
      )
  }
  
  base_theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_pos,
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggplot() +
    { if (nrow(bw) > 0)
      ggplot2::geom_ribbon(
        data = bw,
        ggplot2::aes(x = P, ymin = m_low, ymax = m_high),
        fill = "grey70", alpha = ribbon_alpha
      )
    } +
    ggplot2::geom_line(
      data = dplyr::filter(bf, stable, m < 0.5),
      ggplot2::aes(P, m),
      linewidth = 1.1,
      colour = "grey20"
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(bf, stable, m >= 0.5),
      ggplot2::aes(P, m),
      linewidth = 1.1,
      colour = "grey20"
    ) +
    ggplot2::geom_line(
      data = dplyr::filter(bf, !stable),
      ggplot2::aes(P, m),
      linewidth = 0.9,
      linetype  = "dashed",
      colour    = "grey55"
    ) +
    ggplot2::geom_point(
      data = sim_plot,
      ggplot2::aes(P, m, colour = group),
      alpha = point_alpha,
      size  = point_size
    ) +
    ggplot2::labs(
      title = title,
      subtitle = expression("Equilibrium " * m^"*"*(P) * "; trajectory " * (P[t] * "," * m[t]) * "."),
      x = "Context load  P",
      y = "Mean symptom activation  m",
      colour = NULL
    ) +
    { if (!is.null(cols)) ggplot2::scale_colour_manual(values = cols) } +
    
    # >>> ADD THIS HERE (headroom so points aren't clipped at 1.0) <<<
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    
    base_theme
}



# ------------------------------------------------------------
# Shock-time extractor (pooled or by-group)
# ------------------------------------------------------------
shock_times <- function(sim_df, burn_in = 0, by_group = FALSE) {
  stopifnot("t" %in% names(sim_df))
  
  df <- sim_df |> dplyr::filter(t > burn_in)
  
  # Build a boolean "is_shock" robustly (avoid abs() on non-numeric)
  is_shock <- rep(FALSE, nrow(df))
  
  if ("shock_any" %in% names(df)) {
    is_shock <- is_shock | (df$shock_any == 1L)
  }
  
  if ("shock_exo" %in% names(df) && is.numeric(df$shock_exo)) {
    is_shock <- is_shock | (abs(df$shock_exo) > 0)
  }
  
  if ("shock_endo" %in% names(df) && is.numeric(df$shock_endo)) {
    is_shock <- is_shock | (abs(df$shock_endo) > 0)
  }
  
  if ("shock" %in% names(df) && is.numeric(df$shock)) {
    is_shock <- is_shock | (abs(df$shock) > 0)
  }
  
  df <- df |> dplyr::filter(is_shock)
  
  if (nrow(df) == 0) {
    return(if (by_group) tibble::tibble(group = character(0), t = integer(0))
           else         tibble::tibble(t = integer(0)))
  }
  
  if (by_group) df |> dplyr::distinct(group, t) else df |> dplyr::distinct(t)
}







# ------------------------------------------------------------
# Bifurcation diagram 
# - splits stable into low/high branches to avoid "filled slab"
# - optional fold markers
# ------------------------------------------------------------
plot_bifurcation_diagram <- function(
    bf,
    title    = "Bifurcation diagram of fast Curie–Weiss layer",
    subtitle = "Grey band: bistable region (two stable equilibria)",
    show_ribbon = TRUE,
    show_folds  = TRUE,
    base_size = 12
){
  stopifnot(all(c("P","m","stable") %in% names(bf)))
  
  # stable summaries per P
  stab_by_P <- bf %>%
    dplyr::filter(stable) %>%
    dplyr::group_by(P) %>%
    dplyr::summarise(
      n_stable = dplyr::n(),
      m_low  = min(m),
      m_high = max(m),
      .groups = "drop"
    )
  
  bistable_window <- stab_by_P %>% dplyr::filter(n_stable == 2)
  
  if (nrow(bistable_window) > 0) {
    P_low_fold  <- min(bistable_window$P)
    P_high_fold <- max(bistable_window$P)
  } else {
    P_low_fold  <- NA_real_
    P_high_fold <- NA_real_
  }
  
  fold_df <- tibble::tibble(P_fold = c(P_low_fold, P_high_fold)) %>%
    dplyr::filter(is.finite(P_fold))
  
  p <- ggplot2::ggplot()
  
  # ribbon (optional)
  if (isTRUE(show_ribbon) && nrow(bistable_window) > 0) {
    p <- p + ggplot2::geom_ribbon(
      data = bistable_window,
      ggplot2::aes(x = P, ymin = m_low, ymax = m_high),
      fill = "grey85", alpha = 0.6
    )
  }
  
  # folds (optional)
  if (isTRUE(show_folds) && nrow(fold_df) > 0) {
    p <- p + ggplot2::geom_vline(
      data = fold_df,
      ggplot2::aes(xintercept = P_fold),
      linetype = "dotted", colour = "grey40", linewidth = 0.6
    )
  }
  
  # IMPORTANT: split stable into low/high branches to avoid the “black slab”
  p +
    ggplot2::geom_line(
      data = bf %>% dplyr::filter(!stable),
      ggplot2::aes(P, m),
      linetype = "dashed", colour = "grey60",
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      data = bf %>% dplyr::filter(stable, m < 0.5) %>% dplyr::arrange(P, m),
      ggplot2::aes(P, m),
      linewidth = 1.0,
      color = "skyblue3"
    ) +
    ggplot2::geom_line(
      data = bf %>% dplyr::filter(stable, m >= 0.5) %>% dplyr::arrange(P, m),
      ggplot2::aes(P, m),
      linewidth = 1.0,
      color = "red3"
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x = expression("External context " * P),
      y = expression("Equilibrium mean symptom activation " * m^"*"*(P))
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title.position = "plot",
      panel.grid.minor = ggplot2::element_blank()
    )
}
