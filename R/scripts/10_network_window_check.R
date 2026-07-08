# ============================================================
# 10_network_window_check.R
# ============================================================
# Purpose
# -------
# Compact SUPPLEMENTARY check responding to Vitor's stronger version of
# the network-estimation-check comment: instead of treating slow-context
# variation as an abstract cross-sectional heterogeneity parameter
# (SD_P, as in 09_network_estimation_check.R), this script asks how
# omitted-context bias in estimated symptom networks depends on the
# OBSERVATION-WINDOW LENGTH relative to the autocorrelation timescale of
# the slow context itself -- i.e. the "dynamic" version of the
# equivalence argument (cf. Kruis & Maris, 2016; Marsman et al., 2018).
#
# Design (per-person P is no longer drawn i.i.d. from N(mu, SD_P); it is
# the time-average of an Ornstein-Uhlenbeck path):
#
#   dP = -kappa*(P - mu) dt + sigma_P * dW      (shocks/feedback OFF)
#
# For a given window length W (in steps), each simulated person's
# covariate is
#
#   Pbar_i = (1/W) * sum_{t=1}^{W} P_i(t)        after a long burn-in
#
# and their symptom vector is generated from the SAME pairwise model as
# 09_network_estimation_check.R, conditional on that person's Pbar_i:
#
#   logit Pr(S_i = 1 | S_-i, Pbar) = h_i + sum_j W_true_ij S_j + gamma_i * Pbar
#
# kappa is fixed at the main model's calibrated value (0.20) so that
# tau_P = 1/kappa has the same meaning as elsewhere in the paper. sigma_P
# is chosen so the STATIONARY sd of P matches the largest SD_P already
# tested in the static check (1.0) -- this makes the static check's
# SD_P = 1.0 condition approximately the W -> 0 limit of this design, so
# the two checks are talking about the same underlying confound, just
# sampled differently. W is then varied on a log scale from far below
# tau_P (near-instantaneous / cross-sectional) to far above it (long
# temporal aggregation), and we ask whether the RAW (omit-context)
# estimator's distortion shrinks toward the CONTEXT-ADJUSTED estimator's
# level as the window lengthens.
#
# IMPORTANT interpretive caution (do not overclaim in the writeup): the
# expected result is that the GAP between raw and adjusted narrows with
# window length, not that the raw estimator fully recovers W_true. A
# preliminary Python prototype at this same scale found the adjusted
# estimator itself stays roughly flat across window lengths (consistent
# with 09's finding that the adjusted estimator is already fairly
# insensitive to the confound level) while the raw estimator's excess
# distortion over the adjusted estimator shrinks by roughly 4-5x between
# W/tau_P ~ 0 and W/tau_P ~ 20. Phrase the manuscript text accordingly:
# "distortion decreases as the window becomes long relative to tau_P",
# not "bias vanishes".
#
# Scope note: this is deliberately a SMALLER run than
# 09_network_estimation_check.R (fewer DGP draws / reps / per-group N).
# This is a compact sensitivity check, not a second headline result --
# see network_estimation_check_writeup.md for the scoping discussion.
#
# Outputs
# -------
#   res/network_check/network_window_check_raw.rds
#   res/network_check/network_window_check_summary.csv
#   figs/Figure_S_network_window_check.pdf / .png   (supplementary)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(future)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------
# 0. Design
# ------------------------------------------------------------------------
# --- pairwise DGP (identical spec to 09_network_estimation_check.R, so
#     the two checks are directly comparable) ---
N            <- 12L
edge_density <- 0.5
w_lo         <- 0.15
w_hi         <- 0.35
h_lo         <- -1.6
h_hi         <- -1.0
gamma_lo     <- 0.8
gamma_hi     <- 1.4

mu_low  <- -0.3
mu_high <- 1.0

# --- OU slow-context process (shocks/feedback off; kappa matches the
#     main model's calibrated value so tau_P is paper-consistent) ---
kappa           <- 0.20
dt              <- 0.02
tau_P_steps     <- 1 / (kappa * dt)          # = 250 steps
sd_P_stationary <- 1.0                        # matches static check's SD_P = 1.0
sigma_P         <- sd_P_stationary * sqrt(2 * kappa)
burn_in_steps   <- 4 * tau_P_steps            # let P reach stationarity before the window starts

# window lengths in steps, chosen to span W/tau_P from ~0 to ~20 -- a
# Python prototype at this scale showed the effect is essentially flat
# for W/tau_P < 1 and only becomes clearly visible out to W/tau_P ~ 10-20,
# so the grid is deliberately log-spaced and extends well past tau_P.
window_steps <- c(1, 25, 100, 250, 500, 1000, 2500, 5000)

n_per_group <- 1000L   # smaller than 09's 2000 -- compact supplementary check
n_reps      <- 8L      # smaller than 09's 15
n_dgp       <- 6L      # smaller than 09's 10

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. Data-generating network (duplicated from 09_network_estimation_check.R
#    rather than sourced, so this script stays standalone-runnable; keep
#    the two copies in sync if the DGP spec above ever changes)
# ------------------------------------------------------------------------
build_dgp <- function(dgp_seed) {
  set.seed(dgp_seed)

  W_true <- matrix(0, N, N)
  pairs <- combn(N, 2)
  n_edges <- round(edge_density * ncol(pairs))
  sel <- sample(ncol(pairs), n_edges)
  for (k in sel) {
    i <- pairs[1, k]; j <- pairs[2, k]
    w <- runif(1, w_lo, w_hi)
    W_true[i, j] <- w; W_true[j, i] <- w
  }
  h_vec <- runif(N, h_lo, h_hi)
  gamma_vec <- runif(N, gamma_lo, gamma_hi)

  wtrue_edges <- W_true[iu]
  true_gs <- sum(abs(wtrue_edges))

  states <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
  storage.mode(states) <- "double"
  quad_base <- 0.5 * rowSums((states %*% W_true) * states)

  sample_states <- function(theta_matrix) {
    n_person <- nrow(theta_matrix)
    linear <- states %*% t(theta_matrix)
    logp <- sweep(linear, 1, quad_base, "+")
    logp <- sweep(logp, 2, apply(logp, 2, max), "-")
    p <- exp(logp); p <- sweep(p, 2, colSums(p), "/")
    cdf <- apply(p, 2, cumsum)
    u <- runif(n_person)
    below <- sweep(cdf, 2, u, FUN = "<")
    idx <- colSums(below) + 1L
    states[idx, , drop = FALSE]
  }

  list(W_true = W_true, h_vec = h_vec, gamma_vec = gamma_vec,
       wtrue_edges = wtrue_edges, true_gs = true_gs,
       sample_states = sample_states)
}

fit_edges <- function(S, Pvec = NULL) {
  beta <- matrix(0, N, N)
  for (i in seq_len(N)) {
    y <- S[, i]
    if (sd(y) == 0) next
    others <- setdiff(seq_len(N), i)
    df <- as.data.frame(S[, others, drop = FALSE])
    names(df) <- paste0("S", others)
    df$y <- y
    rhs <- names(df)[names(df) != "y"]
    if (!is.null(Pvec)) { df$P <- Pvec; rhs <- c(rhs, "P") }
    form <- as.formula(paste("y ~", paste(rhs, collapse = " + ")))
    fit <- suppressWarnings(glm(form, data = df, family = binomial()))
    cf <- coef(fit)
    for (j in others) {
      cname <- paste0("S", j)
      if (cname %in% names(cf) && is.finite(cf[[cname]])) beta[i, j] <- cf[[cname]]
    }
  }
  (beta + t(beta)) / 2
}

# ------------------------------------------------------------------------
# 2. OU window-average generator (vectorized across the n_per_group people
#    in one group -- each person is an independent OU realization)
# ------------------------------------------------------------------------
simulate_ou_window_means <- function(mu, n_ind, W_steps, burn_in) {
  P <- rep(mu, n_ind)
  total_steps <- burn_in + W_steps
  running_sum <- numeric(n_ind)
  for (t in seq_len(total_steps)) {
    P <- P - kappa * dt * (P - mu) + sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t > burn_in) running_sum <- running_sum + P
  }
  running_sum / W_steps
}

# ------------------------------------------------------------------------
# 3. One (dgp, W, rep) condition
# ------------------------------------------------------------------------
run_condition_window <- function(dgp_seed, W_steps, rep) {
  net_gen <- build_dgp(dgp_seed)   # see naming note in 09 re: tibble() column shadowing
  seed <- 100000L * dgp_seed + 1000L * rep + W_steps
  set.seed(seed)

  Pbar_low  <- simulate_ou_window_means(mu_low,  n_per_group, W_steps, burn_in_steps)
  Pbar_high <- simulate_ou_window_means(mu_high, n_per_group, W_steps, burn_in_steps)

  theta_low  <- outer(rep(1, n_per_group), net_gen$h_vec) + outer(Pbar_low,  net_gen$gamma_vec)
  theta_high <- outer(rep(1, n_per_group), net_gen$h_vec) + outer(Pbar_high, net_gen$gamma_vec)

  S_low  <- net_gen$sample_states(theta_low)
  S_high <- net_gen$sample_states(theta_high)

  Wr_lo <- fit_edges(S_low)
  Wr_hi <- fit_edges(S_high)
  Wa_lo <- fit_edges(S_low,  Pbar_low)
  Wa_hi <- fit_edges(S_high, Pbar_high)

  tibble(
    dgp = dgp_seed, W_steps = W_steps, rep = rep,
    sd_pbar_lo = sd(Pbar_low), sd_pbar_hi = sd(Pbar_high),
    gs_raw_lo = sum(abs(Wr_lo[iu])), gs_raw_hi = sum(abs(Wr_hi[iu])),
    gs_adj_lo = sum(abs(Wa_lo[iu])), gs_adj_hi = sum(abs(Wa_hi[iu])),
    mae_raw_lo = mean(abs(Wr_lo[iu] - net_gen$wtrue_edges)),
    mae_raw_hi = mean(abs(Wr_hi[iu] - net_gen$wtrue_edges)),
    mae_adj_lo = mean(abs(Wa_lo[iu] - net_gen$wtrue_edges)),
    mae_adj_hi = mean(abs(Wa_hi[iu] - net_gen$wtrue_edges)),
    true_gs = net_gen$true_gs
  )
}

# ------------------------------------------------------------------------
# 4. Run: n_dgp x window_steps x n_reps conditions, parallelized
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), W_steps = window_steps, rep = seq_len(n_reps))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

results <- design |>
  future_pmap_dfr(
    .f = function(dgp_seed, W_steps, rep) run_condition_window(dgp_seed, W_steps, rep),
    .options = furrr_options(seed = TRUE)
  )

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

saveRDS(results, "res/network_check/network_window_check_raw.rds")

# Same two-stage averaging convention as 09: average within each DGP draw
# first, then across draws.
per_dgp <- results |>
  group_by(dgp, W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    mean
  ), .groups = "drop")

summary_tbl <- per_dgp |>
  group_by(W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop") |>
  mutate(W_over_tauP = W_steps / tau_P_steps)

write.csv(summary_tbl, "res/network_check/network_window_check_summary.csv", row.names = FALSE)
print(summary_tbl |> select(W_steps, W_over_tauP, sd_pbar_lo_mean, sd_pbar_hi_mean,
                             gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean))
cat(sprintf("tau_P = %.0f steps; window range covers W/tau_P = %.3f to %.1f\n",
            tau_P_steps, min(window_steps) / tau_P_steps, max(window_steps) / tau_P_steps))

# ------------------------------------------------------------------------
# 5. Analytic validity check: does simulated sd(Pbar) match OU theory?
# ------------------------------------------------------------------------
# For a stationary OU process with stationary variance sigma_s^2 =
# sigma_P^2 / (2*kappa), the variance of the time-average over a window
# of duration T (continuous time) is a standard result, computed here by
# numerical integration rather than a hand-typed closed form to avoid a
# transcription error:
#   Var(Pbar_T) = (2*sigma_s^2 / T^2) * integral_0^T (T - tau) exp(-kappa*tau) dtau
analytic_sd_window <- function(W_steps_val) {
  Tcont <- W_steps_val * dt
  if (Tcont <= 0) return(sd_P_stationary)
  sigma_s2 <- sigma_P^2 / (2 * kappa)
  integrand <- function(tau) (Tcont - tau) * exp(-kappa * tau)
  val <- stats::integrate(integrand, lower = 0, upper = Tcont)$value
  sqrt((2 * sigma_s2 / Tcont^2) * val)
}
summary_tbl$sd_pbar_theory <- vapply(summary_tbl$W_steps, analytic_sd_window, numeric(1))
cat("\nSANITY CHECK: simulated sd(Pbar) (pooled lo/hi) vs analytic OU window-average sd\n")
print(summary_tbl |>
  transmute(W_steps, W_over_tauP,
            sd_pbar_sim = (sd_pbar_lo_mean + sd_pbar_hi_mean) / 2,
            sd_pbar_theory))

# ------------------------------------------------------------------------
# 6. Figure: 3-panel summary, log10(W/tau_P) x-axis
# ------------------------------------------------------------------------
col_green <- "#2E7D32"
col_red   <- "#C62828"
col_raw   <- "grey15"
col_adj   <- "#1565C0"

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

x_scale <- scale_x_log10()
x_lab <- expression(W / tau[P]~"(observation window / context timescale)")

# Panel A: effective sd(Pbar) vs window length, simulated + analytic overlay
panelA <- summary_tbl |>
  transmute(W_over_tauP,
            sd_pbar_sim = (sd_pbar_lo_mean + sd_pbar_hi_mean) / 2,
            sd_pbar_theory)

pA <- ggplot(panelA, aes(W_over_tauP, sd_pbar_sim)) +
  geom_line(aes(y = sd_pbar_theory), colour = "grey50", linetype = "dotted") +
  geom_line(colour = col_raw) + geom_point(colour = col_raw, size = 2) +
  x_scale +
  labs(x = x_lab, y = expression(SD(bar(P)[W])),
       title = expression("(A) Effective "*SD*" of "*bar(P)[W])) +
  theme_pub()

# Panel B: global strength, raw vs adjusted, averaged over groups (this
# figure is about the raw-vs-adjusted GAP as a function of window length,
# not a group comparison -- see header note on interpretation)
panelB <- summary_tbl |>
  transmute(W_over_tauP,
            gs_raw = (gs_raw_lo_mean + gs_raw_hi_mean) / 2,
            gs_adj = (gs_adj_lo_mean + gs_adj_hi_mean) / 2,
            true_gs = true_gs_mean) |>
  pivot_longer(c(gs_raw, gs_adj), names_to = "estimator", values_to = "gs") |>
  mutate(estimator = if_else(estimator == "gs_raw", "Symptom-only", "Context-adjusted"))

pB <- ggplot(panelB, aes(W_over_tauP, gs, colour = estimator)) +
  geom_line(aes(y = true_gs), colour = "grey50", linetype = "dotted") +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  x_scale +
  labs(x = x_lab, y = "Estimated global strength",
       title = "(B) Global strength vs. observation window", colour = NULL) +
  theme_pub()

# Panel C: recovery error, raw vs adjusted
panelC <- summary_tbl |>
  transmute(W_over_tauP,
            mae_raw = (mae_raw_lo_mean + mae_raw_hi_mean) / 2,
            mae_adj = (mae_adj_lo_mean + mae_adj_hi_mean) / 2) |>
  pivot_longer(c(mae_raw, mae_adj), names_to = "estimator", values_to = "mae") |>
  mutate(estimator = if_else(estimator == "mae_raw", "Symptom-only", "Context-adjusted"))

pC <- ggplot(panelC, aes(W_over_tauP, mae, colour = estimator)) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  x_scale +
  labs(x = x_lab, y = expression("Mean "*group("|", hat(W) - W^{true}, "|")),
       title = expression("(C) Recovery error vs. observation window"),
       colour = NULL) +
  theme_pub()

fig_window <- (pA | pB | pC) +
  plot_layout(guides = "keep") &
  theme(legend.position = "bottom")

ggsave("figs/Figure_S_network_window_check.pdf", fig_window, width = 13.5, height = 4.8)
ggsave("figs/Figure_S_network_window_check.png", fig_window, width = 13.5, height = 4.8, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/network_check/network_window_check_raw.rds\n")
cat("  res/network_check/network_window_check_summary.csv\n")
cat("  figs/Figure_S_network_window_check.pdf (+ .png)  <- supplementary\n")
