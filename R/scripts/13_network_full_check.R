# ============================================================
# 13_network_full_check.R
# ============================================================
# Purpose
# -------
# REPLACES 09_network_estimation_check.R as THE network-estimation check
# reported in the main text (Section~sec:network_estimation_check,
# Figures 7-8). Also supersedes 12_network_timescale_check.R, which was
# built as a supplementary appendix check but tested only one axis
# (observation window) while holding the model's diffusion parameter
# sigma_P fixed at its calibrated value -- leaving the two checks (main
# text vs. appendix) run on visibly different scales and inviting the
# "why do these two checks disagree" question. This script folds both
# into one design, so there is a single network-estimation check in the
# paper, entirely inside the real, coupled feedback model.
#
# The old check (09) generated the slow context P as an exogenous,
# per-person covariate with a directly-specified SD -- useful for
# isolating the identification problem cleanly and visibly, but not
# itself a simulation FROM the paper's coupled model. This is what
# Vitor's review specifically flagged:
#
#   "you could simulate from the slow-fast model, fit IsingFit per group
#   (with and without conditioning on P), and characterize what the
#   regression recovers as a function of the timescale separation...
#   For a cross-sectional snapshot, the entire slow-state variance is
#   from between-observation variance. Under long temporal aggregation
#   the field is averaged away and the true coupling is recovered."
#
# Design
# ------
# Same per-person coupled fast/slow simulation as 12
# (simulate_coupled_window()), but sigma_P -- the model's own diffusion
# parameter (Table 2 calibrated value: 0.04) -- is now itself swept, in
# place of directly setting an exogenous SD_P. For a simple
# mean-reverting process with rate kappa, the stationary SD is
# approximately sigma_P / sqrt(2*kappa); at kappa=0.20 this means
# sigma_P=0.632 gives an (near-instantaneous, W=1) stationary spread
# close to the SD_P=1.0 upper end used in the old check. Feedback
# (b=0.15) is on throughout, at every sigma_P level -- not an ablation,
# and not decoupled from the rest of the model.
#
# At each sigma_P level, observation window W is swept exactly as in 12,
# from near-instantaneous (W=1) to long aggregation (W=5000,
# W/tau_P=20). The design is therefore a 2D grid: sigma_P x W, fully
# crossed, at n_dgp independently drawn networks x n_reps replicates.
#
# Outputs
# -------
#   res/network_check/network_full_check_raw.rds
#   res/network_check/network_full_check_summary.csv
#   figs/Figure7_network_estimation_check.pdf / .png   (replaces old Fig 7)
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

if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
suppressPackageStartupMessages(library(matrixStats))

# ------------------------------------------------------------------------
# 0. Design
# ------------------------------------------------------------------------
N            <- 12L
edge_density <- 0.5
w_lo         <- 0.15
w_hi         <- 0.35
h_lo         <- -1.6
h_hi         <- -1.0
gamma_lo     <- 0.8
gamma_hi     <- 1.4

P_base_low  <- -0.3
P_base_high <- 1.0

kappa    <- 0.20
dt       <- 0.02
lambda_m <- 0.001
m_star   <- 0.25
b_real   <- 0.15

tau_P_steps   <- 1 / (kappa * dt)   # = 250 steps
burn_in_steps <- 3000L              # validated sufficient in 11/12

# sigma_P levels: calibrated value (0.04) plus a sweep up to ~0.65,
# chosen so the near-instantaneous (W=1) stationary SD spans roughly the
# same range the old check swept directly (SD_P = 0 to 1.0). Approx.
# stationary SD = sigma_P / sqrt(2*kappa):
#   0.04 -> ~0.06   (the model's real, calibrated value)
#   0.15 -> ~0.24
#   0.35 -> ~0.55
#   0.65 -> ~1.03
sigma_P_levels <- c(0.04, 0.15, 0.35, 0.65)

# window lengths in steps: W/tau_P ~= 0.004, 1, 5, 20
window_steps <- c(1, 250, 1250, 5000)

n_per_group <- 1000L   # validated stable at this scale in 11/12
n_reps      <- 3L
n_dgp       <- 3L

# Quick smoke-test switch: ALWAYS run this first. This design is new
# (sigma_P x W x dgp x rep = up to 144 conditions at full scale) and the
# per-condition cost scales with sigma_P x W in ways not yet timed.
quick_test <- FALSE
if (quick_test) {
  n_reps         <- 1L
  n_dgp          <- 1L
  sigma_P_levels <- c(0.04, 0.65)
  window_steps   <- c(1, 5000)
}

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. Data-generating network (duplicated from 09/11/12, standalone)
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
    logp <- sweep(logp, 2, matrixStats::colMaxs(logp), "-")
    p <- exp(logp); p <- sweep(p, 2, colSums(p), "/")
    cdf <- matrixStats::colCumsums(p)
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
# 2. Coupled per-person simulation with windowed observation. sigma_P is
#    now a function argument (swept), not a fixed constant.
# ------------------------------------------------------------------------
simulate_coupled_window <- function(net_gen, P_base, n_ind, b, sigma_P, burn_in, W_steps) {
  P <- rep(P_base, n_ind)
  m_slow <- rep(m_star, n_ind)
  total_steps <- burn_in + W_steps
  running_sum <- numeric(n_ind)

  for (t in seq_len(total_steps)) {
    theta <- outer(rep(1, n_ind), net_gen$h_vec) + outer(P, net_gen$gamma_vec)
    S <- net_gen$sample_states(theta)
    m_t <- rowMeans(S)
    m_slow <- (1 - lambda_m) * m_slow + lambda_m * m_t
    P <- P + (-kappa * (P - P_base) + b * (m_slow - m_star)) * dt +
      sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t > burn_in) running_sum <- running_sum + P
  }

  Pbar_W <- running_sum / W_steps
  theta_obs <- outer(rep(1, n_ind), net_gen$h_vec) + outer(Pbar_W, net_gen$gamma_vec)
  S_obs <- net_gen$sample_states(theta_obs)
  list(S_obs = S_obs, Pbar_W = Pbar_W)
}

# ------------------------------------------------------------------------
# 3. One (dgp, sigma_P, W, rep) run
# ------------------------------------------------------------------------
run_condition_full <- function(dgp_seed, sigma_P, W_steps, rep) {
  net_gen <- build_dgp(dgp_seed)
  seed <- 900000L * dgp_seed + 5000L * rep + 100L * W_steps + round(sigma_P * 1000)
  set.seed(seed)

  res_low  <- simulate_coupled_window(net_gen, P_base_low,  n_per_group, b_real, sigma_P, burn_in_steps, W_steps)
  res_high <- simulate_coupled_window(net_gen, P_base_high, n_per_group, b_real, sigma_P, burn_in_steps, W_steps)

  Wr_lo <- fit_edges(res_low$S_obs)
  Wr_hi <- fit_edges(res_high$S_obs)
  Wa_lo <- fit_edges(res_low$S_obs,  res_low$Pbar_W)
  Wa_hi <- fit_edges(res_high$S_obs, res_high$Pbar_W)

  tibble(
    dgp = dgp_seed, sigma_P = sigma_P, W_steps = W_steps, rep = rep,
    sd_pbar_lo = sd(res_low$Pbar_W),   sd_pbar_hi = sd(res_high$Pbar_W),
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
# 4. Convergence diagnostic (worst case: largest sigma_P)
# ------------------------------------------------------------------------
check_convergence <- function(net_gen, P_base, b_val, sigma_P, n_ind = 100, steps = burn_in_steps) {
  P <- rep(P_base, n_ind); m_slow <- rep(m_star, n_ind)
  q1 <- floor(steps / 2); q2 <- floor(3 * steps / 4)
  snap_mid <- snap_3q <- NULL
  for (t in seq_len(steps)) {
    theta <- outer(rep(1, n_ind), net_gen$h_vec) + outer(P, net_gen$gamma_vec)
    S <- net_gen$sample_states(theta)
    m_t <- rowMeans(S)
    m_slow <- (1 - lambda_m) * m_slow + lambda_m * m_t
    P <- P + (-kappa * (P - P_base) + b_val * (m_slow - m_star)) * dt +
      sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t == q1) snap_mid <- list(P = mean(P), sd_P = sd(P), m_slow = mean(m_slow))
    if (t == q2) snap_3q  <- list(P = mean(P), sd_P = sd(P), m_slow = mean(m_slow))
  }
  list(midpoint = snap_mid, three_quarter = snap_3q,
       final = list(P = mean(P), sd_P = sd(P), m_slow = mean(m_slow)))
}

net_gen_diag <- build_dgp(1)
cat("CONVERGENCE DIAGNOSTIC (high-baseline, dgp=1, largest sigma_P, burn_in only):\n")
print(check_convergence(net_gen_diag, P_base_high, b_real, sigma_P = max(sigma_P_levels)))
cat("\n")

# ------------------------------------------------------------------------
# 5. Run: n_dgp x sigma_P x window_steps x n_reps
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), sigma_P = sigma_P_levels,
                       W_steps = window_steps, rep = seq_len(n_reps))
cat(sprintf("Design has %d conditions.\n", nrow(design)))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

t0 <- Sys.time()
results <- design |>
  future_pmap_dfr(
    .f = function(dgp_seed, sigma_P, W_steps, rep) run_condition_full(dgp_seed, sigma_P, W_steps, rep),
    .options = furrr_options(seed = TRUE)
  )
t1 <- Sys.time()
cat(sprintf("Design loop wall time: %.1f %s\n", as.numeric(difftime(t1, t0, units = "auto")),
            units(difftime(t1, t0, units = "auto"))))

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
saveRDS(results, "res/network_check/network_full_check_raw.rds")

# ------------------------------------------------------------------------
# 6. Summary + paired gap (raw - adjusted), computed directly from the
#    replicate-level results (paired within dgp/rep, pooled across the
#    low/high-baseline groups since 12 showed both behave similarly)
# ------------------------------------------------------------------------
per_dgp <- results |>
  group_by(dgp, sigma_P, W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    mean
  ), .groups = "drop")

summary_tbl <- per_dgp |>
  group_by(sigma_P, W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop") |>
  mutate(W_over_tauP = W_steps / tau_P_steps)

write.csv(summary_tbl, "res/network_check/network_full_check_summary.csv", row.names = FALSE)
cat("SUMMARY:\n")
print(summary_tbl |> select(sigma_P, W_steps, W_over_tauP, sd_pbar_lo_mean, sd_pbar_hi_mean,
                             gs_raw_hi_mean, gs_adj_hi_mean, true_gs_mean))

gap_tbl <- results |>
  mutate(gap_gs  = ((gs_raw_lo - gs_adj_lo) + (gs_raw_hi - gs_adj_hi)) / 2,
         gap_mae = ((mae_raw_lo - mae_adj_lo) + (mae_raw_hi - mae_adj_hi)) / 2,
         sd_pbar = (sd_pbar_lo + sd_pbar_hi) / 2) |>
  group_by(sigma_P, W_steps) |>
  summarise(
    sd_pbar_mean = mean(sd_pbar), sd_pbar_se = sd(sd_pbar) / sqrt(n()),
    gap_gs_mean  = mean(gap_gs),  gap_gs_se  = sd(gap_gs)  / sqrt(n()),
    gap_mae_mean = mean(gap_mae), gap_mae_se = sd(gap_mae) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(W_over_tauP = W_steps / tau_P_steps,
         sigma_P_lab = factor(sigma_P, levels = sort(unique(sigma_P))))

write.csv(gap_tbl, "res/network_check/network_full_check_gaps.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# 7. Figure (replaces Figure 7): 3 panels, one line per sigma_P level
# ------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 1),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

x_scale <- scale_x_log10()
x_lab <- expression(W / tau[P]~"(observation window / context timescale)")
col_scale <- scale_colour_viridis_d(name = expression(sigma[P]), option = "plasma", end = 0.85)

pA <- ggplot(gap_tbl, aes(W_over_tauP, sd_pbar_mean, colour = sigma_P_lab)) +
  geom_line() + geom_point(size = 2) +
  x_scale + col_scale +
  labs(x = x_lab, y = expression(SD(bar(P)[W])),
       title = expression("(A) Effective "*SD*" of "*bar(P)[W])) +
  theme_pub()

pB <- ggplot(gap_tbl, aes(W_over_tauP, gap_gs_mean, colour = sigma_P_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = gap_gs_mean - gap_gs_se, ymax = gap_gs_mean + gap_gs_se), size = 0.3) +
  x_scale + col_scale +
  labs(x = x_lab, y = "Symptom-only - context-adjusted\n(global strength)",
       title = "(B) Omitted-context bias vs. observation window") +
  theme_pub()

pC <- ggplot(gap_tbl, aes(W_over_tauP, gap_mae_mean, colour = sigma_P_lab)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = gap_mae_mean - gap_mae_se, ymax = gap_mae_mean + gap_mae_se), size = 0.3) +
  x_scale + col_scale +
  labs(x = x_lab,
       y = expression(paste("Symptom-only - context-adjusted  (mean ",
                             group("|", hat(W) - W^{true}, "|"), ")")),
       title = "(C) Recovery-error bias vs. observation window") +
  theme_pub()

fig_full <- (pA | pB | pC) & theme(legend.position = "bottom")
ggsave("figs/Figure7_network_estimation_check.pdf", fig_full, width = 14.5, height = 4.8)
ggsave("figs/Figure7_network_estimation_check.png", fig_full, width = 14.5, height = 4.8, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/network_check/network_full_check_raw.rds\n")
cat("  res/network_check/network_full_check_summary.csv\n")
cat("  res/network_check/network_full_check_gaps.csv\n")
cat("  figs/Figure7_network_estimation_check.pdf (+ .png)\n")
cat("\nquick_test =", quick_test, "-- if TRUE and this ran fast/looks sane,\n")
cat("set quick_test <- FALSE and rerun for the full design.\n")
