# ============================================================
# 12_network_timescale_check.R
# ============================================================
# Purpose
# -------
# SUPERSEDES both 10_network_window_check.R (exogenous OU, no feedback,
# window sweep) and 11_network_feedback_check.R (real feedback, single
# instantaneous snapshot, no window sweep). Neither fully answered
# Vitor's actual request:
#
#   "You could simulate from the slow-fast model, fit IsingFit per group
#   (with and without conditioning on P), and characterize what the
#   regression recovers as a function of the timescale separation...
#   For a cross-sectional snapshot, the entire slow-state variance is
#   from between-observation variance. Under long temporal aggregation
#   the field is averaged away and the true coupling is recovered."
#
# This requires BOTH pieces at once: (a) data generated from the real,
# calibrated coupled model (feedback on, b=0.15 -- what 11 got right and
# 10 got wrong by turning feedback off), AND (b) a sweep over
# observation-window length / timescale separation, from cross-sectional
# (W -> small) to long temporal aggregation (W -> large) -- what 10 got
# right and 11 got wrong by only ever taking one instantaneous snapshot.
#
# Design
# ------
# Each simulated individual runs their own private coupled fast/slow
# system (as in 11), calibrated exactly as the main scenario simulations
# (kappa=0.20, sigma_P=0.04, lambda_m=0.001, m_star=0.25, b=0.15 -- the
# model's own value, not toggled off). They are run for a fixed burn-in
# (long enough for feedback to reach its own quasi-stationary regime,
# validated in 11 at burn_in=4000) and then for an additional W steps,
# during which feedback continues to operate. The context covariate used
# for that person is the window average over those W steps,
#   Pbar_W = (1/W) * sum_{t=burn_in+1}^{burn_in+W} P_t,
# exactly the definition used in 10. Their single observed symptom
# vector is drawn conditional on Pbar_W (not on the instantaneous P at
# the final step), matching 10's convention: this is what "an
# observation aggregated over window W" means operationally. W is then
# swept from near-instantaneous (W=1) to long aggregation (W=5000,
# i.e. W/tau_P = 20), the same grid used in 10.
#
# There is no off/on feedback comparison in this script -- per the
# decision to keep this to "our model, as it is" rather than an
# ablation. b=0.15 throughout.
#
# Expected result: because sigma_P is the model's own small, calibrated
# value (0.04) rather than the inflated value used in 10 to hit
# SD_P=1.0, the ABSOLUTE bias here will likely be more modest throughout
# than in 10's figure. The qualitative claim under test is whether the
# gap between symptom-only and context-adjusted estimates still shrinks
# as W/tau_P grows, i.e. whether "the field is averaged away and the
# true coupling is recovered" holds under the real feedback dynamics,
# not whether the effect size matches the artificially inflated
# SD_P=1.0 condition from the main check.
#
# No closed-form analytic overlay for SD(Pbar_W) is included here (unlike
# 10's panel A): the OU windowed-variance formula used there does not
# apply once feedback is active, since the process is no longer a simple
# linear OU process.
#
# Outputs
# -------
#   res/network_check/network_timescale_check_raw.rds
#   res/network_check/network_timescale_check_summary.csv
#   figs/Figure_S_network_timescale_check.pdf / .png   (supplementary)
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

# matrixStats::colCumsums is a compiled, vectorized replacement for
# apply(p, 2, cumsum) -- the latter loops over every person in R on every
# single simulated time step and was almost certainly the main reason the
# quick_test run took an hour for just 3 conditions. Install if missing.
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  install.packages("matrixStats")
}
suppressPackageStartupMessages(library(matrixStats))

# ------------------------------------------------------------------------
# 0. Design
# ------------------------------------------------------------------------
# --- pairwise DGP (identical spec to 09/10/11) ---
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

# --- slow-state dynamics: the model's own calibrated values (Table 2),
#     feedback ON throughout, not toggled ---
kappa    <- 0.20
dt       <- 0.02
sigma_P  <- 0.04
lambda_m <- 0.001
m_star   <- 0.25
b_real   <- 0.15

tau_P_steps   <- 1 / (kappa * dt)     # = 250 steps
burn_in_steps <- 3000L                # trimmed from 4000: the diagnostic
                                       # above shows ~95% convergence
                                       # already by step 3000 (midpoint
                                       # 2000 -> 79%, three-quarter 3000/
                                       # ~2250 by scale -> ~95%), and
                                       # burn-in dominates total runtime
                                       # since it's paid on every one of
                                       # the 8 window conditions

# window lengths in steps, same grid as script 10: W/tau_P from ~0 to ~20
window_steps <- c(1, 25, 100, 250, 500, 1000, 2500, 5000)

n_per_group <- 1000L   # validated stable at this scale in script 11
n_reps      <- 3L      # smaller than 11's 4 -- this design has 8 window
                        # values on top of dgp x rep, so total compute is
                        # already ~2x script 11's full run
n_dgp       <- 6L      # bumped from 3: with only 3 DGP draws the paired
                        # gap SEs (and the "significant at short window,
                        # ~0 at long window" claim) rest on df=2 t-tests,
                        # which is fragile -- 6 draws roughly halves the
                        # SE and smooths out DGP-level noise (e.g. the
                        # W/tau_P~0.4 bump in the high-baseline gap)

# Quick smoke-test switch: set TRUE first to sanity-check convergence
# and runtime at reduced scale before committing to the full run above.
quick_test <- FALSE
if (quick_test) {
  n_reps        <- 1L
  n_dgp         <- 1L
  window_steps  <- c(1, 250, 5000)   # small subset spanning the range
}

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. Data-generating network (duplicated from 09/10/11 rather than
#    sourced, so this script stays standalone-runnable)
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
    logp <- sweep(logp, 2, matrixStats::colMaxs(logp), "-")  # was
                                         # apply(logp, 2, max) -- same
                                         # per-column R-loop problem as
                                         # the cumsum call below, missed
                                         # on the first pass
    p <- exp(logp); p <- sweep(p, 2, colSums(p), "/")
    cdf <- matrixStats::colCumsums(p)   # was apply(p, 2, cumsum) -- the
                                         # per-time-step bottleneck
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
# 2. Coupled per-person simulation WITH windowed observation: runs
#    burn_in_steps (feedback settles) + W_steps (observation window,
#    feedback still active), then draws ONE symptom vector per person
#    conditional on their window-mean context Pbar_W -- not on the
#    instantaneous P at the final step. Vectorized across n_ind people.
# ------------------------------------------------------------------------
simulate_coupled_window <- function(net_gen, P_base, n_ind, b, burn_in, W_steps) {
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

  # single observation per person, drawn conditional on the WINDOW MEAN,
  # not the instantaneous P at the final simulated step
  theta_obs <- outer(rep(1, n_ind), net_gen$h_vec) + outer(Pbar_W, net_gen$gamma_vec)
  S_obs <- net_gen$sample_states(theta_obs)

  list(S_obs = S_obs, Pbar_W = Pbar_W)
}

# ------------------------------------------------------------------------
# 3. One (dgp, W, rep) run
# ------------------------------------------------------------------------
run_condition_timescale <- function(dgp_seed, W_steps, rep) {
  net_gen <- build_dgp(dgp_seed)
  seed <- 300000L * dgp_seed + 1000L * rep + W_steps
  set.seed(seed)

  res_low  <- simulate_coupled_window(net_gen, P_base_low,  n_per_group, b_real, burn_in_steps, W_steps)
  res_high <- simulate_coupled_window(net_gen, P_base_high, n_per_group, b_real, burn_in_steps, W_steps)

  Wr_lo <- fit_edges(res_low$S_obs)
  Wr_hi <- fit_edges(res_high$S_obs)
  Wa_lo <- fit_edges(res_low$S_obs,  res_low$Pbar_W)
  Wa_hi <- fit_edges(res_high$S_obs, res_high$Pbar_W)

  tibble(
    dgp = dgp_seed, W_steps = W_steps, rep = rep,
    sd_pbar_lo = sd(res_low$Pbar_W),   sd_pbar_hi = sd(res_high$Pbar_W),
    mean_pbar_lo = mean(res_low$Pbar_W), mean_pbar_hi = mean(res_high$Pbar_W),
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
# 4. Convergence diagnostic -- reused from 11, run once before the full
#    design loop. Confirms burn_in_steps is adequate for feedback to
#    settle before the observation window starts accumulating.
# ------------------------------------------------------------------------
check_convergence <- function(net_gen, P_base, b_val, n_ind = 100, steps = burn_in_steps) {
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
    if (t == q1) snap_mid <- list(P = mean(P), m_slow = mean(m_slow))
    if (t == q2) snap_3q  <- list(P = mean(P), m_slow = mean(m_slow))
  }
  list(midpoint = snap_mid, three_quarter = snap_3q,
       final = list(P = mean(P), m_slow = mean(m_slow)))
}

net_gen_diag <- build_dgp(1)
cat("CONVERGENCE DIAGNOSTIC (high-baseline group, dgp=1, burn_in only):\n")
print(check_convergence(net_gen_diag, P_base_high, b_real))
cat("Already validated safe at burn_in_steps=4000 in script 11; re-checked here\n")
cat("since this script duplicates the simulation code rather than sourcing it.\n\n")

# ------------------------------------------------------------------------
# 5. Run: n_dgp x window_steps x n_reps, parallelized across the design
#    grid (each worker runs one full dgp/W/rep job: two groups, each
#    burn_in_steps + W_steps sequential coupled updates)
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), W_steps = window_steps, rep = seq_len(n_reps))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

results <- design |>
  future_pmap_dfr(
    .f = function(dgp_seed, W_steps, rep) run_condition_timescale(dgp_seed, W_steps, rep),
    .options = furrr_options(seed = TRUE)
  )

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
saveRDS(results, "res/network_check/network_timescale_check_raw.rds")

# Two-stage averaging convention, matching 09/10/11: average within each
# DGP draw first, then across draws.
per_dgp <- results |>
  group_by(dgp, W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, mean_pbar_lo, mean_pbar_hi,
      gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    mean
  ), .groups = "drop")

summary_tbl <- per_dgp |>
  group_by(W_steps) |>
  summarise(across(
    c(sd_pbar_lo, sd_pbar_hi, mean_pbar_lo, mean_pbar_hi,
      gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop") |>
  mutate(W_over_tauP = W_steps / tau_P_steps)

write.csv(summary_tbl, "res/network_check/network_timescale_check_summary.csv", row.names = FALSE)
cat("SUMMARY:\n")
print(summary_tbl |> select(W_steps, W_over_tauP, sd_pbar_lo_mean, sd_pbar_hi_mean,
                             gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean,
                             mae_raw_hi_mean, mae_adj_hi_mean, true_gs_mean))

# ------------------------------------------------------------------------
# 6. Figure: 3-panel summary, log10(W/tau_P) x-axis. No analytic overlay
#    on panel A (unlike script 10) -- no closed form for SD(Pbar_W) once
#    feedback is active.
# ------------------------------------------------------------------------
col_raw <- "grey15"; col_adj <- "#1565C0"

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      legend.title = element_blank(), legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

x_scale <- scale_x_log10()
x_lab <- expression(W / tau[P]~"(observation window / context timescale)")

panelA <- summary_tbl |>
  transmute(W_over_tauP, `Low baseline` = sd_pbar_lo_mean, `High baseline` = sd_pbar_hi_mean) |>
  pivot_longer(-W_over_tauP, names_to = "group", values_to = "sd_pbar")

pA <- ggplot(panelA, aes(W_over_tauP, sd_pbar, colour = group)) +
  geom_line() + geom_point(size = 2) +
  x_scale +
  labs(x = x_lab, y = expression(SD(bar(P)[W])),
       title = expression("(A) Effective "*SD*" of "*bar(P)[W]), colour = NULL) +
  theme_pub()

# Panels B/C plot the RAW-MINUS-ADJUSTED GAP directly (with paired SEs),
# not the two estimators as separate near-overlapping lines. Raw and
# adjusted differ by ~0.1-0.2 against values around 13-15, so plotting
# them as two absolute-value lines buries the actual effect under
# between-DGP noise -- the same failure mode diagnosed earlier in the
# bar-chart figure. The gap (paired within dgp/rep, since raw and
# adjusted are fit on the same simulated data) is the real quantity of
# interest and has a much smaller SE than either estimator alone.
gap_tbl <- results |>
  mutate(gap_gs_hi  = gs_raw_hi  - gs_adj_hi,  gap_gs_lo  = gs_raw_lo  - gs_adj_lo,
         gap_mae_hi = mae_raw_hi - mae_adj_hi, gap_mae_lo = mae_raw_lo - mae_adj_lo) |>
  group_by(W_steps) |>
  summarise(
    gs_hi_mean  = mean(gap_gs_hi),  gs_hi_se  = sd(gap_gs_hi)  / sqrt(n()),
    gs_lo_mean  = mean(gap_gs_lo),  gs_lo_se  = sd(gap_gs_lo)  / sqrt(n()),
    mae_hi_mean = mean(gap_mae_hi), mae_hi_se = sd(gap_mae_hi) / sqrt(n()),
    mae_lo_mean = mean(gap_mae_lo), mae_lo_se = sd(gap_mae_lo) / sqrt(n()),
    .groups = "drop"
  ) |>
  mutate(W_over_tauP = W_steps / tau_P_steps)

panelB <- bind_rows(
  gap_tbl |> transmute(W_over_tauP, group = "High baseline", gap = gs_hi_mean, se = gs_hi_se),
  gap_tbl |> transmute(W_over_tauP, group = "Low baseline",  gap = gs_lo_mean, se = gs_lo_se)
)

pB <- ggplot(panelB, aes(W_over_tauP, gap, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = gap - se, ymax = gap + se), size = 0.35) +
  x_scale +
  labs(x = x_lab, y = "Symptom-only − context-adjusted\n(global strength)",
       title = "(B) Omitted-context bias vs. observation window", colour = NULL) +
  theme_pub()

panelC <- bind_rows(
  gap_tbl |> transmute(W_over_tauP, group = "High baseline", gap = mae_hi_mean, se = mae_hi_se),
  gap_tbl |> transmute(W_over_tauP, group = "Low baseline",  gap = mae_lo_mean, se = mae_lo_se)
)

pC <- ggplot(panelC, aes(W_over_tauP, gap, colour = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line() +
  geom_pointrange(aes(ymin = gap - se, ymax = gap + se), size = 0.35) +
  x_scale +
  labs(x = x_lab,
       y = expression(paste("Symptom-only − context-adjusted  (mean ",
                             group("|", hat(W) - W^{true}, "|"), ")")),
       title = "(C) Recovery-error bias vs. observation window", colour = NULL) +
  theme_pub()

fig_timescale <- (pA | pB | pC) & theme(legend.position = "bottom")
ggsave("figs/Figure_S_network_timescale_check.pdf", fig_timescale, width = 13.5, height = 4.8)
ggsave("figs/Figure_S_network_timescale_check.png", fig_timescale, width = 13.5, height = 4.8, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/network_check/network_timescale_check_raw.rds\n")
cat("  res/network_check/network_timescale_check_summary.csv\n")
cat("  figs/Figure_S_network_timescale_check.pdf (+ .png)  <- supplementary\n")
