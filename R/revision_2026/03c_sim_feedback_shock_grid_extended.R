# ============================================================
# R/revision_2026/03c_sim_feedback_shock_grid_extended.R
# ============================================================
# Extends the Simulation 2/3 SHOCK design (baseline -> single shock at
# t=0 -> watch recovery) across the same feedback-strength grid used by
# the history-dependence regime check (05_supp_regime_history_
# dependence.R: b = 0, 0.50, 0.75, 1.00, 1.25, 1.50), instead of just
# comparing b=0 vs the locked b=0.50.
#
# Why this exists: the regime check's b=1.0+ results came from a
# DIFFERENT design -- no shock at all, just two chains started directly
# at a low- vs. high-burden initial state, watched for 1500 steps to see
# if they reconverge. That's a real, valid test of history-dependence,
# but a trajectory panel built from it doesn't share Figure 3 panel A's
# "shock at t=0" visual structure, which made it read as unclear/
# unmotivated on its own (no shock marker, unfamiliar framing).
#
# This script instead asks the more directly comparable question: what
# does the SAME shock-and-recovery experiment from Sim 2/3 look like as
# feedback strength b increases past the locked value? At low b the
# shock should still recover (matches Sim 2/3); at high b the same
# perturbation may instead settle into a persistently elevated state
# instead of decaying back to baseline within the window -- shown with
# the same shock marker and visual language as panel A, just with more
# lines (one per b), not a different experimental design.
#
# post_shock_steps is extended from Sim 2/3's 750 to 1500 (matching the
# history-dependence check's window) specifically so a plateau at high b
# has enough time to become visually unambiguous rather than looking
# like "still slowly recovering, just cut off early."
#
# Outputs
# -------
#   res/revision_2026/sim3c/sim3c_raw.rds
#   res/revision_2026/sim3c/sim3c_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-state parameters -- identical to Sim 2/3
# ------------------------------------------------------------------------
kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02
P_base  <- 0

burn_in_steps    <- 200L
# Extended again 2026-08-25 (1500 -> 2500): Figure 3 panel B is dropping
# b=0/0.50 (now redundant with panel A) in favor of showing values further
# into the transition -- but b=1.25 was ALREADY still visibly rising at
# step 1500 in the previous run, not plateaued. A new b=1.30 point at the
# same 1500-step window would almost certainly have the same problem.
# Rather than add a value we can't yet honestly call "resolved," extend
# the window so the higher-b points have a real chance to actually
# plateau within it. If b=1.30 (or 1.50) is STILL rising at step 2500
# when this is rerun, that's a genuine finding to check, not something to
# paper over -- see the end-of-window printout below.
post_shock_steps <- 2500L
shock_time       <- burn_in_steps + 1L
shock_magnitude  <- 1.0

total_steps <- burn_in_steps + post_shock_steps
n_chains    <- 300L

m_star       <- 0.2766389   # Sim 1 middle-condition (P=0) equilibrium mean_m
alpha_smooth <- 0.05

# Grid extended 2026-08-25: added b=0.90 and b=1.10/1.30 so Figure 3 panel
# B can show points strictly BETWEEN the already-established "recovers"
# (b<=0.75) and "clear elevated plateau" (b=1.00) regimes, and one point
# further into the history-dependent regime, without reusing b=0/0.50
# (now panel A's job) or the previously-unresolved b=1.25. Kept b=1.25
# and the original b=0/0.50 in the grid too (not removed) so the CSV
# stays a superset usable for other checks even though panel B's display
# will only pull a subset of these columns.
b_grid <- c(b000 = 0, b050 = 0.50, b075 = 0.75, b090 = 0.90, b100 = 1.00,
            b110 = 1.10, b125 = 1.25, b130 = 1.30, b150 = 1.50)

# ------------------------------------------------------------------------
# One chain: identical shock-and-recovery design to Sim 2/3
# ------------------------------------------------------------------------
run_chain <- function(b) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  m_smooth <- m_star

  P_trace <- numeric(total_steps)
  M_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {
    if (t == shock_time) {
      P <- P + shock_magnitude
    }

    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    P_trace[t] <- P
    M_trace[t] <- symptom_burden(S)

    m_t <- M_trace[t] / N
    m_smooth <- m_smooth + alpha_smooth * (m_t - m_smooth)

    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1) +
      b * (m_smooth - m_star) * dt
  }
  list(P = P_trace, M = M_trace)
}

# ------------------------------------------------------------------------
# Run all (b, chain) combinations in parallel -- same mclapply pattern as
# Sim 3 / Sim 4 (fork-based, macOS/Linux only)
# ------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260825)

tasks <- expand.grid(cond = names(b_grid), chain = seq_len(n_chains), stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("Simulating %d chains x %d steps x %d feedback levels (b = %s), %d tasks total...\n",
            n_chains, total_steps, length(b_grid), paste(b_grid, collapse = ", "), n_tasks))

run_task <- function(i) {
  cond <- tasks$cond[i]
  b <- b_grid[[cond]]
  out <- run_chain(b)
  tibble(b_name = cond, b = b, chain = tasks$chain[i],
         step = seq_len(total_steps), P = out$P, M = out$M)
}

chain_tbls <- parallel::mclapply(seq_len(n_tasks), run_task, mc.cores = n_cores, mc.set.seed = TRUE)

failed <- vapply(chain_tbls, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d parallel tasks failed. First error: %s",
               sum(failed), n_tasks, attr(chain_tbls[[which(failed)[1]]], "condition")$message))
}

traj <- bind_rows(chain_tbls) |>
  mutate(b_name = factor(b_name, levels = names(b_grid)),
         time_since_shock = step - shock_time, m = M / N)

dir.create("res/revision_2026/sim3c", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(traj = traj,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base,
                            b_grid = b_grid, m_star = m_star, alpha_smooth = alpha_smooth,
                            shock_time = shock_time, shock_magnitude = shock_magnitude,
                            burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim3c/sim3c_raw.rds")

summary_tbl <- traj |>
  group_by(b_name, b, step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(M / 9), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim3c/sim3c_summary.csv", row.names = FALSE)

cat("\n=== END-OF-WINDOW (last 50 steps), by b ===\n")
print(summary_tbl |> filter(time_since_shock >= post_shock_steps - 50) |>
        group_by(b) |>
        summarise(mean_P = mean(mean_P), mean_m = mean(mean_m)))

cat("\nCompare end-of-window mean_m/mean_P across b: values near the\n")
cat("pre-shock baseline indicate recovery; values still clearly elevated\n")
cat("after 1500 steps indicate a persistently elevated (non-recovering)\n")
cat("regime at that b.\n")

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim3c/sim3c_raw.rds\n")
cat("  res/revision_2026/sim3c/sim3c_summary.csv\n")
