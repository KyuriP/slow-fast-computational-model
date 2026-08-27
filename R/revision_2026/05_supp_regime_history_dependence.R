# ============================================================
# R/revision_2026/05_supp_regime_history_dependence.R
# ============================================================
# Supplementary regime check:
# Can the revised 0/1 heterogeneous symptom-level slow-fast model show
# tipping-like or history-dependent behavior under stronger feedback?
#
# This is NOT part of the main four-simulation sequence.
# It is a robustness / regime check replacing the old homogeneous
# mean-field bistability figure.
#
# Main question:
#   If we start the same system from a low-burden state versus a
#   high-burden state, do trajectories converge to the same late state,
#   or do they remain separated for a long time under stronger feedback?
#
# Interpretation:
#   - If late low-start and high-start states overlap: no evidence of
#     history dependence in that parameter regime.
#   - If they remain separated: evidence for tipping-like/metastable
#     behavior in the revised symptom-level model.
#   - If P_t runs away to extreme values: call it runaway/saturation,
#     not bistability.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("R/revision_2026/00_parameters_uncentered01.R")
source("R/revision_2026/utils_uncentered01_model.R")

set.seed(20260825)

# ------------------------------------------------------------------------
# 0. Parameters
# ------------------------------------------------------------------------
N <- length(tau)

# Use the Simulation 1 middle-context mean active-symptom fraction as m_star.
# Prefer reading from the actual Sim 1 summary if available.
sim1_summary_path <- "res/revision_2026/sim1/sim1_summary.csv"
if (file.exists(sim1_summary_path)) {
  sim1_summary <- read_csv(sim1_summary_path, show_col_types = FALSE)
  m_star <- sim1_summary |>
    filter(condition == "middle") |>
    pull(mean_m)
  if (length(m_star) != 1 || is.na(m_star)) {
    stop("Could not extract middle-condition mean_m from sim1_summary.csv")
  }
} else {
  # fallback from earlier locked Simulation 1 result
  m_star <- 0.28
}

# Slow-field parameters.
# Keep kappa, sigma_P, and dt close to Simulations 2-3.
P_base  <- 0
kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02

# Smoothing for symptom burden entering feedback.
alpha_smooth <- 0.05

# Feedback grid.
# Current main simulation uses b = 0.50.
# Here we test stronger regimes, but we do not assume they are clinically calibrated.
# Densified from the original 6-point grid (0, .5, .75, 1, 1.25, 1.5) to a
# 0.1-step grid for a smoother-looking regime-index plot (Figure 3 panel
# C) -- 16 points instead of 6, ~2.7x the compute, but still tractable in
# parallel.
b_grid <- seq(0, 1.5, by = 0.1)

# Simulation horizon.
n_chains <- 200
T_total  <- 1500
late_window <- 1250:1500

# Save every few steps to keep file size reasonable.
save_every <- 5

# ------------------------------------------------------------------------
# 1. One-chain simulator
# ------------------------------------------------------------------------
simulate_history_chain <- function(
    b,
    init = c("low", "high"),
    chain_id = 1,
    T_total = 1500
) {
  init <- match.arg(init)

  if (init == "low") {
    S <- rep(0L, N)
    P <- 0
    m_smooth <- 0
  }
  if (init == "high") {
    S <- rep(1L, N)
    P <- 1
    m_smooth <- 1
  }

  out <- vector("list", length = floor(T_total / save_every) + 1)
  out_i <- 1

  for (t in seq_len(T_total)) {

    # Fast symptom update: one random-order sweep.
    S <- simulate_fast_sweep(S = S, tau = tau, omega = omega, gamma = gamma, P = P)
    M <- sum(S)
    m <- M / N

    # Smooth symptom burden for feedback.
    m_smooth <- (1 - alpha_smooth) * m_smooth + alpha_smooth * m

    # Slow-field update. Signed feedback (2026-08-27, replaces the earlier
    # positive-part clamp) -- see 03_sim_feedback.R for the rationale.
    feedback_term <- b * (m_smooth - m_star) * dt
    P <- P +
      kappa * (P_base - P) * dt +
      sigma_P * sqrt(dt) * rnorm(1) +
      feedback_term

    if (t %% save_every == 0 || t == 1) {
      out[[out_i]] <- tibble(
        b = b,
        init = init,
        chain = chain_id,
        step = t,
        P = P,
        M = M,
        m = m,
        m_smooth = m_smooth
      )
      out_i <- out_i + 1
    }
  }

  bind_rows(out)
}

# ------------------------------------------------------------------------
# 2. Run grid, in parallel
# ------------------------------------------------------------------------
# 6 b-values x 2 inits x 200 chains x 1500 steps = 2400 independent chains,
# substantially more total work than Sim 3's 1000-chain run -- each
# (b, init, chain) combination is an independent, embarrassingly-parallel
# unit of work, so flatten into one task list and fork across cores with
# mclapply() rather than looping serially via pmap_dfr(). Same pattern as
# 03_sim_feedback.R / 04b_sim_network_estimation_replicated.R:
# fork-based mclapply (macOS/Linux only) + "L'Ecuyer-CMRG" + mc.set.seed =
# TRUE so each forked worker gets its own independent, reproducible RNG
# substream instead of silently replaying the parent's stream.
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(20260825)

grid <- expand_grid(
  b = b_grid,
  init = c("low", "high"),
  chain = seq_len(n_chains)
)
n_tasks <- nrow(grid)

message("Running supplementary history-dependence grid...")
message("Total chains: ", n_tasks)

run_task <- function(i) {
  simulate_history_chain(
    b = grid$b[i],
    init = grid$init[i],
    chain_id = grid$chain[i],
    T_total = T_total
  )
}

chain_tbls <- parallel::mclapply(seq_len(n_tasks), run_task,
                                  mc.cores = n_cores, mc.set.seed = TRUE)

# mclapply silently returns try-error objects for tasks that fail in a
# worker rather than stopping the whole run -- check for that before
# trusting the output.
failed <- vapply(chain_tbls, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d parallel tasks failed. First error: %s",
               sum(failed), n_tasks, attr(chain_tbls[[which(failed)[1]]], "condition")$message))
}

sim_hist <- bind_rows(chain_tbls)

dir.create("res/revision_2026/supp_history", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  sim_hist,
  "res/revision_2026/supp_history/history_dependence_raw.rds"
)

# ------------------------------------------------------------------------
# 3. Summaries
# ------------------------------------------------------------------------
late_summary_chain <- sim_hist |>
  filter(step %in% late_window) |>
  group_by(b, init, chain) |>
  summarise(
    late_m = mean(m),
    late_P = mean(P),
    max_abs_P = max(abs(P)),
    .groups = "drop"
  )

late_summary <- late_summary_chain |>
  group_by(b, init) |>
  summarise(
    mean_late_m = mean(late_m),
    se_late_m = sd(late_m) / sqrt(n()),
    mean_late_P = mean(late_P),
    se_late_P = sd(late_P) / sqrt(n()),
    mean_max_abs_P = mean(max_abs_P),
    .groups = "drop"
  )

separation <- late_summary_chain |>
  select(b, init, chain, late_m, late_P, max_abs_P) |>
  pivot_wider(
    names_from = init,
    values_from = c(late_m, late_P, max_abs_P)
  ) |>
  group_by(b) |>
  summarise(
    separation_m = mean(late_m_high - late_m_low, na.rm = TRUE),
    se_separation_m = sd(late_m_high - late_m_low, na.rm = TRUE) / sqrt(n()),
    separation_P = mean(late_P_high - late_P_low, na.rm = TRUE),
    se_separation_P = sd(late_P_high - late_P_low, na.rm = TRUE) / sqrt(n()),
    max_abs_P = max(c(max_abs_P_low, max_abs_P_high), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    regime_flag = case_when(
      max_abs_P > 6 ~ "runaway/saturation",
      separation_m > 0.15 ~ "history-dependent",
      TRUE ~ "convergent"
    )
  )

write_csv(late_summary, "res/revision_2026/supp_history/history_late_summary.csv")
write_csv(separation, "res/revision_2026/supp_history/history_separation_summary.csv")

print(late_summary)
print(separation)

message("Done. Outputs saved in res/revision_2026/supp_history/")
