# ============================================================
# R/revision_2026/02_sim_stress_recovery.R
# ============================================================
# Simulation 2: introduce P_t dynamics (mean reversion + diffusion +
# a single acute shock). Feedback is OFF (b = 0) -- this shows stress and
# recovery cleanly before Simulation 3 adds symptom-to-context feedback.
#
# Model (per Step 7):
#   P_{t+dt} = P_t + kappa*(P_base - P_t)*dt + sigma_P*sqrt(dt)*eps_t
#              + shock_t + b*(m_smooth - m_star)*dt      [b = 0 here]
#
# Fast layer unchanged from Simulation 1:
#   logit Pr(S_i=1 | S_-i, P_t) = tau_i + sum_j!=i omega_ij S_j + gamma_i P_t
#
# At the shock step, the perturbation is applied first. Symptoms are then
# updated conditional on the perturbed P_t. After recording, P_t recovers
# toward baseline for the next step. Burn-in (200 steps, same as the
# validated Sim 1 burn-in) lets the fast layer mix before the shock is
# applied; P_base = 0 so there's no separate slow-state burn-in needed.
#
# Outputs
# -------
#   res/revision_2026/sim2/sim2_raw.rds
#   res/revision_2026/sim2/sim2_summary.csv
#   figs/revision_2026/fig_sim2_stress_recovery.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (with +1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-state parameters
# ------------------------------------------------------------------------
kappa   <- 0.20   # mean-reversion rate
sigma_P <- 0.04   # diffusion scale
dt      <- 0.02   # step size
P_base  <- 0      # "middle" baseline, matching Sim 1's middle condition
b       <- 0      # feedback OFF for this simulation -- Sim 3 turns this on

burn_in_steps    <- 200L   # fast-layer mixing only (validated in Sim 1); P starts at P_base
post_shock_steps <- 750L   # long enough to see substantial (~95%) recovery: exp(-kappa*post_shock_steps*dt) ~ 0.05
shock_time    <- burn_in_steps + 1L
shock_magnitude <- 1.0     # single deterministic jump added to P at shock_time

total_steps <- burn_in_steps + post_shock_steps
n_chains <- 200L

# ------------------------------------------------------------------------
# One chain: P_t trajectory + fast-layer response, before/after one shock
# ------------------------------------------------------------------------
run_chain <- function() {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  P_trace <- numeric(total_steps)
  M_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {

    # Apply the acute perturbation before symptoms respond at this step.
    if (t == shock_time) {
      P <- P + shock_magnitude
    }

    # Symptoms respond to the current slow field.
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    # Record the current state.
    P_trace[t] <- P
    M_trace[t] <- symptom_burden(S)

    # Slow recovery/diffusion for the next step.
    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1)
    # feedback term b*(m_smooth - m_star)*dt omitted entirely since b = 0
  }
  list(P = P_trace, M = M_trace)
}

set.seed(2026L)
cat(sprintf("Simulating %d chains x %d steps (burn-in=%d, shock at step %d, magnitude=%.1f)...\n",
            n_chains, total_steps, burn_in_steps, shock_time, shock_magnitude))

chains <- lapply(seq_len(n_chains), function(c) run_chain())

P_mat <- do.call(rbind, lapply(chains, `[[`, "P"))
M_mat <- do.call(rbind, lapply(chains, `[[`, "M"))

# ------------------------------------------------------------------------
# Save + summarize (time relative to shock, so shock = time 0)
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim2", recursive = TRUE, showWarnings = FALSE)

traj <- tibble(
  step = rep(seq_len(total_steps), each = n_chains),
  chain = rep(seq_len(n_chains), times = total_steps),
  P = as.vector(P_mat), M = as.vector(M_mat)
) |>
  mutate(time_since_shock = step - shock_time, m = M / N)

saveRDS(list(traj = traj, params = list(kappa = kappa, sigma_P = sigma_P, dt = dt,
                                          P_base = P_base, b = b, shock_time = shock_time,
                                          shock_magnitude = shock_magnitude,
                                          burn_in_steps = burn_in_steps,
                                          post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim2/sim2_raw.rds")

summary_tbl <- traj |>
  group_by(step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(m), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim2/sim2_summary.csv", row.names = FALSE)

cat("\n=== PRE-SHOCK BASELINE (last 50 burn-in steps) ===\n")
print(summary_tbl |> filter(time_since_shock < 0, time_since_shock >= -50) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\n=== PEAK RESPONSE (first 20 steps after shock) ===\n")
print(summary_tbl |> filter(time_since_shock >= 0, time_since_shock < 20) |>
        summarise(peak_P = max(mean_P), peak_M = max(mean_M), peak_m = max(mean_m)))

cat("\n=== END-OF-WINDOW RECOVERY (last 50 steps) ===\n")
print(summary_tbl |> filter(time_since_shock >= post_shock_steps - 50) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\nCompare pre-shock baseline vs. end-of-window: if end-of-window mean_P and\n")
cat("mean_M are close to pre-shock baseline, recovery is essentially complete by\n")
cat("the end of the simulated window. If not, extend post_shock_steps.\n")

# ------------------------------------------------------------------------
# Figure: stacked P_t / M_t trajectory around the shock (matches the
# upper-P/lower-symptom-trace convention used in the original manuscript's
# Figure 3)
# ------------------------------------------------------------------------
plot_window <- summary_tbl |> filter(time_since_shock >= -100, time_since_shock <= post_shock_steps)

pP <- ggplot(plot_window, aes(time_since_shock, mean_P)) +
  geom_ribbon(aes(ymin = mean_P - se_P, ymax = mean_P + se_P), alpha = 0.25) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  labs(x = NULL, y = expression(P[t]), title = "Simulation 2: stress and recovery (feedback off, b=0)") +
  theme_classic(base_size = 12)

pM <- ggplot(plot_window, aes(time_since_shock, mean_m)) +
  geom_ribbon(aes(ymin = mean_m - se_M / N, ymax = mean_m + se_M / N), alpha = 0.25) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  labs(x = "Steps since shock (dashed line = shock onset)", y = expression(m[t])) +
  theme_classic(base_size = 12)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
p_combined <- patchwork::wrap_plots(pP, pM, ncol = 1)

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/fig_sim2_stress_recovery.pdf", p_combined, width = 8, height = 6.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim2/sim2_raw.rds\n")
cat("  res/revision_2026/sim2/sim2_summary.csv\n")
cat("  figs/revision_2026/fig_sim2_stress_recovery.pdf\n")
