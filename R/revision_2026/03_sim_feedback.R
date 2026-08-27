# ============================================================
# R/revision_2026/03_sim_feedback.R
# ============================================================
# Simulation 3 (MAIN TEXT): same shock/recovery design as Simulation 2, but
# with symptom-to-context feedback turned on. Compares b = 0 (Sim 2's
# validated baseline) against b = 0.5, isolating the effect of feedback from
# everything else, since fast layer, shock timing/magnitude, kappa, sigma_P,
# and dt are all identical to Sim 2.
#
# b = 0.5 was selected after a pilot grid (b = 0.3/0.5/0.7; see
# 03b_sim_feedback_grid_supplement.R) because it gives a visible feedback-
# related delay in both P and symptom burden, while all trajectories still
# decay toward baseline. The goal is persistence/slower recovery, not
# bistability or runaway dynamics: b=0.3 alone was too subtle in symptom
# burden (M: 2.58 off vs 2.71 on); b=0.7 gave the largest separation but
# risks reading as tipping-like even though the trajectory is still decaying.
#
# Model:
#   P_{t+dt} = P_t + kappa*(P_base - P_t)*dt + sigma_P*sqrt(dt)*eps_t
#              + b*(m_smooth_t - m_star)*dt
#   (shock applied as a one-off jump to P at shock_time, same as Sim 2)
#
# 2026-08-27: switched from the positive-part clamp [m_smooth_t - m_star]_+
# to the signed difference (m_smooth_t - m_star). Feedback is now symmetric:
# elevated burden (m_smooth > m_star) pushes the slow field up, and burden
# below m_star actively pulls P down (in addition to the kappa mean-
# reversion term already doing so). All main-text equations/captions/tables
# should describe this as the signed feedback term b*(m_bar_t - m_star), not
# the clamped [.]_+ form used in the earlier design.
#
# Fast layer unchanged:
#   logit Pr(S_i=1 | S_-i, P_t) = tau_i + sum_j!=i omega_ij S_j + gamma_i P_t
#
# m_smooth_t is an exponential moving average of the single-sweep symptom
# fraction m_t = M_t/N. Feedback needs a smoothed signal because m_t from
# a single fast sweep is noisy (this is the same reason Sim 1 averaged
# over 200 post-burn-in sweeps rather than reading off one sweep); an EMA
# is the cheapest way to get a low-noise, causally-laggable signal without
# storing a full window.
#
# m_star is the reference ("no context effect") burden level: the P=0
# equilibrium mean_m = 0.2766 from Sim 1's middle condition
# (res/revision_2026/sim1/sim1_summary.csv). Choosing m_star this way
# means the system starts already at its own b=0 fixed point (P=P_base=0,
# m_smooth=m_star), so burn-in does not need to be extended for the
# feedback case -- there's nothing to equilibrate away.
#
# b and the EMA smoothing constant (alpha_smooth) are pilot values, flagged
# below, same calibration status as Sim 1's tau shift. b is chosen to be
# clearly visible against the b=0 recovery curve (i.e., end-of-window
# mean_P / mean_M should sit further from the pre-shock baseline than the
# b=0 case) without letting the loop run away or create a second stable
# state -- runaway/bistability/hysteresis is out of scope for this
# simulation and is being deliberately avoided, not tested for.
#
# Outputs
# -------
#   res/revision_2026/sim3/sim3_raw.rds
#   res/revision_2026/sim3/sim3_summary.csv
#   figs/revision_2026/fig_sim3_feedback_comparison.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Slow-state parameters (kappa, sigma_P, dt, shock: identical to Sim 2)
# ------------------------------------------------------------------------
kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02
P_base  <- 0

burn_in_steps    <- 200L
post_shock_steps <- 750L
shock_time       <- burn_in_steps + 1L
shock_magnitude  <- 1.0

total_steps <- burn_in_steps + post_shock_steps
n_chains    <- 1000L   # raised from 200 during calibration -- at 200 chains
                        # the b=0.3 peak (first 20 post-shock steps) came out
                        # lower than the b=0 peak (4.41 vs 4.73); at 1000
                        # chains that gap disappeared (confirmed Monte Carlo
                        # noise, not a real effect), so 1000 is kept as the
                        # locked chain count for the main-text comparison.

# ------------------------------------------------------------------------
# Feedback-specific parameters -- PILOT VALUES, calibrate against output
# ------------------------------------------------------------------------
m_star       <- 0.2766389   # Sim 1 middle-condition (P=0) equilibrium mean_m
alpha_smooth <- 0.05        # EMA rate; half-life = ln(2)/alpha_smooth steps
                             # ~= 13.9 steps ~= 0.28 time units at dt=0.02 --
                             # fast enough to track the post-shock rise/decay,
                             # slow enough to average out single-sweep noise.

b_values <- c(off = 0, on = 0.5)   # LOCKED main-text comparison

# ------------------------------------------------------------------------
# One chain: P_t / m_t trajectory with feedback strength b
# ------------------------------------------------------------------------
run_chain <- function(b) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  m_smooth <- m_star   # start at the reference level (see header note)

  P_trace  <- numeric(total_steps)
  M_trace  <- numeric(total_steps)
  ms_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {

    # Apply the acute perturbation before symptoms respond at this step.
    if (t == shock_time) {
      P <- P + shock_magnitude
    }

    # Symptoms respond to the current slow field.
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    # Record the current state.
    P_trace[t]  <- P
    M_trace[t]  <- symptom_burden(S)
    ms_trace[t] <- m_smooth

    # Update the smoothed burden signal used for feedback.
    m_t <- M_trace[t] / N
    m_smooth <- m_smooth + alpha_smooth * (m_t - m_smooth)

    # Slow recovery/diffusion + feedback for the next step.
    # Signed feedback (2026-08-27, replaces the earlier positive-part
    # clamp): elevated burden (m_smooth > m_star) pushes P up, and burden
    # below m_star pulls P down.
    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1) +
      b * (m_smooth - m_star) * dt
  }
  list(P = P_trace, M = M_trace, m_smooth = ms_trace)
}

# ------------------------------------------------------------------------
# Run all (feedback condition x chain) combinations in parallel
# ------------------------------------------------------------------------
# Each chain is an independent, embarrassingly-parallel unit of work, so we
# flatten (condition, chain) into one task list and fork across cores with
# mclapply() rather than looping serially per condition. mclapply is
# fork-based (macOS/Linux only, which is what this repo is developed on);
# it wouldn't work as-is on Windows.
#
# RNG note: forked workers inherit the parent's RNG state, so without care
# every worker would draw the *same* random stream. Using the "L'Ecuyer-CMRG"
# generator + mc.set.seed=TRUE gives each forked worker its own independent,
# reproducible substream (the standard approach recommended in ?mclapply).

n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

tasks <- expand.grid(cond = names(b_values), chain = seq_len(n_chains),
                      stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("Simulating %d chains x %d steps x %d feedback conditions (b = %s), %d tasks total...\n",
            n_chains, total_steps, length(b_values), paste(b_values, collapse = ", "), n_tasks))

run_task <- function(i) {
  cond <- tasks$cond[i]
  b <- b_values[[cond]]
  out <- run_chain(b)
  tibble(feedback = cond, b = b, chain = tasks$chain[i],
         step = seq_len(total_steps), P = out$P, M = out$M, m_smooth = out$m_smooth)
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

traj <- bind_rows(chain_tbls) |>
  mutate(feedback = factor(feedback, levels = names(b_values)),
         time_since_shock = step - shock_time, m = M / N)

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim3", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(traj = traj,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base,
                            b_values = b_values, m_star = m_star, alpha_smooth = alpha_smooth,
                            shock_time = shock_time, shock_magnitude = shock_magnitude,
                            burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim3/sim3_raw.rds")

summary_tbl <- traj |>
  group_by(feedback, b, step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(m), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim3/sim3_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Pilot checks: does feedback visibly slow recovery relative to b=0?
# ------------------------------------------------------------------------
cat("\n=== PRE-SHOCK BASELINE (last 50 burn-in steps), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock < 0, time_since_shock >= -50) |>
        group_by(feedback) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\n=== PEAK RESPONSE (first 20 steps after shock), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock >= 0, time_since_shock < 20) |>
        group_by(feedback) |>
        summarise(peak_P = max(mean_P), peak_M = max(mean_M), peak_m = max(mean_m)))

cat("\n=== END-OF-WINDOW RECOVERY (last 50 steps), by feedback condition ===\n")
print(summary_tbl |> filter(time_since_shock >= post_shock_steps - 50) |>
        group_by(feedback) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), mean_m = mean(mean_m)))

cat("\nCompare 'on' vs 'off': 'on' should sit further from the pre-shock\n")
cat("baseline than 'off' at end-of-window (feedback slows/incompletes\n")
cat("recovery), while still clearly decaying rather than plateaued or\n")
cat("still rising -- we are not trying to induce bistability/hysteresis\n")
cat("here, just a visibly slower/less complete recovery. This comparison\n")
cat("is LOCKED (b=0.5, chosen via 03b_sim_feedback_grid_supplement.R); this\n")
cat("script is for reproducing/re-checking the final main-text numbers, not\n")
cat("for re-calibrating b.\n")

# ------------------------------------------------------------------------
# Figure: P_t / m_t trajectories, b=off vs b=on (0.5) overlaid -- MAIN FIGURE
# ------------------------------------------------------------------------
plot_window <- summary_tbl |> filter(time_since_shock >= -100, time_since_shock <= post_shock_steps)

pP <- ggplot(plot_window, aes(time_since_shock, mean_P, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_P - se_P, ymax = mean_P + se_P), alpha = 0.2, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(off = "grey50", on = "firebrick"),
                       labels = c(off = "b = 0", on = sprintf("b = %.2f", b_values[["on"]]))) +
  scale_fill_manual(values = c(off = "grey50", on = "firebrick"), guide = "none") +
  labs(x = NULL, y = expression(P[t]), colour = "Feedback",
       title = "Simulation 3: symptom-to-context feedback slows recovery") +
  theme_classic(base_size = 12)

pM <- ggplot(plot_window, aes(time_since_shock, mean_m, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_m - se_M / N, ymax = mean_m + se_M / N), alpha = 0.2, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = c(off = "grey50", on = "firebrick"),
                       labels = c(off = "b = 0", on = sprintf("b = %.2f", b_values[["on"]]))) +
  scale_fill_manual(values = c(off = "grey50", on = "firebrick"), guide = "none") +
  labs(x = "Steps since shock (dashed line = shock onset)", y = expression(m[t]), colour = "Feedback") +
  theme_classic(base_size = 12)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
p_combined <- patchwork::wrap_plots(pP, pM, ncol = 1) + patchwork::plot_layout(guides = "collect")

ggsave("figs/revision_2026/fig_sim3_feedback_comparison.pdf", p_combined, width = 8, height = 6.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim3/sim3_raw.rds\n")
cat("  res/revision_2026/sim3/sim3_summary.csv\n")
cat("  figs/revision_2026/fig_sim3_feedback_comparison.pdf\n")
