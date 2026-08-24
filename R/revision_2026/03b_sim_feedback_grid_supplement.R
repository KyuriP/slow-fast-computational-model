# ============================================================
# R/revision_2026/03b_sim_feedback_grid_supplement.R
# ============================================================
# Supplement to Simulation 3: the b calibration grid used to pick the main-
# text feedback value. Same design as 03_sim_feedback.R (same shock/recovery
# setup as Simulation 2, feedback term b*max(m_smooth - m_star, 0)*dt), but
# sweeps b across off/mild/medium/strong instead of just off vs on.
#
# This is NOT the main-text simulation. b=0.5 ("medium") was selected from
# this grid as the main-text value because it produced a visible feedback
# effect in both P and symptom burden while all four conditions continued to
# decay toward baseline (checked explicitly below and in the original
# calibration run: comparing steps 700-724 vs 725-749 within the post-shock
# window showed every condition, including "strong", still declining, not
# plateaued). b=0.3 alone was too subtle in symptom burden; b=0.7 gave the
# largest separation but risks reading as tipping-like even though it isn't.
# See 03_sim_feedback.R for the locked main-text off-vs-on(b=0.5) comparison.
#
# Outputs
# -------
#   res/revision_2026/sim3/sim3_feedback_grid_raw.rds
#   res/revision_2026/sim3/sim3_feedback_grid_summary.csv
#   figs/revision_2026/fig_sim3_feedback_grid_supplement.pdf
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
# Slow-state parameters (identical to Sim 2 / Sim 3 main)
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
n_chains    <- 1000L

# ------------------------------------------------------------------------
# Feedback-specific parameters
# ------------------------------------------------------------------------
m_star       <- 0.2766389   # Sim 1 middle-condition (P=0) equilibrium mean_m
alpha_smooth <- 0.05

b_grid <- c(off = 0, mild = 0.3, medium = 0.5, strong = 0.7)

# ------------------------------------------------------------------------
# One chain: P_t / m_t trajectory with feedback strength b (positive-part)
# ------------------------------------------------------------------------
run_chain <- function(b) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  m_smooth <- m_star

  P_trace  <- numeric(total_steps)
  M_trace  <- numeric(total_steps)
  ms_trace <- numeric(total_steps)

  for (t in seq_len(total_steps)) {
    if (t == shock_time) {
      P <- P + shock_magnitude
    }
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)

    P_trace[t]  <- P
    M_trace[t]  <- symptom_burden(S)
    ms_trace[t] <- m_smooth

    m_t <- M_trace[t] / N
    m_smooth <- m_smooth + alpha_smooth * (m_t - m_smooth)

    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1) +
      b * max(m_smooth - m_star, 0) * dt
  }
  list(P = P_trace, M = M_trace, m_smooth = ms_trace)
}

# ------------------------------------------------------------------------
# Run all (feedback condition x chain) combinations in parallel
# ------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

tasks <- expand.grid(cond = names(b_grid), chain = seq_len(n_chains),
                      stringsAsFactors = FALSE)
n_tasks <- nrow(tasks)

cat(sprintf("Simulating %d chains x %d steps x %d feedback conditions (b = %s), %d tasks total...\n",
            n_chains, total_steps, length(b_grid), paste(b_grid, collapse = ", "), n_tasks))

run_task <- function(i) {
  cond <- tasks$cond[i]
  b <- b_grid[[cond]]
  out <- run_chain(b)
  tibble(feedback = cond, b = b, chain = tasks$chain[i],
         step = seq_len(total_steps), P = out$P, M = out$M, m_smooth = out$m_smooth)
}

chain_tbls <- parallel::mclapply(seq_len(n_tasks), run_task,
                                  mc.cores = n_cores, mc.set.seed = TRUE)

failed <- vapply(chain_tbls, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d parallel tasks failed. First error: %s",
               sum(failed), n_tasks, attr(chain_tbls[[which(failed)[1]]], "condition")$message))
}

traj <- bind_rows(chain_tbls) |>
  mutate(feedback = factor(feedback, levels = names(b_grid)),
         time_since_shock = step - shock_time, m = M / N)

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim3", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(traj = traj,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base,
                            b_grid = b_grid, m_star = m_star, alpha_smooth = alpha_smooth,
                            shock_time = shock_time, shock_magnitude = shock_magnitude,
                            burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps)),
        "res/revision_2026/sim3/sim3_feedback_grid_raw.rds")

summary_tbl <- traj |>
  group_by(feedback, b, step, time_since_shock) |>
  summarise(mean_P = mean(P), se_P = sd(P) / sqrt(n()),
            mean_M = mean(M), se_M = sd(M) / sqrt(n()),
            mean_m = mean(m), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim3/sim3_feedback_grid_summary.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Pilot checks
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

cat("\n=== STILL DECAYING CHECK: end-of-window first half vs second half ===\n")
print(summary_tbl |>
        filter(time_since_shock >= post_shock_steps - 50) |>
        mutate(half = ifelse(time_since_shock < post_shock_steps - 25, "first_half", "second_half")) |>
        group_by(feedback, half) |>
        summarise(mean_P = mean(mean_P), mean_M = mean(mean_M), .groups = "drop") |>
        pivot_wider(names_from = half, values_from = c(mean_P, mean_M)))
cat("Each condition's second_half should be <= first_half (still decaying, not\n")
cat("plateaued/runaway) -- this is what rules out bistability/hysteresis here.\n")

# ------------------------------------------------------------------------
# Figure: P_t / m_t trajectories across the full b grid
# ------------------------------------------------------------------------
plot_window <- summary_tbl |> filter(time_since_shock >= -100, time_since_shock <= post_shock_steps)

feedback_palette <- setNames(c("grey50", "#fcae91", "#fb6a4a", "#a50f15"), names(b_grid))
feedback_labels  <- setNames(sprintf("%s (b=%.2f)", names(b_grid), b_grid), names(b_grid))

pP <- ggplot(plot_window, aes(time_since_shock, mean_P, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_P - se_P, ymax = mean_P + se_P), alpha = 0.15, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = feedback_palette, labels = feedback_labels) +
  scale_fill_manual(values = feedback_palette, guide = "none") +
  labs(x = NULL, y = expression(P[t]), colour = "Feedback",
       title = "Simulation 3 supplement: feedback calibration grid") +
  theme_classic(base_size = 12)

pM <- ggplot(plot_window, aes(time_since_shock, mean_m, colour = feedback, fill = feedback)) +
  geom_ribbon(aes(ymin = mean_m - se_M / N, ymax = mean_m + se_M / N), alpha = 0.15, colour = NA) +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  scale_colour_manual(values = feedback_palette, labels = feedback_labels) +
  scale_fill_manual(values = feedback_palette, guide = "none") +
  labs(x = "Steps since shock (dashed line = shock onset)", y = expression(m[t]), colour = "Feedback") +
  theme_classic(base_size = 12)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
p_combined <- patchwork::wrap_plots(pP, pM, ncol = 1) + patchwork::plot_layout(guides = "collect")

ggsave("figs/revision_2026/fig_sim3_feedback_grid_supplement.pdf", p_combined, width = 8, height = 6.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim3/sim3_feedback_grid_raw.rds\n")
cat("  res/revision_2026/sim3/sim3_feedback_grid_summary.csv\n")
cat("  figs/revision_2026/fig_sim3_feedback_grid_supplement.pdf\n")
