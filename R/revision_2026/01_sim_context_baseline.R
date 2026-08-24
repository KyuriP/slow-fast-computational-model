# ============================================================
# R/revision_2026/01_sim_context_baseline.R
# ============================================================
# Simulation 1: same tau_i, omega_ij, gamma_i; only P differs.
# No slow diffusion, shocks, feedback, heatmaps, or network estimation --
# those come later (02_sim_stress_recovery.R onward), once this is
# checked clean per Step 6 of the revision plan.
#
# Outputs
# -------
#   res/revision_2026/sim1/sim1_raw.rds
#   res/revision_2026/sim1/sim1_summary.csv
#   figs/revision_2026/fig_sim1_context_baseline.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # defines symptoms, N, tau, gamma, omega

# ------------------------------------------------------------------------
# Design
# ------------------------------------------------------------------------
P_values <- c(low = -0.6, middle = 0, high = 0.6)

# Pilot defaults -- adjust per Step 6 if the burden distribution saturates.
# Per the revision plan: if it saturates, adjust ONLY P values / tau /
# omega strength / gamma strength. Do not change simulate_fast_sweep()
# itself -- that's the one shared model definition every script relies on.
T_burn   <- 200L
T_post   <- 200L
n_chains <- 200L

n_trace_chains <- 10L  # how many chains to keep a full (burn-in + post) trace for,
                        # purely for the convergence check below -- not used in
                        # any summary statistic.

run_condition <- function(P, n_chains, T_burn, T_post) {
  total_sweeps <- T_burn + T_post
  M_post <- matrix(NA_real_, nrow = n_chains, ncol = T_post)
  S_post <- array(NA_real_, dim = c(n_chains, T_post, N))  # per-symptom activation
  M_trace <- matrix(NA_real_, nrow = n_trace_chains, ncol = total_sweeps)  # full trace, burn-in included

  for (c in seq_len(n_chains)) {
    S <- rbinom(N, size = 1, prob = 0.5)  # random init
    for (t in seq_len(total_sweeps)) {
      S <- simulate_fast_sweep(S, tau, omega, gamma, P)
      if (c <= n_trace_chains) M_trace[c, t] <- symptom_burden(S)
      if (t > T_burn) {
        M_post[c, t - T_burn] <- symptom_burden(S)
        S_post[c, t - T_burn, ] <- S
      }
    }
  }
  list(M_post = M_post, S_post = S_post, M_trace = M_trace)
}

# ------------------------------------------------------------------------
# Run all three conditions
# ------------------------------------------------------------------------
set.seed(2026L)
raw_list <- lapply(names(P_values), function(lbl) {
  P <- P_values[[lbl]]
  cat(sprintf("Simulating P = %s (%.1f): %d chains x %d burn-in + %d post-burn-in sweeps...\n",
              lbl, P, n_chains, T_burn, T_post))
  out <- run_condition(P, n_chains, T_burn, T_post)

  burden_tbl <- tibble(condition = lbl, P = P,
                        chain = rep(seq_len(n_chains), T_post),
                        sweep = rep(seq_len(T_post), each = n_chains),
                        M = as.vector(out$M_post))

  # per-symptom activation probability for this condition
  symptom_activation <- tibble(
    condition = lbl, P = P, symptom = symptoms,
    p_active = apply(out$S_post, 3, mean)
  )

  trace_tbl <- tibble(
    condition = lbl, P = P,
    chain = rep(seq_len(n_trace_chains), each = ncol(out$M_trace)),
    sweep = rep(seq_len(ncol(out$M_trace)), times = n_trace_chains),
    M = as.vector(t(out$M_trace))
  )

  list(burden = burden_tbl, symptom_activation = symptom_activation, trace = trace_tbl)
})

burden_all <- bind_rows(lapply(raw_list, `[[`, "burden")) |> mutate(m = M / N)
symptom_activation_all <- bind_rows(lapply(raw_list, `[[`, "symptom_activation"))
trace_all <- bind_rows(lapply(raw_list, `[[`, "trace"))

# ------------------------------------------------------------------------
# Save raw + summary
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim1", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

saveRDS(list(burden = burden_all, symptom_activation = symptom_activation_all,
             params = list(tau = tau, omega = omega, gamma = gamma, symptoms = symptoms),
             design = list(P_values = P_values, T_burn = T_burn, T_post = T_post, n_chains = n_chains)),
        "res/revision_2026/sim1/sim1_raw.rds")

summary_tbl <- burden_all |>
  group_by(condition, P) |>
  summarise(
    mean_M = mean(M), sd_M = sd(M),
    mean_m = mean(m),
    pr_high_burden = mean(M >= N / 2),
    .groups = "drop"
  )
write.csv(summary_tbl, "res/revision_2026/sim1/sim1_summary.csv", row.names = FALSE)

symptom_activation_wide <- symptom_activation_all |>
  select(-P) |>
  pivot_wider(names_from = condition, values_from = p_active) |>
  select(symptom, low, middle, high)
write.csv(symptom_activation_wide, "res/revision_2026/sim1/sim1_symptom_activation.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# Step 6 pilot checks
# ------------------------------------------------------------------------
cat("\n=== MEAN SYMPTOM BURDEN BY P ===\n")
print(summary_tbl |> select(condition, P, mean_M, mean_m, sd_M))

cat("\n=== PROBABILITY OF HIGH BURDEN (M >= N/2) BY P ===\n")
print(summary_tbl |> select(condition, P, pr_high_burden))

cat("\n=== ACTIVATION PROBABILITY BY SYMPTOM AND P ===\n")
print(symptom_activation_all |>
        select(-P) |>
        pivot_wider(names_from = condition, values_from = p_active) |>
        select(symptom, low, middle, high))

cat("\n--- Saturation check ---\n")
cat("Good pattern: low P shows low-but-nonzero activation, high P shows higher\n")
cat("activation but not all nine symptoms pinned near 0 or 1 across the board.\n")
cat("Bad pattern: low P all near 0, high P all near 1 for every symptom -- if\n")
cat("so, adjust P values / tau / omega strength / gamma strength only, per\n")
cat("Step 6 of the revision plan. Do not change simulate_fast_sweep() itself.\n")

# ------------------------------------------------------------------------
# One figure
# ------------------------------------------------------------------------
p1 <- ggplot(burden_all, aes(x = M, fill = condition)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55, colour = "white") +
  scale_fill_viridis_d(name = "Context (P)",
                        breaks = names(P_values),
                        labels = sprintf("%s (P=%.1f)", names(P_values), P_values)) +
  labs(x = "Symptom burden (M = number of active symptoms)", y = "Count",
       title = "Simulation 1: symptom burden distribution under fixed coupling",
       subtitle = sprintf("N=%d symptoms, same tau/omega/gamma across conditions, only P differs", N)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figs/revision_2026/fig_sim1_context_baseline.pdf", p1, width = 8, height = 5.5)

# ------------------------------------------------------------------------
# Burn-in convergence check: does the first half of the post-burn-in
# window (sweeps T_burn+1 .. T_burn+T_post/2) still differ systematically
# from the second half? If mean_M is still trending at T_burn, burn-in
# wasn't long enough and T_burn should be increased.
# ------------------------------------------------------------------------
half <- T_post %/% 2
drift_check <- burden_all |>
  mutate(half = ifelse(sweep <= half, "first_half_post_burnin", "second_half_post_burnin")) |>
  group_by(condition, half) |>
  summarise(mean_M = mean(M), .groups = "drop") |>
  pivot_wider(names_from = half, values_from = mean_M) |>
  mutate(drift = second_half_post_burnin - first_half_post_burnin)

cat("\n=== BURN-IN DRIFT CHECK (want 'drift' close to 0) ===\n")
print(drift_check)
cat("If |drift| is small relative to the between-condition differences in\n")
cat("mean_M above, T_burn = 200 is adequate. If drift is still substantial\n")
cat("(comparable in size to the low-vs-high P gap), increase T_burn and rerun.\n")

p2 <- ggplot(trace_all, aes(x = sweep, y = M, group = chain)) +
  geom_line(alpha = 0.35, linewidth = 0.3) +
  geom_vline(xintercept = T_burn, linetype = "dashed", colour = "red") +
  facet_wrap(~condition, ncol = 1) +
  labs(x = "Sweep", y = "Symptom burden (M)",
       title = sprintf("Simulation 1: burn-in trace (%d example chains per condition)", n_trace_chains),
       subtitle = "Dashed red line = end of burn-in (T_burn). Traces should look flat/stationary by then.") +
  theme_classic(base_size = 12)

ggsave("figs/revision_2026/fig_sim1_burnin_trace.pdf", p2, width = 8, height = 8)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim1/sim1_raw.rds\n")
cat("  res/revision_2026/sim1/sim1_summary.csv\n")
cat("  res/revision_2026/sim1/sim1_symptom_activation.csv\n")
cat("  figs/revision_2026/fig_sim1_context_baseline.pdf\n")
cat("  figs/revision_2026/fig_sim1_burnin_trace.pdf\n")
