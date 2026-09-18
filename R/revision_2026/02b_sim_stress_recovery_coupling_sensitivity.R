# ============================================================
# R/revision_2026/02b_sim_stress_recovery_coupling_sensitivity.R
# ============================================================
# Supplementary sensitivity check for Simulation 2 (Appendix): does
# recovery from the same slow-context perturbation depend on the
# strength of symptom-symptom coupling? Same design as
# 02_sim_stress_recovery.R (single shock, feedback off, b=0), but the
# coupling matrix is scaled by c in c_grid while its topology (which
# pairs are connected) is held fixed. c=1.0 reproduces 02's own omega
# exactly -- that's the main-text Simulation 2 condition, kept here as
# the reference level rather than replaced.
#
# Recovery is reported per condition as normalized excess activation,
#   R(t) = (M(t) - M_pre) / (M_peak - M_pre),
# using each condition's OWN pre-shock baseline and post-shock peak --
# not raw M(t). Scaling omega up also raises baseline activation, so
# comparing raw M(t) across conditions would conflate "recovers slower"
# with "started from a higher baseline".
#
# c in {0.5, 1.0, 1.5} alone (n_chains=1000) showed no coupling
# dependence in recovery timing -- recovery is governed by the
# slow-process relaxation rate (kappa), not by omega, which is what the
# slow-fast separation this model is built on would predict. c_grid now
# extends further (up to 3x the reference omega) to check a different,
# real hypothesis: whether pushing coupling strong enough drives the
# network toward its own critical/bistable regime, where the FAST
# layer's own relaxation could slow down independent of P. That is
# checked directly below (CRITICALITY / BIMODALITY CHECK) rather than
# assumed from the mean trajectory alone -- rising per-chain variance or
# bimodality in pre-shock activation is what that regime would look
# like; a merely higher mean is not evidence of it.
#
# Outputs
# -------
#   res/revision_2026/sim2/sim2_coupling_sensitivity_raw.rds
#   res/revision_2026/sim2/sim2_coupling_sensitivity_summary.csv
#   res/revision_2026/sim2/sim2_coupling_sensitivity_recovery.csv
#   res/revision_2026/sim2/sim2_coupling_sensitivity_criticality_diagnostics.csv
#   figs/revision_2026/fig_sim2_coupling_sensitivity.pdf
#   figs/revision_2026/fig_sim2_coupling_criticality_check.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

source("R/revision_2026/utils_uncentered01_model.R")
source("R/revision_2026/00_parameters_uncentered01.R")  # tau, omega, gamma, symptoms, N -- omega here is c=1.0
source("R/revision_2026/figures/theme_publication.R")   # theme_pub(), main_line_width

kappa   <- 0.20
sigma_P <- 0.04
dt      <- 0.02
P_base  <- 0
b       <- 0

burn_in_steps    <- 200L
post_shock_steps <- 750L
shock_time       <- burn_in_steps + 1L
shock_magnitude  <- 1.0

total_steps <- burn_in_steps + post_shock_steps
n_chains <- 1000L

c_grid <- c(0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)  # 0.3 added for a wider low-end contrast point
stopifnot(length(c_grid) >= 3, 1 %in% c_grid)

run_chain_c <- function(omega_c) {
  P <- P_base
  S <- rbinom(N, size = 1, prob = 0.5)
  P_trace <- numeric(total_steps)
  M_trace <- numeric(total_steps)
  for (t in seq_len(total_steps)) {
    if (t == shock_time) P <- P + shock_magnitude
    S <- simulate_fast_sweep(S, tau, omega_c, gamma, P)
    P_trace[t] <- P
    M_trace[t] <- symptom_burden(S)
    P <- P + kappa * (P_base - P) * dt + sigma_P * sqrt(dt) * rnorm(1)
  }
  list(P = P_trace, M = M_trace)
}

run_condition_c <- function(c_val) {
  omega_c <- c_val * omega
  cat(sprintf("Simulating c = %.2f (%d chains x %d steps)...\n", c_val, n_chains, total_steps))
  chains <- lapply(seq_len(n_chains), function(i) run_chain_c(omega_c))
  P_mat <- do.call(rbind, lapply(chains, `[[`, "P"))
  M_mat <- do.call(rbind, lapply(chains, `[[`, "M"))
  tibble(
    c_coupling = c_val,
    step = rep(seq_len(total_steps), each = n_chains),
    chain = rep(seq_len(n_chains), times = total_steps),
    P = as.vector(P_mat), M = as.vector(M_mat)
  ) |>
    mutate(time_since_shock = step - shock_time, m = M / N)
}

set.seed(2026L)
traj_all <- bind_rows(lapply(c_grid, run_condition_c))

dir.create("res/revision_2026/sim2", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(traj = traj_all,
             params = list(kappa = kappa, sigma_P = sigma_P, dt = dt, P_base = P_base, b = b,
                           shock_time = shock_time, shock_magnitude = shock_magnitude,
                           burn_in_steps = burn_in_steps, post_shock_steps = post_shock_steps,
                           c_grid = c_grid)),
        "res/revision_2026/sim2/sim2_coupling_sensitivity_raw.rds")

summary_tbl <- traj_all |>
  group_by(c_coupling, step, time_since_shock) |>
  summarise(mean_P = mean(P), mean_M = mean(M), se_M = sd(M) / sqrt(n()), mean_m = mean(m), .groups = "drop")
write.csv(summary_tbl, "res/revision_2026/sim2/sim2_coupling_sensitivity_summary.csv", row.names = FALSE)

baseline_tbl <- summary_tbl |>
  filter(time_since_shock < 0, time_since_shock >= -50) |>
  group_by(c_coupling) |>
  summarise(mean_M_pre = mean(mean_M), mean_m_pre = mean(mean_m), .groups = "drop")

cat("\n=== SATURATION CHECK (pre-shock baseline by coupling condition) ===\n")
cat("want mean_m_pre away from 0/1 for every c -- a mean pinned near 0 or 1\n")
cat("means that condition is uninformative on its own, but check the\n")
cat("criticality diagnostic below too: a network can be near a tipping\n")
cat("point without its MEAN baseline looking extreme.\n\n")
print(baseline_tbl)

peak_tbl <- summary_tbl |>
  filter(time_since_shock >= 0, time_since_shock < 20) |>
  group_by(c_coupling) |>
  summarise(mean_M_peak = max(mean_M), .groups = "drop")

cat("\n=== PEAK RESPONSE BY COUPLING CONDITION ===\n")
print(peak_tbl)

recovery_tbl <- summary_tbl |>
  left_join(baseline_tbl, by = "c_coupling") |>
  left_join(peak_tbl, by = "c_coupling") |>
  mutate(R = (mean_M - mean_M_pre) / (mean_M_peak - mean_M_pre))

write.csv(recovery_tbl |> select(c_coupling, step, time_since_shock, mean_M, R),
          "res/revision_2026/sim2/sim2_coupling_sensitivity_recovery.csv", row.names = FALSE)

recovery_time <- function(df, threshold) {
  hit <- df |> filter(time_since_shock >= 0, R <= threshold)
  if (nrow(hit) == 0) return(NA_real_)
  min(hit$time_since_shock)
}
recovery_time_summary <- recovery_tbl |>
  group_by(c_coupling) |>
  group_modify(~ tibble(t_to_R50 = recovery_time(.x, 0.5), t_to_R10 = recovery_time(.x, 0.1))) |>
  ungroup()

cat("\n=== TIME TO 50% / 90% RECOVERY (steps since shock; NA = not reached in window) ===\n")
print(recovery_time_summary)

# --- criticality / bimodality check -------------------------------------
# A network approaching its own tipping point shows rising variance and
# eventually bimodality in its resting-state activation across
# independent chains, even while the MEAN stays well inside (0, N).
# That is a different signature than the saturation check above and has
# to be checked on the per-chain distribution, not the condition mean.
bimodality_coefficient <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  if (s == 0 || n < 4) return(NA_real_)
  skew <- mean((x - m)^3) / s^3
  kurt <- mean((x - m)^4) / s^4
  (skew^2 + 1) / (kurt + 3 * (n - 1)^2 / ((n - 2) * (n - 3)))
}

chain_baseline <- traj_all |>
  filter(time_since_shock >= -50, time_since_shock < 0) |>
  group_by(c_coupling, chain) |>
  summarise(chain_mean_M = mean(M), .groups = "drop")

criticality_diag <- chain_baseline |>
  group_by(c_coupling) |>
  summarise(
    n_chains_used = n(),
    mean_M_pre_chain = mean(chain_mean_M),
    sd_M_pre_chain = sd(chain_mean_M),
    min_M_pre_chain = min(chain_mean_M),
    max_M_pre_chain = max(chain_mean_M),
    bimodality_coef = bimodality_coefficient(chain_mean_M),
    .groups = "drop"
  )

cat("\n=== CRITICALITY / BIMODALITY CHECK (pre-shock M, averaged per chain) ===\n")
cat("Sarle's bimodality coefficient > ~0.555 suggests bimodality (rule of\n")
cat("thumb -- a uniform distribution gives exactly 5/9). A rising\n")
cat("sd_M_pre_chain with c, or bimodality_coef crossing ~0.555, means that\n")
cat("condition should NOT be trusted as a clean 'stronger coupling, same\n")
cat("kind of dynamics' comparison -- the network is behaving qualitatively\n")
cat("differently there, not just more strongly.\n\n")
print(criticality_diag)

write.csv(criticality_diag, "res/revision_2026/sim2/sim2_coupling_sensitivity_criticality_diagnostics.csv", row.names = FALSE)

# --- palette + labels (diverging, anchored at c = 1) ---------------------
c_lab_str <- sprintf("c = %.2g%s", c_grid, ifelse(c_grid == 1, " (reference)", ""))
c_labels <- setNames(c_lab_str, as.character(c_grid))

diverging_ramp <- grDevices::colorRampPalette(c("#3B76AF", "grey40", "#C85C5C"))(201)
rng <- max(abs(c_grid - 1))
get_col <- function(c_val) {
  pos <- if (rng == 0) 0 else (c_val - 1) / rng
  idx <- round((pos + 1) / 2 * 200) + 1
  diverging_ramp[idx]
}
pal_coupling <- setNames(sapply(c_grid, get_col), c_lab_str)

recovery_tbl <- recovery_tbl |>
  mutate(c_label = factor(c_labels[as.character(c_coupling)], levels = c_lab_str))

p_recovery <- ggplot(
  recovery_tbl |> filter(time_since_shock >= -20, time_since_shock <= post_shock_steps),
  aes(time_since_shock, R, colour = c_label)
) +
  geom_hline(yintercept = c(0, 1), linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_line(linewidth = main_line_width) +
  scale_colour_manual(values = pal_coupling, name = NULL) +
  labs(x = "Steps since shock", y = "R(t)  (normalized excess activation)",
       title = "Sensitivity of recovery to symptom-coupling strength",
       subtitle = expression(omega^{(c)} == c %.% omega * ",  topology fixed -- recovery relative to each condition's own baseline")) +
  theme_pub()

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)
ggsave("figs/revision_2026/fig_sim2_coupling_sensitivity.pdf", p_recovery, width = 8, height = 5.5)

p_hist <- ggplot(
  chain_baseline |> mutate(c_label = factor(c_labels[as.character(c_coupling)], levels = c_lab_str)),
  aes(x = chain_mean_M, fill = c_label)
) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.1) +
  facet_wrap(~c_label, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = pal_coupling, guide = "none") +
  labs(x = "Per-chain mean pre-shock symptom activation (M)", y = "Number of chains",
       title = "Pre-shock activation distribution by coupling strength",
       subtitle = "A shift from unimodal toward bimodal signals approach to a critical/bistable regime") +
  theme_pub()

ggsave("figs/revision_2026/fig_sim2_coupling_criticality_check.pdf", p_hist, width = 8, height = 6)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim2/sim2_coupling_sensitivity_raw.rds\n")
cat("  res/revision_2026/sim2/sim2_coupling_sensitivity_summary.csv\n")
cat("  res/revision_2026/sim2/sim2_coupling_sensitivity_recovery.csv\n")
cat("  res/revision_2026/sim2/sim2_coupling_sensitivity_criticality_diagnostics.csv\n")
cat("  figs/revision_2026/fig_sim2_coupling_sensitivity.pdf\n")
cat("  figs/revision_2026/fig_sim2_coupling_criticality_check.pdf\n")
cat("\nCheck the CRITICALITY / BIMODALITY CHECK output above before trusting\n")
cat("any high-c condition as a like-for-like comparison to the reference.\n")
