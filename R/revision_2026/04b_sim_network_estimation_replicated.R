# ============================================================
# R/revision_2026/04b_sim_network_estimation_replicated.R
# ============================================================
# Replicated version of Simulation 4 (see 04_sim_network_estimation.R for
# the single-run pilot that established n_person=10000 was needed to get a
# clean naive/adjusted/baseline separation). This is the version to lock
# for the manuscript: n_reps independent replicates of the same design, so
# we can report mean +/- SE across replicates rather than one run, and check
# that naive - adjusted is consistently positive (not just positive once).
#
# Same design as 04_sim_network_estimation.R, same true network (tau/omega/
# gamma from 00_parameters_uncentered01.R, unchanged across replicates --
# only the person-level randomness, i.e. initial states, sweep updates, and
# P_i draws, differs by replicate). Per replicate:
#   baseline -- n_person people at fixed P=0 (no confounding possible)
#   naive    -- n_person people at P_i ~ Uniform(P_range), estimator omits P
#   adjusted -- identical data to naive, estimator conditions on P_i
#
# n_person per replicate is lower than the single-run pilot's 10000 (see
# n_person note below) -- per Kyuri's suggestion, averaging over replicates
# recovers precision that a single large run would give, without needing
# n_reps x 10000 people.
#
# Outputs
# -------
#   res/revision_2026/sim4/sim4_replicated_by_arm.csv        (rep x arm)
#   res/revision_2026/sim4/sim4_replicated_summary_by_arm.csv (mean +/- SE by arm)
#   res/revision_2026/sim4/sim4_replicated_naive_minus_adjusted.csv (rep-level paired diff)
#   res/revision_2026/sim4/sim4_replicated_summary_diff.csv  (mean +/- SE of the
#     paired diff, plus % of replicates where naive > adjusted)
#   figs/revision_2026/fig_sim4_replicated_diff.pdf          (per-replicate
#     naive-adjusted gap, global strength and mae_true_zero)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(parallel)
})

source("R/revision_2026/utils_uncentered01_model.R")   # simulate_fast_sweep(), fit_edges()
source("R/revision_2026/00_parameters_uncentered01.R")  # tau (+1.3 shift), omega, gamma, symptoms, N

# ------------------------------------------------------------------------
# Design -- PILOT VALUES
# ------------------------------------------------------------------------
T_burn   <- 200L
n_reps   <- 30L
n_person <- 10000L   # per arm, per replicate -- lower than the single-run
                     # pilot's 10000 because averaging over n_reps=30
                     # replicates shrinks the SE of the MEAN by ~sqrt(30)
                     # ~= 5.5x regardless of per-replicate noise; 3000 was
                     # chosen as a middle ground so the whole replicated run
                     # doesn't cost 30x what the 10000-person pilot cost.
                     # If replicate-level SEs below come out too wide, raise
                     # this rather than n_reps (n_reps mainly buys you a
                     # cleaner mean and a %-positive check, not raw precision).

P_fixed <- 0
P_range <- c(-0.6, 0.6)

iu <- which(upper.tri(matrix(0, N, N)))
true_omega_edges <- omega[iu]
pair_names <- outer(symptoms, symptoms, paste, sep = "-")[iu]
sym_i <- symptoms[row(omega)[iu]]
sym_j <- symptoms[col(omega)[iu]]

# ------------------------------------------------------------------------
# One person: simulate forward T_burn sweeps at a given P, return final
# state only (one cross-sectional observation). Same as 04's version.
# ------------------------------------------------------------------------
simulate_person <- function(P) {
  S <- rbinom(N, size = 1, prob = 0.5)
  for (t in seq_len(T_burn)) {
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)
  }
  S
}

# Sequential (not parallel) person-data generator -- used INSIDE each
# replicate, since replicates themselves are the parallel unit below.
# Nesting mclapply inside mclapply is avoided deliberately.
generate_arm_data_seq <- function(Pvec) {
  S <- t(vapply(Pvec, simulate_person, numeric(N)))
  colnames(S) <- symptoms
  list(S = S, Pvec = Pvec)
}

arm_metrics <- function(omega_hat) {
  est <- omega_hat[iu]
  err <- est - true_omega_edges
  nz  <- true_omega_edges != 0
  
  tibble(
    # Keep conventional global strength as a diagnostic,
    # but do not use it as the primary recovery metric
    global_strength_est = sum(abs(est)),
    
    # Primary recovery metric: all generating couplings are positive,
    # so signed sampling error around true zero edges can cancel
    total_coupling_est = sum(est),
    
    # Decompose the absolute-strength statistic
    abs_strength_true_edges = sum(abs(est[nz])),
    spurious_abs_coupling   = sum(abs(est[!nz])),
    
    # Existing recovery diagnostics
    mae_all = mean(abs(err)),
    mae_true_nonzero = mean(abs(err[nz])),
    mae_true_zero = mean(abs(err[!nz])),
    
    # Keep old thresholded metric for supplementary/diagnostic use
    n_phantom_edges_gt_01 = sum(!nz & abs(est) > 0.1)
  )
}

# ------------------------------------------------------------------------
# One replicate: generate baseline + pooled data, fit all 3 arms, return
# per-arm metrics AND the naive-adjusted paired difference for this rep.
# ------------------------------------------------------------------------
run_replicate <- function(rep_id) {
  set.seed(2026L * 1000L + rep_id)

  data_baseline <- generate_arm_data_seq(rep(P_fixed, n_person))
  P_pooled <- runif(n_person, P_range[1], P_range[2])
  data_pooled <- generate_arm_data_seq(P_pooled)

  omega_hat_baseline <- fit_edges(data_baseline$S, Pvec = NULL)
  omega_hat_naive     <- fit_edges(data_pooled$S,   Pvec = NULL)
  omega_hat_adjusted  <- fit_edges(data_pooled$S,   Pvec = data_pooled$Pvec)

  by_arm <- bind_rows(
    arm_metrics(omega_hat_baseline) |> mutate(arm = "baseline", .before = 1),
    arm_metrics(omega_hat_naive)     |> mutate(arm = "naive",     .before = 1),
    arm_metrics(omega_hat_adjusted)  |> mutate(arm = "adjusted",  .before = 1)
  ) |> mutate(rep = rep_id, .before = 1)

  m_naive <- arm_metrics(omega_hat_naive)
  m_adj   <- arm_metrics(omega_hat_adjusted)
  diff_row <- tibble(
    rep = rep_id,
    
    diff_total_coupling =
      m_naive$total_coupling_est - m_adj$total_coupling_est,
    
    diff_spurious_abs_coupling =
      m_naive$spurious_abs_coupling - m_adj$spurious_abs_coupling,
    
    # Keep these for diagnostics
    diff_global_strength =
      m_naive$global_strength_est - m_adj$global_strength_est,
    
    diff_mae_all =
      m_naive$mae_all - m_adj$mae_all,
    
    diff_mae_true_zero =
      m_naive$mae_true_zero - m_adj$mae_true_zero,
    
    diff_mae_true_nonzero =
      m_naive$mae_true_nonzero - m_adj$mae_true_nonzero,
    
    diff_phantom =
      m_naive$n_phantom_edges_gt_01 - m_adj$n_phantom_edges_gt_01
  )

  list(by_arm = by_arm, diff_row = diff_row)
}

# ------------------------------------------------------------------------
# Run all replicates in parallel (replicates are the parallel unit; person
# simulation within each replicate runs sequentially -- see note above).
# ------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores for %d replicates (n_person=%d per arm per replicate)...\n",
            n_cores, n_reps, n_person))

RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

rep_results <- parallel::mclapply(seq_len(n_reps), run_replicate,
                                   mc.cores = n_cores, mc.set.seed = TRUE)

failed <- vapply(rep_results, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) {
  stop(sprintf("%d of %d replicates failed. First error: %s",
               sum(failed), n_reps, attr(rep_results[[which(failed)[1]]], "condition")$message))
}

by_arm_all <- bind_rows(lapply(rep_results, `[[`, "by_arm")) |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")))
diff_all <- bind_rows(lapply(rep_results, `[[`, "diff_row"))

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim4", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

write.csv(by_arm_all, "res/revision_2026/sim4/sim4_replicated_by_arm.csv", row.names = FALSE)
write.csv(diff_all, "res/revision_2026/sim4/sim4_replicated_naive_minus_adjusted.csv", row.names = FALSE)

se <- function(x) sd(x) / sqrt(length(x))

summary_by_arm <- by_arm_all |>
  group_by(arm) |>
  summarise(
    across(
      c(
        global_strength_est,
        total_coupling_est,
        abs_strength_true_edges,
        spurious_abs_coupling,
        mae_all,
        mae_true_nonzero,
        mae_true_zero,
        n_phantom_edges_gt_01
      ),
      list(mean = mean, se = se),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )
write.csv(summary_by_arm, "res/revision_2026/sim4/sim4_replicated_summary_by_arm.csv", row.names = FALSE)

summary_diff <- diff_all |>
  summarise(across(starts_with("diff_"), list(mean = mean, se = se), .names = "{.col}_{.fn}")) |>
  mutate(
    pct_reps_naive_gt_adjusted_total =
      mean(diff_all$diff_total_coupling > 0) * 100,
    
    pct_reps_naive_gt_adjusted_spurious =
      mean(diff_all$diff_spurious_abs_coupling > 0) * 100,
    
    # retain older diagnostics
    pct_reps_naive_gt_adjusted_gs =
      mean(diff_all$diff_global_strength > 0) * 100,
    
    pct_reps_naive_gt_adjusted_mtz =
      mean(diff_all$diff_mae_true_zero > 0) * 100,
    
    n_reps = n_reps
  )
write.csv(summary_diff, "res/revision_2026/sim4/sim4_replicated_summary_diff.csv", row.names = FALSE)

cat("\n=== BY-ARM SUMMARY (mean +/- SE across", n_reps, "replicates) ===\n")
print(summary_by_arm)

cat("\n=== NAIVE - ADJUSTED PAIRED DIFFERENCE (mean +/- SE across", n_reps, "replicates) ===\n")
print(summary_diff)

cat("\nCheck: pct_reps_naive_gt_adjusted_gs and _mtz should be high (ideally\n")
cat("close to 100) -- this is the direct answer to 'is naive > adjusted\n")
cat("consistently, or just in one lucky/unlucky run'. If clearly >50% and the\n")
cat("mean diff is many SEs from zero, Simulation 4 is locked per Kyuri's\n")
cat("criterion. If it's closer to 50%, the confounding effect is not robust\n")
cat("at this n_person/P_range and needs a stronger design (more n_person,\n")
cat("wider P_range), not just more replicates.\n")

# ------------------------------------------------------------------------
# Figure: per-replicate naive-adjusted gap (global strength, mae_true_zero)
# ------------------------------------------------------------------------
diff_long <- diff_all |>
  select(rep, diff_global_strength, diff_mae_true_zero) |>
  pivot_longer(-rep, names_to = "metric", values_to = "value") |>
  mutate(metric = recode(metric,
                          diff_global_strength = "naive - adjusted:\nglobal strength",
                          diff_mae_true_zero  = "naive - adjusted:\nMAE (true-zero edges)"))

p1 <- ggplot(diff_long, aes(x = metric, y = value)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_jitter(width = 0.08, alpha = 0.5, size = 1.6) +
  stat_summary(fun = mean, geom = "point", size = 3, colour = "firebrick") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, colour = "firebrick") +
  labs(x = NULL, y = "naive - adjusted (per replicate)",
       title = sprintf("Simulation 4: confounding gap across %d replicates", n_reps),
       subtitle = "Points above the dashed line = naive more inflated than adjusted, that replicate") +
  theme_classic(base_size = 12)

ggsave("figs/revision_2026/fig_sim4_replicated_diff.pdf", p1, width = 6.5, height = 5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim4/sim4_replicated_by_arm.csv\n")
cat("  res/revision_2026/sim4/sim4_replicated_summary_by_arm.csv\n")
cat("  res/revision_2026/sim4/sim4_replicated_naive_minus_adjusted.csv\n")
cat("  res/revision_2026/sim4/sim4_replicated_summary_diff.csv\n")
cat("  figs/revision_2026/fig_sim4_replicated_diff.pdf\n")
