# ============================================================
# R/revision_2026/04_sim_network_estimation.R
# ============================================================
# Simulation 4: does pooling cross-sectional data across different context
# (P) levels -- without conditioning on P -- make the ESTIMATED symptom
# network look more strongly coupled than the TRUE generative network,
# purely because every symptom shares the same confounder (each moves with
# P via its own gamma_i)? This is the paper's central claim made concrete
# and testable: "when context looks like coupling."
#
# Unlike Sims 1-3, this is a cross-sectional (between-person) design, not a
# trajectory. Each of n_person independent "people" is simulated forward
# from a random start for T_burn sweeps at their own P_i, and only their
# FINAL state is kept as their one observed symptom profile -- analogous to
# a single-timepoint self-report survey, not a repeated-measures chain.
#
# Same fast layer, same tau/omega/gamma/N=9 as Simulations 1-3 (validated,
# +1.3 tau shift). The TRUE network is `omega` from 00_parameters_uncentered01.R.
#
# Three arms, same n_person and T_burn throughout so only the design
# differs between them:
#   baseline -- single context: every person simulated at the same fixed
#               P=0 (Sim 1's middle condition). No context variance, so no
#               confounding is possible here by construction. This isolates
#               ordinary finite-sample estimation error as a reference point
#               for what the naive/adjusted arms should be compared against.
#   naive    -- pooled, no adjustment: each person has their own P_i drawn
#               from Uniform(-0.6, 0.6) (Sim 1's low/high range), but the
#               estimator (fit_edges with Pvec=NULL) does not use P_i at
#               all -- as if the analyst pooled data across contexts without
#               having measured or modeled them.
#   adjusted -- pooled, context-adjusted: identical data to "naive" (same
#               P_i draws, same seed), but fit_edges is given Pvec so each
#               nodewise regression conditions on P_i.
#
# Prediction: naive should show inflated apparent connectivity (edges that
# are truly zero come out systematically nonzero, and true edges may be
# over- or under-estimated) relative to both baseline and adjusted; adjusted
# should look close to baseline, showing the inflation in "naive" is a
# confounding artifact of pooling over unmodeled context, not an estimation
# artifact of the regression itself.
#
# n_person, T_burn, and the P_i spread are PILOT VALUES -- calibrate against
# output like Sim 1's tau shift and Sim 3's b grid.
#
# Outputs
# -------
#   res/revision_2026/sim4/sim4_edge_estimates.csv   (one row per symptom
#     pair x arm: true omega, estimated omega)
#   res/revision_2026/sim4/sim4_summary.csv          (MAE, global strength
#     per arm)
#   figs/revision_2026/fig_sim4_network_estimation.pdf (estimated vs true
#     edge weight, faceted by arm, coloured by whether the true edge is
#     zero or nonzero)
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
T_burn    <- 200L    # same validated burn-in as Sim 1 (one fast sweep per
                      # step; each person's chain starts from a random state)
n_person  <- 10000L   # raised from 2000 -- at 2000, unregularized nodewise
                      # logistic regression (no L1/EBIC shrinkage, several
                      # low-base-rate symptoms like suicidal ~7.5%) had a
                      # high noise floor: even the baseline arm (P fixed,
                      # zero confounding possible) showed global strength
                      # 6.50 vs true 3.05, almost as inflated as naive's
                      # 7.01, burying the naive-vs-adjusted confounding
                      # signal under generic estimation noise. 10000 people
                      # is cheap here (each person is just T_burn sweeps,
                      # no trajectory) and run in parallel below, so raise
                      # n rather than add regularization for now.

P_fixed    <- 0          # baseline arm: everyone simulated at this P
P_range    <- c(-0.6, 0.6)  # naive/adjusted arms: P_i ~ Uniform(P_range),
                             # matching Sim 1's low/high context range

iu <- which(upper.tri(matrix(0, N, N)))   # upper-triangle indices, reused
                                            # for every N x N matrix below
true_omega_edges <- omega[iu]

# ------------------------------------------------------------------------
# One person: simulate forward T_burn sweeps at a given P, return final
# state only (one cross-sectional observation).
# ------------------------------------------------------------------------
simulate_person <- function(P) {
  S <- rbinom(N, size = 1, prob = 0.5)
  for (t in seq_len(T_burn)) {
    S <- simulate_fast_sweep(S, tau, omega, gamma, P)
  }
  S
}

# ------------------------------------------------------------------------
# Generate one arm's data (n_person independent people), in parallel.
# Returns S (n_person x N matrix) and Pvec (n_person, the true P_i used --
# kept even for the "naive" arm since we need it to know what confounding
# was present, just not passed to the naive estimator).
# ------------------------------------------------------------------------
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1, na.rm = TRUE)
cat(sprintf("Using %d cores.\n", n_cores))
RNGkind("L'Ecuyer-CMRG")
set.seed(2026L)

generate_arm_data <- function(Pvec) {
  n <- length(Pvec)
  rows <- parallel::mclapply(seq_len(n), function(i) simulate_person(Pvec[i]),
                              mc.cores = n_cores, mc.set.seed = TRUE)
  failed <- vapply(rows, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) stop(sprintf("%d of %d person-simulations failed.", sum(failed), n))
  S <- do.call(rbind, rows)
  colnames(S) <- symptoms
  list(S = S, Pvec = Pvec)
}

cat("Simulating baseline arm (single context, P=0)...\n")
data_baseline <- generate_arm_data(rep(P_fixed, n_person))

cat("Simulating naive/adjusted arms (pooled, P ~ Uniform(-0.6, 0.6))...\n")
P_pooled <- runif(n_person, P_range[1], P_range[2])
data_pooled <- generate_arm_data(P_pooled)   # same data feeds both naive and adjusted

# ------------------------------------------------------------------------
# Estimate networks
# ------------------------------------------------------------------------
cat("Fitting networks (nodewise logistic regression, symmetrized)...\n")

omega_hat_baseline <- fit_edges(data_baseline$S, Pvec = NULL)
omega_hat_naive     <- fit_edges(data_pooled$S,   Pvec = NULL)                 # P NOT given to the estimator
omega_hat_adjusted  <- fit_edges(data_pooled$S,   Pvec = data_pooled$Pvec)     # P given to the estimator

# ------------------------------------------------------------------------
# Save + summarize
# ------------------------------------------------------------------------
dir.create("res/revision_2026/sim4", recursive = TRUE, showWarnings = FALSE)
dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

pair_names <- outer(symptoms, symptoms, paste, sep = "-")[iu]
sym_i <- symptoms[row(omega)[iu]]
sym_j <- symptoms[col(omega)[iu]]

edge_tbl <- bind_rows(
  tibble(arm = "baseline", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_baseline[iu]),
  tibble(arm = "naive", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_naive[iu]),
  tibble(arm = "adjusted", symptom_i = sym_i, symptom_j = sym_j,
         true_omega = true_omega_edges, est_omega = omega_hat_adjusted[iu])
) |>
  mutate(arm = factor(arm, levels = c("baseline", "naive", "adjusted")),
         true_nonzero = true_omega != 0,
         error = est_omega - true_omega)

write.csv(edge_tbl, "res/revision_2026/sim4/sim4_edge_estimates.csv", row.names = FALSE)

saveRDS(list(edge_tbl = edge_tbl,
             omega_hat = list(baseline = omega_hat_baseline, naive = omega_hat_naive,
                               adjusted = omega_hat_adjusted),
             true_omega = omega,
             data = list(baseline = data_baseline, pooled = data_pooled),
             params = list(T_burn = T_burn, n_person = n_person,
                            P_fixed = P_fixed, P_range = P_range)),
        "res/revision_2026/sim4/sim4_raw.rds")

summary_tbl <- edge_tbl |>
  group_by(arm) |>
  summarise(
    mae_all      = mean(abs(error)),
    mae_true_nonzero = mean(abs(error[true_nonzero])),
    mae_true_zero     = mean(abs(error[!true_nonzero])),
    global_strength_true = sum(abs(true_omega)),
    global_strength_est  = sum(abs(est_omega)),
    n_phantom_edges_gt_01 = sum(!true_nonzero & abs(est_omega) > 0.1),
    .groups = "drop"
  )
write.csv(summary_tbl, "res/revision_2026/sim4/sim4_summary.csv", row.names = FALSE)

cat("\n=== SIMULATION 4 SUMMARY (by arm) ===\n")
print(summary_tbl)

cat("\nCompare 'naive' against 'baseline' and 'adjusted': if pooling across\n")
cat("context without conditioning on P inflates global_strength_est and/or\n")
cat("mae_true_zero (spurious nonzero edges among truly-zero pairs) relative\n")
cat("to both baseline and adjusted, that's the confounding effect this\n")
cat("simulation is designed to show. If 'naive' looks similar to 'baseline'/\n")
cat("'adjusted', the P_range spread may be too small to generate visible\n")
cat("confounding -- widen P_range and rerun. If 'adjusted' does not track\n")
cat("'baseline' closely, something is wrong with the P-conditioning itself\n")
cat("(not just a calibration issue).\n")

# ------------------------------------------------------------------------
# Figure: estimated vs true edge weight, one panel per arm
# ------------------------------------------------------------------------
p1 <- ggplot(edge_tbl, aes(true_omega, est_omega, colour = true_nonzero)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                       labels = c(`FALSE` = "true edge = 0", `TRUE` = "true edge != 0"),
                       name = NULL) +
  facet_wrap(~arm, nrow = 1) +
  labs(x = expression(paste("True ", omega[ij])), y = expression(paste("Estimated ", hat(omega)[ij])),
       title = "Simulation 4: omitted-context confounding in estimated symptom networks",
       subtitle = "Dashed line = perfect recovery. Points off the line among grey (true-zero) pairs are spurious edges.") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figs/revision_2026/fig_sim4_network_estimation.pdf", p1, width = 11, height = 4.5)

cat("\nDone. Files:\n")
cat("  res/revision_2026/sim4/sim4_raw.rds\n")
cat("  res/revision_2026/sim4/sim4_edge_estimates.csv\n")
cat("  res/revision_2026/sim4/sim4_summary.csv\n")
cat("  figs/revision_2026/fig_sim4_network_estimation.pdf\n")
