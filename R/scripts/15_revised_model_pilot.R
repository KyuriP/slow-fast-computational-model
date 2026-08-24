# ============================================================
# 15_revised_model_pilot.R
# ============================================================
# Purpose
# -------
# Pilot grid for the REVISED core model (see Simulation_Spec_Sims1-3_2026-08.md):
#
#   logit Pr(S_i = 1 | S_-i, P) = tau_i + sum_j!=i omega_ij S_j + gamma_i P
#
# No beta. No J(m - 1/2). No shared threshold. This replaces the mean-field
# Curie-Weiss fast layer as the paper's core model, per Denny's feedback.
#
# Unlike 13_network_full_check.R, P is NOT a time-evolving SDE here -- Sims
# 1-3 treat P as a per-person draw (or fixed per-group constant), so there is
# no burn-in loop. This makes the pilot much cheaper than script 13: each
# condition is one direct sample from the joint distribution, not an 8000+
# step simulation.
#
# This script runs the pilot grid only: it exists to (a) check the model
# doesn't saturate at all-symptoms-off/all-symptoms-on across a plausible
# mu_P x sigma_P range, and (b) sanity-check the estimator: at sigma_P = 0,
# the omitted-context gap in GS/MAE should be ~0 (this is a correctness
# check on the code, not just a calibration check -- if it fails, something
# is wrong in fit_edges() or the sampler, not just the parameter choices).
#
# Outputs
# -------
#   res/network_check/revised_model_pilot_raw.rds
#   res/network_check/revised_model_pilot_summary.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
suppressPackageStartupMessages(library(matrixStats))

# ------------------------------------------------------------------------
# 0. Design -- see Simulation_Spec_Sims1-3_2026-08.md "Pilot grid" section
# ------------------------------------------------------------------------
# N = 12 for now (reuses the exact-enumeration approach cheaply: 2^12 =
# 4096 states). Switch to 14 later only if you want an exact match to
# Cramer et al. (2016)'s symptom count -- flagged as an open decision, not
# assumed here.
N            <- 12L
edge_density <- 0.40                 # midpoint of the 0.30-0.50 range in the spec

# tau_i (per-symptom threshold): placeholder range per the spec, pending
# real Cramer et al. (2016) values (task: ask Denny for the fitted VATSPUD
# thresholds -- not machine-readable from the PLOS page, only shown as a
# figure). Uniform(-3.0, -1.0) reflects the qualitative pattern the paper
# describes (rare symptoms like thoughts of death get much more negative
# thresholds than common ones like fatigue), not the literal fitted values.
tau_lo <- -3.0
tau_hi <- -1.0

# omega_ij (nonzero edges): Uniform(0.20, 0.40) per the revised spec
w_lo <- 0.20
w_hi <- 0.40

# gamma_i (context sensitivity)
gamma_lo <- 0.6
gamma_hi <- 1.2

# Pilot grid over mu_P x sigma_P
mu_P_grid    <- c(-1.0, 0.0, 1.0)
sigma_P_grid <- c(0, 0.25, 0.5, 1.0)

n_per_cell <- 1000L   # pilot scale; bump to 2000 for the final run
n_reps     <- 20L     # pilot scale; bump to 100 for the final run
n_dgp      <- 3L      # independently drawn data-generating networks

high_burden_cutoff <- N / 2   # Pr(D >= N/2); flagged as a decision in the spec

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. Data-generating network -- same structure as build_dgp() in
#    13_network_full_check.R, renamed tau/omega to match the revised
#    model's notation, ranges updated per the spec above.
# ------------------------------------------------------------------------
build_dgp <- function(dgp_seed) {
  set.seed(dgp_seed)

  omega_true <- matrix(0, N, N)
  pairs <- combn(N, 2)
  n_edges <- round(edge_density * ncol(pairs))
  sel <- sample(ncol(pairs), n_edges)
  for (k in sel) {
    i <- pairs[1, k]; j <- pairs[2, k]
    w <- runif(1, w_lo, w_hi)
    omega_true[i, j] <- w; omega_true[j, i] <- w
  }
  tau_vec   <- runif(N, tau_lo, tau_hi)
  gamma_vec <- runif(N, gamma_lo, gamma_hi)

  omega_true_edges <- omega_true[iu]
  true_gs <- sum(abs(omega_true_edges))

  states <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
  storage.mode(states) <- "double"
  quad_base <- 0.5 * rowSums((states %*% omega_true) * states)

  # Exact joint sampler: given a per-person linear term theta_i = tau_i +
  # gamma_i * P_n, returns one exact draw from the full N-symptom joint
  # distribution implied by the pairwise model. No time loop, no burn-in --
  # this is a single-shot draw from the model's own stationary/cross-
  # sectional distribution, which is all Sims 1-3 need (P is not dynamic).
  sample_states <- function(theta_matrix) {
    n_person <- nrow(theta_matrix)
    linear <- states %*% t(theta_matrix)
    logp <- linear + quad_base
    logp <- logp - rep(matrixStats::colMaxs(logp), each = nrow(states))
    p <- exp(logp)
    p <- p / rep(colSums(p), each = nrow(states))
    cdf <- matrixStats::colCumsums(p)
    u <- runif(n_person)
    below <- cdf < rep(u, each = nrow(states))
    idx <- colSums(below) + 1L
    states[idx, , drop = FALSE]
  }

  list(omega_true = omega_true, tau_vec = tau_vec, gamma_vec = gamma_vec,
       omega_true_edges = omega_true_edges, true_gs = true_gs,
       sample_states = sample_states)
}

# ------------------------------------------------------------------------
# 2. Estimator -- identical to fit_edges() in 13_network_full_check.R.
#    Nodewise logistic regression (base R glm(family = binomial())),
#    symmetrized by averaging directed coefficients. Named explicitly here
#    (and in the Methods text) per Denny's "what estimator/software"
#    comment.
# ------------------------------------------------------------------------
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
# 3. One (dgp, mu_P, sigma_P, rep) pilot condition
# ------------------------------------------------------------------------
run_pilot_condition <- function(dgp_seed, mu_P, sigma_P, rep) {
  net_gen <- build_dgp(dgp_seed)
  seed <- 700000L * dgp_seed + 5000L * rep + round((mu_P + 10) * 1000) + round(sigma_P * 10000)
  set.seed(seed)

  P_n <- if (sigma_P == 0) rep(mu_P, n_per_cell) else rnorm(n_per_cell, mu_P, sigma_P)
  theta <- outer(rep(1, n_per_cell), net_gen$tau_vec) + outer(P_n, net_gen$gamma_vec)
  S <- net_gen$sample_states(theta)

  D <- rowSums(S)
  W_omit <- fit_edges(S)
  W_adj  <- fit_edges(S, P_n)

  # saturation check: proportion of individuals with all symptoms off/on
  prop_all_off <- mean(D == 0)
  prop_all_on  <- mean(D == N)

  tibble(
    dgp = dgp_seed, mu_P = mu_P, sigma_P = sigma_P, rep = rep,
    mean_D = mean(D), var_D = var(D),
    pr_high_burden = mean(D >= high_burden_cutoff),
    prop_all_off = prop_all_off, prop_all_on = prop_all_on,
    mean_p_active = mean(S),   # overall mean activation prob across nodes/individuals
    gs_omit = sum(abs(W_omit[iu])), gs_adj = sum(abs(W_adj[iu])),
    mae_omit = mean(abs(W_omit[iu] - net_gen$omega_true_edges)),
    mae_adj  = mean(abs(W_adj[iu]  - net_gen$omega_true_edges)),
    true_gs = net_gen$true_gs
  )
}

# ------------------------------------------------------------------------
# 4. Run the pilot grid
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), mu_P = mu_P_grid,
                       sigma_P = sigma_P_grid, rep = seq_len(n_reps))
cat(sprintf("Pilot design has %d conditions (no time loop -- should run in seconds to low minutes).\n",
            nrow(design)))

t0 <- Sys.time()
results <- design |>
  pmap_dfr(function(dgp_seed, mu_P, sigma_P, rep) run_pilot_condition(dgp_seed, mu_P, sigma_P, rep))
t1 <- Sys.time()
cat(sprintf("Pilot wall time: %.1f %s\n", as.numeric(difftime(t1, t0, units = "auto")),
            units(difftime(t1, t0, units = "auto"))))

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
saveRDS(results, "res/network_check/revised_model_pilot_raw.rds")

# ------------------------------------------------------------------------
# 5. Summary + the two pilot checks called out in the spec
# ------------------------------------------------------------------------
summary_tbl <- results |>
  group_by(mu_P, sigma_P) |>
  summarise(
    mean_D = mean(mean_D), var_D = mean(var_D), pr_high_burden = mean(pr_high_burden),
    prop_all_off = mean(prop_all_off), prop_all_on = mean(prop_all_on),
    mean_p_active = mean(mean_p_active),
    gs_omit_mean = mean(gs_omit), gs_adj_mean = mean(gs_adj),
    gap_gs_mean = mean(gs_omit - gs_adj), gap_gs_se = sd(gs_omit - gs_adj) / sqrt(n()),
    mae_omit_mean = mean(mae_omit), mae_adj_mean = mean(mae_adj),
    gap_mae_mean = mean(mae_omit - mae_adj), gap_mae_se = sd(mae_omit - mae_adj) / sqrt(n()),
    true_gs_mean = mean(true_gs),
    .groups = "drop"
  )

write.csv(summary_tbl, "res/network_check/revised_model_pilot_summary.csv", row.names = FALSE)

cat("\n=== SATURATION CHECK (want prop_all_off / prop_all_on small, and mean_p_active away from 0/1) ===\n")
print(summary_tbl |> select(mu_P, sigma_P, mean_D, mean_p_active, prop_all_off, prop_all_on))

cat("\n=== SIGMA_P = 0 CORRECTNESS CHECK (want gap_gs_mean and gap_mae_mean ~ 0 here) ===\n")
print(summary_tbl |> filter(sigma_P == 0) |> select(mu_P, sigma_P, gap_gs_mean, gap_gs_se, gap_mae_mean, gap_mae_se))

cat("\n=== OMITTED-CONTEXT GAP GROWTH (want gap_gs_mean and gap_mae_mean to increase with sigma_P) ===\n")
print(summary_tbl |> select(mu_P, sigma_P, gap_gs_mean, gap_gs_se, gap_mae_mean, gap_mae_se))

cat("\nDone. Files:\n")
cat("  res/network_check/revised_model_pilot_raw.rds\n")
cat("  res/network_check/revised_model_pilot_summary.csv\n")
cat("\nIf saturation and sigma_P=0 checks pass, next step is building the full\n")
cat("Simulation 1 / 2 / 3 scripts from Simulation_Spec_Sims1-3_2026-08.md at\n")
cat("final scale (n=2000/cell, 100 reps).\n")
