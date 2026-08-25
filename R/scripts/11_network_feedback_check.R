# ============================================================
# 11_network_feedback_check.R
# ============================================================
# Purpose
# -------
# REPLACES 10_network_window_check.R. That script generated the slow
# context P as an EXOGENOUS Ornstein-Uhlenbeck process (feedback and
# shocks explicitly turned off), so it could only speak to the
# already-known "omitted, threshold-shifting field gets absorbed into
# apparent coupling" story (Kruis & Maris, 2016; Marsman et al., 2018;
# Mukherjee et al., 2024) -- not to this paper's own novel claim, that
# FEEDBACK from sustained symptom elevation onto the slow context
# creates and can amplify that same identification problem.
#
# This script fixes that by giving each simulated individual their own
# private coupled fast/slow system, using the SAME feedback dynamics
# and calibrated parameters as the main scenario simulations (Table 2:
# kappa=0.20, b=0.15, sigma_P=0.04, lambda_m=0.001, m_star=0.25),
# rather than a purpose-built or rescaled process. At every time step,
# a person's full 12-symptom configuration S_t is drawn via the SAME
# exact-enumeration equilibrium sampler used in 09/10, conditional on
# their INSTANTANEOUS P_t (the fast layer is assumed to equilibrate
# each step given the current field, exactly as in 09). Their smoothed
# symptom signal m_slow,t (EWMA of m_t = mean(S_t)) then feeds back
# into their own P via the standard slow-state update (main-text
# Eq. slow_update).
#
# Two conditions, everything else matched:
#   b = 0     : feedback off  (nested reference -- P is then just an
#               ordinary mean-reverting process driven by that person's
#               own path, using the model's OWN calibrated sigma_P,
#               not a rescaled one as in the old script 10)
#   b = 0.15  : feedback on   (the model's actual calibrated value)
#
# Groups differ only in P_base (low: -0.3, high: 1.0), as in 09/10.
# W_true, h, gamma are held fixed across groups AND across the two b
# conditions -- coupling never changes; only the context-generating
# process does.
#
# Timing / "when do we estimate the network"
# --------------------------------------------
# Each person is run for burn_in_steps from a common starting point
# (P_0 = P_base, m_slow_0 = m_star). This burn-in must be long enough
# for TWO things to settle, not just one: (i) P's own mean-reversion
# timescale, tau_P = 1/kappa = 250 steps, and (ii) the EWMA smoothing
# window, 1/lambda_m = 1000 steps, which is what actually lets feedback
# accumulate. We use burn_in_steps = 4000 (4x the slower of the two
# timescales) and print a convergence diagnostic (mean P and m_slow at
# the midpoint, three-quarter point, and end of burn-in) so this can be
# checked empirically rather than assumed -- inspect this before
# trusting the full run. The network-estimation "observation" for each
# person is a SINGLE cross-sectional snapshot: S and P at the FINAL
# step, t = burn_in_steps. This mirrors 09's single-snapshot convention
# exactly and deliberately avoids ALSO varying an observation window
# here -- that is a separate, already-explored question (the retired
# script 10) and mixing it with the feedback question would make it
# impossible to tell which mechanism is doing the work.
#
# Analytic stability check (why this shouldn't blow up): the drift is
# -kappa*(P-P_base) + b*(m_slow(P) - m_star). m_slow is bounded in
# [0,1], so the feedback term is bounded by b*(1-m_star) = 0.1125,
# while the mean-reversion term grows without bound as P moves away
# from P_base. A stable fixed point must therefore exist, with the
# feedback-induced shift in equilibrium P bounded above by
# b*(1-m_star)/kappa = 0.5625 -- a real but bounded effect, not a
# runaway.
#
# Expected asymmetry: m_star = 0.25 sits close to the LOW-baseline
# group's natural resting symptom rate but well below the HIGH-baseline
# group's (cf. Figure 7C: ~0.21-0.28 vs. ~0.61-0.66), so feedback
# should be roughly neutral for the low group and systematically
# push the high group's P upward and apart -- i.e. feedback-driven
# heterogeneity, and hence any resulting network-estimation distortion,
# should show up mainly in the high-baseline group. Worth checking
# against the actual run rather than assuming.
#
# Outputs
# -------
#   res/network_check/network_feedback_check_raw.rds
#   res/network_check/network_feedback_check_summary.csv
#   figs/Figure_S_network_feedback_check.pdf / .png   (supplementary)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(future)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------------------
# 0. Design
# ------------------------------------------------------------------------
# --- pairwise DGP (identical spec to 09/10, so all three checks are
#     directly comparable) ---
N            <- 12L
edge_density <- 0.5
w_lo         <- 0.15
w_hi         <- 0.35
h_lo         <- -1.6
h_hi         <- -1.0
gamma_lo     <- 0.8
gamma_hi     <- 1.4

P_base_low  <- -0.3
P_base_high <- 1.0

# --- slow-state dynamics: SAME calibrated values as the main scenario
#     simulations (Table 2), not rescaled for this check ---
kappa    <- 0.20
dt       <- 0.02
sigma_P  <- 0.04
lambda_m <- 0.001
m_star   <- 0.25
b_conditions <- c(off = 0, on = 0.15)

tau_P_steps   <- 1 / (kappa * dt)     # = 250 steps
burn_in_steps <- 4000L                # ~4x the EWMA memory (1/lambda_m = 1000)

n_per_group <- 1000L  # matches 10's n_per_group -- the quick_test smoke run
                       # at n=100 showed wildly inflated, near-identical
                       # raw/adjusted global strength (~50-80 vs true ~8),
                       # consistent with small-sample logistic-regression
                       # instability (< 4 events/parameter for a 12-
                       # predictor nodewise fit at n=100) rather than a
                       # real effect; 1000 is the smallest n known to be
                       # stable for this exact pairwise DGP (cf. script 10)
n_reps      <- 4L
n_dgp       <- 4L

# Quick smoke-test switch: set TRUE first to sanity-check convergence
# and runtime before committing to the full run above. Deliberately
# does NOT shrink n_per_group -- the first quick_test pass at n=100
# showed inflated, near-identical raw/adjusted global strength (~50-80
# vs true ~8), consistent with small-sample logistic-regression
# instability rather than a real effect. Keeping n_per_group at its
# full value here (while still shrinking reps/dgp/burn-in for speed)
# isolates whether n=1000 alone brings gs back to a sane range
# (comparable to checks 09/10's true~8, raw up to ~26) before
# committing to the full 4-dgp x 4-rep run.
quick_test <- FALSE
if (quick_test) {
  n_reps        <- 1L
  n_dgp         <- 1L
  burn_in_steps <- 1000L
}

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. Data-generating network (duplicated from 09/10 rather than sourced,
#    so this script stays standalone-runnable; keep in sync if the DGP
#    spec ever changes)
# ------------------------------------------------------------------------
build_dgp <- function(dgp_seed) {
  set.seed(dgp_seed)

  W_true <- matrix(0, N, N)
  pairs <- combn(N, 2)
  n_edges <- round(edge_density * ncol(pairs))
  sel <- sample(ncol(pairs), n_edges)
  for (k in sel) {
    i <- pairs[1, k]; j <- pairs[2, k]
    w <- runif(1, w_lo, w_hi)
    W_true[i, j] <- w; W_true[j, i] <- w
  }
  h_vec <- runif(N, h_lo, h_hi)
  gamma_vec <- runif(N, gamma_lo, gamma_hi)

  wtrue_edges <- W_true[iu]
  true_gs <- sum(abs(wtrue_edges))

  states <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
  storage.mode(states) <- "double"
  quad_base <- 0.5 * rowSums((states %*% W_true) * states)

  sample_states <- function(theta_matrix) {
    n_person <- nrow(theta_matrix)
    linear <- states %*% t(theta_matrix)
    logp <- sweep(linear, 1, quad_base, "+")
    logp <- sweep(logp, 2, apply(logp, 2, max), "-")
    p <- exp(logp); p <- sweep(p, 2, colSums(p), "/")
    cdf <- apply(p, 2, cumsum)
    u <- runif(n_person)
    below <- sweep(cdf, 2, u, FUN = "<")
    idx <- colSums(below) + 1L
    states[idx, , drop = FALSE]
  }

  list(W_true = W_true, h_vec = h_vec, gamma_vec = gamma_vec,
       wtrue_edges = wtrue_edges, true_gs = true_gs,
       sample_states = sample_states)
}

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
# 2. Coupled per-person simulation: draws S_t every step conditional on
#    that person's own instantaneous P_t, updates their own m_slow,t via
#    EWMA, updates P_{t+1} via the slow-state feedback SDE. Vectorized
#    across all n_ind people in a group simultaneously.
# ------------------------------------------------------------------------
simulate_coupled_group <- function(net_gen, P_base, n_ind, b, steps) {
  P <- rep(P_base, n_ind)
  m_slow <- rep(m_star, n_ind)   # neutral start

  for (t in seq_len(steps)) {
    theta <- outer(rep(1, n_ind), net_gen$h_vec) + outer(P, net_gen$gamma_vec)
    S <- net_gen$sample_states(theta)
    m_t <- rowMeans(S)
    m_slow <- (1 - lambda_m) * m_slow + lambda_m * m_t
    P <- P + (-kappa * (P - P_base) + b * (m_slow - m_star)) * dt +
      sigma_P * sqrt(dt) * rnorm(n_ind)
  }

  list(S_final = S, P_final = P, m_slow_final = m_slow)
}

# ------------------------------------------------------------------------
# 3. Convergence diagnostic -- RUN AND INSPECT THIS BEFORE THE FULL
#    DESIGN LOOP. Compares mean P and m_slow at the midpoint,
#    three-quarter point, and end of burn-in for the high-baseline,
#    feedback-on case (the condition most likely to still be drifting,
#    per the expected-asymmetry note above). If these are still moving
#    noticeably between the three checkpoints, increase burn_in_steps.
# ------------------------------------------------------------------------
check_convergence <- function(net_gen, P_base, b_val, n_ind = 100, steps = burn_in_steps) {
  P <- rep(P_base, n_ind); m_slow <- rep(m_star, n_ind)
  q1 <- floor(steps / 2); q2 <- floor(3 * steps / 4)
  snap_mid <- snap_3q <- NULL
  for (t in seq_len(steps)) {
    theta <- outer(rep(1, n_ind), net_gen$h_vec) + outer(P, net_gen$gamma_vec)
    S <- net_gen$sample_states(theta)
    m_t <- rowMeans(S)
    m_slow <- (1 - lambda_m) * m_slow + lambda_m * m_t
    P <- P + (-kappa * (P - P_base) + b_val * (m_slow - m_star)) * dt +
      sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t == q1) snap_mid <- list(P = mean(P), m_slow = mean(m_slow))
    if (t == q2) snap_3q  <- list(P = mean(P), m_slow = mean(m_slow))
  }
  list(midpoint = snap_mid, three_quarter = snap_3q,
       final = list(P = mean(P), m_slow = mean(m_slow)))
}

net_gen_diag <- build_dgp(1)
cat("CONVERGENCE DIAGNOSTIC (high-baseline group, feedback ON, dgp=1):\n")
print(check_convergence(net_gen_diag, P_base_high, b_conditions[["on"]]))
cat("If P and m_slow are still drifting noticeably between midpoint / three-quarter\n")
cat("/ final, increase burn_in_steps before trusting the full run below.\n\n")

# ------------------------------------------------------------------------
# 4. One (dgp, b-condition, rep) run
# ------------------------------------------------------------------------
run_condition_feedback <- function(dgp_seed, b_label, rep) {
  net_gen <- build_dgp(dgp_seed)
  b_val <- b_conditions[[b_label]]
  seed <- 200000L * dgp_seed + 1000L * rep + as.integer(b_val * 1000)
  set.seed(seed)

  res_low  <- simulate_coupled_group(net_gen, P_base_low,  n_per_group, b_val, burn_in_steps)
  res_high <- simulate_coupled_group(net_gen, P_base_high, n_per_group, b_val, burn_in_steps)

  Wr_lo <- fit_edges(res_low$S_final)
  Wr_hi <- fit_edges(res_high$S_final)
  Wa_lo <- fit_edges(res_low$S_final,  res_low$P_final)
  Wa_hi <- fit_edges(res_high$S_final, res_high$P_final)

  tibble(
    dgp = dgp_seed, b_label = b_label, b_val = b_val, rep = rep,
    sd_P_lo = sd(res_low$P_final),   sd_P_hi = sd(res_high$P_final),
    mean_P_lo = mean(res_low$P_final), mean_P_hi = mean(res_high$P_final),
    mean_m_lo = mean(rowMeans(res_low$S_final)), mean_m_hi = mean(rowMeans(res_high$S_final)),
    gs_raw_lo = sum(abs(Wr_lo[iu])), gs_raw_hi = sum(abs(Wr_hi[iu])),
    gs_adj_lo = sum(abs(Wa_lo[iu])), gs_adj_hi = sum(abs(Wa_hi[iu])),
    mae_raw_lo = mean(abs(Wr_lo[iu] - net_gen$wtrue_edges)),
    mae_raw_hi = mean(abs(Wr_hi[iu] - net_gen$wtrue_edges)),
    mae_adj_lo = mean(abs(Wa_lo[iu] - net_gen$wtrue_edges)),
    mae_adj_hi = mean(abs(Wa_hi[iu] - net_gen$wtrue_edges)),
    true_gs = net_gen$true_gs
  )
}

# ------------------------------------------------------------------------
# 5. Run: n_dgp x b_conditions x n_reps, parallelized across the design
#    grid (each worker runs one full dgp/b/rep job, i.e. two groups x
#    burn_in_steps sequential coupled updates)
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), b_label = names(b_conditions), rep = seq_len(n_reps))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

results <- design |>
  future_pmap_dfr(
    .f = function(dgp_seed, b_label, rep) run_condition_feedback(dgp_seed, b_label, rep),
    .options = furrr_options(seed = TRUE)
  )

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
saveRDS(results, "res/network_check/network_feedback_check_raw.rds")

# Two-stage averaging convention, matching 09/10: average within each
# DGP draw first, then across draws.
per_dgp <- results |>
  group_by(dgp, b_label) |>
  summarise(across(
    c(sd_P_lo, sd_P_hi, mean_P_lo, mean_P_hi, mean_m_lo, mean_m_hi,
      gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    mean
  ), .groups = "drop")

summary_tbl <- per_dgp |>
  group_by(b_label) |>
  summarise(across(
    c(sd_P_lo, sd_P_hi, mean_P_lo, mean_P_hi, mean_m_lo, mean_m_hi,
      gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop")

write.csv(summary_tbl, "res/network_check/network_feedback_check_summary.csv", row.names = FALSE)
cat("SUMMARY:\n")
print(summary_tbl)

# ------------------------------------------------------------------------
# 6. Figure: feedback off vs on, three panels
# ------------------------------------------------------------------------
col_raw <- "grey15"; col_adj <- "#1565C0"

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      legend.title = element_blank(), legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25)
    )
}

b_labeller <- function(x) factor(x, levels = c("off", "on"), labels = c("Feedback off", "Feedback on"))

# Panel A: realized context heterogeneity by group and feedback condition
panelA <- summary_tbl |>
  transmute(b_label, `Low baseline` = sd_P_lo_mean, `High baseline` = sd_P_hi_mean) |>
  pivot_longer(-b_label, names_to = "group", values_to = "sd_P_g") |>
  mutate(b_label = b_labeller(b_label))

pA <- ggplot(panelA, aes(b_label, sd_P_g, fill = group)) +
  geom_col(position = "dodge") +
  labs(x = NULL, y = expression(SD(P[final])),
       title = "(A) Realized context heterogeneity", fill = NULL) +
  theme_pub()

# Panel B: global strength, symptom-only vs context-adjusted, by feedback condition
panelB <- summary_tbl |>
  transmute(b_label,
            `Symptom-only` = (gs_raw_lo_mean + gs_raw_hi_mean) / 2,
            `Context-adjusted` = (gs_adj_lo_mean + gs_adj_hi_mean) / 2,
            true_gs = true_gs_mean) |>
  pivot_longer(c(`Symptom-only`, `Context-adjusted`), names_to = "estimator", values_to = "gs") |>
  mutate(b_label = b_labeller(b_label))

pB <- ggplot(panelB, aes(b_label, gs, fill = estimator)) +
  geom_col(position = "dodge") +
  geom_hline(aes(yintercept = true_gs), linetype = "dotted", colour = "grey40") +
  scale_fill_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = NULL, y = "Estimated global strength", title = "(B) Global strength", fill = NULL) +
  theme_pub()

# Panel C: recovery error, symptom-only vs context-adjusted, by feedback condition
panelC <- summary_tbl |>
  transmute(b_label,
            `Symptom-only` = (mae_raw_lo_mean + mae_raw_hi_mean) / 2,
            `Context-adjusted` = (mae_adj_lo_mean + mae_adj_hi_mean) / 2) |>
  pivot_longer(c(`Symptom-only`, `Context-adjusted`), names_to = "estimator", values_to = "mae") |>
  mutate(b_label = b_labeller(b_label))

pC <- ggplot(panelC, aes(b_label, mae, fill = estimator)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = NULL, y = expression("Mean "*group("|", hat(W) - W^{true}, "|")),
       title = "(C) Recovery error", fill = NULL) +
  theme_pub()

fig_feedback <- (pA | pB | pC) & theme(legend.position = "bottom")
ggsave("figs/Figure_S_network_feedback_check.pdf", fig_feedback, width = 13.5, height = 4.8)
ggsave("figs/Figure_S_network_feedback_check.png", fig_feedback, width = 13.5, height = 4.8, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/network_check/network_feedback_check_raw.rds\n")
cat("  res/network_check/network_feedback_check_summary.csv\n")
cat("  figs/Figure_S_network_feedback_check.pdf (+ .png)  <- supplementary\n")

# ------------------------------------------------------------------------
# 7. Companion network-graph figure -- mirrors Figure 8's true/raw/
#    adjusted graph comparison, which is the visual the panel-B/C summary
#    numbers above don't give you on their own. One representative DGP
#    (dgp_seed = 1), high-baseline group (where feedback's effect
#    concentrates -- see the mean_P_hi shift in the summary above),
#    comparing feedback off vs on for both estimators. Re-simulated here
#    with fixed seeds rather than pulled from the design loop above,
#    since that loop only kept scalar summaries (gs/mae), not the full
#    estimated edge matrices.
# ------------------------------------------------------------------------
suppressPackageStartupMessages(library(qgraph))

dgp_graph <- build_dgp(1)
set.seed(999001)
res_off_g <- simulate_coupled_group(dgp_graph, P_base_high, n_per_group, b_conditions[["off"]], burn_in_steps)
set.seed(999002)
res_on_g  <- simulate_coupled_group(dgp_graph, P_base_high, n_per_group, b_conditions[["on"]],  burn_in_steps)

Wr_off_g <- fit_edges(res_off_g$S_final)
Wr_on_g  <- fit_edges(res_on_g$S_final)
Wa_off_g <- fit_edges(res_off_g$S_final, res_off_g$P_final)
Wa_on_g  <- fit_edges(res_on_g$S_final,  res_on_g$P_final)

# Save all four estimated matrices + true W, even though the figure below
# only plots three -- keeps the off-condition numbers available for a
# quick control-check sentence in text without re-simulating, and means
# re-plotting later (e.g. restyling) never requires re-running this block.
saveRDS(list(W_true = dgp_graph$W_true, Wr_off = Wr_off_g, Wr_on = Wr_on_g,
             Wa_off = Wa_off_g, Wa_on = Wa_on_g),
        "res/network_check/network_feedback_check_graph_mats.rds")

# Lead figure: feedback-on only (the real, calibrated model), mirroring
# Figure 8's original 3-panel format exactly -- True / Symptom-only /
# Context-adjusted. The off-condition is reported as a control-check
# number in text (Appendix D), not as extra panels here.
mats_g   <- list(dgp_graph$W_true, Wr_on_g, Wa_on_g)
titles_g <- c("True", "Symptom-only", "Context-adjusted")

max_edge_g <- max(sapply(mats_g, function(m) max(abs(m))))
L <- qgraph(dgp_graph$W_true, layout = "spring", DoNotPlot = TRUE)$layout

plot_graph_panel <- function() {
  layout(matrix(1:3, nrow = 1))
  for (k in seq_along(mats_g)) {
    qgraph(mats_g[[k]], layout = L, maximum = max_edge_g, threshold = 0.02,
           posCol = "#1565C0", negCol = "#C62828",
           labels = paste0("S", 1:N), label.cex = 1.1, vsize = 8,
           title = titles_g[k], title.cex = 1.1, mar = c(2, 2, 4, 2))
  }
}

pdf("figs/Figure_S_network_feedback_check_graphs.pdf", width = 10, height = 4)
plot_graph_panel()
dev.off()

png("figs/Figure_S_network_feedback_check_graphs.png", width = 10, height = 4, units = "in", res = 200)
plot_graph_panel()
dev.off()

cat("  figs/Figure_S_network_feedback_check_graphs.pdf (+ .png)  <- companion network-graph figure (3-panel, feedback-on)\n")
cat("  res/network_check/network_feedback_check_graph_mats.rds   <- all 5 matrices, for re-plotting without re-simulating\n")

# ------------------------------------------------------------------------
# 8. RE-PLOT ONLY -- run this block alone (after loading qgraph and the
#    RDS above) to restyle the figure without re-simulating anything:
#
#   library(qgraph)
#   m <- readRDS("res/network_check/network_feedback_check_graph_mats.rds")
#   mats_g <- list(m$W_true, m$Wr_on, m$Wa_on)
#   titles_g <- c("True", "Symptom-only", "Context-adjusted")
#   max_edge_g <- max(sapply(mats_g, function(x) max(abs(x))))
#   L <- qgraph(m$W_true, layout = "spring", DoNotPlot = TRUE)$layout
#   # ... then the same plot_graph_panel()/pdf()/png() calls as above
# ------------------------------------------------------------------------

