# ============================================================
# 09_network_estimation_check.R
# ============================================================
# Purpose
# -------
# Compact network-estimation check for the "Network-estimation check"
# subsection of the manuscript (Simulation design / Results).
#
# Data-generating model (pairwise extension of the fast layer, used only
# for this check):
#
#   logit Pr(S_i = 1 | S_-i, P) = h_i + sum_j W_true_ij S_j + gamma_i * P
#
# W_true, h_i, gamma_i are held FIXED across two groups. The groups differ
# only in the distribution of the slow context P (low-baseline group:
# mean mu_low; high-baseline group: mean mu_high). We vary the amount of
# slow-context variation present in the sampled data (SD_P) and ask
# whether a symptom-only ("raw") network estimator absorbs the shared,
# unmodeled slow-context effect into apparent group differences in
# estimated coupling, relative to a "context-adjusted" estimator that
# includes P as a covariate.
#
# NOTE ON PROVENANCE
# -------------------
# This script was prototyped and numerically validated in Python
# (no R available in that sandbox), then ported here by hand to match
# this repo's conventions (base R + qgraph, nodewise logistic regression
# instead of IsingSampler/IsingFit, since IsingSampler does not support a
# distinct threshold vector per observation, which we need because each
# simulated person has their own P_i). The data-generating network below
# (W_true, h_vec, gamma_vec) is HARD-CODED to the exact values used in the
# validated Python run, so the results here should reproduce the numbers
# already reported in the manuscript draft. This script has NOT been
# executed in this environment (R was not available) -- please run it
# locally and sanity-check against network_estimation_check_summary.csv
# (Python output) before trusting the numbers for the paper.
#
# Sanity check to look for when you run this: at SD_P = 0 (first row),
# the raw and context-adjusted estimates should be virtually identical
# for every quantity (P is ~constant within each group there, so it is
# collinear with the intercept and "adjusting" for it does ~nothing).
#
# Outputs
# -------
#   res/network_check/network_check_raw.rds       (replicate-level)
#   res/network_check/network_check_summary.csv    (condition-level means/SEs)
#   img/Figure7_network_estimation_check.pdf       (3-panel summary figure)
#   img/Figure8_network_estimation_check_graphs.pdf (qgraph network comparison)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(future)
  library(ggplot2)
  library(qgraph)
})

set.seed(2024)

# ------------------------------------------------------------------------
# 0. Fixed data-generating network (identical to the validated Python run)
# ------------------------------------------------------------------------
N <- 8L

W_true <- matrix(c(
  0.0000000, 0.0000000, 0.2036731, 0.3114364, 0.3451244, 0.0000000, 0.0000000, 0.1912688,
  0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.2056083, 0.0000000, 0.0000000,
  0.2036731, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1926315,
  0.3114364, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.3098857, 0.2150699, 0.2693645,
  0.3451244, 0.0000000, 0.0000000, 0.0000000, 0.0000000, 0.1509259, 0.2430238, 0.2048490,
  0.0000000, 0.2056083, 0.0000000, 0.3098857, 0.1509259, 0.0000000, 0.2385451, 0.3249916,
  0.0000000, 0.0000000, 0.0000000, 0.2150699, 0.2430238, 0.2385451, 0.0000000, 0.0000000,
  0.1912688, 0.0000000, 0.1926315, 0.2693645, 0.2048490, 0.3249916, 0.0000000, 0.0000000
), nrow = N, byrow = TRUE)

h_vec <- c(-1.4391623, -1.5574709, -1.3196747, -1.4414767,
           -1.0666348, -1.4282090, -1.1357398, -1.3076531)

gamma_vec <- c(1.0808114, 1.3789581, 1.3389364, 0.8474206,
               0.9471226, 0.9108722, 1.3432849, 1.1322992)

symptom_labels <- paste0("S", seq_len(N))

mu_low  <- -0.3
mu_high <- 1.0
sdP_levels   <- c(0.0, 0.1, 0.3, 0.5, 0.75, 1.0)
n_per_group  <- 2500L
n_reps       <- 30L

iu <- which(upper.tri(matrix(0, N, N)))
wtrue_edges <- W_true[iu]
true_global_strength <- mean(wtrue_edges)

# ------------------------------------------------------------------------
# 1. Exact sampler for the pairwise model (N = 8 -> 256 states; enumerate
#    exactly rather than via MCMC, since each of the n_per_group people
#    has their own field h + gamma*P_i and IsingSampler does not accept a
#    per-observation threshold vector).
# ------------------------------------------------------------------------
states <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
storage.mode(states) <- "double"
quad_base <- 0.5 * rowSums((states %*% W_true) * states)

#' Exact-enumeration sampler for the pairwise Ising-with-covariate model
#'
#' @param theta_matrix n_person x N matrix of per-person effective fields
#'   (h_vec + gamma_vec * P_i for each row).
#' @return n_person x N binary matrix of sampled symptom states.
sample_states <- function(theta_matrix) {
  n_person <- nrow(theta_matrix)
  linear <- states %*% t(theta_matrix)                 # (2^N x n_person)
  logp <- sweep(linear, 1, quad_base, "+")
  logp <- sweep(logp, 2, apply(logp, 2, max), "-")
  p <- exp(logp)
  p <- sweep(p, 2, colSums(p), "/")
  cdf <- apply(p, 2, cumsum)                            # (2^N x n_person)
  u <- runif(n_person)
  below <- sweep(cdf, 2, u, FUN = "<")
  idx <- colSums(below) + 1L                            # first state with cdf >= u
  states[idx, , drop = FALSE]
}

# ------------------------------------------------------------------------
# 2. Nodewise logistic-regression network estimator (symmetrized), with
#    an optional covariate P for the "context-adjusted" estimator.
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
    if (!is.null(Pvec)) {
      df$P <- Pvec
      rhs <- c(rhs, "P")
    }
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
# 3. Main simulation loop
# ------------------------------------------------------------------------
run_condition <- function(sdP, rep) {
  seed <- 10000L + rep + round(1000 * sdP)
  set.seed(seed)

  sd_use <- max(sdP, 1e-6)
  P_low  <- rnorm(n_per_group, mu_low,  sd_use)
  P_high <- rnorm(n_per_group, mu_high, sd_use)

  theta_low  <- outer(rep(1, n_per_group), h_vec) + outer(P_low,  gamma_vec)
  theta_high <- outer(rep(1, n_per_group), h_vec) + outer(P_high, gamma_vec)

  S_low  <- sample_states(theta_low)
  S_high <- sample_states(theta_high)

  Wr_lo <- fit_edges(S_low)
  Wr_hi <- fit_edges(S_high)
  Wa_lo <- fit_edges(S_low,  P_low)
  Wa_hi <- fit_edges(S_high, P_high)

  d_raw <- mean(abs(Wr_hi[iu] - Wr_lo[iu]))
  d_adj <- mean(abs(Wa_hi[iu] - Wa_lo[iu]))

  tibble(
    sdP = sdP, rep = rep,
    d_raw = d_raw, d_adj = d_adj,
    gs_raw_lo = mean(Wr_lo[iu]), gs_raw_hi = mean(Wr_hi[iu]),
    gs_adj_lo = mean(Wa_lo[iu]), gs_adj_hi = mean(Wa_hi[iu]),
    mae_raw_lo = mean(abs(Wr_lo[iu] - wtrue_edges)),
    mae_raw_hi = mean(abs(Wr_hi[iu] - wtrue_edges)),
    mae_adj_lo = mean(abs(Wa_lo[iu] - wtrue_edges)),
    mae_adj_hi = mean(abs(Wa_hi[iu] - wtrue_edges)),
    m_low = mean(S_low), m_high = mean(S_high)
  )
}

design <- expand_grid(sdP = sdP_levels, rep = seq_len(n_reps))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

results <- design |>
  future_pmap_dfr(
    .f = function(sdP, rep) run_condition(sdP, rep),
    .options = furrr_options(seed = TRUE)
  )

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("img", showWarnings = FALSE)

saveRDS(results, "res/network_check/network_check_raw.rds")

summary_tbl <- results |>
  group_by(sdP) |>
  summarise(across(
    c(d_raw, d_adj, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, m_low, m_high),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop")

write.csv(summary_tbl, "res/network_check/network_check_summary.csv", row.names = FALSE)
print(summary_tbl)
cat(sprintf("true_global_strength = %.4f\n", true_global_strength))

# ------------------------------------------------------------------------
# 4. Figure 7: 3-panel summary (apparent group difference / global
#    strength by group / recovery error vs W_true)
# ------------------------------------------------------------------------
col_green <- "#2E7D32"   # low-baseline group (matches rest of paper)
col_red   <- "#C62828"   # high-baseline group
col_raw   <- "grey15"
col_adj   <- "#1565C0"

panelA <- summary_tbl |>
  transmute(
    sdP,
    diff_raw = abs(gs_raw_hi_mean - gs_raw_lo_mean),
    diff_adj = abs(gs_adj_hi_mean - gs_adj_lo_mean),
    se_raw = sqrt(gs_raw_hi_se^2 + gs_raw_lo_se^2),
    se_adj = sqrt(gs_adj_hi_se^2 + gs_adj_lo_se^2)
  ) |>
  pivot_longer(c(diff_raw, diff_adj), names_to = "estimator", values_to = "diff") |>
  mutate(
    se = if_else(estimator == "diff_raw", se_raw, se_adj),
    estimator = if_else(estimator == "diff_raw", "Symptom-only (raw)", "Context-adjusted")
  )

pA <- ggplot(panelA, aes(sdP, diff, colour = estimator, linetype = estimator)) +
  geom_errorbar(aes(ymin = diff - se, ymax = diff + se), width = 0.02) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only (raw)" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = expression(SD[P]~"(slow-context variation)"),
       y = "Apparent group difference\nin estimated network strength",
       title = "(A) Apparent group network difference", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.85))

panelB <- summary_tbl |>
  select(sdP, gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean) |>
  pivot_longer(-sdP, names_to = "series", values_to = "gs") |>
  mutate(
    group = if_else(grepl("_lo_", series), "Low baseline", "High baseline"),
    estimator = if_else(grepl("raw", series), "raw", "adjusted")
  )

pB <- ggplot(panelB, aes(sdP, gs, colour = group, linetype = estimator, shape = estimator)) +
  geom_hline(yintercept = true_global_strength, colour = "grey50", linetype = "dotted") +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Low baseline" = col_green, "High baseline" = col_red)) +
  labs(x = expression(SD[P]~"(slow-context variation)"), y = "Estimated global strength",
       title = "(B) Estimated global strength by group", colour = NULL, linetype = NULL, shape = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.8))

panelC <- summary_tbl |>
  transmute(
    sdP,
    mae_raw = (mae_raw_lo_mean + mae_raw_hi_mean) / 2,
    mae_adj = (mae_adj_lo_mean + mae_adj_hi_mean) / 2,
    se_raw = sqrt(mae_raw_lo_se^2 + mae_raw_hi_se^2) / 2,
    se_adj = sqrt(mae_adj_lo_se^2 + mae_adj_hi_se^2) / 2
  ) |>
  pivot_longer(c(mae_raw, mae_adj), names_to = "estimator", values_to = "mae") |>
  mutate(
    se = if_else(estimator == "mae_raw", se_raw, se_adj),
    estimator = if_else(estimator == "mae_raw", "Symptom-only (raw)", "Context-adjusted")
  )

pC <- ggplot(panelC, aes(sdP, mae, colour = estimator, linetype = estimator)) +
  geom_errorbar(aes(ymin = mae - se, ymax = mae + se), width = 0.02) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only (raw)" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = expression(SD[P]~"(slow-context variation)"), y = expression("Mean "*group("|", hat(W) - W^{true}, "|")),
       title = expression("(C) Recovery error relative to "*W^{true}), colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.85))

library(patchwork)
fig7 <- (pA | pB | pC) +
  plot_annotation(title = "Figure 7. Network-estimation check under omitted slow context")

ggsave("img/Figure7_network_estimation_check.pdf", fig7, width = 13.5, height = 4.2)
ggsave("img/Figure7_network_estimation_check.png", fig7, width = 13.5, height = 4.2, dpi = 200)

# ------------------------------------------------------------------------
# 5. Figure 8: qgraph network comparison at the largest SD_P (where the
#    identification problem is most visible), averaged across replicates
#    to show the systematic pattern rather than a single noisy draw.
# ------------------------------------------------------------------------
sdP_show <- max(sdP_levels)

avg_edges <- function(sdP_val, estimator = c("raw", "adj"), group = c("lo", "hi")) {
  estimator <- match.arg(estimator); group <- match.arg(group)
  Ws <- vector("list", n_reps)
  for (r in seq_len(n_reps)) {
    seed <- 10000L + r + round(1000 * sdP_val)
    set.seed(seed)
    sd_use <- max(sdP_val, 1e-6)
    P_low  <- rnorm(n_per_group, mu_low,  sd_use)
    P_high <- rnorm(n_per_group, mu_high, sd_use)
    theta_low  <- outer(rep(1, n_per_group), h_vec) + outer(P_low,  gamma_vec)
    theta_high <- outer(rep(1, n_per_group), h_vec) + outer(P_high, gamma_vec)
    S_low  <- sample_states(theta_low)
    S_high <- sample_states(theta_high)
    S  <- if (group == "lo") S_low else S_high
    Pv <- if (group == "lo") P_low else P_high
    Ws[[r]] <- if (estimator == "raw") fit_edges(S) else fit_edges(S, Pv)
  }
  Reduce("+", Ws) / length(Ws)
}

W_raw_lo_avg <- avg_edges(sdP_show, "raw", "lo")
W_raw_hi_avg <- avg_edges(sdP_show, "raw", "hi")
W_adj_lo_avg <- avg_edges(sdP_show, "adj", "lo")
W_adj_hi_avg <- avg_edges(sdP_show, "adj", "hi")

# shared layout (from W_true) and shared edge-weight scale (from the raw
# networks, which have the largest range) so panels are visually comparable
Lmat <- qgraph(W_true, layout = "spring", DoNotPlot = TRUE)$layout
max_edge <- max(abs(c(W_true, W_raw_lo_avg, W_raw_hi_avg, W_adj_lo_avg, W_adj_hi_avg)))

pdf("img/Figure8_network_estimation_check_graphs.pdf", width = 15, height = 3.4)
layout_mat <- matrix(1:5, nrow = 1)
layout(layout_mat)
qgraph(W_true, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = expression(W^{true}), vsize = 10)
qgraph(W_raw_lo_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Raw, low baseline\n(SD[P]=%.2f)", sdP_show), vsize = 10)
qgraph(W_raw_hi_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Raw, high baseline\n(SD[P]=%.2f)", sdP_show), vsize = 10)
qgraph(W_adj_lo_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Adjusted, low baseline\n(SD[P]=%.2f)", sdP_show), vsize = 10)
qgraph(W_adj_hi_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Adjusted, high baseline\n(SD[P]=%.2f)", sdP_show), vsize = 10)
dev.off()

cat("Done.\n")
cat("  res/network_check/network_check_raw.rds\n")
cat("  res/network_check/network_check_summary.csv\n")
cat("  img/Figure7_network_estimation_check.pdf (+ .png)\n")
cat("  img/Figure8_network_estimation_check_graphs.pdf\n")
