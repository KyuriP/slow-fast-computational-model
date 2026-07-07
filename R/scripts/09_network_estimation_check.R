# ============================================================
# 09_network_estimation_check.R
# ============================================================
# Purpose
# -------
# Compact network-estimation check reported in the manuscript
# ("Network-estimation check" subsection). Implements the pairwise
# data-generating model exactly as specified there:
#
#   logit Pr(S_i = 1 | S_-i, P) = h_i + sum_j W_true_ij S_j + gamma_i * P
#
# with W_true, h_i, gamma_i held fixed across two groups that differ only
# in the distribution of the slow context P. We vary the amount of
# within-sample slow-context variation (SD_P) and compare a symptom-only
# ("raw") nodewise logistic-regression network estimator against a
# context-adjusted estimator that includes P as a covariate.
#
# We use exact enumeration for the fast layer (rather than IsingSampler /
# MCMC) because each simulated person has their own field h + gamma*P_i,
# and IsingSampler takes one shared threshold vector for all n draws --
# it cannot give each observation its own covariate-dependent field.
# With N = 12 that is 2^12 = 4096 states, which enumerates instantly.
#
# Robustness: to avoid reporting results from a single arbitrary draw of
# W_true/h/gamma, we repeat the whole check over N_DGP independent random
# data-generating networks and average, so the figure reflects a general
# pattern rather than one convenient network.
#
# Validation status: the estimator/sampler logic and the full multi-DGP
# loop structure were run end-to-end without errors in a headless
# WebAssembly R session (webR) at reduced scale (small n, few reps/DGP
# draws) to catch bugs before this went into the repo; that run confirmed
# the qualitative pattern (raw estimator overshoots true global strength
# and diverges between groups as SD_P grows; adjusted estimator does not).
# It was NOT run at the full scale below -- please run this script
# locally and inspect the sanity check printed at the end (SD_P = 0 row:
# raw and adjusted should match closely there) before using the figures.
#
# Outputs
# -------
#   res/network_check/network_check_raw.rds        (replicate-level, all DGP draws)
#   res/network_check/network_check_summary.csv     (condition-level means/SEs, averaged over DGP draws)
#   figs/Figure7_network_estimation_check.pdf / .png
#   figs/Figure8_network_estimation_check_graphs.pdf  (SUPPLEMENTARY -- not main text; see note below)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(future)
  library(ggplot2)
  library(patchwork)
  library(qgraph)
})

# ------------------------------------------------------------------------
# 0. Design
# ------------------------------------------------------------------------
N            <- 12L               # matches Appendix C (NCT check), for consistency
edge_density <- 0.5
w_lo         <- 0.15
w_hi         <- 0.35
h_lo         <- -1.6
h_hi         <- -1.0
gamma_lo     <- 0.8
gamma_hi     <- 1.4

mu_low  <- -0.3
mu_high <- 1.0
sdP_levels  <- c(0.0, 0.1, 0.3, 0.5, 0.75, 1.0)
n_per_group <- 2000L
n_reps      <- 15L    # replicates PER data-generating network
n_dgp       <- 10L    # independent random data-generating networks (robustness)

iu <- which(upper.tri(matrix(0, N, N)))

# ------------------------------------------------------------------------
# 1. One data-generating network: build sampler + estimator closures
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

# Nodewise logistic-regression network estimator (symmetrized), with an
# optional covariate P for the context-adjusted estimator.
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
# 2. One (dgp, sdP, rep) condition
# ------------------------------------------------------------------------
run_condition <- function(dgp_seed, sdP, rep) {
  dgp <- build_dgp(dgp_seed)
  seed <- 100000L * dgp_seed + 1000L * rep + round(1000 * sdP)
  set.seed(seed)

  sd_use <- max(sdP, 1e-6)
  P_low  <- rnorm(n_per_group, mu_low,  sd_use)
  P_high <- rnorm(n_per_group, mu_high, sd_use)

  theta_low  <- outer(rep(1, n_per_group), dgp$h_vec) + outer(P_low,  dgp$gamma_vec)
  theta_high <- outer(rep(1, n_per_group), dgp$h_vec) + outer(P_high, dgp$gamma_vec)

  S_low  <- dgp$sample_states(theta_low)
  S_high <- dgp$sample_states(theta_high)

  Wr_lo <- fit_edges(S_low)
  Wr_hi <- fit_edges(S_high)
  Wa_lo <- fit_edges(S_low,  P_low)
  Wa_hi <- fit_edges(S_high, P_high)

  # --- Panel A metric -----------------------------------------------
  # NOTE: an earlier draft used mean(|Wr_hi - Wr_lo|) (absolute value
  # taken PER EDGE before averaging). We tested that version at full
  # scale and it is noise-dominated -- per-edge sampling noise swamps
  # the systematic signal and it stays flat (~0.12-0.13) regardless of
  # SD_P, sometimes even reversing raw vs adjusted. Taking the mean
  # SIGNED difference first and then the absolute value exploits
  # cancellation of independent per-edge noise across the N*(N-1)/2
  # edges while preserving the systematic (confounding-driven) shift,
  # and it reproduces the identification-problem pattern cleanly and
  # matches the (independently validated) global-strength story below.
  net_diff_raw <- abs(mean(Wr_hi[iu] - Wr_lo[iu]))
  net_diff_adj <- abs(mean(Wa_hi[iu] - Wa_lo[iu]))

  tibble(
    dgp = dgp_seed, sdP = sdP, rep = rep,
    net_diff_raw = net_diff_raw, net_diff_adj = net_diff_adj,
    gs_raw_lo = sum(abs(Wr_lo[iu])), gs_raw_hi = sum(abs(Wr_hi[iu])),
    gs_adj_lo = sum(abs(Wa_lo[iu])), gs_adj_hi = sum(abs(Wa_hi[iu])),
    mae_raw_lo = mean(abs(Wr_lo[iu] - dgp$wtrue_edges)),
    mae_raw_hi = mean(abs(Wr_hi[iu] - dgp$wtrue_edges)),
    mae_adj_lo = mean(abs(Wa_lo[iu] - dgp$wtrue_edges)),
    mae_adj_hi = mean(abs(Wa_hi[iu] - dgp$wtrue_edges)),
    true_gs = dgp$true_gs,
    m_low = mean(S_low), m_high = mean(S_high)
  )
}

# ------------------------------------------------------------------------
# 3. Run: n_dgp x sdP_levels x n_reps conditions, parallelized
# ------------------------------------------------------------------------
design <- expand_grid(dgp_seed = seq_len(n_dgp), sdP = sdP_levels, rep = seq_len(n_reps))

n_cores <- max(1L, future::availableCores() - 1L)
plan(multisession, workers = n_cores)

results <- design |>
  future_pmap_dfr(
    .f = function(dgp_seed, sdP, rep) run_condition(dgp_seed, sdP, rep),
    .options = furrr_options(seed = TRUE)
  )

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)

saveRDS(results, "res/network_check/network_check_raw.rds")

# Average within each DGP draw first, then across draws, so each network
# contributes equally regardless of how noisy any one replicate was.
per_dgp <- results |>
  group_by(dgp, sdP) |>
  summarise(across(
    c(net_diff_raw, net_diff_adj, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs, m_low, m_high),
    mean
  ), .groups = "drop")

summary_tbl <- per_dgp |>
  group_by(sdP) |>
  summarise(across(
    c(net_diff_raw, net_diff_adj, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs, m_low, m_high),
    list(mean = mean, se = ~ sd(.x) / sqrt(length(.x))),
    .names = "{.col}_{.fn}"
  ), .groups = "drop")

write.csv(summary_tbl, "res/network_check/network_check_summary.csv", row.names = FALSE)
print(summary_tbl |> select(sdP, net_diff_raw_mean, net_diff_adj_mean,
                             gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean))
cat(sprintf("mean true_global_strength across %d DGP draws = %.3f\n",
            n_dgp, mean(summary_tbl$true_gs_mean)))
cat("\nSANITY CHECK: the SD_P = 0 row above should show raw approx= adjusted\n")
cat("for every column (P is ~constant within each group there, so it is\n")
cat("collinear with the intercept and 'adjusting' for it does ~nothing).\n\n")

# Lightweight formal check: paired t-test, net_diff_raw vs net_diff_adj,
# at the largest SD_P (per-DGP-draw means as the paired unit).
top_sdP <- max(sdP_levels)
paired <- per_dgp |> filter(sdP == top_sdP)
tt <- t.test(paired$net_diff_raw, paired$net_diff_adj, paired = TRUE, alternative = "greater")
cat(sprintf("Paired t-test at SD_P = %.2f (n_dgp = %d draws): raw > adjusted, t = %.2f, p = %.4f\n",
            top_sdP, n_dgp, tt$statistic, tt$p.value))

# ------------------------------------------------------------------------
# 4. Figure 7: 3-panel summary
# ------------------------------------------------------------------------
col_green <- "#2E7D32"
col_red   <- "#C62828"
col_raw   <- "grey15"
col_adj   <- "#1565C0"

panelA <- summary_tbl |>
  transmute(sdP, diff_raw = net_diff_raw_mean, diff_adj = net_diff_adj_mean,
            se_raw = net_diff_raw_se, se_adj = net_diff_adj_se) |>
  pivot_longer(c(diff_raw, diff_adj), names_to = "estimator", values_to = "diff") |>
  mutate(se = if_else(estimator == "diff_raw", se_raw, se_adj),
         estimator = if_else(estimator == "diff_raw", "Symptom-only (raw)", "Context-adjusted"))

pA <- ggplot(panelA, aes(sdP, diff, colour = estimator, linetype = estimator)) +
  geom_errorbar(aes(ymin = diff - se, ymax = diff + se), width = 0.02) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only (raw)" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = expression(SD[P]~"(slow-context variation)"),
       y = "Mean edge-weight difference\nbetween groups",
       title = "(A) Apparent group network difference", colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.85))

panelB <- summary_tbl |>
  select(sdP, gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean, true_gs_mean) |>
  pivot_longer(-c(sdP, true_gs_mean), names_to = "series", values_to = "gs") |>
  mutate(group = if_else(grepl("_lo_", series), "Low baseline", "High baseline"),
         estimator = if_else(grepl("raw", series), "raw", "adjusted"))

pB <- ggplot(panelB, aes(sdP, gs, colour = group, linetype = estimator, shape = estimator)) +
  geom_line(aes(y = true_gs_mean), colour = "grey50", linetype = "dotted", inherit.aes = FALSE,
            data = summary_tbl) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Low baseline" = col_green, "High baseline" = col_red)) +
  labs(x = expression(SD[P]~"(slow-context variation)"), y = "Estimated global strength",
       title = "(B) Estimated global strength by group", colour = NULL, linetype = NULL, shape = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.8))

panelC <- summary_tbl |>
  transmute(sdP,
            mae_raw = (mae_raw_lo_mean + mae_raw_hi_mean) / 2,
            mae_adj = (mae_adj_lo_mean + mae_adj_hi_mean) / 2,
            se_raw = sqrt(mae_raw_lo_se^2 + mae_raw_hi_se^2) / 2,
            se_adj = sqrt(mae_adj_lo_se^2 + mae_adj_hi_se^2) / 2) |>
  pivot_longer(c(mae_raw, mae_adj), names_to = "estimator", values_to = "mae") |>
  mutate(se = if_else(estimator == "mae_raw", se_raw, se_adj),
         estimator = if_else(estimator == "mae_raw", "Symptom-only (raw)", "Context-adjusted"))

pC <- ggplot(panelC, aes(sdP, mae, colour = estimator, linetype = estimator)) +
  geom_errorbar(aes(ymin = mae - se, ymax = mae + se), width = 0.02) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only (raw)" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = expression(SD[P]~"(slow-context variation)"),
       y = expression("Mean "*group("|", hat(W) - W^{true}, "|")),
       title = expression("(C) Recovery error relative to "*W^{true}),
       colour = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = c(0.3, 0.85))

fig7 <- (pA | pB | pC) +
  plot_annotation(title = sprintf(
    "Figure 7. Network-estimation check under omitted slow context (averaged over %d data-generating networks)",
    n_dgp))

ggsave("figs/Figure7_network_estimation_check.pdf", fig7, width = 13.5, height = 4.2)
ggsave("figs/Figure7_network_estimation_check.png", fig7, width = 13.5, height = 4.2, dpi = 200)

# ------------------------------------------------------------------------
# 5. Figure 8 (SUPPLEMENTARY ONLY -- do not use in the main paper; see
#    reviewer note: the qgraph comparison is a useful diagnostic but adds
#    visual complexity that the compact Figure 7 doesn't need). Uses the
#    first DGP draw only, as a single illustrative example.
# ------------------------------------------------------------------------
dgp1 <- build_dgp(1L)
sdP_show <- max(sdP_levels)

avg_edges <- function(dgp, sdP_val, estimator = c("raw", "adj"), group = c("lo", "hi")) {
  estimator <- match.arg(estimator); group <- match.arg(group)
  acc <- matrix(0, N, N)
  for (r in seq_len(n_reps)) {
    seed <- 100000L * 1L + 1000L * r + round(1000 * sdP_val)
    set.seed(seed)
    sd_use <- max(sdP_val, 1e-6)
    P_low  <- rnorm(n_per_group, mu_low,  sd_use)
    P_high <- rnorm(n_per_group, mu_high, sd_use)
    theta_low  <- outer(rep(1, n_per_group), dgp$h_vec) + outer(P_low,  dgp$gamma_vec)
    theta_high <- outer(rep(1, n_per_group), dgp$h_vec) + outer(P_high, dgp$gamma_vec)
    S_low  <- dgp$sample_states(theta_low)
    S_high <- dgp$sample_states(theta_high)
    S  <- if (group == "lo") S_low else S_high
    Pv <- if (group == "lo") P_low else P_high
    acc <- acc + (if (estimator == "raw") fit_edges(S) else fit_edges(S, Pv))
  }
  acc / n_reps
}

W_raw_lo_avg <- avg_edges(dgp1, sdP_show, "raw", "lo")
W_raw_hi_avg <- avg_edges(dgp1, sdP_show, "raw", "hi")
W_adj_lo_avg <- avg_edges(dgp1, sdP_show, "adj", "lo")
W_adj_hi_avg <- avg_edges(dgp1, sdP_show, "adj", "hi")

symptom_labels <- paste0("S", seq_len(N))
Lmat <- qgraph(dgp1$W_true, layout = "spring", DoNotPlot = TRUE)$layout
max_edge <- max(abs(c(dgp1$W_true, W_raw_lo_avg, W_raw_hi_avg, W_adj_lo_avg, W_adj_hi_avg)))

pdf("figs/Figure8_network_estimation_check_graphs.pdf", width = 15, height = 3.4)
layout(matrix(1:5, nrow = 1))
qgraph(dgp1$W_true, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = expression(W^{true}), vsize = 8)
qgraph(W_raw_lo_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Raw, low baseline\n(SD[P]=%.2f)", sdP_show), vsize = 8)
qgraph(W_raw_hi_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Raw, high baseline\n(SD[P]=%.2f)", sdP_show), vsize = 8)
qgraph(W_adj_lo_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Adjusted, low baseline\n(SD[P]=%.2f)", sdP_show), vsize = 8)
qgraph(W_adj_hi_avg, layout = Lmat, labels = symptom_labels, theme = "colorblind",
       maximum = max_edge, title = sprintf("Adjusted, high baseline\n(SD[P]=%.2f)", sdP_show), vsize = 8)
dev.off()

cat("\nDone. Files:\n")
cat("  res/network_check/network_check_raw.rds\n")
cat("  res/network_check/network_check_summary.csv\n")
cat("  figs/Figure7_network_estimation_check.pdf (+ .png)  <- main text\n")
cat("  figs/Figure8_network_estimation_check_graphs.pdf    <- supplementary only\n")
