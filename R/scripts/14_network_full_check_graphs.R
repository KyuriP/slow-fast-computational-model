# ============================================================
# 14_network_full_check_graphs.R
# ============================================================
# Purpose
# -------
# Regenerates Figure 8 (the true/symptom-only/context-adjusted network
# diagrams) to match the new main-text network-estimation check (script
# 13: real coupled feedback model, sigma_P x window design), replacing
# the old panel generated under the abandoned exogenous SD_P=1.0 design.
#
# We pick ONE representative condition -- the most dramatic in the new
# design, sigma_P=0.65 (the largest diffusion level swept) at W=1 (a
# near-instantaneous, cross-sectional snapshot) -- and ONE representative
# data-generating network (dgp_seed=1), high-baseline group. This
# mirrors the role the old SD_P=1.0 condition played: the setting where
# the omitted-context distortion is largest and therefore most visible
# in a static network diagram.
#
# Outputs
# -------
#   res/network_check/network_full_check_graph_mats.rds  (W_true, W_raw, W_adj)
#   figs/Figure8_network_estimation_check_graphs.pdf / .png  (replaces old Fig 8)
# ============================================================

suppressPackageStartupMessages({
  library(qgraph)
  library(parallel)
})
if (!requireNamespace("matrixStats", quietly = TRUE)) install.packages("matrixStats")
suppressPackageStartupMessages(library(matrixStats))

# ------------------------------------------------------------------------
# 0. Same design constants and functions as script 13 (duplicated so this
#    script stays standalone-runnable)
# ------------------------------------------------------------------------
N            <- 12L
edge_density <- 0.5
w_lo         <- 0.15
w_hi         <- 0.35
h_lo         <- -1.6
h_hi         <- -1.0
gamma_lo     <- 0.8
gamma_hi     <- 1.4

P_base_high <- 1.0

kappa    <- 0.20
dt       <- 0.02
lambda_m <- 0.001
m_star   <- 0.25
b_real   <- 0.15

burn_in_steps <- 3000L
n_per_group   <- 1000L

# The representative condition: largest sigma_P, shortest window --
# where the omitted-context distortion is most visible, mirroring the
# role SD_P=1.0 played in the old design.
sigma_P_rep <- 0.65
W_rep       <- 1L
dgp_rep     <- 1L

iu <- which(upper.tri(matrix(0, N, N)))

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

  states <- as.matrix(expand.grid(rep(list(c(0, 1)), N)))
  storage.mode(states) <- "double"
  quad_base <- 0.5 * rowSums((states %*% W_true) * states)

  # NOTE on speed: sample_states() is called once per timestep, and for a
  # single replicate here that means ~3001 calls (burn_in + W). Each call
  # does elementwise arithmetic on a 4096 (2^N states) x n_ind matrix.
  # sweep() is convenient but has real per-call overhead (S3 dispatch +
  # a fresh allocation every time) that adds up over tens of thousands of
  # calls. Replaced with direct vector-recycling arithmetic below, which
  # is numerically identical to the sweep() version (same floating-point
  # result) but avoids that overhead -- this is a pure speed fix, not a
  # change to the model or its randomness (rnorm/runif calls untouched).
  sample_states <- function(theta_matrix) {
    n_person <- nrow(theta_matrix)
    n_states <- nrow(states)
    linear <- states %*% t(theta_matrix)
    logp <- linear + quad_base                                  # recycles by row (length == nrow)
    logp <- logp - rep(matrixStats::colMaxs(logp), each = n_states)
    p <- exp(logp)
    p <- p / rep(colSums(p), each = n_states)
    cdf <- matrixStats::colCumsums(p)
    u <- runif(n_person)
    below <- cdf < rep(u, each = n_states)
    idx <- colSums(below) + 1L
    states[idx, , drop = FALSE]
  }

  list(W_true = W_true, h_vec = h_vec, gamma_vec = gamma_vec, sample_states = sample_states)
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

simulate_coupled_window <- function(net_gen, P_base, n_ind, b, sigma_P, burn_in, W_steps) {
  P <- rep(P_base, n_ind)
  m_slow <- rep(m_star, n_ind)
  total_steps <- burn_in + W_steps
  running_sum <- numeric(n_ind)
  for (t in seq_len(total_steps)) {
    theta <- outer(rep(1, n_ind), net_gen$h_vec) + outer(P, net_gen$gamma_vec)
    S <- net_gen$sample_states(theta)
    m_t <- rowMeans(S)
    m_slow <- (1 - lambda_m) * m_slow + lambda_m * m_t
    P <- P + (-kappa * (P - P_base) + b * (m_slow - m_star)) * dt +
      sigma_P * sqrt(dt) * rnorm(n_ind)
    if (t > burn_in) running_sum <- running_sum + P
  }
  Pbar_W <- running_sum / W_steps
  theta_obs <- outer(rep(1, n_ind), net_gen$h_vec) + outer(Pbar_W, net_gen$gamma_vec)
  S_obs <- net_gen$sample_states(theta_obs)
  list(S_obs = S_obs, Pbar_W = Pbar_W)
}

# ------------------------------------------------------------------------
# 1. Simulate the representative condition and fit both estimators,
#    averaged across n_reps_illustration independent replicates.
#
#    The quantitative results (Figure 7) are averaged across many
#    replicates specifically to wash out single-draw sampling noise in
#    individual edges. The first version of this script plotted only ONE
#    replicate's raw estimate, so it was fully exposed to that noise --
#    e.g. an edge that is genuinely near-zero in expectation can still
#    show up as a sizable, wrong-signed edge in any one draw. Averaging
#    the estimated matrices over replicates here (same condition, same
#    data-generating network, independent draws) puts Figure 8 on the
#    same footing as Figure 7: noise-driven edges shrink toward zero on
#    average, while any edge that survives averaging reflects a
#    systematic pattern rather than single-draw noise.
# ------------------------------------------------------------------------
n_reps_illustration <- 20L

# The 20 replicates are fully independent (different RNG seed, same dgp/
# condition) -- there is no reason to run them one after another. Each
# replicate is the expensive part (~3001-timestep simulation), while
# combining the 20 fitted matrices at the end is trivial. We therefore
# fan the replicate loop out across cores: on Mac/Linux this uses
# mclapply (fork-based, no data copying needed since child processes
# inherit the parent's memory); on Windows (no fork) we fall back to a
# PSOCK cluster with parLapply, exporting only what each worker needs.
n_cores <- max(1L, parallel::detectCores() - 1L)
cat(sprintf("Simulating representative condition: dgp=%d, sigma_P=%.2f, W=%d, high-baseline group, %d replicates on %d cores...\n",
            dgp_rep, sigma_P_rep, W_rep, n_reps_illustration, n_cores))

net_gen <- build_dgp(dgp_rep)
W_true <- net_gen$W_true

run_one_rep <- function(r) {
  set.seed(900000L * dgp_rep + 1000L * r + 100L * W_rep + round(sigma_P_rep * 1000))
  res <- simulate_coupled_window(net_gen, P_base_high, n_per_group, b_real, sigma_P_rep, burn_in_steps, W_rep)
  list(W_raw = fit_edges(res$S_obs), W_adj = fit_edges(res$S_obs, res$Pbar_W))
}

t0 <- Sys.time()
if (.Platform$OS.type == "unix") {
  rep_results <- parallel::mclapply(seq_len(n_reps_illustration), run_one_rep, mc.cores = n_cores)
} else {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("net_gen", "P_base_high", "n_per_group", "b_real", "sigma_P_rep",
                                 "burn_in_steps", "W_rep", "dgp_rep", "simulate_coupled_window",
                                 "fit_edges", "N", "lambda_m", "m_star", "kappa", "dt"),
                           envir = environment())
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(matrixStats)))
  rep_results <- parallel::parLapply(cl, seq_len(n_reps_illustration), run_one_rep)
}
cat(sprintf("Replicates done in %.1f sec.\n", as.numeric(Sys.time() - t0, units = "secs")))

# Guard against a failed replicate silently corrupting the average (e.g. a
# worker error would otherwise show up as a non-list/NULL entry).
ok <- vapply(rep_results, function(x) is.list(x) && !is.null(x$W_raw), logical(1))
if (!all(ok)) stop(sprintf("%d of %d replicates failed -- check errors above.", sum(!ok), n_reps_illustration))

W_raw_acc <- Reduce(`+`, lapply(rep_results, `[[`, "W_raw"))
W_adj_acc <- Reduce(`+`, lapply(rep_results, `[[`, "W_adj"))

W_raw <- W_raw_acc / n_reps_illustration
W_adj <- W_adj_acc / n_reps_illustration

dir.create("res/network_check", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(W_true = W_true, W_raw = W_raw, W_adj = W_adj,
             sigma_P = sigma_P_rep, W_steps = W_rep, dgp_seed = dgp_rep,
             n_reps_illustration = n_reps_illustration),
        "res/network_check/network_full_check_graph_mats.rds")

cat(sprintf("True global strength: %.2f | Symptom-only: %.2f | Context-adjusted: %.2f\n",
            sum(abs(W_true[iu])), sum(abs(W_raw[iu])), sum(abs(W_adj[iu]))))

# ------------------------------------------------------------------------
# 2. Three-panel qgraph figure: True / Symptom-only / Context-adjusted,
#    same spring layout (from the true network) and edge-scaling across
#    all three panels, matching the established style from Figure 8.
# ------------------------------------------------------------------------
# Plotting threshold: hide edges below this magnitude so the panels aren't
# cluttered with near-zero sampling noise (matches the old Figure 8's
# convention -- "edges below a common plotting threshold are omitted for
# clarity"). Applied by zeroing matrix entries directly, rather than
# relying on qgraph's own minimum/threshold argument (behavior of that
# argument varies across qgraph versions and didn't visibly change
# anything when tried).
#
# Applied ONLY to the two ESTIMATED matrices (W_raw, W_adj), not to
# W_true. W_true has no noise to threshold away -- every true edge is
# generated to be >=0.15 by construction (Uniform(0.15, 0.35)), and every
# non-edge is exactly 0. Thresholding it too only cuts real edges sitting
# near the 0.15 floor, artificially sparsifying the one panel that should
# be shown exactly as generated -- which unfairly narrows the visual gap
# to the other two panels instead of clarifying it (this is what
# happened in the previous version of this script).
plot_min <- 0.15
zero_small <- function(W) { W[abs(W) < plot_min] <- 0; W }

W_true_plot <- W_true
W_raw_plot  <- zero_small(W_raw)
W_adj_plot  <- zero_small(W_adj)

mats <- list(True = W_true_plot, `Symptom-only` = W_raw_plot, `Context-adjusted` = W_adj_plot)
max_edge <- max(sapply(mats, function(m) max(abs(m))))

L <- qgraph(W_true_plot, layout = "spring", DoNotPlot = TRUE)$layout

plot_one <- function(W, title) {
  qgraph(W, layout = L, maximum = max_edge, fade = TRUE,
         labels = paste0("S", seq_len(N)),
         posCol = "#1565C0", negCol = "#C62828", edge.width = 1.3,
         vsize = 7, label.cex = 1.1, title = title, title.cex = 1.3,
         theme = "classic", DoNotPlot = FALSE)
}

dir.create("figs", showWarnings = FALSE)

pdf("figs/Figure8_network_estimation_check_graphs.pdf", width = 13, height = 4.6)
layout(t(1:3))
plot_one(W_true_plot, "True")
plot_one(W_raw_plot,  "Symptom-only")
plot_one(W_adj_plot,  "Context-adjusted")
dev.off()

png("figs/Figure8_network_estimation_check_graphs.png", width = 13, height = 4.6, units = "in", res = 200)
layout(t(1:3))
plot_one(W_true_plot, "True")
plot_one(W_raw_plot,  "Symptom-only")
plot_one(W_adj_plot,  "Context-adjusted")
dev.off()

cat("\nDone. Files:\n")
cat("  res/network_check/network_full_check_graph_mats.rds\n")
cat("  figs/Figure8_network_estimation_check_graphs.pdf (+ .png)\n")
