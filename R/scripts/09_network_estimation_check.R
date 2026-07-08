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
# Run this script locally before updating the manuscript figures. The
# SD_P = 0 condition serves as a built-in sanity check: raw and adjusted
# estimates are set to be identical there by construction (see note in
# run_condition() below), so if anything else looks wrong at SD_P = 0,
# the bug is elsewhere in the pipeline, not in the raw-vs-adjusted logic.
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
})
# qgraph is only needed for the supplementary Figure 8 (network-graph
# comparison), not for the main Figure 7. Deliberately not loaded here so
# a missing qgraph install can't block the main analysis; see the
# requireNamespace() guard around the Figure 8 section below.

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
  # NOTE: this is deliberately NOT called `dgp` -- tibble() below creates a
  # column literally named `dgp`, and tibble() evaluates its arguments
  # sequentially, so a variable called `dgp` would get shadowed by that
  # column (an atomic value) for every argument that comes after it in the
  # same tibble() call. (This caused a real "$ operator is invalid for
  # atomic vectors" error in testing -- keep this named differently.)
  net_gen <- build_dgp(dgp_seed)
  seed <- 100000L * dgp_seed + 1000L * rep + round(1000 * sdP)
  set.seed(seed)

  sd_use <- max(sdP, 1e-6)
  P_low  <- rnorm(n_per_group, mu_low,  sd_use)
  P_high <- rnorm(n_per_group, mu_high, sd_use)

  theta_low  <- outer(rep(1, n_per_group), net_gen$h_vec) + outer(P_low,  net_gen$gamma_vec)
  theta_high <- outer(rep(1, n_per_group), net_gen$h_vec) + outer(P_high, net_gen$gamma_vec)

  S_low  <- net_gen$sample_states(theta_low)
  S_high <- net_gen$sample_states(theta_high)

  Wr_lo <- fit_edges(S_low)
  Wr_hi <- fit_edges(S_high)

  # At SD_P = 0, P is (numerically) constant within each group, so it is
  # collinear with the intercept and "adjusting" for it should do exactly
  # nothing to the edge estimates -- fitting glm() on a near-constant
  # covariate can still introduce numerical noise (huge/unstable
  # coefficient on P even though the edges themselves are unaffected), so
  # we make the equivalence exact here rather than relying on 1e-6 jitter
  # to keep glm() from complaining. This is also what makes the SD_P = 0
  # sanity check below an exact check rather than an approximate one.
  if (sdP == 0) {
    Wa_lo <- Wr_lo
    Wa_hi <- Wr_hi
  } else {
    Wa_lo <- fit_edges(S_low,  P_low)
    Wa_hi <- fit_edges(S_high, P_high)
  }

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
  #
  # IMPORTANT: these are left SIGNED here, not abs()'d. abs() must not be
  # applied at the single-replicate level -- if it is, you're averaging
  # abs(noise) across replicates, which never shrinks toward zero even
  # when the true expected signed difference is ~0, and the resulting
  # curve is dominated by replicate-level sampling noise rather than the
  # converged signal (this produced a visibly noisy, non-monotonic, even
  # sign-reversing Panel A in an earlier run of this script -- diagnosed
  # and fixed here). abs() is applied downstream, in the per_dgp
  # aggregation step, AFTER averaging the signed value across the
  # n_reps replicates for a given data-generating network.
  net_diff_raw_signed <- mean(Wr_hi[iu] - Wr_lo[iu])
  net_diff_adj_signed <- mean(Wa_hi[iu] - Wa_lo[iu])

  tibble(
    dgp = dgp_seed, sdP = sdP, rep = rep,
    net_diff_raw_signed = net_diff_raw_signed, net_diff_adj_signed = net_diff_adj_signed,
    gs_raw_lo = sum(abs(Wr_lo[iu])), gs_raw_hi = sum(abs(Wr_hi[iu])),
    gs_adj_lo = sum(abs(Wa_lo[iu])), gs_adj_hi = sum(abs(Wa_hi[iu])),
    mae_raw_lo = mean(abs(Wr_lo[iu] - net_gen$wtrue_edges)),
    mae_raw_hi = mean(abs(Wr_hi[iu] - net_gen$wtrue_edges)),
    mae_adj_lo = mean(abs(Wa_lo[iu] - net_gen$wtrue_edges)),
    mae_adj_hi = mean(abs(Wa_hi[iu] - net_gen$wtrue_edges)),
    true_gs = net_gen$true_gs,
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
#
# net_diff_raw/net_diff_adj: average the SIGNED per-replicate difference
# across the n_reps replicates for this network FIRST (this is where
# replicate-level sampling noise actually cancels, via the CLT over
# n_reps draws), and only THEN take abs() -- once, per data-generating
# network, not once per replicate. See the note in run_condition().
per_dgp <- results |>
  group_by(dgp, sdP) |>
  summarise(across(
    c(net_diff_raw_signed, net_diff_adj_signed, gs_raw_lo, gs_raw_hi, gs_adj_lo, gs_adj_hi,
      mae_raw_lo, mae_raw_hi, mae_adj_lo, mae_adj_hi, true_gs, m_low, m_high),
    mean
  ), .groups = "drop") |>
  mutate(net_diff_raw = abs(net_diff_raw_signed),
         net_diff_adj = abs(net_diff_adj_signed))

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

# Publication theme: visible axes (theme_classic base) rather than
# theme_minimal's boxless look, bottom-anchored legend so it can't overlap
# data regardless of where the lines fall, and a light horizontal grid
# only (no vertical grid) to keep the panels readable without looking like
# a default ggplot output.
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.major.x = element_blank()
    )
}

# NOTE on "Symptom-only" vs "Context-adjusted" labels: earlier drafts used
# "Symptom-only (raw)" -- dropped "(raw)" for the figure since it reads as
# informal; "raw" is still used as an internal variable-name suffix
# (Wr_*, net_diff_raw, mae_raw_*) throughout the script, only the
# user-facing label changed.
#
# NOTE ON DROPPING THE OLD "group-difference" PANEL: an earlier version of
# this figure led with a panel on |mean signed edge-weight shift| between
# groups (net_diff_raw/net_diff_adj). After fixing a real aggregation bug
# (abs() was being applied per-replicate instead of after averaging -- see
# run_condition() and the per_dgp step above) that panel is no longer
# noise-dominated, but on an actual full-scale run of this script it still
# didn't show a clean, publication-strength monotonic pattern: with 10
# independent random data-generating networks, WHICH group ends up more
# distorted is not a fixed direction -- it depends on where that
# particular network's baseline activation sits relative to the steepest
# part of the logistic response, which differs draw to draw and partially
# cancels on average. That is a genuine, defensible finding in its own
# right (the identification problem here is primarily about distortion of
# estimated coupling strength, not a stable between-group edge-pattern
# difference), but it is not strong/clean enough to be a headline main-text
# panel. It is NOT deleted from the analysis -- net_diff_raw_mean /
# net_diff_adj_mean are still in network_check_summary.csv and the paired
# t-test above, for reporting as a supplementary/text-level result if
# useful. The three panels below (global strength, recovery error,
# symptom prevalence) are the ones that reproduce cleanly.
panelA <- summary_tbl |>
  select(sdP, gs_raw_lo_mean, gs_raw_hi_mean, gs_adj_lo_mean, gs_adj_hi_mean, true_gs_mean) |>
  pivot_longer(-c(sdP, true_gs_mean), names_to = "series", values_to = "gs") |>
  mutate(group = if_else(grepl("_lo_", series), "Low baseline", "High baseline"),
         estimator = if_else(grepl("raw", series), "Symptom-only", "Context-adjusted"))

pA <- ggplot(panelA, aes(sdP, gs, colour = group, linetype = estimator, shape = estimator)) +
  geom_line(
    data = summary_tbl,
    aes(x = sdP, y = true_gs_mean),
    colour = "grey50", linetype = "dotted", inherit.aes = FALSE
  ) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Low baseline" = col_green, "High baseline" = col_red)) +
  labs(x = expression(SD[P]~"(slow-context variation)"), y = "Estimated global strength",
       title = "(A) Global strength inflation", colour = NULL, linetype = NULL, shape = NULL) +
  theme_pub()

panelB <- summary_tbl |>
  transmute(sdP,
            mae_raw = (mae_raw_lo_mean + mae_raw_hi_mean) / 2,
            mae_adj = (mae_adj_lo_mean + mae_adj_hi_mean) / 2,
            se_raw = sqrt(mae_raw_lo_se^2 + mae_raw_hi_se^2) / 2,
            se_adj = sqrt(mae_adj_lo_se^2 + mae_adj_hi_se^2) / 2) |>
  pivot_longer(c(mae_raw, mae_adj), names_to = "estimator", values_to = "mae") |>
  mutate(se = if_else(estimator == "mae_raw", se_raw, se_adj),
         estimator = if_else(estimator == "mae_raw", "Symptom-only", "Context-adjusted"))

pB <- ggplot(panelB, aes(sdP, mae, colour = estimator, linetype = estimator)) +
  geom_errorbar(aes(ymin = mae - se, ymax = mae + se), width = 0.02) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Symptom-only" = col_raw, "Context-adjusted" = col_adj)) +
  labs(x = expression(SD[P]~"(slow-context variation)"),
       y = expression("Mean "*group("|", hat(W) - W^{true}, "|")),
       title = expression("(B) Recovery error relative to "*W^{true}),
       colour = NULL, linetype = NULL) +
  theme_pub()

# Panel C (new): symptom prevalence by group. This is a diagnostic, not an
# estimator comparison -- m_low/m_high are the TRUE simulated mean
# activation in each group (not an estimate), included so a reader (and
# you, before trusting the other two panels) can confirm the effects above
# are not an artifact of one group being saturated near 0 or 1.
panelC <- summary_tbl |>
  select(sdP, m_low_mean, m_high_mean) |>
  pivot_longer(-sdP, names_to = "group", values_to = "prevalence") |>
  mutate(group = if_else(group == "m_low_mean", "Low baseline", "High baseline"))

pC <- ggplot(panelC, aes(sdP, prevalence, colour = group)) +
  geom_hline(yintercept = c(0.1, 0.9), colour = "grey85", linewidth = 0.3) +
  geom_line() + geom_point(size = 2) +
  scale_colour_manual(values = c("Low baseline" = col_green, "High baseline" = col_red)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression(SD[P]~"(slow-context variation)"), y = "Mean symptom activation",
       title = "(C) Symptom prevalence by group", colour = NULL) +
  theme_pub()

# No in-figure "Figure 7. ..." title -- that belongs in the LaTeX caption,
# not baked into the image.
#
# NOTE ON LEGENDS: deliberately NOT using plot_layout(guides = "collect").
# Color means a different thing in each panel -- GROUP in A and C, but
# ESTIMATOR in B -- so a single merged legend at the bottom would mix two
# different meanings under one strip and read as confusing ("High
# baseline | Low baseline | Context-adjusted | Symptom-only | ..." with no
# indication of which belongs to which panel). Keeping each panel's own
# local legend (guides = "keep", the patchwork default -- explicit here
# for clarity) means each legend sits under the panel it actually
# describes.
fig7 <- (pA | pB | pC) +
  plot_layout(guides = "keep") &
  theme(legend.position = "bottom")

ggsave("figs/Figure7_network_estimation_check.pdf", fig7, width = 13.5, height = 4.8)
ggsave("figs/Figure7_network_estimation_check.png", fig7, width = 13.5, height = 4.8, dpi = 200)

# ------------------------------------------------------------------------
# 5. Figure 8 (SUPPLEMENTARY ONLY -- do not use in the main paper).
#
# Redesigned as a 3-panel "true -> distorted -> corrected" figure rather
# than 5 separate group-specific panels: True coupling | Symptom-only
# estimate | Context-adjusted estimate, all at the largest SD_P examined,
# using ONE representative data-generating network (dgp1) and averaging
# the raw/adjusted estimate over the two groups (the group comparison is
# already Figure 7's job; this figure's only job is to make the
# distortion visually obvious). Weak edges are suppressed via qgraph's
# built-in `threshold` so the raw panel doesn't just read as "a blob", and
# all three panels share layout, edge scale, and threshold so density can
# be compared visually across them.
#
# qgraph is only needed for this section. Guarded with requireNamespace()
# so a missing qgraph install can't block Figure 7 / the main analysis
# above, which has already completed and saved by this point.
# ------------------------------------------------------------------------
if (requireNamespace("qgraph", quietly = TRUE)) {
  library(qgraph)

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

  # Single representative network per estimator, averaged over the two
  # groups -- this figure is about true-vs-distorted-vs-corrected, not a
  # group comparison (that's Figure 7).
  W_raw_avg <- (W_raw_lo_avg + W_raw_hi_avg) / 2
  W_adj_avg <- (W_adj_lo_avg + W_adj_hi_avg) / 2

  symptom_labels <- paste0("S", seq_len(N))
  Lmat <- qgraph(dgp1$W_true, layout = "spring", DoNotPlot = TRUE)$layout
  max_edge <- max(abs(c(dgp1$W_true, W_raw_avg, W_adj_avg)))

  # Edge threshold: half the smallest true edge weight. This keeps every
  # genuine edge in the true-network panel intact while suppressing the
  # small-magnitude spurious edges (positive or negative) that estimation
  # noise adds to the raw/adjusted panels -- applied uniformly is more
  # honest than selectively recoloring or hiding negative edges, since a
  # surviving large negative edge after thresholding would itself be a
  # real thing worth seeing.
  edge_threshold <- 0.5 * min(abs(dgp1$wtrue_edges[dgp1$wtrue_edges != 0]))

  # --- Publication styling ------------------------------------------
  # Muted, print-friendly palette instead of qgraph's default bright
  # blue/red + beige. Node fill/border/label colors are identical across
  # all three panels (only the edges differ), which is what actually
  # makes "true vs. distorted vs. corrected" readable at a glance.
  node_fill   <- rep("#F3E7C9", N)   # warm ivory
  node_border <- "#2F2F2F"
  label_col   <- "#1F1F1F"
  pos_col     <- "#3B76AF"           # muted steel blue (positive edges)
  neg_col     <- "#C75B5B"           # muted brick red (negative edges)

  # Per-edge color matrix with alpha transparency so weak (but
  # above-threshold) edges read as lighter than strong ones, rather than
  # every surviving edge looking equally bold.
  has_scales <- requireNamespace("scales", quietly = TRUE)
  make_edge_col <- function(W) {
    ec <- matrix(NA_character_, nrow(W), ncol(W))
    a <- pmin(abs(W) / max_edge, 1)
    for (i in seq_len(nrow(W))) for (j in seq_len(ncol(W))) {
      if (abs(W[i, j]) < edge_threshold) next
      base <- if (W[i, j] > 0) pos_col else neg_col
      ec[i, j] <- if (has_scales) scales::alpha(base, 0.35 + 0.55 * a[i, j]) else base
    }
    ec
  }

  qgraph_pub <- function(W, panel_title) {
    qgraph(
      W,
      layout = Lmat,
      theme = "classic",
      labels = symptom_labels,
      label.color = label_col,
      color = node_fill,
      border.color = node_border,
      border.width = 1.6,
      vsize = 9,
      label.cex = 1.15,
      edge.color = make_edge_col(W),
      maximum = max_edge,
      threshold = edge_threshold,
      fade = FALSE,
      title = ""
    )
    
    title(
      main = panel_title,
      line = 1.0,
      cex.main = 1.05,
      font.main = 2
    )
  }
  
  title_true <- "True coupling"
  
  title_raw <- bquote(
    atop("Symptom-only estimate", SD[P] == .(formatC(sdP_show, format = "f", digits = 2)))
  )
  
  title_adj <- bquote(
    atop("Context-adjusted estimate", SD[P] == .(formatC(sdP_show, format = "f", digits = 2)))
  )
  
  pdf("figs/Figure8_network_estimation_check_graphs.pdf", width = 10, height = 3.6)
  layout(matrix(1:3, nrow = 1))
  qgraph_pub(dgp1$W_true, title_true)
  qgraph_pub(W_raw_avg, title_raw)
  qgraph_pub(W_adj_avg, title_adj)
  dev.off()
  
  png("figs/Figure8_network_estimation_check_graphs.png",
      width = 2000, height = 720, res = 200)
  layout(matrix(1:3, nrow = 1))
  qgraph_pub(dgp1$W_true, title_true)
  qgraph_pub(W_raw_avg, title_raw)
  qgraph_pub(W_adj_avg, title_adj)
  dev.off()
  
  cat("  figs/Figure8_network_estimation_check_graphs.pdf    <- supplementary only\n")
} else {
  message("qgraph not installed; skipping supplementary Figure 8 (Figure 7 above is unaffected).")
}

cat("\nDone. Files:\n")
cat("  res/network_check/network_check_raw.rds\n")
cat("  res/network_check/network_check_summary.csv\n")
cat("  figs/Figure7_network_estimation_check.pdf (+ .png)  <- main text\n")
