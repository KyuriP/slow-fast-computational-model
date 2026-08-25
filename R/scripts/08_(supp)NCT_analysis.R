# =============================================================================
# Monte Carlo NCT simulation (corrected for binary Ising data)
# Do threshold shifts alone create apparent network differences?
#
# Key corrections:
#  - binary.data = TRUE in NCT()
#  - weighted_mode switch:
#       TRUE  = weighted edge comparison (default; closest to your current plots)
#       FALSE = dichotomized network comparison (closer to adjacency/presence-absence)
#  - explicit AND = TRUE for binary data
#  - explicit p.adjust.methods = "none" for raw edge-level false-positive rates
# =============================================================================

suppressPackageStartupMessages({
  library(IsingSampler)
  library(IsingFit)
  library(NetworkComparisonTest)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
  library(progressr)
  library(parallelly)
})

# -----------------------------------------------------------------------------
# 0. User settings
# -----------------------------------------------------------------------------

# Main analysis:
# TRUE  = weighted networks (NCT structure = max edge-weight difference)
# FALSE = dichotomized networks (closer to adjacency-style comparison)
weighted_mode <- TRUE

# Monte Carlo settings
n_perm_main <- 250      # for quick runs
# For final results, consider 1000+ permutations
# n_perm_main <- 1000

# Ising / design settings
N_nodes        <- 12
w_within_main  <- 0.25
w_between_main <- 0.05
thresh_low_val <- -2.0
gamma_main     <- 0.25


analysis_tag <- if (weighted_mode) "weighted" else "dichotomized"

# -----------------------------------------------------------------------------
# 1. Parallelization setup
# -----------------------------------------------------------------------------
n_cores <- max(1L, availableCores() - 1L)
plan(multisession, workers = n_cores)
cat(sprintf("Using %d workers\n", n_cores))
cat(sprintf("Analysis mode: %s\n", analysis_tag))

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------
make_network <- function(type = c("dense", "modular", "ring"),
                         N = 12,
                         w_within = 0.25,
                         w_between = 0.05) {
  type <- match.arg(type)
  J <- matrix(0, N, N)
  
  if (type == "dense") {
    J[] <- w_within
    diag(J) <- 0
  } else if (type == "modular") {
    stopifnot(N %% 2 == 0)
    g1 <- 1:(N / 2)
    g2 <- (N / 2 + 1):N
    J[g1, g1] <- w_within
    J[g2, g2] <- w_within
    J[g1, g2] <- w_between
    J[g2, g1] <- w_between
    diag(J) <- 0
  } else if (type == "ring") {
    for (i in 1:N) {
      j <- if (i == N) 1L else i + 1L
      J[i, j] <- w_within
      J[j, i] <- w_within
    }
  }
  
  J
}

get_edge_pvals <- function(nct_obj) {
  p <- nct_obj$einv.pvals
  
  if (is.null(p)) {
    return(numeric(0))
  }
  
  if (is.vector(p)) {
    return(as.numeric(p))
  }
  
  cn <- colnames(p)
  col_idx <- if (!is.null(cn)) {
    grep("^p$|pval|p\\.value", cn, ignore.case = TRUE)
  } else {
    integer(0)
  }
  
  as.numeric(p[, if (length(col_idx)) col_idx[1] else min(3L, ncol(p))])
}

safe_mean_upper <- function(W) {
  vals <- W[upper.tri(W)]
  if (!length(vals)) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 3. One Monte Carlo replicate
# -----------------------------------------------------------------------------
run_rep <- function(topology,
                    n_per_grp,
                    threshold_gap,
                    rep,
                    thresh_low_val = -2.0,
                    N = 12,
                    w_within = 0.25,
                    w_between = 0.05,
                    n_perm = 250,
                    gamma = 0.25,
                    weighted_mode = TRUE) {
  
  seed <- 100000L +
    rep +
    1000L * match(topology, c("dense", "modular", "ring")) +
    10L   * n_per_grp +
    round(100 * threshold_gap)
  
  set.seed(seed)
  
  J_true <- make_network(
    type      = topology,
    N         = N,
    w_within  = w_within,
    w_between = w_between
  )
  
  thr_lo <- rep(thresh_low_val, N)
  thr_hi <- rep(thresh_low_val + threshold_gap, N)
  
  # Simulate two groups with identical true coupling, different thresholds
  dat_lo <- IsingSampler(
    n         = n_per_grp,
    graph     = J_true,
    thresholds = thr_lo,
    responses = c(0L, 1L),
    method    = "MH"
  )
  
  dat_hi <- IsingSampler(
    n         = n_per_grp,
    graph     = J_true,
    thresholds = thr_hi,
    responses = c(0L, 1L),
    method    = "MH"
  )
  
  # Separate estimation summaries (not used by NCT directly, but useful diagnostics)
  fit_lo <- IsingFit(
    dat_lo,
    family = "binomial",
    plot = FALSE,
    progressbar = FALSE
  )
  
  fit_hi <- IsingFit(
    dat_hi,
    family = "binomial",
    plot = FALSE,
    progressbar = FALSE
  )
  
  W_lo <- fit_lo$weiadj
  W_hi <- fit_hi$weiadj
  
  # IMPORTANT: binary.data = TRUE because dat_lo/dat_hi are binary Ising samples
  nct <- NCT(
    data1 = dat_lo,
    data2 = dat_hi,
    gamma = gamma,
    it = n_perm,
    binary.data = TRUE,
    paired = FALSE,
    weighted = weighted_mode,
    AND = TRUE,
    abs = TRUE,
    test.edges = TRUE,
    edges = "all",
    progressbar = FALSE,
    make.positive.definite = TRUE,
    p.adjust.methods = "none",
    test.centrality = FALSE,
    verbose = FALSE
  )
  
  ep <- get_edge_pvals(nct)
  ne <- length(ep)
  
  tibble(
    topology        = topology,
    n_per_grp       = n_per_grp,
    threshold_gap   = threshold_gap,
    rep             = rep,
    analysis_mode   = if (weighted_mode) "weighted" else "dichotomized",
    
    act_low         = mean(colMeans(dat_lo)),
    act_high        = mean(colMeans(dat_hi)),
    
    edge_low        = safe_mean_upper(W_lo),
    edge_high       = safe_mean_upper(W_hi),
    
    thresh_est_low  = mean(fit_lo$thresholds, na.rm = TRUE),
    thresh_est_high = mean(fit_hi$thresholds, na.rm = TRUE),
    
    nwinv_stat      = nct$nwinv.real,
    nwinv_p         = nct$nwinv.pval,
    
    glstr_stat      = nct$glstrinv.real,
    glstr_p         = nct$glstrinv.pval,
    
    edge_fp_05      = if (ne > 0) mean(ep < 0.05, na.rm = TRUE) else NA_real_,
    edge_fp_bonf    = if (ne > 0) mean(ep < (0.05 / ne), na.rm = TRUE) else NA_real_,
    n_edges         = ne
  )
}

# -----------------------------------------------------------------------------
# 4. Design grid
# -----------------------------------------------------------------------------
design <- expand_grid(
  topology      = c("dense", "modular", "ring"),
  n_per_grp     = c(200L, 500L, 1000L),
  threshold_gap = c(0.4, 0.8, 1.2),
  rep           = 1:100
)

cat(sprintf("Total replicates: %d\n", nrow(design)))

# -----------------------------------------------------------------------------
# 5. Run Monte Carlo
# -----------------------------------------------------------------------------
handlers(global = TRUE)

with_progress({
  p <- progressor(nrow(design))
  
  results <- design |>
    future_pmap_dfr(
      .f = function(topology, n_per_grp, threshold_gap, rep) {
        p()
        
        run_rep(
          topology      = topology,
          n_per_grp     = n_per_grp,
          threshold_gap = threshold_gap,
          rep           = rep,
          thresh_low_val = thresh_low_val,
          N             = N_nodes,
          w_within      = w_within_main,
          w_between     = w_between_main,
          n_perm        = n_perm_main,
          gamma         = gamma_main,
          weighted_mode = weighted_mode
        )
      },
      .options = furrr_options(seed = TRUE)
    )
})

raw_file <- sprintf("nct_mc_raw_%s.rds", analysis_tag)
saveRDS(results, raw_file)

# -----------------------------------------------------------------------------
# 6. Summary table
# -----------------------------------------------------------------------------
summary_tbl <- results |>
  group_by(analysis_mode, topology, n_per_grp, threshold_gap) |>
  summarise(
    reject_nwinv = mean(nwinv_p < 0.05, na.rm = TRUE),
    reject_glstr = mean(glstr_p < 0.05, na.rm = TRUE),
    
    fp_edge_05   = mean(edge_fp_05,   na.rm = TRUE),
    fp_edge_bonf = mean(edge_fp_bonf, na.rm = TRUE),
    
    act_low      = mean(act_low,      na.rm = TRUE),
    act_high     = mean(act_high,     na.rm = TRUE),
    
    edge_low     = mean(edge_low,     na.rm = TRUE),
    edge_high    = mean(edge_high,    na.rm = TRUE),
    
    thresh_est_low  = mean(thresh_est_low,  na.rm = TRUE),
    thresh_est_high = mean(thresh_est_high, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  mutate(
    nwinv_flag = case_when(
      reject_nwinv <= 0.08 ~ "nominal",
      reject_nwinv <= 0.15 ~ "mild inflation",
      TRUE                 ~ "substantial inflation"
    ),
    glstr_flag = case_when(
      reject_glstr <= 0.08 ~ "nominal",
      reject_glstr <= 0.15 ~ "mild inflation",
      TRUE                 ~ "substantial inflation"
    ),
    edge_flag = case_when(
      fp_edge_05 <= 0.08 ~ "nominal",
      fp_edge_05 <= 0.15 ~ "mild inflation",
      TRUE               ~ "substantial inflation"
    )
  )

print(summary_tbl)

summary_file <- sprintf("nct_mc_summary_%s.csv", analysis_tag)
write.csv(summary_tbl, summary_file, row.names = FALSE)

# -----------------------------------------------------------------------------
# 7. Plot helper
# -----------------------------------------------------------------------------
base_plot <- function(data, y_var, title, subtitle, y_lab) {
  ggplot(
    data,
    aes(
      x      = factor(n_per_grp),
      y      = .data[[y_var]],
      colour = factor(threshold_gap),
      group  = factor(threshold_gap)
    )
  ) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line() +
    geom_point(size = 2) +
    facet_wrap(~ topology) +
    scale_colour_brewer(palette = "Dark2") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Sample size per group",
      y = y_lab,
      colour = "Threshold gap"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# 8. Plots
# -----------------------------------------------------------------------------
plot_structure <- base_plot(
  summary_tbl,
  "reject_nwinv",
  sprintf("NCT structure invariance rejection rate (%s)", analysis_tag),
  "Groups differ only in threshold; dashed = nominal 5%",
  "Rejection rate"
)

plot_strength <- base_plot(
  summary_tbl,
  "reject_glstr",
  sprintf("NCT global strength rejection rate (%s)", analysis_tag),
  "Groups differ only in threshold; dashed = nominal 5%",
  "Rejection rate"
)

plot_edge_fp <- base_plot(
  summary_tbl,
  "fp_edge_05",
  sprintf("Edge-level false-positive rate (%s)", analysis_tag),
  "Proportion of edges flagged at p < .05; dashed = nominal 5%",
  "False-positive rate"
)

ggsave(
  sprintf("nct_structure_%s.pdf", analysis_tag),
  plot_structure,
  width = 8.5,
  height = 4.8
)

ggsave(
  sprintf("nct_strength_%s.pdf", analysis_tag),
  plot_strength,
  width = 8.5,
  height = 4.8
)

ggsave(
  sprintf("nct_edge_fp_%s.pdf", analysis_tag),
  plot_edge_fp,
  width = 8.5,
  height = 4.8
)

# -----------------------------------------------------------------------------
# 9. Console summary
# -----------------------------------------------------------------------------
cat("\nDone. Files saved:\n")
cat(sprintf("  %s\n", raw_file))
cat(sprintf("  %s\n", summary_file))
cat(sprintf("  nct_structure_%s.pdf\n", analysis_tag))
cat(sprintf("  nct_strength_%s.pdf\n", analysis_tag))
cat(sprintf("  nct_edge_fp_%s.pdf\n", analysis_tag))

# -----------------------------------------------------------------------------
# 10. Optional: run both weighted and dichotomized versions
# -----------------------------------------------------------------------------
# To compare the standard weighted analysis with an adjacency-style version,
# rerun the script twice:
#
#   weighted_mode <- TRUE
#   weighted_mode <- FALSE
#
# The NCT docs define:
#   - nwinv.pval  = network structure invariance test
#   - glstrinv.pval = global strength invariance test
#   - einv.pvals = edge invariance p-values
# and note that weighted = FALSE dichotomizes the estimated networks.

