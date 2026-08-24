# ============================================================
# R/revision_2026/utils_uncentered01_model.R
# ============================================================
# Shared fast-layer functions for the revised uncentered 0/1 model.
# logit Pr(S_i = 1 | S_-i, P) = tau_i + sum_j!=i omega_ij S_j + gamma_i P
#
# Reusable functions only -- every simulation script in R/revision_2026/
# should call simulate_fast_sweep() rather than reimplementing the update
# rule, so we don't accidentally end up with different models in
# different scripts.
# ============================================================

simulate_fast_sweep <- function(S, tau, omega, gamma, P) {
  N <- length(S)
  update_order <- sample(seq_len(N), size = N, replace = FALSE)

  for (i in update_order) {
    eta_i <- tau[i] + sum(omega[i, ] * S) + gamma[i] * P
    p_i <- plogis(eta_i)
    S[i] <- rbinom(1, size = 1, prob = p_i)
  }

  S
}

symptom_burden <- function(S) {
  sum(S)
}

active_fraction <- function(S) {
  mean(S)
}

# ------------------------------------------------------------------------
# Network estimator: nodewise logistic regression, symmetrized.
# Adapted unchanged (logic-wise) from fit_edges() in
# R/scripts/13_network_full_check.R -- reproduced here rather than sourced
# from that file, since 13_network_full_check.R is a protected/untouched
# legacy script tied to the pre-revision manuscript. If Pvec is supplied,
# each nodewise regression conditions on it (context-adjusted estimate);
# if NULL, P is omitted (naive/pooled estimate, susceptible to omitted-
# context confounding when the sample pools across different P values).
# ------------------------------------------------------------------------
fit_edges <- function(S, Pvec = NULL) {
  N <- ncol(S)
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
