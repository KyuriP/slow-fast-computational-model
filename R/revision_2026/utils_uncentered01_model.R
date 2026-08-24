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
