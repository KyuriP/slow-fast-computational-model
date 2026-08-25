# ============================================================
# Fast layer utilities: 0/1 Curie–Weiss (mean-field Ising)
# ============================================================
#
# Model ingredients
# -----------------
# Spins: s_i ∈ {0,1}, mean burden m = mean(s)
#
# Mean-field input (field) at fixed context P:
#   H(m, P) = J (m - 1/2) + h0 + gammaP * P
#
# Conditional activation rule:
#   Pr(s_i <- 1 | m, P) = sigmoid( beta * H(m, P) )
#
# Deterministic mean-field fixed points satisfy:
#   m = sigmoid( beta * ( J (m - 1/2) + h0 + gammaP * P ) )
#
# Simulation convention used here
# -------------------------------
# Within one slow-time step, the fast layer is updated by several Monte Carlo
# sweeps. In each sweep, the current mean burden m is computed once, the field
# H(m, P) is held fixed during that sweep, and then n asynchronous single-site
# updates are drawn with replacement. This is a sweep-level mean-field
# approximation, not a fully sequential site-by-site update with field
# recomputed after every elementary spin update.
# ============================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(purrr)
})

# ------------------------------------------------------------
# Numerically stable logistic sigmoid
# ------------------------------------------------------------
sigmoid <- function(x) {
  # Clamp to avoid overflow/underflow in exp()
  x <- pmax(pmin(x, 35), -35)
  1 / (1 + exp(-x))
}

# ------------------------------------------------------------
# One slow-step fast-layer update (multiple sweeps)
# ------------------------------------------------------------
#' Stochastic update for the 0/1 Curie-Weiss fast layer.
#'
#' @param P       Scalar slow context at the current slow step.
#' @param s       Numeric/integer vector of 0/1 spins.
#' @param beta    Inverse-temperature / sensitivity parameter.
#' @param J       Mean-field coupling strength.
#' @param h0      Baseline field.
#' @param gammaP  Coupling from slow context P to the fast layer.
#' @param sweeps  Number of Monte Carlo sweeps within one slow-time step.
#'
#' @details
#' In each sweep, the current mean burden m = mean(s) is computed once, the
#' field H(m, P) is held constant during that sweep, and then n single-site
#' updates are sampled with replacement. Thus the implementation is asynchronous
#' at the site-update level, but uses a sweep-level constant-field
#' mean-field approximation.
#'
#' @return Updated spin vector s.
ising_fast_step01 <- function(
    P,
    s,
    beta   = 1.0,
    J      = 1.0,
    h0     = 0.0,
    gammaP = 1.0,
    sweeps = 40
){
  n <- length(s)
  if (n <= 0) stop("ising_fast_step01(): s must have positive length.")
  if (!is.finite(P)) stop("ising_fast_step01(): non-finite P.")
  
  # Defensive check: require binary spins
  if (any(!s %in% c(0, 1))) {
    stop("ising_fast_step01(): s must contain only 0/1 values.")
  }
  
  for (sw in seq_len(sweeps)) {
    # Sweep-level mean-field approximation:
    # m and therefore H are frozen within each sweep.
    m <- mean(s)
    H <- J * (m - 0.5) + h0 + gammaP * P
    x <- beta * H
    
    if (!is.finite(x)) {
      stop("ising_fast_step01(): beta * H is not finite; check parameter values.")
    }
    
    p <- sigmoid(x)
    if (!is.finite(p) || p < 0 || p > 1) {
      stop("ising_fast_step01(): invalid Bernoulli probability.")
    }
    
    # n asynchronous single-site updates with replacement
    idx <- sample.int(n, size = n, replace = TRUE)
    s[idx] <- stats::rbinom(n = n, size = 1, prob = p)
  }
  
  s
}

# ------------------------------------------------------------
# Mean-field map and derivative
# ------------------------------------------------------------
#' Mean-field update map F(m; P).
F_map <- function(m, beta, J, h_eff) {
  sigmoid(beta * (J * (m - 0.5) + h_eff))
}

#' Derivative dF/dm.
#'
#' For J > 0, F'(m) >= 0, so discrete-time local stability of a fixed point
#' reduces from |F'(m*)| < 1 to F'(m*) < 1.
Fprime <- function(m, beta, J, h_eff) {
  p <- F_map(m, beta, J, h_eff)
  beta * J * p * (1 - p)
}

# ------------------------------------------------------------
# Fixed points m*(P) in [0,1] via bracketing + uniroot
# ------------------------------------------------------------
#' Compute fixed points for a single context value P.
#'
#' The fixed-point equation is
#'   m = F(m; P) = sigmoid(beta * (J (m - 1/2) + h0 + gammaP * P)).
#'
#' We solve f(m) = F(m; P) - m = 0 on [0, 1] by combining:
#'   (i) detection of near-zero grid values, and
#'   (ii) sign-change bracketing followed by uniroot().
#'
#' A fixed point is labeled stable when |F'(m*)| < 1. Because F'(m) >= 0
#' in this model for J > 0, this reduces to F'(m*) < 1.
#'
#' @return Tibble with columns P, m, stable.
cw01_fixed_points <- function(P, beta, J, h0, gammaP, grid_n = 801, tol = 1e-10) {
  
  h_eff <- h0 + gammaP * P
  f <- function(m) F_map(m, beta, J, h_eff) - m
  
  grid <- seq(0, 1, length.out = grid_n)
  vals <- f(grid)
  roots <- numeric(0)
  
  # (1) Catch near-exact roots on the grid
  near0 <- which(abs(vals) < tol)
  if (length(near0) > 0) {
    roots <- c(roots, grid[near0])
  }
  
  # (2) Bracket sign changes
  sgn <- sign(vals)
  sgn[sgn == 0] <- NA
  idx <- which(
    !is.na(sgn[-1]) &
      !is.na(sgn[-length(sgn)]) &
      (sgn[-1] != sgn[-length(sgn)])
  )
  
  for (k in idx) {
    a <- grid[k]
    b <- grid[k + 1]
    r <- try(stats::uniroot(f, c(a, b))$root, silent = TRUE)
    if (!inherits(r, "try-error")) {
      roots <- c(roots, r)
    }
  }
  
  roots <- sort(unique(round(roots, 10)))
  if (length(roots) == 0) return(NULL)
  
  stable_flag <- Fprime(roots, beta, J, h_eff) < 1
  
  tibble(P = P, m = roots, stable = stable_flag)
}

# ------------------------------------------------------------
# Bifurcation table across a grid of P values
# ------------------------------------------------------------
#' Compute equilibria across a grid of context values P.
#'
#' @return Tibble with columns P, m, stable.
cw01_bifurcation_curve <- function(P_grid, beta, J, h0, gammaP, grid_n = 801) {
  purrr::map_dfr(
    P_grid,
    ~ cw01_fixed_points(.x, beta, J, h0, gammaP, grid_n = grid_n)
  )
}