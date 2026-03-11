# ============================================================
# Slow-fast simulation utilities
# ============================================================
# Active simulation functions:
#   - simulate_slowfast_cw01_v3()
#       final simplified slow-fast simulator with optional exogenous /
#       endogenous jump shocks
#   - simulate_slowfast_cw01()
#       no-shock baseline version of the coupled slow-fast model
#   - simulate_slowfast_cw01_timevaryingPbase()
#       time-varying baseline version for hysteresis experiments
#   - make_ramp_Pbase()
#       helper for up-then-down baseline forcing
#
# Note on notation:
#   In the manuscript, the baseline shock hazard is written as lambda_0.
#   In this code, that baseline hazard is implemented as lambda_exo.
#   The endogenous hazard increment is implemented as
#     lambda1 * max(0, m_slow - m_crit).
# ============================================================

plan(multisession, workers = max(1, parallel::detectCores() - 1))

# ------------------------------------------------------------
# Final slow-fast simulator (v3)
# ------------------------------------------------------------
# Per slow step:
#   1) update the fast Curie-Weiss symptom layer at fixed P_t
#   2) update a slow symptom trace m_slow via EWMA
#   3) update the slow context P via:
#        - mean reversion toward P_base
#        - feedback from m_slow
#        - Gaussian diffusion
#        - optional jump shocks
#
# Shock modes:
#   - "none" : no shocks
#   - "exo"  : exogenous shocks with constant hazard lambda_exo
#   - "endo" : endogenous shocks with hazard
#                lambda1 * max(0, m_slow[t] - m_crit)
#   - "both" : one shock maximum per step using total hazard
#                lambda_exo + lambda_endo(t),
#              with type assigned by hazard share
simulate_slowfast_cw01_v3 <- function(
    T_steps  = 4000,
    dt       = 0.02,
    P_base   = 0,
    P0       = 0,
    m0       = 0.05,
    kappa    = 0.2,
    b        = 0.15,
    m_star   = 0.25,
    sigmaP   = 0.04,
    lambda_m = 0.001,
    
    # fast-layer parameters
    n_nodes = 100,
    sweeps  = 40,
    beta    = 1.0,
    J       = 1.0,
    h0      = 0.0,
    gammaP  = 1.0,
    seed    = 1,
    
    # shock mode
    shock_mode = c("none", "exo", "endo", "both"),
    
    # exogenous shocks
    # manuscript notation: this corresponds to baseline hazard lambda_0
    lambda_exo   = 0.0,
    shock_mu_exo = 0.0,
    shock_sd_exo = 0.0,
    
    # endogenous shocks
    lambda1       = 0.0,
    m_crit        = 0.6,
    shock_mu_endo = 0.0,
    shock_sd_endo = 0.0
){
  set.seed(seed)
  shock_mode <- match.arg(shock_mode)
  
  P      <- numeric(T_steps)
  m_fast <- numeric(T_steps)
  m_slow <- numeric(T_steps)
  
  shock_exo  <- numeric(T_steps)
  shock_endo <- numeric(T_steps)
  shock_any  <- integer(T_steps)
  
  P[1]      <- P0
  m_fast[1] <- m0
  m_slow[1] <- m0
  
  s <- stats::rbinom(n_nodes, 1, m0)
  
  hazard_to_p <- function(lambda){
    lambda <- max(0, lambda)
    p <- 1 - exp(-lambda * dt)
    if (!is.finite(p)) p <- 0
    pmin(pmax(p, 0), 1)
  }
  
  for (t in 1:(T_steps - 1)) {
    
    # 1) fast-layer update at fixed P_t
    s <- ising_fast_step01(
      P = P[t], s = s,
      beta = beta, J = J, h0 = h0, gammaP = gammaP,
      sweeps = sweeps
    )
    m_fast[t] <- mean(s)
    
    # 2) slow symptom memory
    m_slow[t+1] <- m_slow[t] + lambda_m * (m_fast[t] - m_slow[t])
    
    # 3) continuous slow drift + diffusion
    drift <- -kappa * (P[t] - P_base) + b * (m_slow[t] - m_star)
    dW    <- stats::rnorm(1, 0, sigmaP * sqrt(dt))
    
    # 4) optional shock process
    if (shock_mode != "none") {
      
      lambda_endo_t <- lambda1 * max(0, m_slow[t] - m_crit)
      lambda_endo_t <- max(0, lambda_endo_t)
      
      if (shock_mode == "exo") {
        
        p <- hazard_to_p(lambda_exo)
        if (stats::rbinom(1, 1, p) == 1) {
          A <- if (shock_sd_exo > 0) stats::rnorm(1, shock_mu_exo, shock_sd_exo) else shock_mu_exo
          if (!is.finite(A)) A <- 0
          shock_exo[t] <- A
          shock_any[t] <- 1L
        }
        
      } else if (shock_mode == "endo") {
        
        p <- hazard_to_p(lambda_endo_t)
        if (stats::rbinom(1, 1, p) == 1) {
          A <- if (shock_sd_endo > 0) stats::rnorm(1, shock_mu_endo, shock_sd_endo) else shock_mu_endo
          if (!is.finite(A)) A <- 0
          shock_endo[t] <- A
          shock_any[t]  <- 1L
        }
        
      } else if (shock_mode == "both") {
        
        lambda_tot <- lambda_exo + lambda_endo_t
        p <- hazard_to_p(lambda_tot)
        
        if (stats::rbinom(1, 1, p) == 1) {
          
          w_endo  <- if (lambda_tot > 0) lambda_endo_t / lambda_tot else 0
          is_endo <- (stats::runif(1) < w_endo)
          
          if (is_endo) {
            A <- if (shock_sd_endo > 0) stats::rnorm(1, shock_mu_endo, shock_sd_endo) else shock_mu_endo
            if (!is.finite(A)) A <- 0
            shock_endo[t] <- A
          } else {
            A <- if (shock_sd_exo > 0) stats::rnorm(1, shock_mu_exo, shock_sd_exo) else shock_mu_exo
            if (!is.finite(A)) A <- 0
            shock_exo[t] <- A
          }
          
          shock_any[t] <- 1L
        }
      }
    }
    
    # 5) update slow context
    P[t+1] <- P[t] + dt * drift + dW + shock_exo[t] + shock_endo[t]
    
    if (!is.finite(P[t+1])) {
      stop(sprintf(
        "Non-finite P at t=%d (P=%g drift=%g dW=%g exo=%g endo=%g)",
        t, P[t], drift, dW, shock_exo[t], shock_endo[t]
      ))
    }
  }
  
  m_fast[T_steps] <- m_fast[T_steps - 1]
  
  tibble::tibble(
    t = 1:T_steps,
    P = P,
    m = m_fast,
    m_slow = m_slow,
    shock_exo  = shock_exo,
    shock_endo = shock_endo,
    shock_any  = shock_any
  )
}

# ------------------------------------------------------------
# No-shock baseline simulator
# ------------------------------------------------------------
simulate_slowfast_cw01 <- function(
    T_steps  = 4000,
    dt       = 0.02,
    P_base   = 0,
    P0       = 0,
    m0       = 0.05,
    kappa    = 0.2,
    b        = 0.15,
    m_star   = 0.25,
    sigmaP   = 0.04,
    lambda_m = 0.001,
    
    # fast-layer parameters
    n_nodes = 100,
    sweeps  = 40,
    beta    = 1.0,
    J       = 1.0,
    h0      = 0.0,
    gammaP  = 1.0,
    seed    = 1
){
  set.seed(seed)
  
  P      <- numeric(T_steps)
  m_fast <- numeric(T_steps)
  m_slow <- numeric(T_steps)
  
  P[1]      <- P0
  m_fast[1] <- m0
  m_slow[1] <- m0
  
  s <- stats::rbinom(n_nodes, 1, m0)
  
  for (t in 1:(T_steps - 1)) {
    
    s <- ising_fast_step01(
      P = P[t], s = s,
      beta = beta, J = J, h0 = h0, gammaP = gammaP,
      sweeps = sweeps
    )
    m_fast[t] <- mean(s)
    
    m_slow[t+1] <- m_slow[t] + lambda_m * (m_fast[t] - m_slow[t])
    
    drift <- -kappa * (P[t] - P_base) + b * (m_slow[t] - m_star)
    
    P[t+1] <- P[t] + dt * drift + stats::rnorm(1, 0, sigmaP * sqrt(dt))
  }
  
  m_fast[T_steps] <- m_fast[T_steps - 1]
  
  tibble::tibble(
    t = 1:T_steps,
    P = P,
    m = m_fast,
    m_slow = m_slow
  )
}

# ------------------------------------------------------------
# Time-varying baseline simulator for hysteresis experiments
# ------------------------------------------------------------
simulate_slowfast_cw01_timevaryingPbase <- function(
    T_steps    = 6000,
    dt         = 0.02,
    P_base_fun = function(t) 0,
    P0         = 0,
    m0         = 0.05,
    kappa      = 0.2,
    b          = 0.0,
    m_star     = 0.25,
    sigmaP     = 0.04,
    lambda_m   = 0.001,
    
    # fast-layer parameters
    n_nodes = 100,
    sweeps  = 40,
    beta    = 2.0,
    J       = 3.0,
    h0      = 0.0,
    gammaP  = 1.0,
    seed    = 1
){
  set.seed(seed)
  
  P      <- numeric(T_steps)
  m_fast <- numeric(T_steps)
  m_slow <- numeric(T_steps)
  P_base <- numeric(T_steps)
  
  P[1]      <- P0
  m_fast[1] <- m0
  m_slow[1] <- m0
  
  s <- stats::rbinom(n_nodes, 1, m0)
  
  for (t in 1:(T_steps - 1)) {
    P_base[t] <- P_base_fun(t)
    
    s <- ising_fast_step01(
      P = P[t], s = s,
      beta = beta, J = J, h0 = h0, gammaP = gammaP,
      sweeps = sweeps
    )
    m_fast[t] <- mean(s)
    
    m_slow[t+1] <- m_slow[t] + lambda_m * (m_fast[t] - m_slow[t])
    
    drift <- -kappa * (P[t] - P_base[t]) + b * (m_slow[t] - m_star)
    
    P[t+1] <- P[t] + dt * drift + stats::rnorm(1, 0, sigmaP * sqrt(dt))
  }
  
  P_base[T_steps] <- P_base_fun(T_steps)
  m_fast[T_steps] <- m_fast[T_steps - 1]
  
  tibble::tibble(
    t = 1:T_steps,
    P = P,
    P_base = P_base,
    m = m_fast,
    m_slow = m_slow
  )
}

# ------------------------------------------------------------
# Helper for ramped baseline forcing
# ------------------------------------------------------------
make_ramp_Pbase <- function(T_steps, Pmin, Pmax, frac_up = 0.5){
  Tup <- floor(T_steps * frac_up)
  
  function(t){
    if (t <= Tup){
      Pmin + (Pmax - Pmin) * (t - 1) / (Tup - 1)
    } else {
      td    <- t - Tup
      Tdown <- T_steps - Tup
      Pmax - (Pmax - Pmin) * (td - 1) / (Tdown - 1)
    }
  }
}