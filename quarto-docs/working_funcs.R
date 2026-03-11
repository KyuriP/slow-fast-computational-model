#' Make group-level parameters for hybrid symptom–precarity model
#'
#' This function estimates empirical parameters for a given subgroup of the data 
#' (e.g., Dutch vs. NonDutch) that will be used in the hybrid Ising–GGM simulation model. 
#' It extracts symptom and precarity networks, calibrates cross-domain links, and 
#' stores empirical residuals for noise resampling.
#'
#' @param df_group A data frame containing data for one subgroup (e.g., Dutch or NonDutch).
#' @param precarity_vars Character vector of column names representing precarity variables.
#' @param symptom_vars Character vector of column names representing symptom variables.
#' @param glasso_gamma Numeric, EBICglasso tuning parameter for estimating the precarity network (default = 0.5).
#' @param shrink_W Numeric scaling factor applied to symptom network weights after spectral radius normalization (default = 0.8).
#' @param shrink_Ap Numeric scaling factor applied to precarity network weights after spectral radius normalization (default = 0.5).
#' @param scale_Gamma Numeric multiplier applied to regression-derived precarity → symptom effects (default = 1.0).
#' @param scale_B Numeric multiplier applied to regression-derived symptom → precarity effects (default = 0.2).
#'
#' @details
#' The function proceeds in several steps:
#' \itemize{
#'   \item Standardizes precarity variables (mean 0, SD 1).
#'   \item Binarizes symptom variables (0 = inactive, 1 = active).
#'   \item Estimates a symptom Ising model using \code{IsingFit}, returning 
#'   weight matrix \eqn{W} and thresholds \eqn{\theta}.
#'   \item Estimates a precarity network (\eqn{A_p}) with EBICglasso.
#'   \item Fits logistic regressions of each symptom on precarity factors to form 
#'   the precarity → symptom matrix (\eqn{\Gamma}).
#'   \item Fits linear regressions of each precarity variable on symptoms to form 
#'   the symptom → precarity matrix (\eqn{B}), and stores regression residuals 
#'   as empirical noise pools.
#' }
#'
#' Spectral radius normalization ensures that both the symptom network (W) and 
#' precarity network (Ap) remain stable during simulation. Residual pools allow 
#' simulation to resample empirical noise (skewed/heavy-tailed) instead of assuming Gaussian innovations.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{mu_p}{Means of precarity variables.}
#'   \item{sd_p}{Standard deviations of precarity variables.}
#'   \item{W}{Symptom network weight matrix.}
#'   \item{theta}{Symptom thresholds from Ising model.}
#'   \item{Ap}{Precarity network adjacency matrix.}
#'   \item{Gamma}{Matrix of precarity → symptom effects.}
#'   \item{B}{Matrix of symptom → precarity effects.}
#'   \item{resid_pools}{List of empirical residuals per precarity variable for resampling.}
#'   \item{precarity_vars}{Names of precarity variables.}
#'   \item{symptom_vars}{Names of symptom variables.}
#' }
#'
#' @examples
#' \dontrun{
#' params_dutch <- make_group_params(df_dutch, 
#'                                    precarity_vars = c("income", "discrimination"),
#'                                    symptom_vars = c("depression", "sleep"))
#' }
#'
#' @seealso \code{\link[IsingFit]{IsingFit}}, \code{\link[qgraph]{EBICglasso}}
#' 
#' @export

make_group_params <- function(df_group,
                              precarity_vars, symptom_vars,
                              glasso_gamma = 0.5,
                              shrink_W = 0.8,
                              shrink_Ap = 0.5,
                              scale_Gamma = 1.0,
                              scale_B = 0.2) {
  
  # Empirical distribution of precarity (raw units)
  mu_p  <- sapply(precarity_vars, function(v) mean(df_group[[v]], na.rm = TRUE))
  sd_p  <- sapply(precarity_vars, function(v) sd(df_group[[v]],   na.rm = TRUE))
  # Replace any zero SD with a tiny positive number to avoid division by zero later
  sd_p[sd_p == 0] <- 1e-6
  
  # Z-score precarity (for modeling Γ and A_p)
  # We use z-scores to (a) estimate the precarity–precarity network (A_p) on a common scale and 
  # (b) get Γ (precarity→symptoms) as per-SD effects.
  # We still simulate precarity in raw units, but standardize internally when needed.
  Pz <- as.data.frame(mapply(function(x, m, s) (x - m)/s, df_group[precarity_vars], mu_p, sd_p,
                             SIMPLIFY = FALSE))
  colnames(Pz) <- precarity_vars
  
  # Binarize symptoms for Ising and GLMs
  Sbin <- df_group[, symptom_vars, drop = FALSE]
  Sbin[] <- lapply(Sbin, function(x) ifelse(x < 1, 0, 1))
  
  # Symptom Ising model (W & θ) + spectral shrink
  ising <- IsingFit::IsingFit(as.data.frame(Sbin), plot = FALSE, progressbar = FALSE)
  W <- ising$wei
  theta <- ising$thresholds
  rownames(W) <- colnames(W) <- symptom_vars
  names(theta) <- symptom_vars
  # Shrink W so its spectral radius becomes shrink_W (default 0.8). 
  # This is a stability knob: big W can make the network too “sticky” or explosive; 
  # scaling it keeps dynamics well-behaved without changing its pattern.
  spW <- max(Mod(eigen(W, only.values = TRUE)$values))
  if (spW > 0) W <- shrink_W * W / spW
  
  # Precarity GGM (A_p) on z-scores + spectral shrink 
  # Compute the correlation matrix of z-scored precarity.
  corP <- suppressWarnings(cor(Pz, use = "pairwise.complete.obs"))
  # Run EBICglasso to get a sparse precarity–precarity network A_p.
  Ap   <- qgraph::EBICglasso(corP, n = nrow(Pz), gamma = glasso_gamma)
  # Zero the diagonal (no self-loops).
  diag(Ap) <- 0
  # Shrink by spectral radius to shrink_Ap (default 0.5) for stability of the slow precarity dynamics.
  spA <- max(Mod(eigen(Ap, only.values = TRUE)$values))
  if (spA > 0) Ap <- shrink_Ap * Ap / spA
  rownames(Ap) <- colnames(Ap) <- precarity_vars
  
  # Γ: precarity → symptoms (node-specific external fields)
  # At simulation time, the external field for symptom i is h_i = theta_i + (Gamma %*% p_z)_i,   # so each precarity variable can shift the symptom’s effective threshold.
  Gamma <- matrix(0, nrow = length(symptom_vars), ncol = length(precarity_vars),
                  dimnames = list(symptom_vars, precarity_vars))
  # For each symptom sy, fit a logistic regression sy ~ all z-scored precarity.
  for (sy in symptom_vars) {
    dat <- cbind(y = Sbin[[sy]], Pz)
    fit <- try(suppressWarnings(glm(y ~ ., data = as.data.frame(dat), family = binomial())), silent = TRUE)
    # Store the non-intercept coefficients as row sy in Γ (dimensions: symptoms × precarity).
    co  <- if (inherits(fit, "try-error")) rep(0, length(precarity_vars)) else coef(fit)[-1]
    Gamma[sy, colnames(Pz)] <- co
  }
  # Fill NAs with 0 (e.g., separation or small samples).
  Gamma[is.na(Gamma)] <- 0
  # Optionally multiply by scale_Gamma to globally strengthen/soften precarity’s push on symptoms.
  Gamma <- Gamma * scale_Gamma
  
  # B: symptoms → precarity (feedback in raw units) + residual pools
  B <- matrix(0, nrow = length(precarity_vars), ncol = length(symptom_vars),
              dimnames = list(precarity_vars, symptom_vars))
  resid_pools <- list()
  # For each precarity variable pv, regress it in raw units on all binary symptoms.
  for (pv in precarity_vars) {
    dat <- cbind(y = df_group[[pv]], Sbin)
    fit <- try(lm(y ~ ., data = as.data.frame(dat)), silent = TRUE)
    # Store non-intercept coefficients as the row of B (dimensions: precarity × symptoms).
    # These are the feedback effects: when symptoms are on, they push precarity up/down in raw units.
    if (inherits(fit, "try-error")) {
      B[pv, ] <- 0
      # Build a residual pool for each precarity variable from the linear model residuals:
      ## Later, the simulator can resample from these real residuals instead of using Gaussian noise.
      ## That preserves skew/heavy tails of shocks (more realistic).
      ## If too few residuals, fall back to a small synthetic pool based on the variable’s SD.
      resid_pools[[pv]] <- rep(0, sum(complete.cases(dat)))
    } else {
      co <- coef(fit)[-1]
      B[pv, colnames(Sbin)] <- co
      r <- residuals(fit)
      r <- r[is.finite(r)]
      if (length(r) < 5) r <- rep(sd(df_group[[pv]], na.rm = TRUE), 20)
      resid_pools[[pv]] <- r - mean(r, na.rm = TRUE)
    }
  }
  # Replace NAs with 0 and optionally multiply B by scale_B to control feedback strength.
  B[is.na(B)] <- 0
  B <- B * scale_B
  
  # Return a tidy parameter bundle
  list(
    # mu_p, sd_p: to convert between raw and z-scores during simulation.
    mu_p = mu_p, sd_p = sd_p,
    # W, theta: Ising symptom network.
    W = W, theta = theta,
    # Ap: precarity network on z-scores (used to compute drift, then mapped back to raw units).
    Ap = Ap,
    # Gamma: how precarity (z) shifts symptom thresholds.
    Gamma = Gamma,
    # how symptoms push precarity (raw).
    B = B,
    # empirical shock distributions for each precarity variable.
    resid_pools = resid_pools,
    # Variable name vectors to keep dimensions aligned downstream.
    precarity_vars = precarity_vars,
    symptom_vars = symptom_vars)
}




#' Simulate hybrid symptom–precarity dynamics
#'
#' This function simulates coupled dynamics of symptoms (fast, Ising network) 
#' and precarity factors (slow, linear Gaussian/empirical process with feedback), 
#' using parameter packs created by \code{make_group_params}.
#'
#' @param par A parameter list produced by \code{make_group_params}, containing
#'   symptom network (\code{W}, \code{theta}), precarity network (\code{Ap}), 
#'   cross-domain matrices (\code{Gamma}, \code{B}), empirical means/SDs, 
#'   and residual pools.
#' @param T Integer. Number of simulation time steps (default = 300).
#' @param delta_p Numeric. Step size controlling how fast precarity evolves. 
#'   Smaller values slow precarity relative to symptoms (default = 0.005).
#' @param kappa Numeric. Mean reversion strength toward empirical precarity means. 
#'   Larger values keep precarity anchored near its mean (default = 0.5).
#' @param beta Numeric. Inverse temperature for the Ising update. Lower values 
#'   add randomness (symptoms flip more often), higher values make dynamics more 
#'   deterministic (default = 0.8).
#' @param noise_mode Character. Either \code{"resample"} to draw precarity shocks 
#'   from empirical regression residual pools (preserves skew/heavy tails), or 
#'   \code{"gauss"} to draw Gaussian noise (smoother shocks). 
#'   Default = c("resample","gauss") (argument matching).
#' @param noise_scale Numeric. Scale of Gaussian shocks relative to empirical SD 
#'   (used if \code{noise_mode="gauss"}, default = 0.02).
#' @param noise_rescale Numeric. Rescaling factor applied to empirical residuals 
#'   when \code{noise_mode="resample"} (default = 0.2). Use values < 1 to shrink 
#'   heavy-tailed shocks, > 1 to amplify them.
#' @param clamp_sd Numeric. If >0, clamp precarity values within mean ± clamp_sd * SD.
#'        If NULL, no clamping is applied (default = 2.5).
#' @param seed Integer. Random seed for reproducibility (default = 123).
#'
#' @details
#' The simulation evolves two coupled layers:
#' \itemize{
#'   \item \strong{Precarity layer (slow):} Each factor drifts toward its group mean 
#'         (\eqn{\mu_p}) with strength \code{kappa}, interacts with other factors via 
#'         \eqn{A_p} (estimated on z-scores), is nudged by current symptoms via \eqn{B}, 
#'         and receives stochastic shocks. Evolution is slowed by \code{delta_p}.
#'   \item \strong{Symptom layer (fast Ising):} At each step, symptoms update 
#'         according to their network couplings \eqn{W}, thresholds \eqn{\theta}, 
#'         and external fields from precarity (\eqn{\Gamma p_z}). Updates use 
#'         Glauber dynamics with inverse temperature \code{beta}.
#' }
#'
#' Precarity states are stored in raw units for interpretability, 
#' while z-scores are used internally for model-based updates. 
#' Residual pools (if used) allow resampling of shocks with the empirical 
#' distribution’s skewness and heavy tails; Gaussian mode gives smoother variability.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{precarity}{A \eqn{T \times n_p} matrix of precarity trajectories (raw units).}
#'   \item{symptoms}{A \eqn{T \times n_s} binary matrix of symptom activation states.}
#' }
#'
#' @examples
#' \dontrun{
#' sim_dutch <- simulate_hybrid_empirical(par_dutch, T = 500, 
#'                                        delta_p = 0.005, kappa = 0.5, beta = 0.7,
#'                                        noise_mode = "resample", noise_rescale = 0.2)
#' sim_nondutch <- simulate_hybrid_empirical(par_nondutch, T = 500,
#'                                           delta_p = 0.005, kappa = 0.5, beta = 0.7,
#'                                           noise_mode = "resample", noise_rescale = 0.2)
#' matplot(sim_dutch$precarity, type = "l")
#' image(t(sim_dutch$symptoms), col = c("white","black"))
#' }
#'
#' @seealso \code{\link{make_group_params}}
#' 
#' @export
simulate_hybrid_empirical <- function(par,
                                      T = 300,                  
                                      delta_p = 0.005,          
                                      kappa = 0.5,              
                                      beta = 0.8,               
                                      noise_mode = c("resample","gauss"),
                                      noise_scale = 0.2,       
                                      noise_rescale = 0.2,   
                                      clamp_sd = 2.5,
                                      seed = 123) {
  set.seed(seed)                            # reproducible randomness
  noise_mode <- match.arg(noise_mode)       # validate/choose noise mode
  
  pvars <- par$precarity_vars               # names of precarity variables
  svars <- par$symptom_vars                 # names of symptom variables
  np <- length(pvars); ns <- length(svars)  # counts
  
  # Precarity: start near group means (raw units) with mild random spread
  p <- sapply(pvars, function(v) sample(na.omit(par$mu_p[[v]] + par$sd_p[[v]] *
                                                  rnorm(1,0,0.5)), 1))
  p <- as.numeric(p)   # ensure vector
  
  # Symptoms: start as independent Bernoulli with base prevalence (here 0.2: ~20% prevalence)
  s_prev <- rep(0.2, ns); names(s_prev) <- svars
  s <- rbinom(ns, 1, s_prev); names(s) <- svars
  
  # Baseline precarity reference (empirical mean -> z=0)
  pz_baseline <- rep(0, np)
  
  # Storage for full trajectories over T steps
  P <- matrix(NA_real_, T, np, dimnames = list(NULL, pvars))
  S <- matrix(NA_real_, T, ns, dimnames = list(NULL, svars))
  
  # Random shock generator for precarity (one draw per variable)
  draw_noise <- function(v) {
    if (noise_mode == "resample") {
      # Empirical shocks: resample regression residuals (skew/heavy tails preserved)
      pool <- par$resid_pools[[v]]
      if (length(pool) == 0) return(0) # fallback if no residuals
      return(as.numeric(noise_rescale * sample(pool, 1, replace = TRUE))) # rescale factor applied
    } else {
      # Gaussian shocks: synthetic noise with SD scaled to empirical variability
      return(rnorm(1, 0, noise_scale * par$sd_p[[v]]))
    }
  }
  
  # ---- main simulation loop ------------------------------------------------
  for (t in 1:T) {
    # ---- Precarity (slow layer) -----------------------------------------------
    # Update rule:
    #   p_{t+1} = p_t + δp * [ (A_p z_t)           # network interactions
    #                        + (B s_t)             # feedback from symptoms
    #                        - κ (p_t - μ_p) ]     # mean reversion toward group mean
    #                        + sqrt(δp) * ε_t      # stochastic shocks

    #
    # All drivers (network, symptoms, shocks, and mean reversion) evolve on the
    # same slow timescale δp, keeping dynamics balanced and preventing runaway drift.
    pz <- (p - par$mu_p) / par$sd_p                       # z-scores for model terms
    
    # network and symptom drifts
    drift_A_raw <- as.numeric(par$Ap %*% pz) * par$sd_p   # precarity network drift (back to raw units)
    drift_B_raw <- as.numeric(par$B %*% s)                # symptom → precarity feedback (raw)
    eps <- vapply(pvars, draw_noise, numeric(1))          # stochastic shocks
    
    # apply update
    p <- as.numeric(p + delta_p * (drift_A_raw + drift_B_raw - kappa * (p - par$mu_p)) + sqrt(delta_p) * eps)
    
    
    # Optional clamp
    if (!is.null(clamp_sd) && clamp_sd > 0) {
      p <- pmax(pmin(p, par$mu_p + clamp_sd * par$sd_p),
                par$mu_p - clamp_sd * par$sd_p)
    }
    
    # ---- Symptoms (fast Ising layer)
    # External field update:
    #   h_i = theta_i + (Gamma %*% pz)_i
    #
    # Interpretation:
    # - More negative h_i → symptom harder to activate.
    # - With positive Gamma, if precarity increases (pz > 0),
    #   then h_i increases (less negative),
    #   making activation more likely (symptoms more active).
    h <- par$theta + as.numeric(par$Gamma %*% pz)   # external field (thresholds + precarity influence)
    eta <- as.numeric(par$W %*% s) + h                              # net input = neighbors + field
    prob <- inv_logit(2 * beta * eta)                               # activation probability
    s <- rbinom(ns, 1, pmin(pmax(prob, 1e-4), 1 - 1e-4))            # sample new states; avoid 0/1 saturation
    
    # Save trajectories at time t
    P[t, ] <- p
    S[t, ] <- s
  }
  list(precarity = P, symptoms = S)
}



## ---------------------------------
## plotting function
## ---------------------------------
plot_sim <- function(sim, group_name) {
  par(mfrow = c(2,1), mar = c(3,4,2,1))
  Z <- t(sim$symptoms)
  image(Z, zlim = c(0,1), breaks = c(-0.5,0.5,1.5),
        col = c("white","black"), useRaster = TRUE, axes = FALSE,
        main = paste(group_name, "Symptoms"))
  axis(1, at = seq(0,1,len=5), labels = round(seq(0, nrow(sim$symptoms), len=5)))
  axis(2, at = seq(0,1,len=nrow(Z)), labels = rownames(Z), las = 1)
  
  matplot(sim$precarity, type = "l", lty = 1, lwd = 2,
          ylab = "Precarity (raw units)", xlab = "Time",
          main = paste(group_name, "Precarity"))
  legend("topright", legend = colnames(sim$precarity), lty = 1,
         col = seq_len(ncol(sim$precarity)), cex = 0.7, bty = "n")
}




#' Overlay simulated precarity trajectories with empirical mean ± SD
#'
#' @param sim Simulation output from simulate_hybrid_empirical
#' @param par Parameter pack (from make_group_params)
#' @param group_name Label for plot title
plot_precarity_diagnostics <- function(sim, par, group_name) {
  df <- as.data.frame(sim$precarity)
  df$time <- 1:nrow(df)
  
  df_long <- tidyr::pivot_longer(df, -time, names_to = "variable", values_to = "value")
  
  emp_means <- par$mu_p
  emp_sds   <- par$sd_p
  
  ggplot(df_long, aes(x = time, y = value, color = variable)) +
    geom_line(alpha = 0.6) +
    geom_hline(aes(yintercept = emp_means[variable]), color = "black", linetype = "dashed") +
    geom_hline(aes(yintercept = emp_means[variable] + emp_sds[variable]), 
               color = "grey40", linetype = "dotted") +
    geom_hline(aes(yintercept = emp_means[variable] - emp_sds[variable]), 
               color = "grey40", linetype = "dotted") +
    facet_wrap(~ variable, scales = "free_y") +
    labs(title = paste("Precarity trajectories vs empirical ±1SD -", group_name),
         y = "Precarity (raw units)", x = "Time") +
    theme_minimal()
}



# Plot symptom network W (Ising)
plot_symptom_net <- function(par, group_name, ...) {
  qgraph(par$W,
         layout = "spring",
         labels = par$symptom_vars,
         theme = "colorblind",
         vsize = 6,
         title = paste("Symptom Network -", group_name))
}

# Plot precarity network A_p (GGM)
plot_precarity_net <- function(par, group_name, ...) {
  qgraph(par$Ap,
         layout = "spring",
         labels = par$precarity_vars,
         theme = "colorblind",
         vsize = 6,
         title = paste("Precarity Network -", group_name))
}




#' Diagnostic: Compare simulated vs empirical variability of precarity
#'
#' @param par Parameter pack from make_group_params().
#' @param sim Simulation output from simulate_hybrid_empirical().
#' @param group_name Label for printing/plotting.
#'
#' @return A data.frame with empirical SD, simulated SD, and ratio (sim/emp).
#' @export
check_precarity_sd <- function(par, sim, group_name = "Group") {
  sim_sd <- apply(sim$precarity, 2, sd, na.rm = TRUE)
  emp_sd <- par$sd_p
  ratio  <- sim_sd / emp_sd
  
  df <- data.frame(
    variable   = names(sim_sd),
    empirical  = round(emp_sd, 3),
    simulated  = round(sim_sd, 3),
    ratio      = round(ratio, 2),
    group      = group_name
  )
  
  print(df)
  invisible(df)
}



#' Plot diagnostic: Sim vs Empirical SD ratios
#'
#' @param results A data.frame produced by check_precarity_sd, possibly rbind-ed across runs
#' @export
plot_precarity_sd_ratios <- function(df) {
  ggplot(df, aes(x = variable, y = ratio_mean, color = group)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_point(position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = (sim_mean - 1.96*sim_sd) / empirical,
                      ymax = (sim_mean + 1.96*sim_sd) / empirical),
                  width = 0.15, position = position_dodge(width = 0.4)) +
    labs(title = "Precarity SD: simulated / empirical (mean ± 1.96 SE across runs)",
         y = "SD ratio", x = NULL) +
    theme_minimal() + scale_color_discrete()
}


#' Grid search over noise_rescale values for resample mode
#'
#' @param par Parameter pack from make_group_params().
#' @param rescale_grid Numeric vector of noise_rescale values to try.
#' @param T Integer. Simulation length (default 500 for speed).
#' @param seed Random seed (default 123).
#'
#' @return A combined data.frame of sim vs empirical SD ratios across grid.
#' @export
grid_search_rescale <- function(par, rescale_grid = c(0.1, 0.15, 0.2),
                                T = 500, seed = 123) {
  results <- list()
  for (r in rescale_grid) {
    sim <- simulate_hybrid_empirical(par, T = T, noise_mode = "resample",
                                     noise_rescale = r, seed = seed)
    df <- check_precarity_sd(par, sim, group_name = paste0("rescale=", r))
    results[[as.character(r)]] <- df
  }
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Plot results of grid search
#'
#' @param results Data.frame returned by grid_search_rescale().
#' @export
plot_rescale_grid <- function(results) {
  ggplot(results, aes(x = variable, y = ratio, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Noise rescale grid search: Simulated vs Empirical SD",
         y = "Simulated / Empirical SD", x = "Precarity Variable") +
    theme_minimal()
}


#' Joint diagnostics: symptoms + precarity
#'
#' @param par Parameter pack from make_group_params()
#' @param sim Simulation output from simulate_hybrid_empirical()
#' @param group_name Label for group
#' @param empirical_sym Optional: empirical prevalence of symptoms (named vector)
#' 
#' @return List of two data.frames: prevalence + precarity SD ratios
#' @export
check_joint_diagnostics <- function(par, sim, group_name = "Group", empirical_sym = NULL) {
  ## ---- Symptom prevalence ----
  sim_prev <- colMeans(sim$symptoms, na.rm = TRUE)
  df_prev <- data.frame(
    symptom    = names(sim_prev),
    simulated  = round(sim_prev, 3),
    group      = group_name,
    stringsAsFactors = FALSE
  )
  if (!is.null(empirical_sym)) {
    df_prev$empirical <- empirical_sym[names(sim_prev)]
    df_prev$ratio <- round(df_prev$simulated / df_prev$empirical, 2)
  }
  
  ## ---- Precarity variability ----
  df_prec <- check_precarity_sd(par, sim, group_name = group_name)
  
  list(prevalence = df_prev, precarity_sd = df_prec)
}


#' Plot simulated vs empirical symptom prevalence
#'
#' @param df Data.frame from check_joint_diagnostics()$prevalence
#' @export
# Symptom prevalence plot
plot_symptom_prevalence <- function(df) {
  ggplot(df, aes(x = variable)) +
    geom_point(aes(y = empirical), shape = 1, size = 2) +
    geom_point(aes(y = sim_mean, color = group), position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = pmax(sim_mean - 1.96*sim_sd, 0),
                      ymax = pmin(sim_mean + 1.96*sim_sd, 1),
                      color = group),
                  width = 0.15, position = position_dodge(width = 0.4)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Symptom prevalence: empirical vs simulated (mean ± 1.96 SE across runs)",
         y = "Prevalence", x = NULL) +
    theme_minimal() + scale_color_discrete()
}




plot_sim_summary <- function(sims, label) {
  df <- lapply(seq_along(sims), function(r) {
    d <- as.data.frame(sims[[r]]$precarity)
    d$time <- 1:nrow(d)
    d$rep  <- r
    d
  }) |> dplyr::bind_rows()
  
  df_long <- tidyr::pivot_longer(df, -c(time, rep),
                                 names_to = "variable", values_to = "value")
  df_summary <- df_long |> 
    group_by(time, variable) |> 
    summarise(mean = mean(value), sd = sd(value), .groups = "drop")
  
  ggplot(df_summary, aes(x = time, y = mean, color = variable, fill = variable)) +
    geom_line(size=1) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    labs(title = paste("Precarity dynamics -", label),
         y = "Precarity (raw units)", x = "Time") +
    theme_minimal()
}


plot_symptom_summary <- function(sims, label) {
  df <- lapply(seq_along(sims), function(r) {
    d <- as.data.frame(sims[[r]]$symptoms)
    d$time <- 1:nrow(d)
    d$rep  <- r
    d
  }) |> dplyr::bind_rows()
  
  df_long <- tidyr::pivot_longer(df, -c(time, rep),
                                 names_to = "symptom", values_to = "active")
  df_summary <- df_long |> 
    group_by(time, symptom) |> 
    summarise(mean = mean(active), sd = sd(active), .groups = "drop")
  
  ggplot(df_summary, aes(x = time, y = mean, color = symptom, fill = symptom)) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2, color = NA) +
    labs(title = paste("Symptom prevalence trajectories -", label),
         y = "Prevalence (proportion active)", x = "Time") +
    theme_minimal()
}








# ------------------------------------------------------------------
# frozen_stats_emp(): Compute fast-layer steady-state branches
#                     for empirical HELIUS parameter packs
# ------------------------------------------------------------------

#' Estimate frozen fast-layer statistics (empirical HELIUS model)
#'
#' @description
#' For each fixed precarity level P, simulate the empirical fast symptom Ising
#' network (given parameter pack `par`) to measure its mean activation `m(P)`.
#' Detect unimodal vs. bistable behavior and estimate upward/downward
#' switching rates.
#'
#' @param P_grid Numeric vector of precarity values to scan.
#' @param par Parameter pack from `make_group_params()`.
#' @param v Numeric unit vector defining the direction in precarity z-space.
#' @param steps Number of Ising update sweeps (default 4000).
#' @param thin Keep one observation every `thin` sweeps.
#' @param seed Integer for reproducibility.
#' @param beta Inverse temperature of Ising layer (default 1.0).
#'
#' @return Tibble with columns:
#'   * `P` precarity level  
#'   * `unimodal` TRUE/FALSE  
#'   * `m_low`, `m_high` estimated mean activation on each stable branch  
#'   * `k_up`, `k_down` estimated upward/downward switching rates
#'
#' @export
frozen_stats_emp <- function(P_grid, par, v,
                             steps = 4000, thin = 2, seed = 1, beta = 1.0) {
  
  # inner fast-layer simulator for fixed P
  simulate_fast_empirical <- function(P, par, v, steps = 4000,
                                      burn = 1000, thin = 2, beta = 1.0,
                                      seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    ns <- length(par$symptom_vars)
    s  <- rbinom(ns, 1, 0.2)
    zP <- P * v
    h  <- as.numeric(par$theta + par$Gamma %*% zP)
    rec <- numeric(0)
    for (sw in seq_len(steps)) {
      i <- sample.int(ns, 1)
      eta <- sum(par$W[i, ] * s) + h[i]
      p_i <- 1 / (1 + exp(-2 * beta * eta))
      s[i] <- rbinom(1, 1, pmin(pmax(p_i, 1e-6), 1 - 1e-6))
      if (sw > burn && (sw - burn) %% thin == 0) rec <- c(rec, mean(s))
    }
    rec
  }
  
  # loop over grid of P values
  out <- lapply(P_grid, function(P) {
    x <- tryCatch(simulate_fast_empirical(P, par, v,
                                          steps = steps, thin = thin,
                                          seed = seed + round(1000 * P),
                                          beta = beta),
                  error = function(e) NA_real_)
    if (length(x) < 30 || any(!is.finite(x)))
      return(data.frame(P = P, unimodal = NA,
                        m_low = NA, m_high = NA,
                        k_up = 0, k_down = 0))
    
    # detect unimodal vs. bimodal with dip test or clustering
    mod <- suppressWarnings(diptest::dip.test(x)$p.value)
    unimodal <- (mod > 0.05)
    
    if (unimodal) {
      data.frame(P = P, unimodal = TRUE,
                 m_low = mean(x), m_high = NA,
                 k_up = 0, k_down = 0)
    } else {
      # split high vs. low branch by midpoint
      m_low  <- quantile(x, 0.25, na.rm = TRUE)
      m_high <- quantile(x, 0.75, na.rm = TRUE)
      thr <- mean(c(m_low, m_high))
      labs <- as.integer(x >= thr)
      r <- rle(labs)
      low_dw  <- r$lengths[r$values == 0]
      high_dw <- r$lengths[r$values == 1]
      k_up   <- if (length(low_dw))  1 / mean(low_dw)  else 0
      k_down <- if (length(high_dw)) 1 / mean(high_dw) else 0
      data.frame(P = P, unimodal = FALSE,
                 m_low = m_low, m_high = m_high,
                 k_up = k_up, k_down = k_down)
    }
  })
  
  tibble::as_tibble(do.call(rbind, out))
}
