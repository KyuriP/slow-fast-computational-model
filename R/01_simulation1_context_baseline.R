# ============================================================
# 01_simulation1_context_baseline.R
# ============================================================
# Purpose
# -------
# Simulation 1 from Simulation_Spec_Sims1-3_2026-08.md / the manuscript's
# new "Simulation 1: Slow context shifts symptom burden under fixed
# coupling" subsection. Five things only, per instruction:
#   1. Define symptom names
#   2. Define tau, omega, gamma
#   3. Simulate trajectories for P = -0.6, 0, 0.6
#   4. Save summary data
#   5. Make one pilot figure
# No shocks, no feedback, no network estimation here.
#
# Model (matches the manuscript formula exactly):
#   logit Pr(S_{t+1,i} = 1 | S_{t,-i}, P) = tau_i + sum_j!=i omega_ij S_{t,j} + gamma_i P
#
# This is sequential single-site (Gibbs/heat-bath) updating, matching
# Cramer et al. (2016)'s own model structure -- NOT the exact-joint-draw
# approach used in 15_revised_model_pilot.R. Deliberate choice for this
# subsection: ties the paper's core model directly to the literature
# Denny wants it connected to.
#
# T_burn / T_post are left as explicit, easy-to-change constants below
# (per the manuscript table, "to be set in code") -- current values are
# pilot defaults, not final. Chain-convergence behavior should be checked
# before locking these into the Table~\ref{tab:sim1_parameters} value.
#
# Outputs
# -------
#   res/sim1/sim1_trajectories_raw.rds
#   res/sim1/sim1_summary.csv
#   figs/Sim1_pilot_burden_distribution.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ------------------------------------------------------------------------
# 1. Symptom names -- standard PHQ-9 items, matching N = 9 in the
#    manuscript's Simulation 1 subsection (Table~\ref{tab:sim1_parameters}).
# ------------------------------------------------------------------------
symptom_names <- c(
  "anhedonia",        # little interest or pleasure
  "depressed_mood",   # feeling down, depressed, hopeless
  "sleep",             # sleep trouble (insomnia/hypersomnia)
  "fatigue",           # fatigue or low energy
  "appetite",          # appetite change
  "worthlessness",     # feelings of worthlessness / guilt
  "concentration",     # concentration problems
  "psychomotor",       # psychomotor agitation/retardation
  "suicidal_ideation"  # thoughts of death or self-harm
)
N <- length(symptom_names)   # 9
stopifnot(N == 9L)

# ------------------------------------------------------------------------
# 2. tau, omega, gamma -- ONE fixed network used across all three P
#    conditions (per the manuscript text: "The same threshold vector tau,
#    coupling matrix omega, and context-loading vector gamma were used in
#    all conditions"). This is a placeholder network, not yet calibrated
#    to real Cramer et al. (2016) values (tracked separately: get the
#    actual fitted VATSPUD thresholds/weights from Denny). Everything
#    below this block should be a one-line swap once those are in hand.
# ------------------------------------------------------------------------
set.seed(20260824L)  # fixed seed -- this IS "the" Simulation 1 network

tau_lo <- -3.0; tau_hi <- -1.0
w_lo   <-  0.20; w_hi  <-  0.40
gamma_lo <- 0.6; gamma_hi <- 1.2
edge_density <- 0.40

tau <- setNames(runif(N, tau_lo, tau_hi), symptom_names)

omega <- matrix(0, N, N, dimnames = list(symptom_names, symptom_names))
pairs <- combn(N, 2)
n_edges <- round(edge_density * ncol(pairs))
sel <- sample(ncol(pairs), n_edges)
for (k in sel) {
  i <- pairs[1, k]; j <- pairs[2, k]
  w <- runif(1, w_lo, w_hi)
  omega[i, j] <- w; omega[j, i] <- w
}

gamma <- setNames(runif(N, gamma_lo, gamma_hi), symptom_names)

cat("tau:\n");   print(round(tau, 2))
cat("\ngamma:\n"); print(round(gamma, 2))
cat(sprintf("\nomega: %d nonzero edges out of %d possible (density = %.2f)\n",
            n_edges, ncol(pairs), n_edges / ncol(pairs)))

# ------------------------------------------------------------------------
# 3. Simulate trajectories for P = -0.6, 0, 0.6
# ------------------------------------------------------------------------
P_values <- c(-0.6, 0, 0.6)
names(P_values) <- c("low", "middle", "high")

# Pilot defaults -- NOT final. Re-examine chain traces before locking
# these into the manuscript table.
T_burn    <- 200L
T_post    <- 200L
n_chains  <- 200L   # independent trajectories per P condition

# One sweep = N single-site (Gibbs/heat-bath) updates, vectorized across
# all n_chains simultaneously. plogis() is the inverse-logit function.
run_chains <- function(P, tau, omega, gamma, n_chains, T_burn, T_post) {
  S <- matrix(rbinom(n_chains * N, 1, 0.5), nrow = n_chains, ncol = N)  # random init
  total_sweeps <- T_burn + T_post
  M_post <- matrix(NA_real_, nrow = n_chains, ncol = T_post)

  for (t in seq_len(total_sweeps)) {
    site_order <- sample.int(N)  # random sweep order each step
    for (i in site_order) {
      others <- setdiff(seq_len(N), i)
      input_i <- tau[i] + S[, others, drop = FALSE] %*% omega[i, others] + gamma[i] * P
      p_active <- plogis(as.numeric(input_i))
      S[, i] <- rbinom(n_chains, 1, p_active)
    }
    if (t > T_burn) M_post[, t - T_burn] <- rowSums(S)
  }
  M_post  # n_chains x T_post matrix of symptom burden M_t
}

set.seed(2026L)
traj_list <- lapply(names(P_values), function(lbl) {
  P <- P_values[[lbl]]
  cat(sprintf("Simulating P = %s (%.1f): %d chains x %d burn-in + %d post-burn-in sweeps...\n",
              lbl, P, n_chains, T_burn, T_post))
  M_post <- run_chains(P, tau, omega, gamma, n_chains, T_burn, T_post)
  tibble(condition = lbl, P = P, chain = rep(seq_len(n_chains), T_post),
         sweep = rep(seq_len(T_post), each = n_chains), M = as.vector(M_post))
})
traj <- bind_rows(traj_list) |> mutate(m = M / N)

# ------------------------------------------------------------------------
# 4. Save summary data
# ------------------------------------------------------------------------
dir.create("res/sim1", recursive = TRUE, showWarnings = FALSE)
dir.create("figs", showWarnings = FALSE)
saveRDS(traj, "res/sim1/sim1_trajectories_raw.rds")

summary_tbl <- traj |>
  group_by(condition, P) |>
  summarise(
    mean_M = mean(M), sd_M = sd(M),
    mean_m = mean(m),
    pr_high_burden = mean(M >= N / 2),
    .groups = "drop"
  )
write.csv(summary_tbl, "res/sim1/sim1_summary.csv", row.names = FALSE)

cat("\n=== SIMULATION 1 SUMMARY ===\n")
print(summary_tbl)
cat("\nSame omega throughout (symptom-symptom coupling unchanged) -- any\n")
cat("difference in mean_M / pr_high_burden across rows above is generated\n")
cat("entirely by the slow contextual field P.\n")

# ------------------------------------------------------------------------
# 5. One pilot figure: burden distribution by condition
# ------------------------------------------------------------------------
p1 <- ggplot(traj, aes(x = M, fill = condition)) +
  geom_histogram(binwidth = 1, position = "identity", alpha = 0.55, colour = "white") +
  scale_fill_viridis_d(name = "Context (P)",
                        breaks = names(P_values),
                        labels = sprintf("%s (P=%.1f)", names(P_values), P_values)) +
  labs(x = "Symptom burden (M = number of active symptoms)", y = "Count",
       title = "Simulation 1 (pilot): symptom burden distribution under fixed coupling",
       subtitle = sprintf("N=%d symptoms, same tau/omega/gamma across conditions, only P differs", N)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figs/Sim1_pilot_burden_distribution.pdf", p1, width = 8, height = 5.5)
ggsave("figs/Sim1_pilot_burden_distribution.png", p1, width = 8, height = 5.5, dpi = 200)

cat("\nDone. Files:\n")
cat("  res/sim1/sim1_trajectories_raw.rds\n")
cat("  res/sim1/sim1_summary.csv\n")
cat("  figs/Sim1_pilot_burden_distribution.pdf (+ .png)\n")
cat("\nBefore trusting the numbers: check chain traces for burn-in adequacy\n")
cat("(not yet plotted here -- add a trace plot for a few chains if pr_high_burden\n")
cat("or mean_M look unstable / still drifting at T_burn).\n")
