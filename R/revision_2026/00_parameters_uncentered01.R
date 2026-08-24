# ============================================================
# R/revision_2026/00_parameters_uncentered01.R
# ============================================================
# Symptom names, tau_i, omega_ij, gamma_i for the revised uncentered 0/1
# model. Replaces 01_create_finalparams.R for the revision track.
#
# NOTE: these tau/gamma/omega values are hand-set placeholders, not yet
# calibrated to real data. Tracked separately: get the actual IsingFit
# thresholds/weights Cramer et al. (2016) estimated from the VATSPUD data
# (ask Denny -- he's a co-author and likely has the fitted values; they
# aren't in machine-readable form on the PLOS page, only shown as a
# figure). Swap-in point once obtained: replace the tau/omega values
# below directly, keep the rest of the pipeline unchanged.
# ============================================================

symptoms <- c(
  "anhedonia",
  "depressed",
  "sleep",
  "energy",
  "appetite",
  "guilt",
  "concentration",
  "psychomotor",
  "suicidal"
)

N <- length(symptoms)

tau <- c(
  anhedonia     = -2.40,
  depressed     = -2.30,
  sleep         = -2.00,
  energy        = -1.90,
  appetite      = -2.50,
  guilt         = -2.80,
  concentration = -2.30,
  psychomotor   = -3.00,
  suicidal      = -4.00
) + 1.3  # pilot adjustment (was too floor-heavy at the original values --
         # mean_m stayed under 0.14 even at the high-context condition;
         # +1.3 is a first attempt at getting the middle condition closer
         # to mean_m ~ 0.3-0.4). Revert to the un-shifted values above if
         # this overshoots into ceiling saturation instead.

gamma <- c(
  anhedonia     = 0.90,
  depressed     = 1.00,
  sleep         = 0.70,
  energy        = 0.80,
  appetite      = 0.70,
  guilt         = 1.00,
  concentration = 0.80,
  psychomotor   = 0.70,
  suicidal      = 0.60
)

omega <- matrix(0, nrow = N, ncol = N)
rownames(omega) <- symptoms
colnames(omega) <- symptoms

add_edge <- function(mat, a, b, w) {
  mat[a, b] <- w
  mat[b, a] <- w
  mat
}

omega <- add_edge(omega, "anhedonia", "depressed", 0.45)
omega <- add_edge(omega, "depressed", "guilt", 0.35)
omega <- add_edge(omega, "depressed", "suicidal", 0.30)
omega <- add_edge(omega, "sleep", "energy", 0.40)
omega <- add_edge(omega, "energy", "concentration", 0.35)
omega <- add_edge(omega, "concentration", "psychomotor", 0.25)
omega <- add_edge(omega, "appetite", "energy", 0.25)
omega <- add_edge(omega, "guilt", "suicidal", 0.25)
omega <- add_edge(omega, "anhedonia", "energy", 0.25)
omega <- add_edge(omega, "anhedonia", "concentration", 0.20)

params_uncentered01 <- list(
  symptoms = symptoms,
  N = N,
  tau = tau,
  gamma = gamma,
  omega = omega
)

dir.create("res/revision_2026", recursive = TRUE, showWarnings = FALSE)
saveRDS(params_uncentered01, "res/revision_2026/params_uncentered01.rds")
