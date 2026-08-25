# ============================================================
# R/revision_2026/figures/fig4_network_trio.R
# ============================================================
# The "context looks like coupling" figure -- true / symptom-only /
# context-adjusted network diagrams side by side. This is the new Figure
# 4's network-trio panel, styled to match the old Figure 8 exactly: same
# spring layout (fixed from the true network) reused across all three
# panels so edge differences are visually comparable, same qgraph recipe
# (posCol/negCol, edge.width, vsize, label.cex, classic theme). Reproduced
# from R/scripts/14_network_full_check_graphs.R's plotting recipe rather
# than sourcing that script, since it's a protected legacy script tied to
# the pre-revision manuscript and uses a different (N=12, exogenous SD_P)
# design entirely.
#
# Data: res/revision_2026/sim4/sim4_raw.rds, from the SINGLE-RUN pilot
# (04_sim_network_estimation.R, n_person=10000) -- not the replicated
# version (04b). This mirrors how the old Figure 8 also used one
# representative estimate (averaged over 20 replicates there, a single
# n=10000 draw here) rather than the 30 x n=3000 replicated design used
# for the quantitative claim in Figure 4's summary panels. The replicated
# run (04b) never saved full edge matrices, only summary metrics, so it
# can't drive this figure directly -- if a replicate-averaged version is
# wanted instead, 04b would need to additionally accumulate and average
# omega_hat matrices across replicates the way 14 did.
#
# Uses only the "baseline" and "adjusted"/"naive" arms from sim4_raw.rds,
# relabeled "Symptom-only" / "Context-adjusted" to match old Figure 8's
# exact panel titles. The "baseline" (fixed-context, no-confounding) arm
# is not shown here -- it belongs in Figure 4's quantitative summary
# panels (by-arm bar/point plot), not this qualitative three-network
# comparison, matching the old figure's structure (which never had a
# fixed-context baseline panel either).
#
# Outputs
# -------
#   figs/revision_2026/Figure4_network_trio.pdf / .png
# ============================================================

suppressPackageStartupMessages({
  library(qgraph)
})

source("R/revision_2026/figures/theme_publication.R")   # qgraph_pos_col, qgraph_neg_col
source("R/revision_2026/00_parameters_uncentered01.R")   # symptoms, N

raw <- readRDS("res/revision_2026/sim4/sim4_raw.rds")

W_true <- raw$true_omega
W_naive <- raw$omega_hat$naive
W_adjusted <- raw$omega_hat$adjusted

# ------------------------------------------------------------------------
# Plotting threshold: hide edges below this magnitude so the naive panel
# isn't cluttered with near-zero sampling noise. Matches the old Figure
# 8's convention exactly. True edges in this model are Uniform(0.20,
# 0.45), so plot_min=0.15 keeps every real edge while suppressing the
# bulk of confounding-driven noise. Applied only to the two ESTIMATED
# matrices, never to W_true -- W_true has no noise to threshold and every
# non-edge is exactly 0 by construction.
# ------------------------------------------------------------------------
plot_min <- 0.15
zero_small <- function(W) { W[abs(W) < plot_min] <- 0; W }

W_true_plot <- W_true
W_naive_plot <- zero_small(W_naive)
W_adjusted_plot <- zero_small(W_adjusted)

mats <- list(True = W_true_plot, `Symptom-only` = W_naive_plot, `Context-adjusted` = W_adjusted_plot)
max_edge <- max(sapply(mats, function(m) max(abs(m))))

# Fixed layout from the true network -- same node positions in all three
# panels, so a reader can track a given symptom pair across panels.
L <- qgraph(W_true_plot, layout = "spring", DoNotPlot = TRUE)$layout

# PHQ-9-style abbreviations, keyed by symptom NAME (not position) so this
# stays correct even if `symptoms` is ever reordered upstream. Replaces
# the earlier auto-generated 4-char truncation (ANHE/DEPR/SLEE/ENER/
# APPE/GUIL/CONC/PSYC/SUIC), which read as slightly awkward/inconsistent
# -- these read as more standard, PHQ-9-recognizable abbreviations.
phq_abbrev <- c(
  anhedonia     = "ANH",
  depressed     = "DEP",
  sleep         = "SLP",
  energy        = "ENE",
  appetite      = "APP",
  guilt         = "GLT",
  concentration = "CON",
  psychomotor   = "MOT",
  suicidal      = "SUI"
)
short_labels <- unname(phq_abbrev[symptoms])
stopifnot(!anyNA(short_labels))  # catch silently if `symptoms` ever changes

# ------------------------------------------------------------------------
# Phantom-edge highlighting: this is the main content upgrade of this
# design pass. Previously every edge in every panel was the same blue
# (posCol), so the reader had to compare panels side by side to notice
# that the symptom-only (P omitted) panel fabricates extra coupling.
# Now: edges that are real in the true network stay blue; edges that
# appear in an ESTIMATED panel but do NOT exist in the true network
# ("phantom" edges, i.e. abs(W_true) below the plotting threshold while
# the estimate clears it) are drawn in orange -- the same colour already
# used for the naive/symptom-only estimator elsewhere in this figure set
# (col_naive from theme_publication.R), so the color itself now carries
# the "this is what naive estimation gets wrong" meaning across the whole
# figure, not just this one panel. Genuine negative edges (not phantom)
# keep the locked negCol, matching qgraph's normal sign-based coloring.
# Applied to BOTH estimated panels (naive and adjusted) using the same
# rule -- the adjusted panel should show few or no phantom edges, and
# seeing that directly (not just inferring it from the summary panel) is
# the point.
# ------------------------------------------------------------------------
phantom_col <- col_naive

build_edge_colors <- function(W_est, W_true_ref) {
  cmat <- matrix(qgraph_pos_col, nrow(W_est), ncol(W_est))
  cmat[W_est < 0] <- qgraph_neg_col
  is_phantom <- (abs(W_true_ref) < plot_min) & (abs(W_est) >= plot_min)
  cmat[is_phantom] <- phantom_col
  cmat
}

edge_col_naive    <- build_edge_colors(W_naive_plot, W_true_plot)
edge_col_adjusted <- build_edge_colors(W_adjusted_plot, W_true_plot)

# Panel titles, node labels, and margins all tightened per the
# design-pass review -- smaller margins and a smaller export canvas
# (12x4.3in -> 8.8x2.9in) remove the excess white space around the trio,
# and larger relative node size (vsize) fills more of the freed-up room.
#
# 2026-08-25 readability pass: vsize raised 9.5 -> 13 and label.cex/
# title.cex nudged up slightly (0.95->1.05, 1.0->1.08) -- nodes/labels
# were reported hard to read at the original size once the figure was
# scaled down to \textwidth in the compiled PDF. title.cex bumped to stay
# visually consistent with the larger panel-D/E titles in
# fig4_metric_strip.R once that script's own titles were enlarged to
# compensate for its narrower (0.82\textwidth) placement -- see that
# script's header note for the actual print-size arithmetic.
plot_one <- function(W, title, edge_color_mat = NULL) {
  # qgraph()'s first formal argument is literally named `input` (not `x`
  # or the first positional slot in the usual S3-plot-method sense) --
  # naming this list element "x" meant do.call() couldn't match it to
  # anything, `input` stayed unfilled, and qgraph errored on its own
  # required argument. Renamed to match qgraph's actual signature.
  args <- list(
    input = W, layout = L, maximum = max_edge, fade = TRUE,
    labels = short_labels,
    edge.width = 1.3, vsize = 12, label.cex = 1.05,
    title = title, title.cex = 1.08,
    theme = "classic", DoNotPlot = FALSE,
    # Margins widened (was c(1,1,2,1)) to compensate for the larger vsize:
    # the fixed spring layout's node positions were computed back when
    # vsize=9.5, so peripheral nodes (APP at top, SLP on the right) sat
    # close to the plot boundary already -- bigger circles at the same
    # positions with the same margin pushed past the edge, most visibly
    # in panel C. More margin pulls the drawn plotting region inward
    # without moving/shrinking the nodes themselves.
    mar = c(2, 2.2, 3.2, 2.2)
  )
  if (is.null(edge_color_mat)) {
    args$posCol <- qgraph_pos_col
    args$negCol <- qgraph_neg_col
  } else {
    args$edge.color <- edge_color_mat
  }
  do.call(qgraph, args)
}

# 2026-08-25: the 3-colour edge-legend row was dropped (redundant with the
# caption, and the reason it was removed), but the PHQ-9 abbreviation key
# was explicitly requested back -- it's genuinely useful in-figure (readers
# shouldn't have to hold ANH/DEP/SLP/... in their head or flip to the
# caption) and isn't the cluttered part. Kept as a single slim text line,
# not a full legend() row.
abbrev_key <- paste(sprintf("%s = %s", phq_abbrev, names(phq_abbrev)), collapse = "; ")

draw_abbrev_strip <- function() {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  text(0.5, 0.5, abbrev_key, cex = 0.75, col = "grey30")
}

draw_trio <- function() {
  layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE), heights = c(3.1, 0.32))
  plot_one(W_true_plot, "(A) True coupling")
  plot_one(W_naive_plot, "(B) Estimated network, P omitted", edge_col_naive)
  plot_one(W_adjusted_plot, "(C) Estimated network, P included", edge_col_adjusted)
  draw_abbrev_strip()
}

dir.create("figs/revision_2026", recursive = TRUE, showWarnings = FALSE)

# Canvas height nudged back up just enough for the slim abbreviation row
# (9.2x3.1 -> 9.2x3.42) -- much less than the earlier full legend row
# needed (3.65), since this is one line of text, not a legend + key.
pdf("figs/revision_2026/Figure4_network_trio.pdf", width = 9.2, height = 3.42)
draw_trio()
dev.off()

png("figs/revision_2026/Figure4_network_trio.png", width = 9.2, height = 3.42, units = "in", res = 220)
draw_trio()
dev.off()

cat(sprintf("True global strength: %.2f | Symptom-only: %.2f | Context-adjusted: %.2f\n",
            sum(abs(W_true[upper.tri(W_true)])), sum(abs(W_naive[upper.tri(W_naive)])),
            sum(abs(W_adjusted[upper.tri(W_adjusted)]))))

cat("\nDone. Files:\n")
cat("  figs/revision_2026/Figure4_network_trio.pdf (+ .png)\n")

cat("\nNote: symptom labels are PHQ-9-style abbreviations (", paste(short_labels, collapse=", "),
    ").\nOrange edges = phantom couplings (absent in the true network, but\n",
    "cross the plotting threshold in that panel's estimate).\n", sep = "")

cat("\nCheck the rendered PNG: qgraph's edge.color argument is being passed a\n")
cat("full NxN matrix (build_edge_colors()) rather than a vector -- this is a\n")
cat("documented qgraph input format, but if edges render in the wrong color\n")
cat("or qgraph errors on the matrix form, the fallback is to get the edge\n")
cat("list via qgraph(...)$Edgelist and pass edge.color as a vector aligned\n")
cat("to that edge order instead.\n")
