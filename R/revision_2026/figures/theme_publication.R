# ============================================================
# R/revision_2026/figures/theme_publication.R
# ============================================================
# Shared theme + palette for the four publication figures. Every
# fig*_publication.R script in this folder should source this file and
# use theme_pub() rather than theme_classic()/defaults directly, so the
# four figures read as one coherent set rather than four one-off plots.
#
# Modeled on the styling actually used in the old Figure 7
# (R/scripts/13_network_full_check.R's theme_pub(), panel letters, bottom
# legend, log-scale axis, minimal gridlines) -- reproduced/extended here
# rather than sourcing that script directly, since it's a protected
# legacy script tied to the pre-revision manuscript.
#
# qgraph figures (the network-trio panels) are base-R graphics, not
# ggplot, so this theme doesn't apply to them directly -- their styling
# recipe (spring layout fixed from the true network, posCol/negCol,
# edge.width, vsize, label.cex, theme="classic") is documented inline in
# fig4_network_trio.R instead, matching the old Figure 8 recipe exactly.
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# ------------------------------------------------------------------------
# Global final-polish style tokens (locked 2026-08-25 final figure-design
# pass). Every fig*.R script should reference these rather than hardcoding
# its own sizes/widths/alphas, so the whole set reads as one consistent
# journal-style figure family instead of four separately-tuned plots.
# Individual panels MAY still deviate slightly (e.g. a smaller point size
# for a very dense panel) but should start from these values.
# ------------------------------------------------------------------------
base_size_panel_title <- 10.5
base_size_axis_title  <- 10
base_size_axis_text   <- 9
base_size_label       <- 8.5

main_line_width       <- 0.65   # summary/mean lines
secondary_line_width  <- 0.55   # secondary comparison lines (e.g. dashed reference)
spaghetti_width        <- 0.18   # individual-chain background trajectories
spaghetti_alpha         <- 0.10
ribbon_alpha           <- 0.08
point_size_main        <- 2.2
point_size_metric      <- 2.4
grid_colour            <- "grey90"
axis_colour            <- "grey15"

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "plain", size = base_size_panel_title, hjust = 0),
      plot.subtitle = element_text(size = base_size - 2, colour = "grey40", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(face = "plain", size = base_size, hjust = 0),
      legend.title = element_text(size = base_size_axis_text),
      legend.text = element_text(size = base_size_axis_text),
      legend.position = "bottom",
      axis.title = element_text(size = base_size_axis_title),
      axis.text = element_text(size = base_size_axis_text, colour = axis_colour),
      panel.grid.major.y = element_line(colour = grid_colour, linewidth = 0.25),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 14, 10, 10)
    )
}

# Panel letters (A), (B), (C)... added as a plot title prefix rather than
# a separate annotation layer -- keeps each panel self-contained so
# patchwork::wrap_plots() doesn't need tag_levels bookkeeping to match.
panel_title <- function(letter, text) {
  bquote(bold(.(paste0("(", letter, ") "))) * .(text))
}

# ------------------------------------------------------------------------
# Palette -- locked, one set of colors reused across all four figures so
# the paper reads as one visual identity rather than four separate plots.
# Restrained (no rainbow/viridis-per-figure), each color tied to one
# recurring quantity across the whole simulation set.
# ------------------------------------------------------------------------
col_P        <- "#3B76AF"   # slow context / P_t
col_M        <- "#D97A2B"   # symptom burden / M_t or m_t
col_low      <- "#4C9A6A"   # low-context condition (muted green)
col_mid      <- "#7A8DA6"   # middle-context condition (muted blue-grey)
col_high     <- "#C85C5C"   # high-context condition (muted red)
col_off      <- "#7F7F7F"   # feedback off (gray)
col_on       <- "#2F5AA8"   # feedback on (darker blue)
col_true     <- "#5E8CC0"   # true network / true value reference
col_naive    <- "#D95F02"   # naive / symptom-only estimator
col_adjusted <- "#1B9E77"   # context-adjusted estimator

# Named vectors for direct use in scale_colour_manual(values = ...):
pal_context  <- c(low = col_low, middle = col_mid, high = col_high)
pal_feedback <- c(off = col_off, on = col_on)
# Estimation-arm palette: baseline reuses col_off (both are "reference,
# nothing being corrected" conditions); naive/adjusted use the locked
# naive/adjusted colors so Figure 4's summary panels and network trio
# read as the same two conditions everywhere they appear.
pal_arm      <- c(baseline = col_off, naive = col_naive, adjusted = col_adjusted)

# Network edge colors (qgraph posCol/negCol) -- kept as the original blue
# used in the already-approved network-trio figure (not tied to col_true,
# which is a separate "true value reference line" color used in Figure
# 4's summary panels -- no reason to risk changing a figure that already
# works to force palette purity here).
qgraph_pos_col <- "#1565C0"
qgraph_neg_col <- "#C62828"
