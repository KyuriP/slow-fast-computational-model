# Slow–Fast Coupling: Curie–Weiss Symptoms Under Context Dynamics

This repository contains the simulation code, parameter-tuning workflow, and figure-generation
scripts for a slow–fast computational model in which a fast Curie–Weiss symptom layer is coupled
to a slower contextual process with feedback, diffusion, and optional jump shocks.

The manuscript itself ("When Context Looks Like Coupling: A Slow–Fast Theory of Psychological
Dynamics") is written and maintained in Overleaf and is not tracked in this repository. This
repo holds the simulation code, cached results, and the figures generated from them.

## Repository structure

```
slow-fast-computational-model/
├── README.md
├── .gitignore
├── LICENSE
├── slow-fast-computational-model.Rproj
│
├── R/
│   ├── scripts/
│   │   ├── 01_create_finalparams.R          Recreates tuned parameter objects
│   │   ├── 02_main_sim.R                    Runs the four main scenarios (M, M+S, B, B+S)
│   │   ├── 03_create_heatmaps.R             Regime heatmaps over (βJ, P_base)
│   │   ├── 04_figure_bistability_geometry.R Bistability / fixed-point geometry figure
│   │   ├── 05_figures_make_all.R            Collects and exports the main figures
│   │   ├── 06_create_hysteresis_gif.R       Hysteresis animation
│   │   ├── 07_inspect_finalparams.R         Inspection / validation helper
│   │   ├── 08_(supp)NCT_analysis.R          Supplementary NCT check (threshold differences only)
│   │   ├── 09_network_estimation_check.R    Early network-estimation check (superseded)
│   │   ├── 10_network_window_check.R        Early window-based check (superseded)
│   │   ├── 11_network_feedback_check.R      Early feedback-only check (superseded)
│   │   ├── 12_network_timescale_check.R     Merged timescale+feedback check (superseded by 13)
│   │   ├── 13_network_full_check.R          Main-text network-estimation check (σ_P × window sweep)
│   │   └── 14_network_full_check_graphs.R   Generates the companion network-diagram figure
│   └── utils/
│       ├── utils_fastlayer.R
│       ├── utils_slowfast.R
│       ├── utils_diagnostics.R
│       ├── utils_plotting.R
│       ├── utils_legacy.R
│       └── utils_gridsearch.R
│
├── res/
│   ├── tuned/            Tuned parameter objects (final_params.rds, etc.)
│   ├── main/              Cached results for the four main scenarios
│   ├── heatmaps/          Cached (βJ, P_base) regime-map results
│   ├── hysteresis/        Cached hysteresis-ramp simulation
│   ├── network_check/     Cached results for all network-estimation check versions (09–14)
│   ├── NCT_analysis/      Cached results for the supplementary NCT check
│   └── archive/           Historical/intermediate results, kept for reference only
│
├── figs/
│   └── ... all figures currently used in the manuscript (see below) ...
│
├── img/
│   └── ... legacy/exploratory figures not used in the current manuscript ...
│
└── quarto-docs/
    └── ... exploratory .qmd notebooks and rendered .html versions ...
```

## Figures (`figs/`)

`figs/` contains exactly the figures referenced by the current manuscript, under the filenames
used in the LaTeX source:

| File | Manuscript figure |
|---|---|
| `slow-fast-model-illustration_v2.png` | Fig. 1 — conceptual schematic |
| `Figure_AB_fixedpoints_bifurcation_v2.pdf` | Fig. 2 — bistability / bifurcation geometry |
| `Figure_timeseries_4panel_M.pdf` | Fig. 3 — representative trajectories (M/M+S/B/B+S) |
| `Figure_overlay_4panel_M.pdf` | Fig. 4 — trajectories on the equilibrium curve |
| `Figure_ensemble_ribbons_M.pdf` | Fig. 5 — across-replicate ensemble summaries |
| `Figure_heatmaps_2panel_plasma_clean.pdf` | Fig. 6 — regime heatmaps |
| `Figure7_network_estimation_check.pdf` | Fig. 7 — network-estimation check (σ_P × window) |
| `Figure8_network_estimation_check_graphs.pdf` | Fig. 8 — example estimated networks |
| `Figure_hysteresis.pdf` | Appendix B — hysteresis under an exogenous ramp |
| `nct_structure_weighted.pdf`, `nct_strength_weighted.pdf`, `nct_edge_fp_weighted.pdf` | Appendix D — supplementary NCT check |

Each `.pdf` has a matching `.png` where one was generated, for quick preview.

`img/` holds earlier drafts, exploratory plots, and superseded figure versions (including the
abandoned exogenous-context network-check figures, `Figure_S_network_timescale_check.*` and
`Figure_S_network_window_check.*`, both replaced by Fig. 7/8). Nothing in `img/` is referenced by
the current manuscript; it is kept for reference only.

## Main scenario scripts

* `01_create_finalparams.R` — Recreates tuned parameter objects and writes `final_params.rds`.
* `02_main_sim.R` — Runs the four main manuscript scenarios (M, M+S, B, B+S).
* `03_create_heatmaps.R` — Generates heatmap-based summaries across parameter ranges.
* `04_figure_bistability_geometry.R` — Produces the bistability/geometry figure.
* `05_figures_make_all.R` — Collects and exports the main-text figures.
* `06_create_hysteresis_gif.R` — Generates the hysteresis animation.
* `07_inspect_finalparams.R` — Quick inspection/validation helper for tuned parameters.
* `08_(supp)NCT_analysis.R` — Supplementary NCT check under threshold differences only.
* `13_network_full_check.R` — Main-text network-estimation check: simulates the coupled
  slow–fast model directly, sweeping the diffusion parameter σ_P and observation window W, and
  compares a symptom-only versus context-adjusted network estimator (Fig. 7).
* `14_network_full_check_graphs.R` — Regenerates the companion example-network figure (Fig. 8)
  from the same design, averaged across independent replicates.

Scripts `09`–`12` are earlier iterations of the network-estimation check, kept for reference;
`13` and `14` are the versions used in the current manuscript.

## License

This project is licensed under the GNU General Public License v3.0.

You are free to use, modify, and redistribute the code, provided that
derivative works are also distributed under the same license.
