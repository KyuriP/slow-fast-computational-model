# Slow–Fast Coupling: A Heterogeneous-Threshold Symptom Network Under Context Dynamics

This repository contains the simulation code, cached results, and figure-generation scripts for
a slow–fast computational model in which a fast binary symptom network (heterogeneous thresholds
`tau_i`, symptom–symptom coupling `omega_ij`, context sensitivity `gamma_i`) is coupled to a
slower Ornstein–Uhlenbeck contextual process `P_t`, with optional feedback from symptoms back to
context and acute shocks.

The manuscript itself ("When Context Looks Like Coupling: A Slow–Fast Theory of Psychological
Dynamics") is written and maintained in Overleaf and is not tracked in this repository. This
repo holds the simulation code, cached results, and the figures generated from them.

**Everything used by the current manuscript lives under `R/revision_2026/`, `res/revision_2026/`,
and `figs/revision_2026/`**, plus a small number of files kept outside those folders because the
manuscript still references them directly (see below). Everything else in the repo (`R/scripts/`,
`R/utils/`, `res/archive/` and related cache folders, loose files in `figs/`/`img/`, and
`quarto-docs/`) is from an earlier, superseded version of the model (a mean-field Curie–Weiss
formulation) and is kept only for historical reference — none of it is reproducible against, or
cited by, the current manuscript.

## Repository structure

```
slow-fast-computational-model/
├── README.md
├── .gitignore
├── LICENSE
├── slow-fast-computational-model.Rproj
│
├── R/
│   ├── revision_2026/                              CURRENT model — everything below is used
│   │   ├── 00_parameters_uncentered01.R             tau_i, omega_ij, gamma_i for all 9 symptoms
│   │   ├── 01_sim_context_baseline.R                Simulation 1: context shifts symptom activation
│   │   ├── 02_sim_stress_recovery.R                 Simulation 2: P_t dynamics + single shock, feedback off
│   │   ├── 03_sim_feedback.R                        Simulation 3 (main text): feedback on vs. off
│   │   ├── 03b_sim_feedback_grid_supplement.R        Supplement: b calibration grid
│   │   ├── 03c_sim_feedback_shock_grid_extended.R    Extends Sim 2/3 shock design across the b grid
│   │   ├── 04_sim_network_estimation.R              Simulation 4: single-run network-estimation pilot
│   │   ├── 04b_sim_network_estimation_replicated.R  Simulation 4: replicated design (locked results)
│   │   ├── 05_supp_regime_history_dependence.R      Supplement: history-dependence under feedback
│   │   ├── 06_supp_window_check_revised.R           Appendix B: observation-window omitted-context check
│   │   ├── utils_uncentered01_model.R               Shared fast-layer update + nodewise-logistic network estimator
│   │   └── figures/                                 One script per manuscript figure (see below)
│   │
│   ├── scripts/    legacy (Curie–Weiss model) — superseded, not used by the current manuscript
│   ├── utils/       legacy — superseded
│   └── 01_simulation1_context_baseline.R  legacy — superseded
│
├── res/
│   ├── revision_2026/     CURRENT cached results, one subfolder per simulation (sim1–sim4, window_check, supp_history)
│   ├── NCT_analysis/       still used — Appendix C's supplementary NCT robustness check (see Figures below)
│   └── archive/, main/, heatmaps/, hysteresis/, network_check/, tuned/    legacy — superseded
│
├── figs/
│   ├── revision_2026/                              CURRENT figures — everything referenced by the manuscript's
│   │                                                 main text (Figures 1–4) and Appendix B
│   ├── nct_structure_weighted.pdf,
│   │   nct_strength_weighted.pdf,
│   │   nct_edge_fp_weighted.pdf                     still used — Appendix C's NCT check
│   ├── slow-fast-model-illustration_v2.png / .af    still used — Figure 1 schematic + editable source
│   └── ... other loose files ...                    legacy — superseded
│
├── img/
│   ├── model-illustration/    editable design-iteration source files for the Figure 1 schematic
│   └── ...                    legacy/exploratory — superseded
│
└── quarto-docs/
    └── ...                    legacy exploratory .qmd notebooks — superseded
```

## Current model pipeline (`R/revision_2026/`)

Run in order; each simulation script writes its cached output to the matching `res/revision_2026/`
subfolder, and each figure script (in `R/revision_2026/figures/`) reads from `res/revision_2026/`
and writes to `figs/revision_2026/`.

* `00_parameters_uncentered01.R` — defines the 9 PHQ-9-style symptoms and their `tau_i`, `omega_ij`,
  `gamma_i` values, shared by every simulation script.
* `01_sim_context_baseline.R` — **Simulation 1**: holds symptom coupling fixed and varies only the
  slow context level `P`, showing that context alone shifts the distribution of active symptoms
  (Figure 2).
* `02_sim_stress_recovery.R` — **Simulation 2**: adds the full `P_t` OU process (mean reversion +
  diffusion) plus a single acute shock, feedback off (`b = 0`).
* `03_sim_feedback.R` — **Simulation 3 (main text)**: same shock/recovery design as Simulation 2,
  comparing feedback off vs. on (`b = 0.5`) to isolate the effect of symptom-to-context feedback
  (Figure 3).
* `03b_sim_feedback_grid_supplement.R` / `03c_sim_feedback_shock_grid_extended.R` — supplementary
  grids used to calibrate `b` and to extend the shock/recovery design across feedback strengths.
* `04_sim_network_estimation.R` / `04b_sim_network_estimation_replicated.R` — **Simulation 4**:
  tests whether pooling cross-sectional data across context levels without conditioning on `P`
  inflates the estimated symptom network relative to the true generative coupling (Figure 4). `04b`
  is the replicated (30-run) design behind the locked quantitative results; `04` is the earlier
  single-run pilot used for the qualitative network-diagram panels.
* `05_supp_regime_history_dependence.R` — supplementary check for tipping-like/history-dependent
  behavior under stronger feedback.
* `06_supp_window_check_revised.R` — **Appendix B**: reruns the observation-window
  omitted-context check under the current model (superseding the old `res/network_check/` version,
  which used the pre-revision N=12 model).
* `utils_uncentered01_model.R` — shared `simulate_fast_sweep()` (fast-layer update) and
  `fit_edges()` (nodewise-logistic network estimator), reused by every simulation script above so
  the fast-layer dynamics and network-estimation procedure are identical across simulations.

## Figures (`figs/revision_2026/` + kept exceptions)

| File | Manuscript figure |
|---|---|
| `figs/slow-fast-model-illustration_v2.png` | Figure 1 — model schematic |
| `Figure2_context_baseline.pdf` | Figure 2 — context shifts symptom activation (Simulation 1) |
| `Figure3_recovery_feedback.pdf` / `Figure3_regimes.pdf` | Figure 3 — stress, recovery, and feedback (Simulations 2–3) |
| `Figure4_network_trio.pdf` + `Figure4_metric_strip.pdf` | Figure 4 — true vs. estimated symptom networks (Simulation 4) |
| `FigureS_history_dependence.pdf` | Supplement — history-dependence check |
| `Figure_S_window_check_revised.pdf` | Appendix B — observation-window check |
| `figs/nct_structure_weighted.pdf`, `nct_strength_weighted.pdf`, `nct_edge_fp_weighted.pdf` | Appendix C — supplementary NCT robustness check |

Each `.pdf` has a matching `.png` for quick preview. The NCT check (Appendix C) and the Figure 1
schematic source files live outside `figs/revision_2026/` for historical reasons but are still
directly referenced by the manuscript, so they're kept alongside it rather than archived.

## Legacy content

`R/scripts/`, `R/utils/`, `res/archive/`, `res/main/`, `res/heatmaps/`, `res/hysteresis/`,
`res/network_check/`, `res/tuned/`, most of `figs/` and `img/`, and `quarto-docs/` all belong to
an earlier mean-field Curie–Weiss version of the model (different parameterization, N=12 symptoms,
`beta`/`J` notation) that predates the current revision. They're kept for provenance but are not
reproducible against, and are not cited by, the current manuscript.

## License

This project is licensed under the GNU General Public License v3.0.

You are free to use, modify, and redistribute the code, provided that
derivative works are also distributed under the same license.
