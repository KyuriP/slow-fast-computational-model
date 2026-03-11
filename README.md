# Slow–Fast Coupling: Curie–Weiss Symptoms Under Context Dynamics

This repository contains the simulation code, parameter-tuning workflow, and figure-generation scripts for a slow–fast computational model in which a fast Curie–Weiss symptom layer is coupled to a slower contextual process with feedback, diffusion, and optional jump shocks.

## Repository structure

```
slow-fast-linking/
├── README.md
├── .gitignore
├── slow-fast-computational-model.Rproj
│
├── R/
│   ├── scripts/
│   │   ├── 01_create_finalparams.R
│   │   ├── 02_main_sim.R
│   │   ├── 03_create_heatmaps.R
│   │   ├── 04_figure_bistability_geometry.R
│   │   ├── 05_figures_make_all.R
│   │   ├── 06_create_hysteresis_gif.R
│   │   └── 07_inspect_finalparams.R
│   └── utils/
│       ├── utils_fastlayer.R
│       ├── utils_slowfast.R
│       ├── utils_diagnostics.R
│       ├── utils_plotting.R
│       ├── utils_legacy.R
│       └── utils_gridsearch.R
│
├── res/
│   ├── tuned/
│   │   ├── b0_table.rds
│   │   ├── bb_table.rds
│   │   └── final_params.rds
│   ├── main/
│   │   ├── res_A_clean.rds
│   │   ├── res_AS_clean.rds
│   │   ├── res_B0_clean.rds
│   │   └── res_BS_clean.rds
│   ├── heatmaps/
│   │   ├── heat_noshock.rds
│   │   ├── heat_zoom_noshock_hyst.rds
│   │   └── heat_zoom_noshock_hyst2.rds
│   ├── hysteresis/
│   │   └── sim_B2.rds
│   └── archive/
│       └── ... historical intermediate RDS files ...
│
├── figures/
│   └── ... exported figures, gifs, and illustrations ...
└── quarto-docs/
    ├── theoretical_slowfast2.qmd
    ├── theoretical_slowfast3.qmd
    ├── theoretical_slowfast4.qmd
    ├── theoretical_slowfast5_shock.qmd
    ├── timeseparation.qmd
    ├── theoretical_slowfast2.html
    ├── theoretical_slowfast3.html
    ├── theoretical_slowfast4.html
    ├── theoretical_slowfast5_shock.html
    ├── timeseparation.html
    └── working_funcs.R
```

## Main scenario scripts

* `01_create_finalparams.R`
  Recreates tuned parameter objects and writes `final_params.rds`.

* `02_main_sim.R`
  Runs the four main manuscript scenarios:

  * A: βJ < 4, no shocks
  * A+S: βJ < 4, shocks
  * B0: βJ > 4, no shocks
  * B+S: βJ > 4, shocks

* `03_create_heatmaps.R`
  Generates heatmap-based summaries across parameter ranges.

* `04_figure_bistability_geometry.R`
  Produces bistability / geometry figures.

* `05_figures_make_all.R`
  Collects and exports final figures.

* `06_create_hysteresis_gif.R`
  Generates hysteresis animation output.

* `07_inspect_finalparams.R`
  Quick inspection / validation helper for tuned parameters.


## License

This project is licensed under the GNU General Public License v3.0.

You are free to use, modify, and redistribute the code, provided that
derivative works are also distributed under the same license.
