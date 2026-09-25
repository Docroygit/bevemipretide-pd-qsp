# Bevemipretide QSP Model Analyses

Independent Python rebuild of the 28-ODE QSP model and R code for publication figures.

## Code

| File | Language | Purpose |
|------|----------|---------|
| `code/model.py` | Python | 28-ODE model: 5 PK states (analytical) + 23 PD states (numerical) |
| `code/fast_kernel.py` | Python | Numba-compiled ODE kernel (algebraically identical to model.py) |
| `code/bjp_corrected.py` | Python | All analyses: mouse validation, PK, human translation, virtual clinical trial, parameter envelope, Sobol SA |
| `code/make_bjp_figures.py` | Python | Publication figures 4-7 and S1-S3 from saved results |
| `code/Integrated Model.R` | R | Full 28-ODE coupled system with calibration, validation, and `run_integrated()` function |
| `code/Generate_New_Plots.R` | R | Publication figures 1-3 (model architecture, integrated dynamics, calibration/validation) |
| `data/parameters.csv` | Data | All model parameters (matches manuscript Table ST1) |

## Running

### Python analyses and figures 4-7, S1-S3

```bash
cd code
python -X utf8 bjp_corrected.py --workers 22 --patients 250 --candidates 1200 --boots 10000 --sobol-base 256
python -X utf8 make_bjp_figures.py
```

Requires: Python 3.9+, numpy, scipy, numba, matplotlib, pandas.

Results are written to `results/` (JSON, CSV) and figures to `figures/` (PNG, SVG).

### R figures 1-3

```bash
cd code
Rscript --vanilla Generate_New_Plots.R
```

Requires: R with deSolve, ggplot2, tidyr, dplyr, gridExtra.

Figures are written to `../New publication plots/`.

## What each analysis produces

| Function in `bjp_corrected.py` | Output files | Manuscript element |
|---|---|---|
| `mouse()` | `mouse.json`, `mouse_timeseries.csv` | Table 2 |
| `pk_table()` | `pk.json`, `pk_accumulation.csv` | Table 4 |
| `human_nominal()` | `human_nominal.json`, `human_20y_trajectories.csv` | Fig 5, Table S4 |
| `trial()` | `vct_summary.json`, `vct_endpoints.csv` | Table 5, Fig 6 |
| `envelope()` | `envelope.json`, `envelope.csv` | Table S3, Fig 7 |
| `sobol()` | `sobol_pm50.json` | Table 3, Fig S1 |

| R script | Figures |
|---|---|
| `Generate_New_Plots.R` | Fig 1 (model architecture), Fig 2 (integrated dynamics), Fig 3 (calibration/validation) |
