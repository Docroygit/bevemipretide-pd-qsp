# Bevemipretide (SBT-272) QSP Model for Parkinson's Disease

Quantitative Systems Pharmacology (QSP) model simulating bevemipretide's mechanism of action through cardiolipin stabilisation in Parkinson's disease. Developed at the Department of Pharmacology, AIIMS Bhubaneswar.

## Model Overview

A 28-ODE, 8-module system capturing the self-amplifying vicious cycle linking α-synuclein aggregation, cardiolipin oxidation, Complex I impairment, and mitochondrial ROS overproduction:

```
M0 (PK) → M1 (Cardiolipin) → M3 (ETC Bioenergetics) → M4 (Oxidative Stress)
                ↑                                              ↓
           M2 (α-Synuclein) ←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←←
                ↓                    ↓
           M6 (Neuroinflammation)  M5 (Cell Death) → M7 (Motor Endpoint)
```

| Module | States | ODEs | Key Outputs |
|--------|--------|------|-------------|
| M0 Pharmacokinetics | Depot, A_plasma, A_periph, A_brain, A_mito | 5 | Drug exposure → M1 |
| M1 Cardiolipin | CL_n, CL_ox, CL_ext, ALCAT1, TAZ | 5 | CL_ratio → M3 |
| M2 α-Synuclein | aSyn_mono, aSyn_olig, aSyn_ext | 3 | aSyn_olig → M1, M3 |
| M3 ETC Bioenergetics | CI_activity, SC_integrity, Δψ_m, ATP, mPTP_open | 5 | CI, SC → M4 |
| M4 Oxidative Stress | mROS, SOD2_act, GPx4_act | 3 | mROS → M1, M2, M3 |
| M5 Cell Death | PINK1_act, mito_damage, DA_neuron | 3 | DA_neuron → M7 |
| M6 Neuroinflammation | MG_active, TNF, IL1b | 3 | TNF, IL1b → M5 |
| M7 Motor Endpoint | Motor_score | 1 | Clinical endpoint |

## Repository Structure

```
Bevemipretide QSP model analyses/
├── code/
│   ├── model.py              # 28-ODE model: PK (analytical) + 23 PD states (numerical)
│   ├── fast_kernel.py         # Numba-compiled ODE kernel (algebraically identical to model.py)
│   ├── bjp_corrected.py       # All analyses: mouse validation, PK, human translation,
│   │                          # virtual clinical trial, parameter envelope, Sobol SA
│   ├── make_bjp_figures.py    # Publication figures from saved results (Figs 5-7, S1-S3)
│   ├── Integrated Model.R    # Full 28-ODE coupled system (R implementation)
│   └── Generate_New_Plots.R   # Publication figures 1-4 (model architecture, dynamics,
│                              # calibration, phase portrait)
├── data/
│   └── parameters.csv         # All model parameters (matches manuscript Table ST1)
└── README.md
```

## Requirements

**Python** (for analyses and Figs 5-7, S1-S3):
```bash
pip install numpy scipy numba matplotlib pandas
```

**R** (for Figs 1-4):
- Packages: `deSolve`, `ggplot2`, `tidyr`, `dplyr`, `gridExtra`

## Running the Code

```bash
cd "Bevemipretide QSP model analyses/code"

# Python: run all analyses (~7 min on 22 cores)
python -X utf8 bjp_corrected.py --workers 22 --patients 250 --candidates 1200 --boots 10000 --sobol-base 256

# Python: generate publication figures 5-7, S1-S3
python -X utf8 make_bjp_figures.py

# R: generate publication figures 1-4
Rscript --vanilla Generate_New_Plots.R
```

## Citation

If you use this model, please cite the accompanying manuscript (details to be added upon publication).

## License

This code is provided for peer review and academic use.
