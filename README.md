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

## Requirements

- Python 3.9+
- Packages: `numpy`, `scipy`, `numba`, `matplotlib`, `pandas`

```bash
pip install numpy scipy numba matplotlib pandas
```

## Repository Structure

```
Bevemipretide QSP model analyses/
├── code/
│   ├── model.py              # 28-ODE model: PK (analytical) + 23 PD states (numerical)
│   ├── fast_kernel.py         # Numba-compiled ODE kernel (algebraically identical to model.py)
│   ├── bjp_corrected.py       # All analyses: mouse validation, PK, human translation,
│   │                          # virtual clinical trial, parameter envelope, Sobol SA
│   └── make_bjp_figures.py    # Publication figures (Figs 5-7, S1-S3) from saved results
├── data/
│   └── parameters.csv         # All model parameters (matches manuscript Table ST1)
└── README.md                  # Corrections applied, corrected results, and run instructions
```

## Running the Code

```bash
cd "Bevemipretide QSP model analyses/code"

# Run all analyses (~7 min on 22 cores)
python -X utf8 bjp_corrected.py --workers 22 --patients 250 --candidates 1200 --boots 10000 --sobol-base 256

# Generate publication figures
python -X utf8 make_bjp_figures.py
```

Results are written to `Bevemipretide QSP model analyses/results/` and figures to `Bevemipretide QSP model analyses/figures/`.

### What the analysis produces

| Analysis | Output files | Manuscript tables/figures |
|----------|-------------|--------------------------|
| Mouse calibration and validation | `mouse.json`, `mouse_timeseries.csv` | Table 2 |
| PK characterisation | `pk.json`, `pk_accumulation.csv` | Table 4 |
| Human 20-year trajectories | `human_nominal.json`, `human_20y_trajectories.csv` | Fig 5, Table S4 |
| Phased simulation (prodromal → treatment) | `ST4_phased.csv` | Table S4 |
| Virtual clinical trial (N=250/arm) | `vct_summary.json`, `vct_endpoints.csv` | Table 5, Fig 6 |
| Parameter uncertainty envelope | `envelope.json`, `envelope.csv` | Table S3, Fig 7 |
| Sobol sensitivity analysis (±50%) | `sobol_pm50.json` | Table 3, Fig S1 |
| Mouse dose-response | `mouse_dose_response_day50.csv` | Fig S3 |
| Figures 5, 6, 7, S1, S2, S3 | `figures/*.png`, `figures/*.svg` | Figs 5-7, S1-S3 |

## Key Parameters

### Calibrated (mouse)

| Parameter | Value | Biological Meaning |
|-----------|-------|--------------------|
| CI_max | 0.74 | Maximum Complex I repair capacity |
| α_clear | 0.25 | α-Synuclein clearance capacity |
| k_impair | 0.35 | Microglial impairment of clearance |
| K_death | 0.25 | Death threshold (Hill K_m) |
| k_damage_mPTP | 4.0 | mPTP damage sensitivity |
| n_death | 4 | Death Hill cooperativity |

### Human Translation

| Parameter | Mouse | Human |
|-----------|-------|-------|
| CI_max | 0.74 | 0.90 |
| α_clear | 0.25 | 0.35 |
| death_scale | 1.0 | 0.06 |
| k_e,plasma | 0.231 h⁻¹ | 0.032 h⁻¹ |
| Reference dose | 5 mg/kg IP | 30 mg SC daily |

## Corrections Applied

This code applies 7 corrections to the original R-based analyses (detailed in `Bevemipretide QSP model analyses/README.md`):

1. Drug exposure normalised by exact steady-state mean (not day-120 trough C_ref)
2. k_drug re-anchored to preserve published mouse results (scenario A) with sensitivity check (scenario B)
3. Phased simulation: DA monotonicity enforced (dDA/dt ≤ 0)
4. Virtual trial: symptomatic enrolment (DA < 0.70 screening threshold)
5. Death attribution: cumulative pathway loss (sums to total loss exactly)
6. Response thresholds relabelled (simulation thresholds, not validated MCIDs)
7. Sobol re-run at stated ±50% range with bootstrap intervals

## Validation Summary

The model passes 18 independent validation checks across three tiers:

- **Tier 1 (Pathway Fidelity):** 7 directional checks verifying documented causal links
- **Tier 2 (Preclinical Calibration):** 6 quantitative checks against Gao 2017, Choi 2022, Ivanova 2024
- **Tier 3 (Extrapolation):** 5 behavioural checks beyond the calibration window

## Citation

If you use this model, please cite the accompanying manuscript (details to be added upon publication).

## License

This code is provided for peer review and academic use.
