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
| M0 Pharmacokinetics | Depot, C_plasma, C_periph, C_brain, C_mito | 5 | C_mito → M1 |
| M1 Cardiolipin | CL_n, CL_ox, CL_ext, ALCAT1, TAZ | 5 | CL_ratio → M3 |
| M2 α-Synuclein | aSyn_mono, aSyn_olig, aSyn_ext | 3 | aSyn_olig → M1, M3 |
| M3 ETC Bioenergetics | CI_activity, SC_integrity, Δψ_m, ATP, mPTP_open | 5 | CI, SC → M4 |
| M4 Oxidative Stress | mROS, SOD2_act, GPx4_act | 3 | mROS → M1, M2, M3 |
| M5 Cell Death | PINK1_act, mito_damage, DA_neuron | 3 | DA_neuron → M7 |
| M6 Neuroinflammation | MG_active, TNF, IL1b | 3 | TNF, IL1b → M5 |
| M7 Motor Endpoint | Motor_score | 1 | Clinical endpoint |

## Requirements

- **R** >= 4.3.0
- R packages: `deSolve`, `ggplot2`, `tidyr`, `dplyr`, `gridExtra`, `grid`, `parallel`

Install required packages:
```r
install.packages(c("deSolve", "ggplot2", "tidyr", "dplyr", "gridExtra"))
```

Optional (for manuscript table/text generation):
- Python 3.9+ with `python-docx`

## Repository Structure

```
Bevemipretide QSP/                    # Model source code
├── Integrated Model.R                # Full 28-ODE coupled system + calibration + validation
│                                      # (the single source of truth: every analysis script below
│                                      #  sources ONLY this file)
├── Module 0 PK.R                     # Standalone module-level verification (not part of the
├── Module 1 pass.R                   # analysis pipeline — each file re-implements and
├── Module 2 aSyn.R                   # independently tests one module's ODEs in isolation,
├── Module 3 ETC.R                    # with its own V1-V6 checks, as a development-time
├── Module 4 ROS.R                    # cross-check against the integrated model's behaviour)
├── Module 5 Mitophagy.R
├── Module 6 Neuroinflammation.R
├── Module 7 Motor.R
├── Sensitivity Analysis.R            # OAT + Morris + Sobol (34 parameters) -> Table 3, Fig S1
├── Human Translation.R               # Allometric mouse-to-human scaling -> Table 4, Fig 5
├── run_arm.R                         # Per-arm virtual-population simulation (run 4x, one per dose)
├── Robustness Analysis.R             # Assembles run_arm.R output -> Table 5, Fig 6
├── Uncertainty Analysis.R            # Bootstrap CIs + parameter envelope -> Table S3, Fig 7
├── Phased Simulation.R               # Prodromal -> diagnosis -> treatment -> Table S4
├── Identifiability Analysis.R        # Structural identifiability + collinearity -> Fig S4, S6
├── profile_likelihood_fast.R         # Practical identifiability (profile likelihood) -> Fig S5
├── Generate_New_Plots.R              # Main figure generation (Fig 1, 2, 4, 5, 6, S1-S3, plus
│                                      # first-pass Fig 3 and Fig 7)
├── Fix_Fig3_Fig7.R                   # Corrects Fig 3 and Fig 7 — must be run AFTER
│                                      # Generate_New_Plots.R to produce the final versions
├── export_st1_classification.py      # Exports the parameter classification table to JSON
│                                      # (consumed by Identifiability Analysis.R)
├── ST1_classification.json           # Parameter classification data (output of the script above)
└── Bevemipretide QSP.Rproj           # RStudio project file

New publication plots/                # Generated publication figures and tables
├── Fig1-Fig7 (main figures)
├── FigS1-FigS6 (supplementary figures)
├── Main_Tables.docx
└── Supplementary_Tables.docx
```

## Running the Code

All commands should be run from the `Bevemipretide QSP/` directory. Use `--vanilla` to ensure clean R sessions. Run order matters where noted.

### 1. Integrated Model + Validation

```bash
cd "Bevemipretide QSP"
Rscript --vanilla "Integrated Model.R"
```

Runs the full 28-ODE system with 4-stage calibration and 18 independent validation checks against literature targets (Gao 2017, Choi 2022, Ivanova 2024, Bernheimer 1973). This is the only file every other script depends on.

### 2. (Optional) Standalone Module Verification

```bash
Rscript --vanilla "Module 0 PK.R"
Rscript --vanilla "Module 1 pass.R"
Rscript --vanilla "Module 2 aSyn.R"
Rscript --vanilla "Module 3 ETC.R"
Rscript --vanilla "Module 4 ROS.R"
Rscript --vanilla "Module 5 Mitophagy.R"
Rscript --vanilla "Module 6 Neuroinflammation.R"
Rscript --vanilla "Module 7 Motor.R"
```

Each module independently re-implements and tests its own ODEs in isolation (V1-V6 checks, 48 tests total). These are development-time cross-checks, not inputs to any downstream analysis — safe to skip if you only want to reproduce the manuscript's reported figures and tables.

### 3. Analyses

```bash
Rscript --vanilla "Sensitivity Analysis.R"       # ~20 min (OAT + Morris + Sobol) -> Table 3, Fig S1
Rscript --vanilla "Human Translation.R"          # ~5 min -> Table 4, Fig 5
Rscript --vanilla run_arm.R 1                    # ~15 min each, run once per arm (1-4)
Rscript --vanilla run_arm.R 2
Rscript --vanilla run_arm.R 3
Rscript --vanilla run_arm.R 4
Rscript --vanilla "Robustness Analysis.R"        # assembles the 4 arms above -> Table 5, Fig 6
Rscript --vanilla "Uncertainty Analysis.R"       # ~15 min -> Table S3, Fig 7
Rscript --vanilla "Phased Simulation.R"          # ~10 min -> Table S4
Rscript --vanilla "Identifiability Analysis.R"   # ~10 min -> Fig S4, S6 (reads ST1_classification.json)
Rscript --vanilla "profile_likelihood_fast.R"    # ~35 min -> Fig S5
```

### 4. Generate Publication Figures

```bash
Rscript --vanilla "Generate_New_Plots.R"         # Fig 1, 2, 4, 5, 6, S1-S3, plus first-pass Fig 3, 7
Rscript --vanilla "Fix_Fig3_Fig7.R"              # corrects Fig 3 and Fig 7 — run this AFTER the line above
```

Output appears in `../New publication plots/`.

### 5. (If parameter classifications change) Regenerate ST1_classification.json

```bash
python export_st1_classification.py
```

Only needed if Table ST1's parameter classifications are edited; `ST1_classification.json` is already committed with its current output.

## Key Parameters

### Calibrated (6 parameters, fitted to experimental data)

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

## Parameter Classification (90 total)

| Category | Count | Description |
|----------|-------|-------------|
| Literature-derived | 32 (35.6%) | From published experimental data with DOI/PMID |
| Calibrated | 6 (6.7%) | Fitted to quantitative experimental targets |
| Estimated | 33 (36.7%) | Fitted to qualitative behaviour or steady-state constraints |
| Assumed | 19 (21.1%) | Mechanistically justified with literature support |

## Validation Summary

The model passes 18 independent validation checks across three tiers:

- **Tier 1 (Pathway Fidelity):** 7 directional checks verifying documented causal links
- **Tier 2 (Preclinical Calibration):** 6 quantitative checks against Gao 2017, Choi 2022, Ivanova 2024
- **Tier 3 (Extrapolation):** 5 behavioural checks beyond the calibration window

## Citation

If you use this model, please cite the accompanying manuscript (details to be added upon publication).

## License

This code is provided for peer review and academic use.
