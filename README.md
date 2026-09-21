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
Bevemipretide QSP/              # Model source code
├── Integrated Model.R          # Full 28-ODE coupled system + validation
├── Module 0 PK.R               # Pharmacokinetics (5 ODEs)
├── Module 1 pass.R             # Cardiolipin dynamics (5 ODEs)
├── Module 2 aSyn.R             # α-Synuclein pathology (3 ODEs)
├── Module 3 ETC.R              # ETC bioenergetics (5 ODEs)
├── Module 4 ROS.R              # Oxidative stress (3 ODEs)
├── Module 5 Mitophagy.R        # Cell death & mitophagy (3 ODEs)
├── Module 6 Neuroinflammation.R# Neuroinflammation (3 ODEs)
├── Module 7 Motor.R            # Motor endpoint (1 ODE)
├── Sensitivity Analysis.R      # OAT + Morris + Sobol (34 parameters)
├── Human Translation.R         # Allometric mouse-to-human scaling
├── Robustness Analysis.R       # Virtual population (N=250/arm)
├── run_arm.R                   # Per-arm VPop simulation worker
├── Uncertainty Analysis.R      # Bootstrap CIs + parameter envelope
├── Phased Simulation.R         # Prodromal → diagnosis → treatment
├── Identifiability Analysis.R  # Structural + practical identifiability
├── profile_likelihood_fast.R   # Profile likelihood (Figure S5)
├── Generate_New_Plots.R        # Publication figures (Fig 1-7, S1-S3)
├── Novel Figures.R             # Phase portrait, death decomposition, radar
├── Publication Figures.R       # Earlier figure generation script
├── Fix_Fig3_Fig7.R             # Targeted Fig 3/7 regeneration
├── Model Audit.R               # Bifurcation & solver stability diagnostics
├── build_tables.py             # Generate Main/Supplementary Tables (DOCX)
├── build_methods_results.py    # Generate Methods & Results (DOCX)
├── build_manuscript_new.py     # Manuscript assembly
├── build_narrative_review.py   # Narrative review section
├── ST1_classification.json     # Parameter classification data
└── Bevemipretide QSP.Rproj    # RStudio project file

New publication plots/          # Generated publication figures
├── Fig1-Fig7 (main figures)
├── FigS1-FigS6 (supplementary figures)
├── Main_Tables.docx
└── Supplementary_Tables.docx
```

## Running the Code

All commands should be run from the `Bevemipretide QSP/` directory. Use `--vanilla` to ensure clean R sessions.

### 1. Individual Module Verification (each has 6 built-in tests)

```bash
cd "Bevemipretide QSP"
Rscript --vanilla "Module 0 PK.R"
Rscript --vanilla "Module 1 pass.R"
Rscript --vanilla "Module 2 aSyn.R"
Rscript --vanilla "Module 3 ETC.R"
Rscript --vanilla "Module 4 ROS.R"
Rscript --vanilla "Module 5 Mitophagy.R"
Rscript --vanilla "Module 6 Neuroinflammation.R"
Rscript --vanilla "Module 7 Motor.R"
```

Each module runs V1-V6 verification tests (48 tests total across all modules).

### 2. Integrated Model + Validation

```bash
Rscript --vanilla "Integrated Model.R"
```

Runs the full 28-ODE system with 4-stage calibration verification and 18 independent validation checks against literature targets (Gao 2017, Choi 2022, Ivanova 2024, Bernheimer 1973).

### 3. Analyses

```bash
Rscript --vanilla "Sensitivity Analysis.R"      # ~20 min (OAT + Morris + Sobol)
Rscript --vanilla "Human Translation.R"          # ~5 min
Rscript --vanilla "Robustness Analysis.R"        # ~60 min (N=250/arm, 4 arms)
Rscript --vanilla "Uncertainty Analysis.R"       # ~15 min
Rscript --vanilla "Phased Simulation.R"          # ~10 min
Rscript --vanilla "Identifiability Analysis.R"   # ~10 min
Rscript --vanilla "profile_likelihood_fast.R"    # ~35 min
```

### 4. Generate Publication Figures

```bash
Rscript --vanilla "Generate_New_Plots.R"         # Generates Fig 1-7, S1-S3
```

Output appears in `../New publication plots/`.

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
