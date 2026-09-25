# GitHub R implementation vs Python implementation

Compared 24 September 2026:
- **R:** https://github.com/Docroygit/bevemipretide-pd-qsp at commit 2e259db (21 Sep 2026).
- **Python:** the rebuild in `code/`, with the corrections in `bjp_corrected.py`.

The R scripts were run on this machine with R 4.4.3 and deSolve:
- `Integrated Model.R`
- `Phased Simulation.R`
- the placebo arm of `run_arm.R`
- one timed patient per arm

The three R drug arms were stopped after 30+ minutes because they were projected to take hours.
Their design was replicated exactly with the Python engine instead (`code/r_design_replication.py`, `results/r_design_replication.json`).

## 1. Same model

The ODE right-hand side in `Integrated Model.R` (lines 200–356) matches `code/model.py` term for term, and every parameter value matches.

Running both codes gives identical results:

| Output | R (this machine) | Python |
|---|---|---|
| Healthy DA, day 50 | 0.9887 | 0.9886 |
| Disease CI / ROS ratio / CL drop, day 35 | 0.4592 / 161.64% / 15.75% | 0.4592 / 161.64% / 15.75% |
| Disease DA loss, day 35 | 36.02% | 36.04% |
| TLR2 KO oligomer reduction | 20.3% | 20.3% |
| DA gap rescued, 5 mg/kg (R C_ref 4.43) | 26.3% | 26.0% |
| Human onset / clinical PD / Bernheimer | 5.6 / 10.8 / 14.3 y | 5.55 / 10.78 / 14.26 y |
| DA at diagnosis (onset + 12 mo) | 0.6564 | 0.6564 |
| Human C_ref (R definition) | 2.151082 | 2.151082 |

The small rescue difference (26.3% vs 26.0%) comes from starting states:
- R starts from a hard-coded vector that is not an equilibrium of the integrated model (e.g. mPTP_open 0.002 vs 0.019 at equilibrium).
- Python equilibrates the healthy biochemistry first.

## 2. Implementation differences

| Aspect | R | Python | Which is sounder |
|---|---|---|---|
| Starting state | Fixed vector from standalone modules; not an equilibrium, so an early transient runs | 5,000 h burn-in to equilibrium | Python |
| State bounds | Every state clamped to [0,1] inside the ODE function, and output clamped again. This makes the right-hand side non-smooth and can hide negative or runaway states | No clamping; negative or out-of-range values fail tests | Python |
| PK | Numerical, with 547 dose events (a 30 mg patient takes ~23 s vs 0.5 s for placebo) | Exact analytic solution | Python (same answer, verifiable) |
| Mouse C_ref | 4.43, described as "trough at 5 mg/kg SS". The stated PK gives 11.25 at steady state; 4.43 is reached around day 40 | Exact steady-state mean (11.27); k_drug re-anchored (A) or kept (B) | Python |
| Human C_ref | Day-120 trough of 30 mg (2.151), only 52% of steady state, so exposure keeps rising to about 1.9 × C_ref. Effective human potency at steady state 0.195 vs mouse 0.254 | Exact steady-state mean; human 30 mg = mouse 5 mg/kg (A), or literal k_drug (B) | Python: explicit, and consistent across species |
| Tests | Print-based PASS checks with loosened windows (below) | 15 automated numerical tests | Python |

## 3. Problems in the R analyses themselves

1. **The virtual "Phase 2 trial" is a prevention simulation.**
   - `run_arm.R` starts every patient at `integrated_init()`, with DA = 1 and no disease yet, and doses from day 0 for 547 days.
   - No patient has symptoms at entry. The nominal patient's placebo DA at 18 months is 0.909.
   - The manuscript describes symptomatic enrolment and cites Phase 2 trial design.
2. **The sample depends on the machine.**
   - Each patient has a 60-second wall-clock timeout (`PATIENT_TIMEOUT <- 60`).
   - Placebo patients take about 0.5 s but drug patients about 23 s on an idle machine, so only drug-arm patients time out, and more of them on slower or busier hardware.
   - Patients who time out in any arm are then dropped from all arms. This is the source of 250 → 244.
3. **The published ST4 table was not produced by this code.**
   - `Phased Simulation.R` gives a 30 mg benefit of +0.0221 from diagnosis, with treated DA of 0.618, below the starting 0.656.
   - The manuscript reports +0.088 with DA rising to 0.744.
   - The code measures delay from diagnosis; the table labels it as delay from symptoms.
   - Only the −12% relative loss (87.9% retained) matches.
4. **Calibration is narrower than reported.**
   - Only CI_max, α_clear and k_impair are searched. K_death, k_damage_mPTP and n_death are returned as constants (`return(list(K_death = 0.25, ...))`), yet they are listed as calibrated.
   - Stage 4 accepts a TLR2 reduction of ≥ 20% rather than 23–45%.
5. **Validation windows are looser than the stated targets.**
   - The healthy check prints "NEEDS ADJUSTMENT" (0.9887), but the paper reports 0.994 PASS.
   - Validation V1 passes DA in the range [0.97, 1.01].
   - The TLR2 checks pass anything in 15–50% while stating a 23–45% target.
   - The code counts 24/24 and 19/19 checks; the paper reports 18/18.
6. **"6.7-fold amplification vs ~1.5-fold standalone" (T3.5).**
   - 6.7 is disease ÷ healthy oligomer at day 35.
   - "~1.5x" is a literal string in the print statement, never computed.
   - T3.4 compares one-year disease and drug values with the day-50 healthy value (0.989; the healthy value at one year is 0.920).
7. **The profile likelihood is not a likelihood.**
   - Its "data" are the model's own nominal outputs, so the objective is zero by construction.
   - Its σ values are half-widths of the acceptance ranges (e.g. 50 percentage points for ROS). Flat profiles below χ² = 3.84 are nearly guaranteed.
   - With 7 outputs and 39 parameters, the FIM rank of 7 is the maximum possible, not a demonstration of identifiability.
8. **Sobol details.**
   - Plain Monte Carlo sampling.
   - S1 and S_T are clipped at zero, which inflates the reported "interaction" share.
   - CI_max is sampled up to 1.11, which is non-physical and hidden by the clamp.
   - A 13-parameter ±50% re-run in Python ranks k_SOD2 (0.32) ≈ k_ROS_basal (0.29) above K_mPTP (0.23), not k_ROS_basal 0.50 and K_mPTP 0.33.

## 4. What drives the difference in clinical results

All rows use 250 virtual patients sampled from the R distributions:

| Design | Exposure normalisation | 30 mg d (DA) | 30 mg d (motor) | N/arm for 80% |
|---|---|---|---|---|
| Manuscript (R, as published) | R C_ref | 0.221 | 0.221 | 322 |
| R design replicated in Python (healthy start, dose from day 0) | R C_ref | 0.191 | 0.216 | 432 |
| Symptomatic enrolment + 12-month delay | R C_ref | 0.104 | 0.122 | ~1,440 |
| Symptomatic enrolment (corrected, scenario A) | Steady state, mouse-anchored | 0.117 | 0.136 | 1,153 |
| Symptomatic enrolment (scenario B) | Steady state, literal k_drug | 0.080 | 0.093 | 2,433 |

The replicated R design falls within the published bootstrap interval (0.177–0.273); it differs from 0.221 only because the random draws differ.
Moving to symptomatic enrolment halves the effect whichever normalisation is used.

## Verdict

- **Same model.** Both implementations encode the same biology and give the same numbers when set up the same way.
- **The Python implementation is the scientifically sounder basis for the paper's claims**:
  - equilibrated starting states;
  - no clamping;
  - exact, species-consistent exposure;
  - a symptomatic trial population with complete accounting;
  - honest identifiability statements;
  - automated tests.
- **The R code is a faithful source for** the mouse results and the untreated human milestones.
- **The R code does not support** the manuscript's virtual-trial claims (effect size, responder rates, sample size), the ST4 table, or the identifiability and validation claims as worded.
- **The corrected manuscript already reflects this.** The authors' own R phased simulation independently reproduces the corrected ST4 benefit (+0.022 at 30 mg from diagnosis).
