# Corrected analyses for the BJP bevemipretide QSP manuscript

Regenerated 23 September 2026 with the independent Python rebuild of the 28-ODE model.
Nothing outside this folder was modified.

## Provenance

`code/model.py`, `code/fast_kernel.py`, `code/test_model.py` and `data/parameters.csv` are unmodified copies
from `G:\My Drive\Bevemipretide QSP\Revisions_2026-09-18`:

| File | SHA-256 |
|---|---|
| code/model.py | a3e6f288659ae38fee9edab6d25aac36c74412d7c187c5875d6a43459dc992a2 |
| code/fast_kernel.py | 98649944acec7be224c1550e248b02325811634f1dcdb22a64ab3150a16125fe |
| code/test_model.py | 6100a1187ce56f42eb8ce65ae51522772c1cd11fc6f5ae05b1b2c22871827fc7 |
| data/parameters.csv | aadbd18dbf7b68ebd45f11caad461d47e334ba8eed50611c71736bf8e555cb1b |

The rebuild's equations match manuscript Table ST5 term by term, and its parameters match Table ST1.
On this machine it reproduces its own saved outputs to about 1e-10 and passes its 11 tests.
New code: `code/bjp_corrected.py` (analyses), `code/make_bjp_figures.py` (figures), `code/test_bjp_corrected.py` (4 tests).

## Run

```
cd code
python -X utf8 -m pytest test_model.py test_bjp_corrected.py -q --junitxml=../results/verification.xml
python -X utf8 bjp_corrected.py --workers 22 --patients 250 --candidates 1200 --boots 10000 --sobol-base 256
python -X utf8 make_bjp_figures.py
```

Tested with Python 3.9.13, NumPy 1.26.4, SciPy 1.11.4 and numba 0.60.0. The full run took 7.3 minutes on 22 processes.
Run the tests before and after the analysis: the correction tests read the saved results.

## What was corrected

| # | Problem in the BJP package | Correction |
|---|---|---|
| 1 | C_ref = 4.43 (mouse) and 2.15 (human) are described as 120-day steady-state values. With the stated PK rates they are not. The slowest disposition half-life is 56 days (mouse) and 112 days (human). At day 120 the mitochondrial trough is 77% and 52% of steady state. 4.43 is only 39% of the true mouse steady-state mean. | Exposure = mitochondrial amount / exact steady-state mean at the reference dose of each species (5 mg/kg; 30 mg). Doses scale linearly: 10/30/60 mg are 1/3, 1 and 2 reference doses. |
| 2 | k_drug = 0.1 was tuned together with C_ref = 4.43. Correcting C_ref alone would change the published mouse drug response. | **Scenario A (primary):** k_drug = 0.1 × 11.273/4.43 = 0.2545. This reproduces every published mouse drug result exactly (max difference 4e-9) and carries the mouse steady-state effect to humans at 30 mg. **Scenario B (sensitivity):** literal k_drug = 0.1. |
| 3 | ST4 reports survival *rising* under treatment (0.700 → 0.744). The model forbids this, because dDA/dt ≤ 0. The "without treatment" column (0.656) is the value 12 months later, not 18. | Phased simulation recomputed. Treated ≤ start and treated ≥ placebo are both enforced by test. |
| 4 | The virtual trial is described as symptomatic, but in Fig 6A the placebo median is about 0.85 after 18 months. 250 were generated but 244 analysed, without explanation. | Candidates are screened until untreated survival falls below 0.70 (30-year horizon), then 12 months untreated, then 18 months in 4 matched arms. All 271 screened candidates are recorded; 250 were enrolled and 250 analysed. |
| 5 | Death attribution used instantaneous rate shares and a fixed "37% unreachable" ceiling. | Cumulative loss per pathway, which adds up exactly to total loss. The rescue fraction is the reduction in cumulative pathway loss relative to placebo. |
| 6 | "MCID" responder thresholds. | Same cut-offs (0.05 DA, 0.03 motor, absolute), relabelled as simulation response thresholds. |
| 7 | Sobol ranges not stated. | Re-run at the manuscript's ±50% OAT range, N = 256, with bootstrap intervals, on the 13 distinct parameters that appear in Table 3's Sobol columns. |

The VPop distributions follow Table ST2 exactly: 13 log-normal parameters, plus truncated-normal CI_max and α_clear. The exposure reference is fixed at the nominal human value, so between-patient PK variability is kept.

## Corrected results

### Unchanged by the corrections (scenario A mouse = published model)

| Quantity | BJP | Corrected |
|---|---|---|
| Complex I activity, day 35 | 0.459 | 0.459 absolute; **0.484 relative to healthy** (the value to compare with "49% of WT") |
| ROS / healthy | 161.6% | 161.6% |
| Native CL depletion | 15.7% | 15.7% |
| DA loss, day 35 | 36.0% | 36.0% |
| Motor score, day 35 | 15.6% | 15.7% |
| TLR2 KO (k_act_MG = 0) α-syn reduction | 20.3% | 20.3%, below the 23–45% literature range |
| **Healthy DA, day 50** | **0.994 (PASS > 0.99)** | **0.9886, which fails > 0.99**. Use 0.989 (the text value) and restate the criterion |
| Drug rescue of DA gap, 5 mg/kg day 35 | 26.3% | 26.0% (A) / 15.1% (B) |
| CL gap restored | 94.6% | 94.6% (A) / 60.2% (B) |
| α-syn oligomer reduction | 33.5% | 33.5% (A) / 22.8% (B) |
| Motor improvement | 26.6% | 26.6% (A) / 17.5% (B) |
| Dose ordering 0.5–5 mg/kg | monotonic | monotonic (A and B) |
| k_drug = 0 ablation | no rescue | ΔDA = 1e-10 |
| 1-year DA untreated / treated | 0.006 / 0.079 | 0.006 / 0.079 (A); 0.046 (B) |
| DA checkpoints 30/90/180/365 d | 0.687/0.292/0.081/0.006 | 0.687/0.291/0.080/0.006 |
| "6.7-fold amplification" | loop gain | 6.66 is simply disease/healthy oligomer at day 35, not a gain. Remove the claim |
| Death shares, day 35 | 63/12/25% | 63/12/25% cumulative (65/12/23% instantaneous) |
| Human prodrome / clinical PD / Bernheimer | 5.6 / 10.8 / 14.3 y | 5.55 / 10.78 / 14.26 y |
| Healthy human DA at 20 y | 90.4% | 90.4% |

All mouse states stayed non-negative, fractions stayed ≤ 1 and DA was monotone across healthy, disease, 1 and 10 mg/kg scenarios.
Pool indices peaked at 0.956, so the "all 28 states in [0,1]" statement holds for the PD states but not for the PK amounts.

### Table 4 (PK), corrected columns

| Quantity | Mouse | Human |
|---|---|---|
| ln2 / k_e (central elimination half-time) | 3.0 h | 21.7 h |
| Slowest disposition-mode half-life | 56.5 d | 112.4 d |
| Days to 50 / 90 / 95% of steady-state trough | 57 / 188 / 245 | 115 / 376 / 488 |
| Trough at day 120, % of steady state | 77.0% | 51.7% |
| Steady-state mean mitochondrial amount per daily reference dose | 11.27 | 48.83 |
| BJP C_ref | 4.43 | 2.15 |

The stated BW^−0.25 rule gives human k_brain_out = 0.0022 h⁻¹, not 0.010. Either change the value or relabel it as assumed.
The FDA K_m HED of 5 mg/kg for a 70 kg adult is 28.4 mg.

### Human course, 30 mg from year 0 (Fig 5B; BJP: onset delay 4.1 y, clinical PD delay 8.1 y)

| | Scenario A | Scenario B |
|---|---|---|
| Delay in crossing 70% (symptoms) | 4.5 y | 3.0 y |
| Delay in crossing 50% (clinical PD) | 8.9 y | 5.9 y |
| DA benefit at year 1 / 5 / 10 / 20 | +2.3 / +11.1 / +17.6 / +21.9 points | +1.7 / +8.4 / +13.3 / +16.1 points |

BJP's "+73.8% at year 20" is not an absolute gain. In relative terms, scenario A gives +79% and B gives +58%.
BJP's own human onset delays correspond to an effective k_drug of about 0.19–0.20, which lies between scenarios A and B.

### ST4, phased simulation (scenario A; B in `results/ST4_phased.csv`)

| Delay | DA at start | Placebo 18 mo | 10 mg ΔDA | 30 mg ΔDA | 60 mg ΔDA | 30 mg relative to 0-mo |
|---|---|---|---|---|---|---|
| 0 mo | 0.700 | 0.636 | +0.017 | +0.025 | +0.029 | reference |
| 6 mo | 0.678 | 0.616 | +0.017 | +0.025 | +0.028 | −3% |
| 12 mo | 0.656 | 0.596 | +0.016 | +0.024 | +0.027 | −6% |
| 24 mo | 0.616 | 0.559 | +0.015 | +0.022 | +0.025 | −12% |

BJP reported 30 mg ΔDA of +0.088. The corrected value is +0.025 (A) / +0.019 (B). The −12% relative loss at 24 months is unchanged.

### Table 5, virtual trial at 18 months (n = 250 per arm, matched; 95% bootstrap intervals, 10,000 resamples)

Enrolment flow: 271 screened, 21 never reached 70% survival within 30 years, 250 enrolled and 250 analysed.
Survival at diagnosis was 0.27–0.69 (median 0.64). Median time to symptom onset was 3.7 y (IQR 1.0–8.1). Placebo survival at 18 months had a median of 0.55 (range 0.06–0.68).

**Scenario A (primary)**

| Dose | Endpoint | Mean Δ vs placebo | d (unpaired) | Cliff's δ | Response rate | N/arm, 80% power |
|---|---|---|---|---|---|---|
| 10 mg | DA | +0.014 | 0.075 [0.062, 0.090] | 0.059 | 4.8% | 2,810 |
| 30 mg | DA | +0.021 | **0.117 [0.099, 0.137]** | 0.090 | 9.6% | **1,153** |
| 60 mg | DA | +0.027 | 0.148 [0.127, 0.172] | 0.109 | 14.0% | 720 |
| 10 mg | Motor | −0.018 | 0.086 [0.072, 0.103] | 0.081 | 15.6% | 2,106 |
| 30 mg | Motor | −0.028 | **0.136 [0.116, 0.159]** | 0.124 | 28.4% | **855** |
| 60 mg | Motor | −0.036 | 0.175 [0.151, 0.202] | 0.152 | 38.8% | 515 |

**Scenario B (sensitivity):** 30 mg d = 0.080 (DA) and 0.093 (motor); 60 mg d = 0.107 and 0.124. Full values are in `results/vct_summary.json`.

BJP for comparison: 30 mg d = 0.221 [0.177, 0.273], 25.8% responders, 322 per arm. 60 mg d = 0.273.
In the corrected trial the top quartile of patients accounts for about 60% of the total benefit.

The bootstrap intervals resample virtual patients from one assumed population. They describe simulation precision, not uncertainty about the real effect.

### Fig 6C, reduction in cumulative neuron loss vs placebo (scenario A)

| Pathway | 10 mg | 30 mg | 60 mg |
|---|---|---|---|
| Mitochondrial | 9% | 14% | 18% |
| Neuroinflammation | 12% | 20% | 26% |
| Proteotoxicity | 19% | 30% | 36% |
| Total | 12% | 19% | 24% |

In symptomatic humans the drug reduces the *non-mitochondrial* pathways proportionally **more** than the mitochondrial one. This reverses BJP Fig 6C (mitochondrial 60% vs 35–38% at 60 mg).
The likely cause is that mitochondrial damage in advanced disease sits on the saturated part of its Hill curve, while oligomer and cytokine hazards sit on the steep part.
The "37% of death is out of reach" monotherapy ceiling is therefore not supported. Placebo cumulative shares were 57% mitochondrial, 17% inflammatory and 26% proteotoxic.

### ST3B, translational parameter grid (125/125 combinations evaluated)

| | Scenario A | Scenario B |
|---|---|---|
| 30 mg DA saved, range (2.5–97.5%) | 0.009–0.051 (0.011–0.045) | 0.007–0.037 (0.008–0.033) |
| Approximate d range (DA saved / pooled SD) | 0.05–0.28 | 0.04–0.21 |
| Fold range of mean DA saved over α_clear / death scale / CI_max | 2.4 / 1.8 / 1.3 | 2.4 / 1.8 / 1.3 |

α_clear remains the largest source of uncertainty (BJP: 2.6-fold), so the GBA-stratification rationale survives as a hypothesis.

### Table 3, Sobol total-order indices at ±50% (N = 256, 3,840 runs; `results/sobol_pm50.json`)

| Rank | DA survival, day 35: S_T (S1) | Motor score, day 35: S_T (S1) |
|---|---|---|
| 1 | k_SOD2 0.32 (0.16) | K_motor 0.44 (0.24) |
| 2 | k_ROS_basal 0.29 (0.17) | k_SOD2 0.25 (0.09) |
| 3 | K_mPTP 0.23 (0.06) | k_ROS_basal 0.23 (0.04) |
| 4 | k_agg 0.20 (0.11) | K_mPTP 0.19 (≈0) |
| 5 | CI_max 0.18 (0.07) | k_agg 0.13 (0.08) |

Compared with BJP:
- ROS-module parameters still lead DA survival, but k_ROS_basal is not uniquely dominant; its interval overlaps k_SOD2's.
- K_mPTP is largely interactive (76% of its S_T for DA, S1 ≈ 0 for motor), which supports the "threshold switch" reading.
- The motor endpoint is led by the motor-mapping threshold K_motor. "ROS parameters dominate both endpoints" should be removed.
- BJP Table 3 lists 13 different parameters across its "top 10" Sobol columns; state which parameters were actually varied.

## Manuscript text that must change

- **Abstract / 3.7 / Discussion:** replace d = 0.221, 25.8% responders and 322 per arm with the Table 5 values above. Do not place the d value alongside ADAGIO or SPARK. At d ≈ 0.12, the corrected projection sits at or below the effect sizes of trials that failed.
- **2.4 / Table 4:** describe C_ref as a steady-state normalising constant. Report modal half-lives and time to steady state. Fix or relabel k_brain_out.
- **2.5:** describe the enrolment rule, the screening flow (271 → 250) and the matched-arm design. Replace "MCID" with "simulation response threshold".
- **2.7 / 3.x / Discussion:** a rank-7 sensitivity matrix over 39 parameters leaves at least 32 unconstrained directions, and Figure S5 shows flat profiles. Say that parameters are not practically identifiable from the calibration targets, and drop "well-determined".
- **3.2 Tier 3:** remove "6.7-fold amplification from coupling".
- **3.3 / Figure 4:** call these trajectory endpoints, not attractors or a saddle.
- **3.4 / Discussion / Conclusion:** remove the 63/37 "ceiling" and the 52% cap. Use the cumulative rescue table.
- **Discussion, long-term:** year-1 +2.3 points to year-20 +21.9 points (A), and state whether values are absolute or relative.
- **Table 2:** healthy DA 0.989, and restate the criterion. Report CI activity relative to healthy (0.484). Note that the TLR2 result is below the literature range.
- **Model_info risk assessment:** lower "model influence" to trial-planning scenarios. The corrected sample size (about 850–1,150 per arm at 30 mg) is itself a planning finding worth stating.

## Corrected submission documents

`documents/` holds copies of the BJP files with the corrections as Word tracked changes (author "Claude (corrected analyses)"). The originals in the project folder are untouched.

| File | Content |
|---|---|
| Manuscript_bjp_tracked.docx | 576 tracked revisions: abstract, Methods 2.1/2.3–2.9, Results 3.1–3.9, Discussion, Limitations, Conclusion, legends of Figs 4–7, Tables 2–5, plus four new table notes and two new Table 4 rows |
| Supplementary_Tables_tracked.docx | 206 tracked revisions in ST1 (C_ref, k_drug, K_CL, six unit labels), ST2 enrolment note, ST3A/B, ST4B/C, and the ST4/ST5 notes |
| *_clean.docx | The same files with all changes accepted (manuscript 40 pages, supplement 12 pages) |
| ../qa/*.pdf | Word renderings of the tracked (with markup) and clean versions |

Checks run:
- Both files pass the OOXML schema validator.
- The redlining validator confirms every change against the original is wrapped in a tracked change.
- Every revised paragraph's visible text was asserted equal to the intended text.

To regenerate: `python -X utf8 code/apply_bjp_edits.py`. It reads the original .docx files and ../results and writes documents/.

Rendering note: Word COM automation worked when the commands ran inline. Running `code/render_documents.ps1` directly from the Google Drive volume hung before Word started, so copy the script to a local folder first.
If a hidden Word instance is killed, its entries under `HKCU:\Software\Microsoft\Office\16.0\Word\Resiliency\DocumentRecovery` must be cleared. Otherwise the next hidden instance waits on an invisible recovery prompt.

Corrected figures under BJP file names are in `figures_bjp/`, with legends in `figures_bjp/LEGENDS.md`.

Not edited:
- `Model_info_bjp.docx`: the risk assessment still describes Medium–High influence on dose, sample size and go/no-go.
- `Bullet_Point_Summary_bjp.docx`: still claims benefit on the other pathways; this now holds, but the size needs stating.
- The Figure 4 image: its "attractor" labels are in the picture.
- `Supplementary_Figures.pdf`.

The GitHub repository cited in Section 2.9 should also receive this Python code, because the manuscript now says the corrected analyses used it.

## Limits of this regeneration

- Scenario A assumes 30 mg in humans reproduces the mouse 5 mg/kg steady-state exposure at the mitochondrion. This is the manuscript's HED assumption, now applied consistently; it is not measured.
- Scenario B shows how much rides on that assumption.
- Mouse calibration, the human recalibration (CI_max 0.90, α_clear 0.35, death × 0.06) and the ST2 distributions are inherited unchanged. They are not refitted.
- The original R code was not available, so BJP numbers not reproduced here cannot be traced to specific code lines. This applies to the Table 5 values, the ST4 treated column, the 0.994 baseline and the 244/250 attrition.
- Morris screening, the OAT tornado and the time-course identifiability rank of 11 were not re-run.
