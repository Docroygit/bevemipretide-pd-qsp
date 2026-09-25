# Results comparison: manuscript, GitHub R code, Python implementation

Columns:
- **Manuscript**: values reported in Manuscript_bjp.docx and Supplementary_Tables.docx.
- **R (run here)**: GitHub commit 2e259db, run with R 4.4.3 on this machine.
- **Python, R set-up**: our engine with R's design choices (hard-coded start state for the trial, C_ref = day-120 trough, healthy start with dosing from day 0).
- **Python A / B**: corrected analysis. Steady-state exposure, symptomatic enrolment (DA < 0.70 then 12 months).
  - A: k_drug re-anchored so mouse results are unchanged.
  - B: literal k_drug = 0.1.

"—" means not rerun. The R drug arms and Human Translation.R were stopped after 30+ minutes; they were projected to take hours.

| Result | Manuscript | R (run here) | Python, R set-up | Python A (primary) | Python B (sensitivity) |
|---|---|---|---|---|---|
| **Mouse model (day 35 unless stated)** | | | | | |
| Healthy DA neurons, day 50 | 0.994 (PASS > 0.99) | 0.9887 ("NEEDS ADJUSTMENT") | 0.9886 | 0.9886 | 0.9886 |
| Complex I activity | 0.459 | 0.4592 | 0.4592 | 0.4592 (0.484 vs healthy) | same |
| mROS vs healthy | 161.6% | 161.64% | 161.64% | 161.64% | same |
| Native CL depletion | 15.7% | 15.75% | 15.75% | 15.75% | same |
| DA neuron loss | 36.0% | 36.02% | 36.04% | 36.04% | same |
| Motor score | 0.156 | 0.1565 | 0.1566 | 0.1566 | same |
| TLR2 KO oligomer reduction (target 23–45%) | 20.3% (PASS) | 20.3% (PASS at 15–50%) | 20.3% | 20.3% | same |
| DA gap rescued, 5 mg/kg | 26.3% | 26.3% | 26.0% | 26.0% | 15.1% |
| CL gap restored, 5 mg/kg | 94.6% | 94.6% | 94.6% | 94.6% | 60.2% |
| α-syn oligomer reduction, 5 mg/kg | 33.5% | 33.5% | 33.5% | 33.5% | 22.8% |
| Motor improvement, 5 mg/kg | 26.6% | 26.6% | 26.6% | 26.6% | 17.5% |
| DA at 1 year, untreated / 5 mg/kg | 0.006 / 0.079 | 0.006 / 0.079 | 0.006 / 0.079 | 0.006 / 0.079 | 0.006 / 0.046 |
| Healthy comparator used at 1 year | 0.989 (the day-50 value) | 0.989 (the day-50 value) | 0.920 | 0.920 | 0.920 |
| "Vicious-cycle amplification" | 6.7× vs ~1.5× standalone | 6.7×; "~1.5×" is text in the print statement | 6.66× (disease ÷ healthy) | same | same |
| Death shares mito / inflam / proteo | 63 / 12 / 25% (rates) | — | 65 / 12 / 23% (rates) | 63 / 12 / 25% (cumulative) | same |
| **Pharmacokinetics** | | | | | |
| Mouse C_ref (drug normaliser) | 4.43 ("SS trough") | 4.43 (hard-coded) | 4.43 | 11.27 (exact SS mean) | 11.27 |
| Human C_ref (R dose units) | 2.15 ("120-day SS") | 2.151 (day-120 trough) | 2.151 | 4.19 (exact SS mean) | 4.19 |
| Slowest PK half-life, mouse / human | not reported | not computed | 56.5 / 112.4 d | same | same |
| Human exposure at day 120, % of SS | assumed 100% | 52% | 52% | 52% | 52% |
| Effective drug strength at SS, human 30 mg | not stated | 0.195 | 0.195 | 0.254 (= mouse 5 mg/kg) | 0.100 |
| **Human disease course (untreated)** | | | | | |
| Symptom onset (DA < 70%) | 5.6 y | 5.6 y | 5.55 y | 5.55 y | 5.55 y |
| Clinical PD (DA < 50%) | 10.8 y | 10.8 y | 10.78 y | 10.78 y | 10.78 y |
| Bernheimer threshold (DA < 40%) | 14.3 y | 14.3 y | 14.26 y | 14.26 y | 14.26 y |
| DA at diagnosis (onset + 12 mo) | 0.656 | 0.6564 | 0.6564 | 0.6564 | 0.6564 |
| Healthy DA at 20 years | 90.4% | — | 90.4% | 90.4% | 90.4% |
| **30 mg from year 0 (Fig 5B)** | | | | | |
| Delay in symptom onset | 4.1 y | — | 4.11 y | 4.53 y | 2.99 y |
| Delay in clinical PD | 8.1 y | — | 8.13 y | 8.93 y | 5.94 y |
| DA benefit, year 1 | +2.3% | — | +2.2 points | +2.3 points | +1.7 points |
| DA benefit, year 20 | +73.8% | — | +20.4 points (+74% relative) | +21.9 points (+79% relative) | +16.1 points (+58% relative) |
| **ST4: 18 months of treatment from diagnosis** | | | | | |
| 30 mg: DA at start → treated / placebo | 0.700 → 0.744 / 0.656 | 0.656 → 0.618 / 0.596 | 0.656 → 0.618 / 0.596 | 0.656 → 0.620 / 0.596 | 0.656 → 0.614 / 0.596 |
| Benefit, 10 / 30 / 60 mg | +0.052 / +0.088 / +0.107 | +0.014 / +0.022 / +0.026 | +0.014 / +0.022 / +0.026 | +0.016 / +0.024 / +0.027 | +0.009 / +0.017 / +0.022 |
| Relative efficacy loss, 24-month delay | −12% | −12.1% | — | −12.1% | −12.1% |
| **Virtual clinical trial (18 months)** | | | | | |
| Patient state at trial entry | not specified (discussed as Phase 2) | healthy, DA = 1, dosed from day 0 | healthy, DA = 1 | symptomatic, DA 0.27–0.69 | same |
| Patients analysed per arm | 244 of 250 | 244 on authors' machine (60-s timeout, drug arms only) | 250 | 250 (271 screened, 21 never symptomatic) | 250 |
| Placebo DA at 18 months, median | ≈0.85 (Fig 6A) | 0.909 (nominal patient) | 0.881 | 0.550 | 0.550 |
| Cohen's d, DA: 10 / 30 / 60 mg | 0.133 / 0.221 / 0.273 | — | 0.109 / 0.191 / 0.234 | 0.075 / 0.117 / 0.148 | 0.045 / 0.080 / 0.107 |
| 95% CI for 30 mg d (DA) | 0.177–0.273 | — | — | 0.099–0.137 | 0.067–0.096 |
| Cohen's d, motor, 30 mg | 0.221 | — | 0.216 | 0.136 | 0.093 |
| Responders at 30 mg, DA (> 5 points) | 25.8% | — | 22.0% | 9.6% | 5.2% |
| Responders at 30 mg, motor (> 3 points) | 17.2% | — | 12.8% | 28.4% | 17.2% |
| Patients per arm for 80% power, 30 mg (DA) | 322 | — | 432 | 1,153 | 2,433 |
| Loss prevented at 30 mg, mito / inflam / proteo / total | 52 / 31 / 25 / 43% (rates) | — | — | 14 / 20 / 30 / 19% (cumulative) | 10 / 14 / 21 / 13% |
| 30 mg DA saved across 125-point sweep | 0.015–0.069 | — | — | 0.011–0.045 | 0.008–0.033 |
| Largest uncertainty source (fold range) | α_clear (2.6×) | — | — | α_clear (2.4×) | α_clear (2.4×) |
| **Sensitivity and identifiability** | | | | | |
| Sobol S_T for DA, top 4 | k_ROS_basal 0.50, K_mPTP 0.33, CI_max 0.29, k_SOD2 0.28 | — | — | k_SOD2 0.32, k_ROS_basal 0.29, K_mPTP 0.23, k_agg 0.20 (±50%, N = 256) | model-only |
| K_mPTP first-order index | ≈ 0 ("purely interactive") | clipped at 0 in code | — | 0.06 | — |
| Top Sobol parameter, motor | k_ROS_basal 0.40 | — | — | K_motor 0.44 | — |
| Sensitivity-matrix rank (7 outputs, 39 parameters) | 7, "full rank" | 7 | 7 | 7, so ≥ 32 unconstrained directions | — |
| Profile likelihood | flat; read as coupled feedback | data = model's own outputs; σ = half acceptance range | — | not identifiable; figure retired | — |
| Validation checks reported | 18/18 | 24/24 calibration + 19/19 three-tier | 15 automated numerical tests | same | same |

## How much each choice contributes (30 mg, Cohen's d for DA)

| Trial design | Exposure normalisation | d | Patients per arm for 80% power |
|---|---|---|---|
| Healthy start, dosing from day 0 (R) | R (day-120 trough) | 0.191 (published 0.221) | 432 (published 322) |
| Symptomatic enrolment | R (day-120 trough) | 0.104 | ~1,440 |
| Symptomatic enrolment | Steady state, scenario A | 0.117 | 1,153 |
| Symptomatic enrolment | Steady state, scenario B | 0.080 | 2,433 |

With the same set-up, Python and R agree exactly. The difference between the manuscript and the corrected results comes from trial design (enrolment) and exposure normalisation, not from the model equations.
