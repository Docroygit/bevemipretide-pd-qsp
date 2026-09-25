# Corrected figures for the BJP submission

Replace these five files in `All_Figures.zip` and `Supplementary_Figures.pdf`. Keep every other original figure:

- Figures 1–4 and S3–S6 are unchanged. Scenario A leaves all mouse simulations identical.
- The Figure 4 legend has been reworded in the manuscript. Its "terminal attractor" labels are in the image itself, so they should be removed or relabelled "day-50 end states" in the original plotting code.

The main-figure legends below are already in `Manuscript_bjp_tracked.docx`.
The supplementary legends need to be updated in `Supplementary_Figures.pdf`.

## Main figures (legends as in the tracked manuscript)

**Figure 5.** Human translation composite. (A) Accumulation of mitochondrial drug exposure under once-daily dosing, shown as the trough relative to the exact steady-state trough, for mouse (orange) and human (blue) PK parameters; the dotted line marks day 120. (B) Twenty-year DA neuron survival for healthy (green), untreated disease (red), and 30 mg daily from year 0 under scenario A (solid blue) and scenario B (dashed blue); horizontal lines mark symptom onset (70% survival), clinical PD (50%), and the Bernheimer threshold (40%). (C) Benefit in DA neuron survival versus placebo after 18 months of treatment at 10, 30, and 60 mg, by diagnostic delay after symptom onset (scenario A).

**Figure 6.** Virtual clinical trial composite (n = 250 symptomatic virtual patients, four matched arms, scenario A). (A) Violin and box plots of DA neuron survival at 18 months; the dotted line marks the 70% enrolment threshold. (B) Waterfall plot of individual benefit at 30 mg relative to each patient’s own placebo simulation; dark bars exceed the 5-percentage-point simulation threshold. (C) Reduction in cumulative 18-month neuron loss relative to placebo, by death pathway and dose.

**Figure 7.** Uncertainty quantification composite. (A) DA neuron survival saved at 18 months by 30 mg across the 125-combination sweep of CI_max_human, α_clear_human, and death_scale_human (scenario A dark, scenario B light). (B) Patient-level benefit (matched ΔDA) by dose, scenario A. (C) Unpaired Cohen’s d with 95% bootstrap intervals (B = 10,000) for DA neuron survival and motor score, by dose and scenario.

## Supplementary figures (new legends)

**Figure S1.** Sobol sensitivity indices for DA neuron survival (A) and motor score (B) at day 35 in the mouse disease model. Thirteen parameters were varied independently over uniform ±50% ranges (N = 256 base samples; 3,840 evaluations; CI_max capped at 1). Bars show total-order indices with 95% bootstrap intervals (500 resamples); points show first-order indices. One-at-a-time and Morris rankings are reported in Table 3.

**Figure S2.** Cumulative neuron loss attributed to the mitochondrial, neuroinflammatory, and proteotoxic pathways over 50 days in the mouse model, untreated (A) and at 5 mg kg⁻¹ daily (B). Each band integrates its pathway hazard multiplied by surviving neurons, so the bands sum to total loss. The drug lowers all three components, the non-mitochondrial ones indirectly through lower oligomer and cytokine levels.

## Also available

`../figures/` holds SVG versions of every figure. It also has a scenario A vs B mouse dose-response (`FigS3_dose_response_corrected`), which can go in the supplement as a sensitivity figure if wanted.
