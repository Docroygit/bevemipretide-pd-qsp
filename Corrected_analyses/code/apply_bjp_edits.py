"""Write the corrected analyses into copies of the BJP manuscript and supplementary tables as tracked changes.

Sources are read-only: ../../Manuscript_bjp.docx and ../../Supplementary_Tables.docx.
Outputs: ../documents/Manuscript_bjp_tracked.docx and ../documents/Supplementary_Tables_tracked.docx.
Numbers are taken from ../results where they are computed; body indices refer to the original documents.
"""
import json
import zipfile
from pathlib import Path
from lxml import etree as E
from redline import Tracker, cell_paragraph, write_docx, q

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROJECT = ROOT.parent
RES = ROOT / 'results'
OUT = ROOT / 'documents'
OUT.mkdir(exist_ok=True)
T = Tracker('Claude (corrected analyses)', '2026-09-23T12:00:00Z')

vct = json.loads((RES / 'vct_summary.json').read_text(encoding='utf-8'))
sob = json.loads((RES / 'sobol_pm50.json').read_text(encoding='utf-8'))


def load(name):
    with zipfile.ZipFile(PROJECT / name) as z:
        root = E.fromstring(z.read('word/document.xml'))
    return root, list(root.find(q('body')))


def sobol_column(j):
    st = [row[j] for row in sob['ST']]
    order = sorted(range(len(st)), key=lambda i: -st[i])[:10]
    name = lambda n: 'α_clear' if n == 'alpha_clear' else n
    return [f'{name(sob["parameters"][i])} ({st[i]:.2f})' for i in order]


def arm(key, ep): return vct['arms'][key][ep]


def fmt_n(n): return f'{n:,}'


# ====================================================================== manuscript
M, mb = load('Manuscript_bjp.docx')
P = {}  # new full text for body paragraphs (original indices)

P[5] = ('An eight-module, 28-ODE model integrated pharmacokinetics, cardiolipin dynamics, α-synuclein pathology, mitochondrial bioenergetics, '
        'oxidative stress, cell death, and neuroinflammation, calibrated against five published experimental benchmarks, checked through an '
        '18-check pathway framework, tested for identifiability, and translated allometrically to a human virtual clinical trial enrolling '
        '250 symptomatic virtual patients in four matched arms.')
P[7] = ('The calibrated model reproduced the experimental benchmarks and passed the directional pathway checks. With drug exposure normalised to '
        'steady state, allometric translation predicted a multi-year delay in symptom onset when treatment began before symptoms. In symptomatic '
        'virtual patients, 30 mg daily produced a small effect on neuron preservation and motor score (Cohen’s d of 0.12 and 0.14), implying '
        'roughly 850 to 1,150 patients per arm for 80% power, with benefit maintained despite realistic diagnostic delays. Through shared-cascade '
        'cross-talk, treatment reduced neuroinflammatory and proteotoxic neuron loss as well as mitochondrial loss.')
P[9] = ('This QSP model provides a mechanistic rationale for cardiolipin as a druggable node in Parkinson’s disease and supplies the dose, timing, '
        'and sample-size estimates needed to plan an informative Phase 2 trial. Because the predicted effect in symptomatic patients is small and '
        'bevemipretide targets only one of several convergent causes of neuronal loss, confirming its effectiveness will require a large trial in '
        'a heterogeneous population assessed against validated motor scales.')
P[20] = None  # targeted replacement below
P[30] = ('Allometric PK scaling from mouse (25 g) to human (70 kg) employed the BW^(−0.25) rule for elimination and distribution rate constants '
         '(Table 4). BBB penetration was scaled using preclinical monkey CSF-to-plasma ratios of Bevemipretide [34], while intracellular kinetics '
         '(mitochondrial uptake and efflux) were left unscaled as physiology-independent processes; the human brain efflux constant (0.010 h⁻¹) '
         'was retained as an assumption, because the BW^(−0.25) rule would give 0.0022 h⁻¹. The FDA K_m body-surface-area method converted the '
         'mouse reference dose (5 mg kg⁻¹ IP) to a human equivalent of 28.4 mg, rounded to 30 mg SC daily. Because the slowest disposition mode '
         'of the linear PK system has a half-life of 56 days in mouse and 112 days in human, steady state is not reached within 120 days; '
         'mitochondrial drug exposure was therefore normalised to the exact steady-state mean mitochondrial amount at the reference dose of each '
         'species, obtained from the periodic solution of the PK system (C_ref = 11.27 in mouse and 48.83 in human, in model amount units). In the '
         'primary analysis (scenario A), k_drug was re-expressed for this normalisation (0.2545 h⁻¹), which leaves every mouse simulation unchanged '
         'and assumes that 30 mg in humans reproduces the mouse steady-state mitochondrial exposure at 5 mg kg⁻¹. A sensitivity analysis '
         '(scenario B) retained k_drug = 0.100 h⁻¹ with the same normalisation. Human disease timescale parameters were recalibrated: CI_max was '
         'set to 0.90, α_clear to 0.35, and all three death rate constants were scaled to 6% of mouse values, reflecting the substantially slower '
         'nigral degeneration observed in humans compared with toxin-based mouse models.')
P[32] = ('A virtual population was generated by sampling 15 parameters from log-normal (mechanistic and PK parameters) or truncated normal '
         '(disease modifiers) distributions (Table ST2). Each candidate was simulated without treatment until DA neuron survival fell below 70% '
         '(symptom onset; 30-year horizon) and then for a further 12-month diagnostic delay; the first 250 eligible candidates were enrolled, and '
         'every screened candidate and exclusion was recorded. Each enrolled patient was simulated in four matched arms: placebo, 10 mg, 30 mg, and '
         '60 mg SC daily for 547 days (18 months), with daily dosing events and drug exposure normalised to the nominal human reference so that '
         'between-patient PK variability was retained. Effect sizes were reported as unpaired Cohen’s d and Cliff’s delta; response rates were '
         'defined as benefit relative to each patient’s own placebo simulation exceeding 5 percentage points of DA neuron survival or 3 points of '
         'motor score, which are simulation thresholds rather than validated minimal clinically important differences. Sample size for 80% power '
         'was calculated from unpaired d using the exact noncentral t distribution (two-sided α = 0.05).')
P[34] = ('A three-tier sensitivity analysis pipeline was applied to 34 parameters across all modules. One-at-a-time (OAT) perturbation of ±50% '
         'generated tornado diagrams (68 model evaluations). Morris elementary effects screening with 20 trajectories (700 evaluations) identified '
         'parameters with strong nonlinear or interactive behaviour. Sobol variance-based decomposition on 13 parameters selected from the OAT and '
         'Morris screens (independent uniform ±50% ranges; N = 256 base samples; 3,840 evaluations; Jansen estimators with 500 bootstrap '
         'resamples) partitioned output variance into first-order and total-order indices (Figure S1, Table 3). Uncertainty was quantified through '
         '10,000 nonparametric bootstrap resamples of the virtual population and a 125-combination joint sweep of the three translational '
         'parameters (Table ST3).')
P[37] = ('For structural identifiability, a sensitivity-based Fisher Information Matrix (FIM) was constructed from the endpoint sensitivity matrix '
         'assuming independent, unit observation errors, and its rank was determined via eigenvalue decomposition. The FIM reached rank 7, the '
         'maximum possible with seven calibration outputs, so the calibration dataset can constrain at most seven combinations of the 39 estimated '
         'parameters and leaves at least 32 directions locally unconstrained (Figure S4B). When the analysis was extended to time-course '
         'sensitivity data sampled at multiple time points along the simulation, the identifiable subspace expanded to 11 directions, indicating '
         'that longitudinal measurements, if available, would resolve approximately 1.6-fold more parameter combinations. The leading eigenvalue '
         'captured 97.5% of the total information content (Figure S4B), indicating that the calibration outputs chiefly inform a single dominant '
         'parameter combination. The normalised sensitivity heatmap (Figure S4C) revealed a well-defined parameter hierarchy, with disease '
         'modifiers (CI_max, α_clear) and ROS pathway parameters (k_ROS_basal, k_shunt) contributing the largest sensitivity coefficients across '
         'multiple outputs, consistent with the Sobol analysis results.')
P[38] = ('Practical identifiability was assessed via profile likelihood analysis of the six calibrated parameters (CI_max, α_clear, k_impair, '
         'K_death, k_damage_mPTP, n_death), in which each parameter was fixed at values spanning its plausible range while all remaining calibrated '
         'parameters were re-optimised (Figure S5). None of the six profiles crossed the 95% χ² threshold within the explored ranges, so these '
         'parameters are not practically identifiable from the calibration targets alone; this reflects the coupled feedback architecture of the '
         'vicious cycle, in which compensatory parameter adjustments maintain similar outputs across wide ranges of individual parameter values, a '
         'property characteristic of tightly coupled biological feedback systems [26]. Collinearity diagnostics (Figure S6) confirmed that the '
         'highest collinearity indices occurred between parameters that participate in the same regulatory switch (e.g., k_SOD2_on and k_SOD2_off, '
         'the activation and deactivation rates of the SOD2 antioxidant response, and k_mPTP_open and k_mPTP_close, the mPTP gating rates), further '
         'supporting the mechanistic interpretation of parameter interdependence. The 32 literature-derived parameters (35.6% of the total '
         'parameter set) are independently constrained by published experimental data and do not depend on the calibration procedure for their '
         'values, reducing, but not removing, the dimensionality of the estimation problem.')
P[42] = None  # append below
P[45] = ('The four-stage calibration identified CI_max = 0.74 and α_clear = 0.25 as the optimal disease modifier combination (Figure 2A). All five '
         'quantitative calibration targets, including DA neuron loss and motor score at 35 days, fell within their pre-specified acceptable ranges '
         '(Figure 2B, Table 2); Complex I activity was 0.484 relative to the healthy simulation, and the TLR2 knockout reduction (20.3%) fell within '
         'the acceptance range but below the 23–45% literature range. The healthy baseline lost 1.1% of DA neurons over 50 days (DA_neuron = 0.989), '
         'narrowly missing the >0.99 criterion, with all pharmacodynamic state variables within 5% of their standalone module steady-state values.')
P[56] = ('Projections of the trajectories onto CL ratio–α-synuclein state space (Figure 4) separated a healthy end state (high CL, low α-synuclein), '
         'a disease end state (low CL, high α-synuclein), and an intermediate drug-treated state, with the drug trajectory initially tracking the '
         'disease path before diverging. Because these runs use different parameters or dosing, and neuron loss continues in all of them, the end '
         'states are not attractors of a single dynamical system. The cascade propagation timeline (Figure S3B) showed that CI dysfunction and '
         'α-synuclein aggregation were the earliest events (onset within 2–3 days), followed by downstream ATP depletion, membrane depolarisation, '
         'and motor impairment.')
P[58] = ('Cumulative decomposition of neuronal loss into its three constituent pathways (Figure S2A) attributed approximately 63% of disease-state '
         'loss at 35 days to mitochondrial damage, with neuroinflammation (12%) and proteotoxicity (25%) contributing the remaining 37%. At '
         '5 mg kg⁻¹, cumulative mitochondrial loss at 35 days fell by 33%, and neuroinflammatory and proteotoxic loss fell by 16% and 10%, because '
         'lower oligomer and cytokine levels reduce those hazards indirectly (Figure S2B).')
dsob, msob = sobol_column(0), sobol_column(1)
P[60] = ('The OAT tornado diagram identified k_agg, k_ROS_basal, CI_max, and k_SOD2 as the top drivers of DA neuron survival under ±50% '
         'perturbation (Table 3). Sobol decomposition over the same ±50% ranges (Figure S1, Table 3) ranked k_SOD2 (S_T = 0.32) and k_ROS_basal '
         '(S_T = 0.29) highest for DA neuron survival, with overlapping bootstrap intervals and roughly half of their influence mediated through '
         'parameter interactions, followed by K_mPTP (S_T = 0.23), whose influence was largely interactive (S_1 = 0.06), consistent with its role '
         'as a threshold switch. For motor score, the threshold of the motor mapping (K_motor, S_T = 0.44) dominated, followed by k_SOD2 and '
         'k_ROS_basal. ROS module parameters therefore dominated DA neuron survival, confirming oxidative stress as the primary amplifier within '
         'the vicious cycle, whereas the motor endpoint also depended strongly on its mapping from neuron loss.')
P[64] = ('Allometric scaling produced human plasma and mitochondrial pharmacokinetic parameters consistent with the mouse-to-human scaling '
         'relationships (Table 4). The slowest disposition mode of the linear PK system had a half-life of 56.5 days in mouse and 112.4 days in '
         'human, so mitochondrial exposure reached 90% of steady state only after 188 and 376 days of daily dosing, respectively (Figure 5A, '
         'Table 4). Twenty-year disease trajectory simulations (Figure 5B) showed a prodromal phase of 5.6 years (consistent with the 5–10-year '
         'range estimated by Postuma and Berg [3]), clinical PD at 10.8 years, and crossing of the Bernheimer threshold (60% DA loss) at 14.3 years. '
         'Treatment with 30 mg bevemipretide from year 0 delayed symptomatic onset by 4.5 years and clinical PD by 8.9 years (3.0 and 5.9 years '
         'under scenario B). The healthy baseline showed 9.6% DA loss at 20 years, consistent with the approximately 5%/decade age-related decline '
         'reported by Fearnley and Lees [35]. Figure 5C summarises the benefit of 18 months of treatment by dose and diagnostic delay (Section 3.9).')
a30, m30, a60 = arm('A_30', 'DA'), arm('A_30', 'motor'), arm('A_60', 'DA')
flow = vct['flow']; excluded = flow['candidates_screened'] - flow['eligible_enrolled']
on = vct['onset_years_median_IQR']; dr = vct['diagnosis_DA_range']
P[68] = (f'Of {flow["candidates_screened"]} candidates screened, {excluded} did not reach symptom onset within 30 years; '
         f'{flow["eligible_enrolled"]} were enrolled a median of {on[0]:.1f} years (IQR {on[1]:.1f}–{on[2]:.1f}) after simulation start, with DA '
         f'neuron survival of {dr[0]:.2f}–{dr[1]:.2f} at diagnosis (median {vct["diagnosis_DA_median"]:.2f}), and all {flow["included_analysed"]} '
         f'were analysed. Treatment produced dose-dependent separation at 18 months (Figure 6A, Table 5): at 30 mg, d = {a30["d_unpaired"]:.3f} '
         f'(95% CI {a30["d_unpaired_95"][0]:.3f}–{a30["d_unpaired_95"][1]:.3f}) for DA neuron preservation and {m30["d_unpaired"]:.3f} for motor '
         f'score, with {100 * a30["response_rate"]:.1f}% and {100 * m30["response_rate"]:.1f}% of patients exceeding the simulation response '
         f'thresholds; approximately {fmt_n(a30["n_per_arm_80pct_power"])} (DA neuron) and {fmt_n(m30["n_per_arm_80pct_power"])} (motor score) '
         f'patients per arm would be required for 80% power. Effect sizes and response rates increased monotonically from 10 mg up to 60 mg, the '
         f'phase 1 multiple ascending dose (MAD) ceiling, with diminishing returns above 30 mg. Individual patient responses (Figure 6B) showed '
         f'substantial heterogeneity, with the top quartile of responders accounting for about 60% of total DA neuron preservation. Under '
         f'scenario B, 30 mg gave d = {arm("B_30", "DA")["d_unpaired"]:.3f} and {arm("B_30", "motor")["d_unpaired"]:.3f}.')
r30 = vct['arms']['A_30']['rescue_fraction_cumulative_loss']
P[69] = (f'Relative to placebo, 30 mg reduced cumulative 18-month neuron loss by {100 * r30["mitochondrial"]:.0f}% in the mitochondrial pathway, '
         f'{100 * r30["inflammation"]:.0f}% in the neuroinflammatory pathway, and {100 * r30["proteotoxic"]:.0f}% in the proteotoxic pathway '
         f'({100 * r30["total"]:.0f}% overall; Figure 6C). In symptomatic patients, indirect reductions in oligomer and cytokine levels therefore '
         f'contributed at least as much benefit as the direct mitochondrial effect.')
env = json.loads((RES / 'envelope.json').read_text(encoding='utf-8'))['A']
P[73] = (f'Bootstrap 95% confidence intervals for the 30 mg unpaired d were {a30["d_unpaired_95"][0]:.3f}–{a30["d_unpaired_95"][1]:.3f} for DA '
         f'neuron preservation and {m30["d_unpaired_95"][0]:.3f}–{m30["d_unpaired_95"][1]:.3f} for motor score (Figure 7C, Table ST3); because '
         f'they resample virtual patients from one assumed population, these intervals describe simulation precision rather than uncertainty in '
         f'the clinical effect. The joint parameter uncertainty envelope, sweeping CI_max_human, α_clear_human, and death_scale_human across 125 '
         f'combinations, produced a 30 mg DA_saved range of [{env["DA_saved_min_max"][0]:.3f}, {env["DA_saved_min_max"][1]:.3f}] (2.5th–97.5th '
         f'percentile [{env["DA_saved_2.5_97.5"][0]:.3f}, {env["DA_saved_2.5_97.5"][1]:.3f}]) and an approximate d range of '
         f'[{env["approx_d_min_max"][0]:.2f}, {env["approx_d_min_max"][1]:.2f}] (Figure 7A). The parameter α_clear_human was the single largest '
         f'source of uncertainty, producing a {env["fold_range_of_mean_DA_saved_over_alpha_clear"]:.1f}-fold range in mean DA_saved, followed by '
         f'death_scale_human ({env["fold_range_of_mean_DA_saved_over_death_scale"]:.1f}-fold) and CI_max_human '
         f'({env["fold_range_of_mean_DA_saved_over_CI_max"]:.1f}-fold). Patient-level benefit at 30 mg ranged from near zero to more than 12 '
         f'percentage points of DA neuron survival (Figure 7B).')
P[75] = ('The phased simulation showed that drug benefit was preserved across diagnostic delays of 0–24 months, with only 12% relative efficacy loss '
         'at 24 months compared with immediate treatment (Table ST4). The dose-by-delay matrix showed that 30 mg and 60 mg retained most of their '
         'benefit (ΔDA of +0.022 and +0.025 after 18 months) even when treatment was delayed by two years after symptom onset (Figure 5C). The '
         'dose-response curve in the mouse model (Figure S3A) demonstrated a steep initial response at low doses with progressive flattening above '
         '5 mg kg⁻¹, consistent with the diminishing returns observed in the human virtual trial.')
P[77] = ('Our QSP model integrated pharmacokinetics, cardiolipin dynamics, α-synuclein pathology, mitochondrial bioenergetics, oxidative stress, '
         'mitophagy, cell death, and neuroinflammation into a single coupled 28-ODE system, simulating the vicious cycle across mouse and '
         'allometrically translated human physiology. The coupled system behaved as an emergent amplifier, converting modest standalone '
         'perturbations into substantial DA neuron loss and motor impairment, with bevemipretide only partially reversing this trajectory given its '
         'single-node mechanism. In symptomatic virtual patients the predicted effect at the clinically relevant dose was small (d ≈ 0.12–0.14), '
         'and the identifiability analysis showed that the calibration targets constrain only a few parameter combinations, so the quantitative '
         'predictions are conditional on the inherited parameter values.')
P[79] = ('This amplification arises because each node feeds forward into the next: α-synuclein oligomers accelerate cardiolipin oxidation [16], '
         'depleted cardiolipin destabilises respiratory supercomplexes [17], impaired electron transport elevates ROS [39], and ROS promotes further '
         'α-synuclein misfolding [40], a positive feedback loop in which the coupled perturbation exceeds any single-pathway perturbation. The '
         'trajectory projections illustrate this: the drug-treated trajectory ends between the healthy and disease end states in CL ratio–'
         'α-synuclein state space rather than returning to health. This intermediate end state reflects that bevemipretide stabilises cardiolipin '
         'at a single node of a four-node cycle and cannot compensate for upstream impairments in Complex I repair capacity and α-synuclein '
         'clearance.')
P[80] = None  # targeted
P[81] = (f'The three-pathway death architecture, comprising mitochondrial damage, neuroinflammation-mediated toxicity via TNF-α and IL-1β [43], and '
         f'direct α-synuclein proteotoxicity [15], shapes but does not cap monotherapy efficacy. In the mouse disease state, approximately 63% of '
         f'cumulative neuronal loss is attributed to the mitochondrial pathway and 37% to neuroinflammation (12%) and proteotoxicity (25%). These '
         f'pathways are not independent of the drug: in symptomatic virtual patients, 30 mg reduced cumulative loss by '
         f'{100 * r30["mitochondrial"]:.0f}% in the mitochondrial pathway but by {100 * r30["inflammation"]:.0f}% and '
         f'{100 * r30["proteotoxic"]:.0f}% in the neuroinflammatory and proteotoxic pathways, because cardiolipin stabilisation lowers oligomer '
         f'burden and microglial activation (Figure 6C). Monotherapy efficacy is therefore limited less by an unreachable fraction of death than by '
         f'the modest size of the overall effect ({100 * r30["total"]:.0f}% reduction in total loss at 30 mg and '
         f'{100 * vct["arms"]["A_60"]["rescue_fraction_cumulative_loss"]["total"]:.0f}% at 60 mg). A monotherapy trial powered to detect a large '
         f'effect size will likely fail, not because the drug is inactive but because the achievable effect is small; combination strategies '
         f'pairing cardiolipin stabilisation with anti-inflammatory agents or aggregation inhibitors could add to this effect and widen the '
         f'therapeutic window.')
P[82] = ('The predicted unpaired Cohen’s d at 30 mg (0.12 for DA neuron preservation and 0.14 for motor score) is below the values cited for the '
         'ADAGIO rasagiline trial (d ≈ 0.20) [44] and the SPARK cinpanemab trial (d ≈ 0.15) [5]. Although this places a model-simulated projection '
         'alongside empirically measured clinical outcomes and should be read as a trial-planning benchmark rather than evidence of comparable '
         'real-world efficacy, it indicates that a trial powered for d ≈ 0.2 would be underpowered for the effect predicted here. The model’s '
         'prediction of diminishing returns above 30 mg is independently supported by preclinical ALS TDP-43 mouse data for SBT-272 [27], where '
         'upper motor neuron retention and reduced neuroinflammation showed a dose-response plateau, while Phase 1 data confirmed systemic '
         'exposure profiles consistent with the allometric PK scaling used here [28], an independent check on the translational bridge.')
P[83] = None  # targeted
P[84] = ('The three-tier sensitivity analysis converged on ROS module parameters as the main determinants of DA neuron survival: k_SOD2 and '
         'k_ROS_basal had the highest Sobol total-order indices (S_T = 0.32 and 0.29, with overlapping bootstrap intervals) and ranked among the '
         'top parameters in the OAT and Morris screens. Roughly half of their influence was mediated through parameter interactions, consistent '
         'with their role as amplifier gain within the vicious cycle and aligning with post-mortem evidence of elevated '
         '8-hydroxy-2′-deoxyguanosine and decreased glutathione in PD substantia nigra [48, 49]. The motor score was dominated by the threshold of '
         'its mapping from neuron loss (K_motor, S_T = 0.44), so the composite motor endpoint adds structural uncertainty beyond the disease '
         'biology. This ROS dominance implies that interventions upstream of ROS production, cardiolipin stabilisation, supercomplex preservation, '
         'may have a disproportionate effect by attenuating the cycle’s primary amplifier.')
P[85] = ('The identifiability analysis defines how the model’s predictions should be interpreted. The seven calibration outputs can constrain at '
         'most seven combinations of the 39 estimated parameters, and the profile likelihoods of all six calibrated parameters were flat, so '
         'individual parameter values are not practically identifiable from the calibration targets. These interdependencies reflect the tightly '
         'coupled feedback structure of the vicious cycle, a property shared by other published QSP models of comparable complexity [26]. That '
         '35.6% of the parameter set is constrained by independent literature data limits, but does not remove, this uncertainty; the quantitative '
         'predictions are therefore conditional on the inherited parameter values and should be refined when longitudinal data become available.')
P[86] = None  # targeted
P[87] = ('The allometrically scaled human model reproduced several independent clinical observations without additional fitting: a prodromal '
         'phase of 5.6 years (within the 5–10-year range estimated by Postuma and Berg [3]), age-related DA loss of 9.6% at 20 years (consistent '
         'with ~5%/decade, per Fearnley and Lees [35]), and Bernheimer threshold crossing at 14.3 years. The predicted 4.5-year delay in symptomatic '
         'onset when 30 mg is started before symptoms (3.0 years under scenario B), if realised in patients, would represent a clinically '
         'substantial benefit. The phased simulation showed efficacy preserved despite diagnostic delay, with a 24-month delay reducing relative '
         'benefit by only 12%, indicating the vicious cycle, once partially interrupted, can be slowed even after substantial neurodegeneration has '
         'occurred, since bevemipretide acts on cardiolipin, a node accessible regardless of disease stage rather than an upstream trigger that may '
         'have already passed.')
P[88] = ('A prediction of the model is that absolute drug benefit grows over time: with 30 mg started at year 0, DA neuron survival was 2.3 '
         'percentage points higher than untreated at year 1 and 21.9 points higher at year 20 (a 79% relative increase); this 20-year projection '
         'extends well beyond the model’s explicitly validated one-year extrapolation window and should be read as illustrating the compounding '
         'structure of partial cycle interruption rather than a precise long-term efficacy estimate. This accelerating separation arises because '
         'partial interruption of the positive feedback loop compounds over long timescales, implying that short trials (12–18 months) will '
         'underestimate long-term benefit. The model quantifies this gap and could inform long-term extension studies or Bayesian adaptive trials '
         'using early biomarker signals to project long-term clinical benefit.')
P[90] = (f'The virtual clinical trial revealed substantial inter-patient heterogeneity in treatment response, with the top quartile of responders '
         f'accounting for about 60% of total DA neuron preservation, reflecting variability in the 15 sampled parameters. At 30 mg, '
         f'{100 * a30["response_rate"]:.1f}% of patients exceeded the 5-point DA neuron threshold and {100 * m30["response_rate"]:.1f}% the 3-point '
         f'motor threshold; heterogeneity of this kind routinely dilutes population-level effect sizes in PD trials [51]. No single sampled '
         f'parameter, including the disease modifiers, predicted individual benefit strongly (all |Spearman ρ| ≤ 0.12), so identifying likely '
         f'responders in a Phase 2 adaptive design will probably require composite baseline biomarker profiles rather than a single stratifier.')
P[91] = None  # targeted
P[93] = None  # targeted
P[94] = None  # targeted
P[98] = ('This work presents a 28-ODE quantitative systems pharmacology model encoding the α-synuclein, cardiolipin, Complex I, and reactive oxygen '
         'species vicious cycle underlying dopaminergic neurodegeneration in Parkinson’s disease. The model was calibrated against five published '
         'experimental benchmarks, checked through 18 quantitative and directional tests, and stress-tested with a three-tier sensitivity analysis '
         'and a structural and practical identifiability analysis, which together showed that its predictions depend chiefly on ROS-module '
         'parameters and that the calibration targets constrain only a few parameter combinations.')
P[99] = ('The analysis identifies cardiolipin stabilisation as a mechanistically coherent point of intervention: bevemipretide acts at the most '
         'upstream druggable node of the vicious cycle, and even partial interruption there propagates through the coupled system to produce '
         'disease-modifying benefit that grows over treatment rather than saturating early. At the same time, roughly a third of neuronal loss '
         'proceeds through neuroinflammatory and proteotoxic routes that a cardiolipin-directed drug reaches only indirectly. Allometric translation '
         'reproduced several independent clinical benchmarks without additional fitting, and the resulting virtual clinical trial yields concrete '
         'guidance for a Phase 2 programme: a small model-predicted effect size in symptomatic patients (d ≈ 0.12–0.14 at 30 mg) implying roughly '
         '850 to 1,150 patients per arm, a dose of 30 mg beyond which returns diminish, and efficacy preserved across realistic diagnostic delays. '
         'As discussed above, confirming these predictions requires a large trial in a real, aetiologically heterogeneous patient population '
         'assessed against validated clinical motor scales.')
P[164] = ('Figure 4. Vicious cycle phase portrait. Trajectories in CL ratio–α-synuclein oligomer state space. Diamonds mark day-50 end states; '
          'arrows indicate direction of flow. Healthy, disease, and drug-treated trajectories use different parameters or dosing, so the end '
          'states are not attractors of a single dynamical system.')
P[165] = ('Figure 5. Human translation composite. (A) Accumulation of mitochondrial drug exposure under once-daily dosing, shown as the trough '
          'relative to the exact steady-state trough, for mouse (orange) and human (blue) PK parameters; the dotted line marks day 120. '
          '(B) Twenty-year DA neuron survival for healthy (green), untreated disease (red), and 30 mg daily from year 0 under scenario A (solid '
          'blue) and scenario B (dashed blue); horizontal lines mark symptom onset (70% survival), clinical PD (50%), and the Bernheimer threshold '
          '(40%). (C) Benefit in DA neuron survival versus placebo after 18 months of treatment at 10, 30, and 60 mg, by diagnostic delay after '
          'symptom onset (scenario A).')
P[166] = ('Figure 6. Virtual clinical trial composite (n = 250 symptomatic virtual patients, four matched arms, scenario A). (A) Violin and box '
          'plots of DA neuron survival at 18 months; the dotted line marks the 70% enrolment threshold. (B) Waterfall plot of individual benefit at '
          '30 mg relative to each patient’s own placebo simulation; dark bars exceed the 5-percentage-point simulation threshold. (C) Reduction in '
          'cumulative 18-month neuron loss relative to placebo, by death pathway and dose.')
P[167] = ('Figure 7. Uncertainty quantification composite. (A) DA neuron survival saved at 18 months by 30 mg across the 125-combination sweep of '
          'CI_max_human, α_clear_human, and death_scale_human (scenario A dark, scenario B light). (B) Patient-level benefit (matched ΔDA) by dose, '
          'scenario A. (C) Unpaired Cohen’s d with 95% bootstrap intervals (B = 10,000) for DA neuron survival and motor score, by dose and '
          'scenario.')

for i, text in P.items():
    if text is not None: T.revise(mb[i], text)

# Targeted replacements within otherwise unchanged paragraphs.
T.replace(mb[20], 'All 28 state variables are normalised to the interval [0, 1], representing fractional activities or concentrations relative to physiological maxima.',
          'The 23 pharmacodynamic state variables are dimensionless indices on the interval [0, 1], representing fractional activities or pool sizes relative to physiological maxima; the five pharmacokinetic states are drug amounts in model units.')
T.replace(mb[28], ', and emergent vicious cycle amplification exceeding standalone module predictions.', ', and oligomer elevation in the coupled system exceeding that of the α-synuclein module alone.')
T.replace(mb[28], 'All 28 state variables were verified to remain bounded within [0, 1]',
          'All pharmacodynamic state variables were verified to remain within [0, 1], and pharmacokinetic amounts to remain non-negative,')
T.replace(mb[28], 'Tier 3 (Extrapolation)', 'Because k_drug was estimated from the same study, Tier 2 tests internal consistency rather than independent prediction. Tier 3 (Extrapolation)')
T.replace(mb[42], 'https://github.com/Docroygit/bevemipretide-pd-qsp.',
          'https://github.com/Docroygit/bevemipretide-pd-qsp. Corrected analyses (steady-state exposure normalisation, virtual-trial enrolment, '
          'phased simulation, cumulative pathway attribution, translational sweep, and Sobol indices) were regenerated with an independent Python '
          'implementation of the same equations and parameters (Python 3.9; SciPy LSODA, relative tolerance 10⁻⁸, absolute tolerance 10⁻¹⁰; no '
          'state clamping), which reproduced the untreated R outputs to three decimal places.')
T.replace(mb[50], 'from 0.989 to 0.640', 'from 0.992 to 0.640')
T.replace(mb[51], '26.3%', '26.0%')
T.replace(mb[51], 'All six checks passed.', 'Because k_drug was estimated from this study, these checks confirm consistency with the data used for estimation rather than independent prediction. All six checks passed.')
T.replace(mb[52], '90-day (0.292), 180-day (0.081)', '90-day (0.291), 180-day (0.080)')
T.replace(mb[52], 'against a healthy baseline of 0.989', 'against a healthy value of 0.920')
T.replace(mb[52], 'The coupled vicious cycle produced 6.7-fold α-synuclein amplification compared with standalone module predictions of approximately 1.5-fold, confirming emergent nonlinear behaviour from inter-module coupling.',
          'In the coupled system, α-synuclein oligomers rose 6.7-fold over healthy levels, compared with approximately 1.5-fold when the α-synuclein module was simulated in isolation, indicating that inter-module coupling accounts for most of the elevation; this is a disease-to-healthy ratio, not a measured loop gain.')
T.replace(mb[53], 'boundedness verification of all 28 state variables', 'verification that pharmacodynamic states remained within [0, 1] and pharmacokinetic amounts non-negative')
T.replace(mb[80], 'identified by Sobol analysis as almost purely interactive (S_1 ≈ 0, S_T = 0.33) and thus an all-or-nothing gate rather than a graded effector',
          'identified by Sobol analysis as acting largely through interactions (S_1 = 0.06, S_T = 0.23 for DA neuron survival), consistent with a threshold-like gate rather than a graded effector')
T.replace(mb[83], 'The bootstrap 95% confidence interval for the 30 mg effect size [0.177, 0.273] places even the model’s lower-bound projection above the effect sizes observed in SPARK [5], STEADY-PD III [47], and SURE-PD3 [46]; as above, this is a comparison between a simulated prediction and measured clinical outcomes, and cross-trial',
          'The corrected 30 mg projection (d = 0.12; 95% bootstrap interval 0.10–0.14) should therefore not be read as exceeding the effect sizes observed in SPARK [5], STEADY-PD III [47], and SURE-PD3 [46]; as above, cross-trial')
T.replace(mb[86], '(2.6-fold range in DA neuron preservation)', '(2.4-fold range in mean DA neuron preservation across the translational sweep)')
T.replace(mb[91], 'by supplying the effect-size and sample-size estimates needed to power it and',
          'by supplying the effect-size and sample-size estimates needed to power it, which at roughly 850 to 1,150 patients per arm for 30 mg indicate that an adequately powered trial will be large, and')
T.replace(mb[93], 'rather than from first-principles biophysics.',
          'rather than from first-principles biophysics; human drug potency additionally assumes that 30 mg reproduces the mouse steady-state mitochondrial exposure at 5 mg kg⁻¹, and retaining the original k_drug (scenario B) reduced the 30 mg effect size by about one-third.')
T.replace(mb[94], 'how well these map onto slowly progressive human PD remains open.',
          'how well these map onto slowly progressive human PD remains open. Seventh, the calibration targets constrain only a few parameter combinations, so predictions are conditional on literature-derived and estimated parameter values.')

# ---------------- Table 2
t2 = mb[47]
T.revise(cell_paragraph(t2, 1, 4), '0.459 (0.484 vs healthy)')
T.revise(cell_paragraph(t2, 5, 6), 'PASS (below literature range)')
T.revise(cell_paragraph(t2, 7, 4), '+14.3% (26.0% of gap)')
T.revise(cell_paragraph(t2, 13, 1), 'PD states in [0,1]; PK amounts ≥ 0')
T.revise(cell_paragraph(t2, 14, 4), '0.989')
T.revise(cell_paragraph(t2, 14, 6), 'Not met (1.1% loss)')
note = mb[23]
T.new_paragraph(t2, 'Complex I activity is shown as the absolute state and relative to the healthy simulation (0.484), the quantity comparable with '
                '~49% of wild type. The TLR2 knockout reduction lies within the acceptance range but below the 23–45% literature range. Drug rescue: '
                'relative gain in DA neuron survival at 5 mg kg⁻¹ on day 35 (fraction of the disease-to-healthy gap rescued). The healthy baseline '
                'lost 1.1% of DA neurons over 50 days, narrowly missing the >0.99 criterion.', note)

# ---------------- Table 3 (Sobol columns replaced)
t3 = mb[62]
for r in range(10):
    T.revise(cell_paragraph(t3, r + 1, 3), dsob[r])
    T.revise(cell_paragraph(t3, r + 1, 4), msob[r])
T.new_paragraph(t3, f'OAT and Morris columns are from the original analysis. Sobol total-order indices were recomputed on 13 parameters with '
                f'independent uniform ±50% ranges (N = {sob["n_base"]} base samples; {sob["evaluations"]:,} evaluations; CI_max capped at 1) for '
                f'day-35 mouse endpoints; 95% bootstrap intervals and first-order indices are shown in Figure S1.', note)

# ---------------- Table 4
t4 = mb[66]
T.revise(cell_paragraph(t4, 6, 3), 'Assumed (BW⁻⁰·²⁵ gives 0.0022 h⁻¹)')
T.revise(cell_paragraph(t4, 9, 1), '11.27')
T.revise(cell_paragraph(t4, 9, 2), '48.83')
T.revise(cell_paragraph(t4, 9, 3), 'Exact steady-state mean at reference dose')
T.revise(cell_paragraph(t4, 10, 3), 'FDA K_m conversion (28.4 mg, rounded)')
T.insert_row(t4, 10, ['Time to 90% of SS trough', '—', '—', 'Daily dosing from zero', '188 d', '376 d'])
T.insert_row(t4, 10, ['Slowest disposition mode', '—', '—', 'Eigenvalue of linear PK system', '56.5 d', '112.4 d'])
T.new_paragraph(t4, 'Rate-constant half-lives are ln2/k for the individual process, not terminal half-lives of the coupled system. C_ref is the exact '
                'steady-state mean mitochondrial amount for one reference dose daily (previously 4.43 and 2.15, described as 120-day values); '
                'k_drug was re-expressed as 0.2545 h⁻¹ so that mouse simulations are unchanged (scenario A).', note)

# ---------------- Table 5
t5 = mb[70]
T.revise(t5, 'Table 5. Virtual clinical trial simulation results at 18 months (n = 250 per arm, matched; scenario A).')
t5t = mb[71]
T.replace(cell_paragraph(t5t, 0, 5), 'MCID', 'Threshold ')
rows = [('A_10', 'DA'), ('A_30', 'DA'), ('A_60', 'DA'), ('A_10', 'motor'), ('A_30', 'motor'), ('A_60', 'motor')]
for r, (k, ep) in enumerate(rows, start=1):
    s = arm(k, ep); md = s['mean_difference']
    vals = [f'{"+" if md > 0 else "−"}{abs(md):.3f}', f'{s["d_unpaired"]:.3f}', f'+{s["cliffs_delta"]:.3f}',
            f'{100 * s["response_rate"]:.1f}%', fmt_n(s['n_per_arm_80pct_power'])]
    for c, v in enumerate(vals, start=2): T.revise(cell_paragraph(t5t, r, c), v)
b = lambda k, ep: f'{arm(k, ep)["d_unpaired"]:.3f}'
T.new_paragraph(t5t, 'Scenario A exposure normalisation (Section 2.4). Response: benefit relative to each patient’s own placebo simulation exceeding 5 '
                'percentage points (DA neuron survival) or 3 points (motor score); these are simulation thresholds, not validated minimal clinically '
                'important differences. N/arm from the exact noncentral t distribution (two-sided α = 0.05). Scenario B (k_drug = 0.100 h⁻¹): '
                f'd = {b("B_10", "DA")}, {b("B_30", "DA")}, and {b("B_60", "DA")} (DA neuron) and {b("B_10", "motor")}, {b("B_30", "motor")}, and '
                f'{b("B_60", "motor")} (motor score) at 10, 30, and 60 mg. Bootstrap 95% intervals are given in Table ST3.', note)

write_docx(PROJECT / 'Manuscript_bjp.docx', M, OUT / 'Manuscript_bjp_tracked.docx')

# ====================================================================== supplementary tables
S, sb = load('Supplementary_Tables.docx')
st1 = sb[9]
T.revise(cell_paragraph(st1, 6, 1), '11.27')
T.revise(cell_paragraph(st1, 6, 4), 'Reference mitochondrial amount (exact steady-state mean, 5 mg kg⁻¹ daily)')
T.revise(cell_paragraph(st1, 6, 5), 'Computed')
T.revise(cell_paragraph(st1, 6, 6), 'Periodic PK solution')
T.revise(cell_paragraph(st1, 14, 1), '0.2545')
T.revise(cell_paragraph(st1, 14, 4), 'Drug-mediated CL_ox reduction (equivalent to 0.100 with C_ref = 4.43)')
T.revise(cell_paragraph(st1, 17, 4), 'CL synthesis pool scale')
T.revise(cell_paragraph(st1, 17, 6), 'Linear feedback form')
T.revise(cell_paragraph(st1, 18, 2), '—')
for r in (50, 69, 71, 72, 73): T.revise(cell_paragraph(st1, r, 2), 'h⁻¹')

T.replace(sb[3], '(N = 250/arm)', '(n = 250 enrolled; four matched arms)')
T.replace(sb[10], '(N = 250/arm)', '(n = 250 enrolled; four matched arms)')
T.replace(sb[11], '4 arms: placebo, 10 mg, 30 mg, 60 mg SC daily; 547 days; 244/250 patients completed per arm (97.6%).',
          f'Candidates were screened until untreated DA neuron survival fell below 0.70 within 30 years ({flow["candidates_screened"]} screened; '
          f'{excluded} did not reach onset), then progressed untreated for a 12-month diagnostic delay; the {flow["eligible_enrolled"]} enrolled '
          f'patients were each simulated in four matched arms (placebo, 10, 30, and 60 mg SC daily; 547 days), and all {flow["included_analysed"]} '
          f'were analysed. Drug exposure was normalised to the nominal human steady-state reference, preserving between-patient PK variability.')

T.replace(sb[14], '(B = 10,000, percentile method)', '(B = 10,000, percentile method; scenario A)')
st3a = sb[15]
spec = [('A_10', 'DA'), ('A_30', 'DA'), ('A_60', 'DA'), ('A_30', 'motor')]
for g, (k, ep) in enumerate(spec):
    s = arm(k, ep)
    trip = [(f'{s["d_unpaired"]:.3f}', *(f'{x:.3f}' for x in s['d_unpaired_95'])),
            (f'{s["cliffs_delta"]:.3f}', *(f'{x:.3f}' for x in s['cliffs_delta_95'])),
            (f'{100 * s["response_rate"]:.1f}%', *(f'{100 * x:.1f}%' for x in s['response_rate_95']))]
    for j, vals in enumerate(trip):
        r = 1 + 3 * g + j
        if j == 2: T.replace(cell_paragraph(st3a, r, 2), 'MCID ', 'Threshold ')
        for c, v in enumerate(vals, start=3): T.revise(cell_paragraph(st3a, r, c), v)
st3b = sb[17]; envB = json.loads((RES / 'envelope.json').read_text(encoding='utf-8'))['B']
import csv
with (RES / 'envelope.csv').open(encoding='utf-8') as f:
    nominal = next(r for r in csv.DictReader(f) if abs(float(r['CI_max']) - .9) < 1e-9 and abs(float(r['alpha_clear']) - .35) < 1e-9 and abs(float(r['death_scale']) - .06) < 1e-9)
T.revise(cell_paragraph(st3b, 1, 3), f'{env["fold_range_of_mean_DA_saved_over_CI_max"]:.1f}-fold range in mean DA_saved')
T.revise(cell_paragraph(st3b, 2, 3), f'Largest ({env["fold_range_of_mean_DA_saved_over_alpha_clear"]:.1f}-fold range in mean DA_saved)')
T.revise(cell_paragraph(st3b, 3, 3), f'{env["fold_range_of_mean_DA_saved_over_death_scale"]:.1f}-fold range in mean DA_saved')
T.revise(cell_paragraph(st3b, 5, 1), '2.5th–97.5th percentile')
T.revise(cell_paragraph(st3b, 5, 2), f'{float(nominal["DA_saved_A"]):.3f}')
T.revise(cell_paragraph(st3b, 5, 3), f'[{env["DA_saved_2.5_97.5"][0]:.3f}, {env["DA_saved_2.5_97.5"][1]:.3f}]')
T.revise(cell_paragraph(st3b, 6, 2), f'{a30["d_unpaired"]:.3f}')
T.revise(cell_paragraph(st3b, 6, 3), f'[{env["approx_d_min_max"][0]:.2f}, {env["approx_d_min_max"][1]:.2f}]')
T.revise(sb[18], 'Bootstrap: nonparametric resampling of matched virtual patients across all 4 arms; intervals describe simulation precision within '
         'one assumed population, not uncertainty in the clinical effect. Parameter envelope: 5 × 5 × 5 grid sweep of three translational '
         f'parameters (scenario A; all 125 combinations reached symptom onset); approximate d = DA_saved divided by the pooled 18-month SD of the '
         f'virtual trial ({env["pooled_sd_used"]:.3f}). Nominal d_unpaired ({a30["d_unpaired"]:.3f}) sits within this range. Scenario B '
         f'(k_drug = 0.100 h⁻¹): 30 mg d = {b("B_30", "DA")} (DA neuron) and {b("B_30", "motor")} (motor score); DA_saved 2.5th–97.5th '
         f'percentile [{envB["DA_saved_2.5_97.5"][0]:.3f}, {envB["DA_saved_2.5_97.5"][1]:.3f}]. α_clear dominates uncertainty because it directly '
         'gates vicious cycle entry via aSyn clearance.')

with (RES / 'ST4_phased.csv').open(encoding='utf-8') as f:
    st4 = [r for r in csv.DictReader(f) if r['scenario'] == 'A']
get = lambda d, mg: next(r for r in st4 if int(r['delay_months']) == d and int(r['dose_mg']) == mg)
T.replace(sb[22], '(18-month treatment, 30 mg)', '(18-month treatment, 30 mg, scenario A)')
st4b, st4c = sb[23], sb[25]
for r, d in enumerate([0, 6, 12, 24], start=1):
    x = get(d, 30)
    vals = [f'{float(x["DA_at_start"]):.3f}', f'{float(x["DA_after_18mo_treated"]):.3f}', f'{float(x["DA_after_18mo_placebo"]):.3f}',
            f'+{float(x["delta_DA"]):.3f}']
    for c, v in enumerate(vals, start=1): T.revise(cell_paragraph(st4b, r, c), v)
    if d: T.revise(cell_paragraph(st4b, r, 5), f'−{abs(100 * float(x["relative_change_vs_0mo_delay"])):.0f}%')
    for c, mg in enumerate([10, 30, 60], start=1): T.revise(cell_paragraph(st4c, r, c), f'+{float(get(d, mg)["delta_DA"]):.3f}')
T.replace(sb[26], 'Drug benefit is robust across delays',
          'Treated DA neuron survival never exceeds survival at treatment start, and the untreated comparator is evaluated at the same 18-month time point. Drug benefit is robust across delays')
T.replace(sb[30], 'All state variables normalised to [0, 1].',
          'Pharmacodynamic state variables are dimensionless on [0, 1]; pharmacokinetic states are drug amounts in model units, and C_mito,norm = C_mito / C_ref, with C_ref the exact steady-state mean at the reference dose.')
T.replace(sb[30], 'atol/rtol = 10⁻¹⁰.', 'atol/rtol = 10⁻¹⁰; corrected analyses used SciPy LSODA (rtol 10⁻⁸, atol 10⁻¹⁰).')

write_docx(PROJECT / 'Supplementary_Tables.docx', S, OUT / 'Supplementary_Tables_tracked.docx')
print('written', OUT / 'Manuscript_bjp_tracked.docx', OUT / 'Supplementary_Tables_tracked.docx', 'last revision id', T.next_id)
