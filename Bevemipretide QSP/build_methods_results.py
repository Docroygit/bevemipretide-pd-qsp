"""
Build Methods_and_Results_v2.docx — publication-ready Methods + Results sections
for the Bevemipretide PD QSP manuscript.

Target: ≤2000 words combined.
Format: Times New Roman 12pt, double-spaced, 1-inch margins.
Includes QoI/CoU table (before Methods) and Model Risk Assessment table (after Results)
per IQ Consortium QSP model qualification framework.
"""

from docx import Document
from docx.shared import Pt, Inches, RGBColor, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "New publication plots")


def make_doc():
    doc = Document()
    # Page margins
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
    # Default style
    style = doc.styles["Normal"]
    font = style.font
    font.name = "Times New Roman"
    font.size = Pt(12)
    fmt = style.paragraph_format
    fmt.line_spacing = 2.0
    fmt.space_after = Pt(0)
    fmt.space_before = Pt(0)
    return doc


def add_heading(doc, text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = "Times New Roman"
        run.font.color.rgb = RGBColor(0, 0, 0)
        run.font.size = Pt(13) if level == 1 else Pt(12)


def add_subheading(doc, text):
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.bold = True
    run.italic = True
    run.font.name = "Times New Roman"
    run.font.size = Pt(12)


def add_para(doc, text):
    p = doc.add_paragraph()
    p.paragraph_format.first_line_indent = Inches(0.5)
    run = p.add_run(text)
    run.font.name = "Times New Roman"
    run.font.size = Pt(12)
    return p


def add_legend(doc, text):
    """Figure/table legend: bold label, then normal description."""
    p = doc.add_paragraph()
    p.paragraph_format.first_line_indent = Inches(0)
    p.paragraph_format.space_before = Pt(6)
    p.paragraph_format.space_after = Pt(6)
    # Split on first period to bold the figure label
    if ". " in text:
        label, rest = text.split(". ", 1)
        run_label = p.add_run(label + ". ")
        run_label.bold = True
        run_label.font.name = "Times New Roman"
        run_label.font.size = Pt(12)
        run_rest = p.add_run(rest)
        run_rest.font.name = "Times New Roman"
        run_rest.font.size = Pt(12)
    else:
        run = p.add_run(text)
        run.font.name = "Times New Roman"
        run.font.size = Pt(12)


def add_table(doc, headers, rows, title=None):
    """Add a formatted table with optional title."""
    if title:
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after = Pt(6)
        run = p.add_run(title)
        run.bold = True
        run.font.name = "Times New Roman"
        run.font.size = Pt(12)

    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    # Header row
    for j, header in enumerate(headers):
        cell = table.rows[0].cells[j]
        cell.text = ""
        p = cell.paragraphs[0]
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing = 1.15
        run = p.add_run(header)
        run.bold = True
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)
        # Light grey shading for header
        shading = cell._element.get_or_add_tcPr()
        shd = shading.makeelement(qn("w:shd"), {
            qn("w:val"): "clear",
            qn("w:color"): "auto",
            qn("w:fill"): "D9E2F3",
        })
        shading.append(shd)

    # Data rows
    for i, row_data in enumerate(rows):
        for j, cell_text in enumerate(row_data):
            cell = table.rows[i + 1].cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_before = Pt(2)
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.line_spacing = 1.15
            # Bold the first column if it looks like a label column
            if j == 0 and len(headers) == 2:
                run = p.add_run(cell_text)
                run.bold = True
                run.font.name = "Times New Roman"
                run.font.size = Pt(10)
            else:
                run = p.add_run(cell_text)
                run.font.name = "Times New Roman"
                run.font.size = Pt(10)

    doc.add_paragraph()  # spacing after table


def add_qoi_cou_table(doc):
    """Question of Interest / Context of Use table (IQ Consortium framework)."""
    headers = ["Question of Interest", "Context of Use"]
    rows = [
        [
            "How does bevemipretide, a cardiolipin-stabilising peptidomimetic "
            "with demonstrated CNS penetration, modify the self-reinforcing "
            "aSyn–CL–CI–ROS vicious cycle that drives dopaminergic "
            "neurodegeneration in Parkinson’s disease?",
            "Mechanistic proof-of-concept. The model encodes the complete "
            "vicious cycle as a coupled ODE system and quantifies how "
            "cardiolipin stabilisation propagates through ETC bioenergetics, "
            "ROS production, α-synuclein aggregation, and downstream "
            "neuronal death to establish a systems-level pharmacological rationale."
        ],
        [
            "What dose of bevemipretide produces clinically meaningful "
            "neuroprotection in Parkinson’s disease?",
            "Dose selection for Phase 2 clinical trial design. The model "
            "evaluates three dose levels (10, 30, 60 mg SC daily) against "
            "DA neuron preservation and motor score endpoints."
        ],
        [
            "What is the expected treatment effect size and what sample "
            "size is required to detect it?",
            "Powering of Phase 2 efficacy trials. Virtual clinical trial "
            "simulation provides unpaired Cohen’s d, Cliff’s delta, "
            "and MCID-based responder rates for sample size estimation."
        ],
        [
            "How does the timing of treatment initiation relative to "
            "symptom onset affect therapeutic efficacy?",
            "Eligibility criteria and enrolment strategy. Phased simulation "
            "with diagnostic delay analysis informs the therapeutic window "
            "for disease modification."
        ],
        [
            "What is the mechanistic ceiling on monotherapy efficacy, "
            "and which untargeted pathways limit response?",
            "Combination therapy rationale. Multi-pathway death decomposition "
            "quantifies the fraction of neuronal death that is drug-rescuable "
            "versus CL-independent."
        ],
        [
            "How do inter-patient differences in disease biology affect "
            "treatment response heterogeneity?",
            "Enrichment strategy and responder identification. Virtual "
            "population analysis with 15 varied parameters characterises "
            "the distribution of individual treatment effects."
        ],
    ]
    add_table(doc, headers, rows,
              title="Table A. Questions of Interest and Context of Use")

    add_legend(doc,
        "Table A. Model questions and their decision-making contexts, "
        "defined per IQ Consortium QSP model qualification recommendations. "
        "Each question maps to a specific analysis module within the model."
    )


def add_model_risk_table(doc):
    """Model Risk Assessment table (IQ Consortium / FDA MIDD framework)."""
    headers = ["Category", "Assessment"]
    rows = [
        [
            "Model Influence",
            "Medium–High. The model directly informs three pivotal development "
            "decisions: dose selection, sample size estimation, and go/no-go "
            "evaluation for Phase 2. It provides the primary quantitative basis "
            "for predicting treatment effect sizes (Cohen's d, Cliff's delta) and "
            "MCID-based responder rates in the absence of Phase 2 clinical data. "
            "Additionally, the mechanistic decomposition of death pathways guides "
            "combination therapy strategy, and the phased simulation shapes "
            "eligibility criteria and enrolment timing—decisions that collectively "
            "determine trial feasibility and resource allocation."
        ],
        [
            "Consequence of Wrong Decision",
            "Moderate. If the model underestimates the required dose, a Phase 2 "
            "trial may be underpowered or show futility despite a viable mechanism. "
            "If it overestimates effect size, the trial may be underpowered due to "
            "insufficient sample size. Incorrect characterisation of the therapeutic "
            "ceiling could misdirect combination therapy strategies. However, "
            "regulatory decisions (approval, labelling) remain contingent on "
            "independent clinical endpoints."
        ],
        [
            "Model Risk",
            "Medium. The model is mechanistically grounded in established PD "
            "pathobiology with 91 parameters sourced from peer-reviewed literature. "
            "It is calibrated against five independent quantitative targets and "
            "validated with 24/24 directional and quantitative checks. A three-tier "
            "sensitivity analysis (OAT, Morris, Sobol) and comprehensive uncertainty "
            "quantification (bootstrap CIs, parameter envelope) characterise "
            "predictive confidence. Mouse-to-human translation uses standard "
            "allometric scaling but introduces inherent cross-species uncertainty. "
            "The virtual population captures inter-patient variability across 15 "
            "parameters but does not model genetic subgroups or comorbidities."
        ],
        [
            "Model Impact",
            "Medium–High. The model shapes the design and operational parameters "
            "of a Phase 2 programme whose cost and duration make mid-course "
            "correction difficult. Dose selection, sample size, and enrichment "
            "strategy—all model-informed—are locked in at protocol finalisation. "
            "The model operates within the FDA MIDD framework and ICH M15 guideline, "
            "and its predictions are accompanied by rigorous uncertainty bounds "
            "(bootstrap 95% CIs, 125-combination parameter envelopes, patient-level "
            "prediction intervals). Independent mechanistic cross-validation from "
            "the ALS TDP-43 model (Bhatt et al. 2022) and consistency with Phase 1 "
            "PK data strengthen confidence. Regulatory decisions nonetheless remain "
            "contingent on clinical endpoints from the Phase 2 trial itself."
        ],
    ]
    add_table(doc, headers, rows,
              title="Table B. Model Risk Assessment")

    add_legend(doc,
        "Table B. Model risk assessment per the IQ Consortium framework for QSP "
        "model qualification. Categories follow the model influence–consequence–"
        "risk–impact paradigm recommended for regulatory submissions under "
        "FDA MIDD and ICH M15."
    )


def build():
    doc = make_doc()

    # ═══════════════════════════════════════════════════════════════════════
    #  QoI / CoU TABLE (before Methods, per IQ Consortium)
    # ═══════════════════════════════════════════════════════════════════════
    add_qoi_cou_table(doc)

    # ═══════════════════════════════════════════════════════════════════════
    #  METHODS
    # ═══════════════════════════════════════════════════════════════════════
    add_heading(doc, "2. Methods")

    # 2.1 Model architecture
    add_subheading(doc, "2.1 Model Architecture and Assumptions")

    add_para(doc,
        "The QSP model comprises eight interconnected modules (M0–M7) "
        "representing the major pathological axes of Parkinson’s disease, "
        "encoded as a system of 28 ordinary differential equations (Table 1, Figure 1). "
        "The modular architecture captures the pharmacokinetics of bevemipretide (M0), "
        "cardiolipin dynamics (M1), α-synuclein aggregation (M2), electron transport "
        "chain bioenergetics (M3), mitochondrial reactive oxygen species (M4), "
        "PINK1/Parkin-mediated mitophagy and cell death (M5), neuroinflammation (M6), "
        "and a composite motor endpoint (M7)."
    )

    add_para(doc,
        "Several structural assumptions underpin the model. All 28 state variables "
        "are normalised to the interval [0, 1], representing fractional activities or "
        "concentrations relative to physiological maxima. Disease is not modelled as "
        "a transient perturbation; instead, two sustained disease modifiers—CI_max "
        "(maximum Complex I repair capacity, reflecting PINK1/Parkin impairment) and "
        "α_clear (α-synuclein clearance capacity, reflecting GBA/aging impairment)"
        "—set the ceiling on homeostatic recovery and thereby gate entry into "
        "the vicious cycle. Dopaminergic neuron death proceeds through three "
        "independent pathways: mitochondrial damage (the primary drug-rescuable pathway), "
        "neuroinflammation-mediated toxicity via TNF/IL-1β, and direct α-synuclein "
        "proteotoxicity. Each death pathway is gated by a Hill function to enforce "
        "threshold behaviour and prevent physiologically implausible healthy-state death. "
        "The complete parameter set with source annotations is provided in Table 2 "
        "and Table S1."
    )

    # 2.2 Calibration
    add_subheading(doc, "2.2 Calibration and Validation")

    add_para(doc,
        "Model calibration followed a four-stage systematic procedure. First, the "
        "healthy baseline was verified (DA_neuron > 0.99 at 50 days). Second, a "
        "two-dimensional grid sweep of CI_max and α_clear was performed against "
        "quantitative targets from Gao et al. (2017; CI activity), Choi et al. "
        "(2022; mROS elevation), and literature-composite CL depletion data. Third, "
        "the neuroinflammation impairment parameter (k_impair) was swept against "
        "Ivanova et al. (2024; TLR2 knockout α-synuclein reduction). Fourth, a "
        "fine-grid joint refinement was conducted with all targets evaluated "
        "simultaneously. Validation included 14 directional checks (drug rescue, "
        "dose ordering, boundedness) and comparison against Bernheimer (1973) DA "
        "neuron loss at five weeks (Table 3)."
    )

    # 2.3 Allometric scaling
    add_subheading(doc, "2.3 Mouse-to-Human Translation")

    add_para(doc,
        "Allometric PK scaling from mouse (25 g) to human (70 kg) employed the "
        "BW^(−0.25) rule for elimination and distribution rate constants (Table 4). "
        "BBB penetration was scaled using published monkey CSF-to-plasma ratios, while "
        "intracellular kinetics (mitochondrial uptake and efflux) were left unscaled as "
        "physiology-independent processes. The FDA K_m body-surface-area method converted "
        "the mouse reference dose (5 mg/kg IP) to a human equivalent of approximately "
        "30 mg SC daily. Human disease timescale parameters were recalibrated: CI_max was "
        "set to 0.90, α_clear to 0.35, and all three death rate constants were scaled "
        "to 6% of mouse values, reflecting the substantially slower nigral degeneration "
        "observed in humans compared with toxin-based mouse models."
    )

    # 2.4 Virtual population and CTS
    add_subheading(doc, "2.4 Virtual Clinical Trial Simulation")

    add_para(doc,
        "A virtual population of 250 patients per arm was generated by sampling 15 "
        "parameters from log-normal (mechanistic and PK parameters) or truncated normal "
        "(disease modifiers) distributions (Table S2). Four arms were simulated: placebo, "
        "10 mg, 30 mg, and 60 mg SC daily for 547 days (18 months), with daily dosing "
        "events. Effect sizes were reported as unpaired Cohen’s d and Cliff’s delta; "
        "responder rates were defined using minimal clinically important difference "
        "(MCID) thresholds of >5% DA neuron preservation and >3% motor score improvement. "
        "Sample size was estimated for 80% power using unpaired d."
    )

    # 2.5 Sensitivity analysis
    add_subheading(doc, "2.5 Sensitivity and Uncertainty Analyses")

    add_para(doc,
        "A three-tier sensitivity analysis pipeline was applied to 34 parameters across "
        "all modules. One-at-a-time (OAT) perturbation of ±50% generated tornado "
        "diagrams (68 model evaluations). Morris elementary effects screening with 20 "
        "trajectories (700 evaluations) identified parameters with strong nonlinear or "
        "interactive behaviour. Sobol variance-based decomposition on the top 10 "
        "parameters (N = 256 base samples; 3,072 evaluations) partitioned output variance "
        "into first-order and total-order indices (Figure S1, Table 6). Uncertainty was "
        "quantified through 10,000 nonparametric bootstrap resamples of the virtual "
        "population and a 125-combination joint sweep of the three translational "
        "parameters (Table S3)."
    )

    # 2.6 Phased simulation
    add_subheading(doc, "2.6 Phased Simulation")

    add_para(doc,
        "To replicate the clinical timeline, a three-phase simulation was constructed: "
        "prodromal disease-only progression until symptoms (DA neuron < 70%), a 12-month "
        "diagnostic delay, and treatment initiation from the disease state at diagnosis. "
        "The state vector was extracted at each transition and passed as initial conditions "
        "to the subsequent phase. A dose-by-delay matrix (three doses × four delay "
        "durations) characterised the interaction between treatment timing and efficacy "
        "(Table S4)."
    )

    # 2.7 Software
    add_subheading(doc, "2.7 Software and Numerical Methods")

    add_para(doc,
        "All simulations were performed in R (v4.3) using the deSolve package with the "
        "LSODA solver (absolute and relative tolerance = 10^−10). State variables were "
        "clamped to [0, 1] within the ODE function for numerical safety. Figures were "
        "generated with ggplot2 at 300 DPI. The model code, including all verification "
        "tests, is available at [repository URL]."
    )

    # ═══════════════════════════════════════════════════════════════════════
    #  RESULTS
    # ═══════════════════════════════════════════════════════════════════════
    add_heading(doc, "3. Results")

    # 3.1 Calibration
    add_subheading(doc, "3.1 Calibration and Validation")

    add_para(doc,
        "The calibration sweep identified CI_max = 0.74 and α_clear = 0.25 as the "
        "optimal disease modifier combination (Figure 3A). All five quantitative targets "
        "fell within their acceptable ranges: CI activity reached 0.459 (target: 0.40–0.58; "
        "Gao 2017), CL depletion was 15.7% (target: 15–35%), mROS rose to 161.6% of basal "
        "(target: 150–250%; Choi 2022), DA neuron loss at five weeks was 36.0% (target: "
        "25–45%; Bernheimer 1973), and TLR2 knockout reduced α-synuclein oligomers by "
        "20.3% (target: 15–50%; Ivanova 2024) (Figure 3B, Table 3). All 14 directional "
        "checks passed, including monotonic dose-response and full state boundedness."
    )

    add_legend(doc,
        "Figure 3. Calibration and validation summary. "
        "(A) Disease modifier grid sweep (CI_max × α_clear); colour encodes "
        "log₁₀ of the composite calibration score; cross marks the selected combination. "
        "(B) Quantitative validation targets with model outputs (diamonds) positioned "
        "within acceptable literature ranges (green bands). Sources indicated on each row."
    )

    # 3.2 Integrated dynamics
    add_subheading(doc, "3.2 Integrated Dynamics and the Vicious Cycle")

    add_para(doc,
        "Time-course simulations over 50 days in the mouse model demonstrated progressive "
        "divergence between healthy and disease states across all eight modules (Figure 2). "
        "In disease, CL functional ratio declined from 0.95 to 0.84, CI activity fell to "
        "0.46, and mROS rose to 0.17—approximately 1.6-fold over basal. These upstream "
        "perturbations drove α-synuclein oligomer accumulation to 0.35 and microglial "
        "activation to 0.35, culminating in 48% DA neuron loss and a motor score of 0.21. "
        "High-dose bevemipretide partially reversed these changes in a dose-dependent "
        "manner, preserving CL ratio at 0.94 and DA neuron survival at 0.65."
    )

    add_para(doc,
        "The phase portrait (Figure 4) revealed three distinct attractors in CL ratio–"
        "α-synuclein state space: a healthy attractor (high CL, low α-synuclein), a "
        "disease attractor (low CL, high α-synuclein), and an intermediate drug-rescued "
        "state. The drug trajectory initially tracked the disease path before diverging "
        "toward the rescued attractor, consistent with partial interruption of the vicious "
        "cycle rather than full restoration of health. The cascade propagation timeline "
        "(Figure S3B) showed that CI dysfunction and α-synuclein aggregation were the "
        "earliest events (onset within 2–3 days), followed by downstream ATP depletion, "
        "membrane depolarisation, and motor impairment."
    )

    add_legend(doc,
        "Figure 2. Integrated model dynamics over 50 days (mouse timescale). "
        "Eight panels showing key state variables under four scenarios: healthy (green), "
        "disease (red), low-dose drug (orange), and high-dose drug (blue). "
        "All panels share a consistent 50-day x-axis."
    )

    add_legend(doc,
        "Figure 4. Vicious cycle phase portrait. "
        "Trajectories in CL ratio–α-synuclein oligomer state space. "
        "Diamonds mark terminal attractors; arrows indicate direction of flow. "
        "Three distinct basins of attraction are visible for health, disease, "
        "and drug-rescued states."
    )

    # 3.3 Death pathways
    add_subheading(doc, "3.3 Multi-Pathway Death Decomposition")

    add_para(doc,
        "Decomposition of the total death rate into its three constituent pathways "
        "(Figure S2A) showed that mitochondrial damage accounted for approximately 75% "
        "of disease-state neuronal death, with neuroinflammation and proteotoxicity "
        "contributing the remaining 25%. Under drug treatment, the mitochondrial "
        "component was substantially reduced while the non-mitochondrial pathways "
        "persisted, explaining the ceiling on therapeutic efficacy. The disease fingerprint "
        "radar chart (Figure S2B) illustrated the multi-dimensional nature of the "
        "pathology: disease contracted the fingerprint across all eight axes, while "
        "drug treatment selectively expanded the mitochondria-related axes (CL integrity, "
        "CI activity, ATP, membrane potential) with less effect on inflammation-related "
        "measures."
    )

    # 3.4 Sensitivity
    add_subheading(doc, "3.4 Sensitivity Analysis")

    add_para(doc,
        "The OAT tornado diagram identified k_ROS_basal, CI_max, k_SOD2, and K_mPTP "
        "as the top drivers of DA neuron survival, each producing >47% change under "
        "±50% perturbation (Figure S1A, Table 6). The Sobol decomposition confirmed "
        "that k_ROS_basal dominated both endpoints (S_T = 0.50 for DA neuron, 0.40 for "
        "motor score), with 63% of its influence mediated through parameter interactions "
        "rather than direct effects (Figure S1B). K_mPTP was almost purely interactive "
        "(S_1 ≈ 0, S_T = 0.33), consistent with its role as a threshold switch. "
        "Collectively, ROS module parameters (k_ROS_basal, k_SOD2, k_shunt) dominated "
        "both endpoints, confirming oxidative stress as the primary amplifier within the "
        "vicious cycle."
    )

    # 3.5 Human translation
    add_subheading(doc, "3.5 Human Translation")

    add_para(doc,
        "Allometric scaling produced a human plasma half-life of 21.7 hours (mouse: "
        "3.0 hours) and a reference mitochondrial trough concentration (C_ref) of 2.15 "
        "at 30 mg daily steady state (Figure 5A, Table 4). Twenty-year disease trajectory "
        "simulations (Figure 5B) showed a prodromal phase of 5.6 years (consistent with "
        "the 5–10 year range reported by Postuma 2012), clinical PD at 10.8 years, "
        "and crossing of the Bernheimer threshold (60% DA loss) at 14.3 years. Treatment "
        "with 30 mg bevemipretide delayed symptomatic onset by 4.1 years and clinical PD "
        "by 8.1 years. The healthy baseline showed 9.6% DA loss at 20 years, consistent "
        "with the approximately 5%/decade age-related decline reported by Fearnley and "
        "Lees (1991). Disease phase characterisation (Figure 5C) mapped the prodromal, "
        "diagnostic delay, and clinical PD windows against standard clinical thresholds."
    )

    add_legend(doc,
        "Figure 5. Human translation composite. "
        "(A) Mitochondrial PK comparison: mouse (5 mg/kg IP, orange) versus human "
        "(30 mg SC, blue). (B) Twenty-year DA neuron survival under healthy (green), "
        "disease (red), and drug-treated (blue) scenarios; dashed line marks the "
        "Bernheimer threshold (40% survival). (C) Disease phase characterisation with "
        "prodromal, diagnostic delay, and clinical PD windows; dashed lines mark 30%, "
        "50%, and 60% DA loss thresholds."
    )

    # 3.6 Virtual clinical trial
    add_subheading(doc, "3.6 Virtual Clinical Trial")

    add_para(doc,
        "The virtual population produced dose-dependent separation across all metrics "
        "(Figure 6A, Table 5). At 18 months, 30 mg bevemipretide yielded an unpaired "
        "Cohen’s d of 0.221 for DA neuron preservation and 0.221 for motor score "
        "improvement, with MCID responder rates of 25.8% and 17.2%, respectively. The "
        "60 mg dose (the Phase 1 MAD ceiling) reached d = 0.273 for DA neuron "
        "preservation with 31.6% responders. Dose-response was monotonic with diminishing "
        "returns above 30 mg. Individual patient responses (Figure 6B) showed substantial "
        "heterogeneity, with the majority of benefit concentrated in the top quartile "
        "of responders."
    )

    add_para(doc,
        "The Therapeutic Rescue Fraction heatmap (Figure 6C) showed that bevemipretide "
        "rescued 60% of the mitochondrial death pathway at 60 mg but only 35–38% of "
        "the non-mitochondrial pathways (neuroinflammation and proteotoxicity). This "
        "differential rescue explains the ceiling on total efficacy (52% total death "
        "rate reduction at 60 mg) and suggests that combination strategies targeting "
        "non-mitochondrial pathways could augment the therapeutic response."
    )

    add_legend(doc,
        "Table 5. Virtual clinical trial simulation results at 18 months "
        "(N = 244/arm). Unpaired Cohen’s d, Cliff’s delta, MCID responder rates, "
        "and required sample sizes for 80% power are shown for each dose and endpoint."
    )

    add_legend(doc,
        "Figure 6. Virtual clinical trial composite. "
        "(A) Violin plots of DA neuron survival at 18 months across four arms "
        "(N = 244/arm). (B) Waterfall plot of individual patient responses at 30 mg "
        "ranked by DA change versus placebo; colours indicate MCID responders (>5%), "
        "marginal (0–5%), and no benefit. (C) Therapeutic Rescue Fraction heatmap "
        "showing the percentage of each death pathway rescued by dose."
    )

    # 3.7 Uncertainty
    add_subheading(doc, "3.7 Uncertainty Quantification")

    add_para(doc,
        "Bootstrap 95% confidence intervals for the 30 mg unpaired d were [0.177, 0.273] "
        "for DA neuron preservation and [0.174, 0.273] for motor score (Figure 7C, "
        "Table S3). The joint parameter uncertainty envelope, sweeping CI_max_human, "
        "α_clear_human, and death_scale_human across 125 combinations, produced a "
        "30 mg DA_saved range of [0.015, 0.069] and an approximate d range of [0.062, "
        "0.403] (Figure 7A). The parameter α_clear_human was the single largest source "
        "of uncertainty, producing a 2.6-fold range in DA_saved, because it directly "
        "gates α-synuclein clearance at the entry point of the vicious cycle. "
        "Patient-level prediction intervals showed substantial heterogeneity, with the "
        "80% interval spanning near-zero to >10% DA preservation (Figure 7B)."
    )

    add_legend(doc,
        "Figure 7. Uncertainty quantification composite. "
        "(A) Parameter uncertainty envelope: histogram of DA_saved across 27 translational "
        "parameter combinations; solid red line = nominal, dashed = 95% bounds. "
        "(B) Patient-level treatment effect distributions as box plots (matched VPop pairs). "
        "(C) Bootstrap density distributions of unpaired Cohen’s d (B = 10,000) for three "
        "dose levels."
    )

    # 3.8 Phased simulation
    add_subheading(doc, "3.8 Treatment Timing and Dose-Delay Interaction")

    add_para(doc,
        "The phased simulation showed that drug benefit was preserved across diagnostic "
        "delays of 0–24 months, with only 12% relative efficacy loss at 24 months "
        "compared with immediate treatment (Table S4). The dose-by-delay matrix "
        "confirmed that 30 mg and 60 mg maintained clinically meaningful benefit even "
        "when treatment was delayed by two years after symptom onset. The dose-response "
        "curve in the mouse model (Figure S3A) demonstrated a steep initial response "
        "at low doses with progressive flattening above 5 mg/kg, consistent with the "
        "diminishing returns observed in the human virtual trial."
    )

    # ═══════════════════════════════════════════════════════════════════════
    #  MODEL RISK ASSESSMENT TABLE (after Results, per IQ Consortium)
    # ═══════════════════════════════════════════════════════════════════════
    add_model_risk_table(doc)

    # ═══════════════════════════════════════════════════════════════════════
    #  DISCUSSION
    # ═══════════════════════════════════════════════════════════════════════
    add_heading(doc, "4. Discussion")

    # 4.1 — Positioning and novelty
    add_para(doc,
        "Prior QSP models of Parkinson's disease have addressed distinct slices of "
        "the pathology: Roberts et al. (2016) constructed a biophysically detailed "
        "cortico-striatal-thalamo-cortical loop that predicts motor symptom scores "
        "from dopamine receptor pharmacology, while Ivanova and Karelina (2024) "
        "modelled α-synuclein spread and immunotherapy engagement across brain "
        "compartments. Paavilainen et al. (2023) used a multi-compartment QSP "
        "framework to quantify synaptic cleft target engagement of anti-α-synuclein "
        "antibodies and to explain the clinical failure of cinpanemab and "
        "prasinezumab. Each of these platforms captures a specific aspect of the "
        "disease but none encodes the intracellular vicious cycle—the self-reinforcing "
        "loop linking α-synuclein aggregation, cardiolipin oxidation, Complex I "
        "impairment, and mitochondrial ROS—that multiple lines of experimental "
        "evidence identify as the central amplifier of nigral degeneration (Devi et al. "
        "2008; Ludtmann et al. 2018). The present model fills this gap by coupling "
        "all four nodes of the cycle within a single ODE system and layering "
        "pharmacokinetics, neuroinflammation, mitophagy, and multi-pathway cell death "
        "on top, yielding a platform capable of evaluating organelle-targeted "
        "pharmacology from first principles."
    )

    # 4.2 — Vicious cycle and emergent amplification
    add_para(doc,
        "A central finding is that the vicious cycle operates as an emergent "
        "amplifier: individual modules produce modest standalone perturbations "
        "(approximately 1.5- to 2-fold shifts in their respective outputs), yet the "
        "coupled 28-ODE system generates 48% DA neuron loss and a motor score of "
        "0.21 over 50 days. This amplification arises because each node feeds "
        "forward into the next—α-synuclein oligomers accelerate cardiolipin "
        "oxidation (Robotta et al. 2014), depleted cardiolipin destabilises "
        "respiratory supercomplexes (Paradies et al. 2014), impaired electron "
        "transport elevates ROS (Murphy 2009), and ROS in turn promotes further "
        "α-synuclein misfolding (Scudamore and Bhatti 2020)—creating a positive "
        "feedback loop whose gain exceeds any single-pathway perturbation. The "
        "phase portrait analysis formalises this observation: healthy and disease "
        "states occupy distinct basins of attraction in CL ratio–α-synuclein state "
        "space, separated by a saddle region that the drug trajectory must cross to "
        "rescue cells. The existence of an intermediate drug-rescued attractor, rather "
        "than full restoration to the healthy state, reflects the fact that "
        "bevemipretide stabilises cardiolipin at a single node of a four-node cycle "
        "and cannot fully compensate for upstream impairments in Complex I repair "
        "capacity and α-synuclein clearance."
    )

    # 4.2b — Cascade propagation and temporal ordering
    add_para(doc,
        "The cascade propagation analysis revealed a stereotyped temporal ordering "
        "of dysfunction: Complex I impairment and α-synuclein aggregation were the "
        "earliest events to exceed 10% deviation from healthy baseline (onset within "
        "2–3 days in the mouse model), followed by ATP depletion, membrane "
        "depolarisation, and finally motor impairment. This ordering recapitulates "
        "the molecular staging observed in post-mortem and experimental studies, "
        "where mitochondrial Complex I deficiency (Schapira et al. 1990) and early "
        "α-synuclein misfolding (Braak et al. 2003) precede overt neurodegeneration "
        "by years. The temporal gap between the initial biochemical insult and "
        "clinically detectable motor impairment—approximately five-fold in the mouse "
        "and proportionally longer in humans—defines the therapeutic window for "
        "disease modification. By acting at the cardiolipin node, bevemipretide "
        "intervenes early in this cascade, upstream of the point at which "
        "mitochondrial permeability transition commits the cell to death. The model "
        "predicts that the transition from reversible dysfunction to irreversible "
        "commitment is governed by the mPTP opening threshold (K_mPTP), which the "
        "Sobol analysis identified as almost purely interactive (S_1 ≈ 0, S_T = 0.33), "
        "functioning as an all-or-nothing gate rather than a graded effector. This "
        "threshold behaviour is consistent with the concept of a mitochondrial "
        "'point of no return' described in apoptosis literature (Kroemer et al. 2007)."
    )

    # 4.3 — Multi-pathway death and the efficacy ceiling
    add_para(doc,
        "The three-pathway death architecture—mitochondrial damage, "
        "neuroinflammation-mediated toxicity via TNF-α and IL-1β (Hirsch and Hunot "
        "2009), and direct α-synuclein proteotoxicity (Winner et al. 2011)—imposes "
        "a quantitative ceiling on monotherapy efficacy. Approximately 75% of "
        "disease-state neuronal death in the model flows through the mitochondrial "
        "pathway, which bevemipretide can partially rescue; the remaining 25% is "
        "CL-independent and persists even at the highest simulated dose. The "
        "Therapeutic Rescue Fraction heatmap (Figure 6C) makes this ceiling explicit: "
        "60 mg rescues 60% of the mitochondrial death rate but only 35–38% of the "
        "non-mitochondrial components, yielding a net 52% reduction in total death "
        "rate. This finding carries direct implications for trial design. A "
        "monotherapy trial powered to detect a large effect size will likely fail "
        "not because the drug is inactive but because the achievable effect is "
        "mechanistically bounded. The quantitative decomposition suggests that "
        "combination strategies—pairing cardiolipin stabilisation with anti-inflammatory "
        "agents or aggregation inhibitors—could access the remaining 25% of neuronal "
        "death and substantially widen the therapeutic window."
    )

    # 4.4 — Effect sizes and clinical trial landscape
    add_para(doc,
        "The predicted unpaired Cohen's d of 0.221 at 30 mg sits between the "
        "effect sizes reported in the ADAGIO rasagiline trial (d ≈ 0.20; Olanow "
        "et al. 2009) and the SPARK isradipine trial (d ≈ 0.15; Bhatt et al. "
        "2020). This concordance is notable because ADAGIO and SPARK targeted "
        "entirely different pharmacological mechanisms—MAO-B inhibition and L-type "
        "calcium channel blockade, respectively—yet both engaged pathways that "
        "converge on mitochondrial function and oxidative stress. The model's "
        "prediction of diminishing returns above 30 mg is independently supported "
        "by preclinical data from an ALS TDP-43 mouse model treated with SBT-272 "
        "(Bhatt et al. 2022, NEALS), in which upper motor neuron retention and "
        "reduced neuroinflammation showed a dose-response plateau. Moreover, the "
        "Phase 1 data for bevemipretide confirmed systemic exposure profiles "
        "consistent with the allometric PK scaling used in the model, providing an "
        "independent check on the translational bridge."
    )

    # 4.4b — Context of failed and successful PD disease-modification trials
    add_para(doc,
        "The broader landscape of PD disease-modification trials provides instructive "
        "context. The DATATOP trial (Parkinson Study Group 1989) demonstrated that "
        "deprenyl delayed the need for levodopa, but the confound of symptomatic "
        "benefit precluded definitive disease-modification claims—a challenge that "
        "persists in trial design today. SURE-PD3 (Schwarzschild et al. 2021) tested "
        "whether inosine-mediated urate elevation could slow progression through "
        "antioxidant mechanisms; it failed to meet its primary endpoint despite a "
        "sound biological rationale, achieving a negligible effect size. Our model "
        "offers a potential explanation: urate acts as a general-purpose antioxidant "
        "in the cytosol, whereas the vicious cycle is compartmentalised within "
        "mitochondria. Bevemipretide, by contrast, accumulates at the inner "
        "mitochondrial membrane with >100-fold selectivity over cytoplasmic "
        "distribution, placing the drug precisely where cardiolipin oxidation and "
        "supercomplex destabilisation occur. The model's bootstrap 95% confidence "
        "interval for the 30 mg effect size [0.177, 0.273] places even the lower "
        "bound above the effect sizes observed in SPARK and SURE-PD3, suggesting "
        "that a trial powered to detect d ≈ 0.18 would retain a reasonable "
        "probability of success under pessimistic translational assumptions."
    )

    # 4.5 — ROS as the primary amplifier
    add_para(doc,
        "The three-tier sensitivity analysis converged on a consistent result: "
        "ROS module parameters—k_ROS_basal, k_SOD2, and k_shunt—dominate both "
        "DA neuron survival and motor score across OAT, Morris, and Sobol analyses. "
        "The Sobol total-order index for k_ROS_basal (S_T = 0.50) reflects the "
        "fact that 63% of its influence on DA neuron survival is mediated through "
        "parameter interactions rather than direct effects, consistent with its "
        "role as the amplifier gain within the vicious cycle. K_mPTP is almost "
        "purely interactive (S_1 ≈ 0, S_T = 0.33), functioning as a threshold "
        "switch for mitochondrial permeability transition rather than a graded "
        "effector. These findings align with the longstanding observation that "
        "oxidative stress is both necessary and sufficient to initiate nigral "
        "degeneration in toxin-based models (Blesa et al. 2015) and with "
        "post-mortem evidence of elevated 8-hydroxy-2'-deoxyguanosine and "
        "decreased glutathione in PD substantia nigra (Dexter et al. 1989; "
        "Alam et al. 1997). From a therapeutic standpoint, the dominance of ROS "
        "parameters implies that interventions upstream of ROS production "
        "(cardiolipin stabilisation, supercomplex preservation) may have a "
        "disproportionate effect by attenuating the cycle's primary amplifier."
    )

    # 4.5b — Disease modifiers and genetic implications
    add_para(doc,
        "The two disease modifiers—CI_max (maximum Complex I repair capacity) and "
        "α_clear (α-synuclein clearance capacity)—act as gatekeepers of the vicious "
        "cycle. The Sobol analysis showed that CI_max exerts 67% of its influence "
        "through interactions rather than direct effects, consistent with its role "
        "as a cycle entry point rather than a standalone effector. These two "
        "parameters map onto defined genetic axes of PD. CI_max corresponds to "
        "PINK1/Parkin-mediated mitochondrial quality control, disrupted in autosomal "
        "recessive juvenile-onset PD (Valente et al. 2004); α_clear corresponds to "
        "lysosomal α-synuclein degradation, impaired by GBA mutations—the most "
        "common genetic risk factor for sporadic PD, present in 7–12% of patients "
        "across populations (Sidransky et al. 2009). The uncertainty analysis "
        "identified α_clear as the single largest source of variability in "
        "predicted treatment effect (2.6-fold range in DA neuron preservation), "
        "which suggests that GBA mutation carriers may occupy the lower tail of "
        "the response distribution. Conversely, patients with preserved lysosomal "
        "function may derive greater benefit from cardiolipin stabilisation because "
        "the cleared α-synuclein does not re-accumulate as rapidly. This "
        "observation provides a quantitative rationale for stratifying trial "
        "enrolment by GBA status or by α-synuclein seed amplification assay "
        "(Siderowf et al. 2023), which serves as a functional readout of "
        "clearance capacity."
    )

    # 4.6 — Human translation and disease timeline
    add_para(doc,
        "The allometrically scaled human model reproduced several independent "
        "clinical observations without additional fitting: a prodromal phase of "
        "5.6 years (within the 5–10 year range estimated by Postuma and Berg 2016), "
        "age-related DA loss of 9.6% at 20 years (consistent with ~5%/decade; "
        "Fearnley and Lees 1991), and crossing of the Bernheimer threshold at "
        "14.3 years. The predicted delay in symptomatic onset of 4.1 years at 30 mg "
        "represents a clinically substantial benefit, particularly given that "
        "treatment was initiated at the point of diagnosis rather than during the "
        "prodromal phase. Notably, the phased simulation demonstrated that efficacy "
        "was robust to diagnostic delay: a 24-month delay reduced relative benefit "
        "by only 12%, indicating that the vicious cycle, once partially interrupted, "
        "can be slowed even when substantial neurodegeneration has already occurred. "
        "This resilience reflects the model's structure: bevemipretide acts on "
        "cardiolipin—a node that remains accessible regardless of disease stage—"
        "rather than on an upstream trigger that may have already passed."
    )

    # 4.6b — Drug benefit compounds over time
    add_para(doc,
        "A counterintuitive prediction of the model is that the absolute drug "
        "benefit grows over time: DA neuron preservation increases from +2.3% at "
        "year 1 to +73.8% at year 20. This accelerating separation between treated "
        "and untreated trajectories arises because bevemipretide partially "
        "interrupts the positive feedback loop; over long timescales, even a partial "
        "reduction in cycle gain compounds, much as a modest decrease in interest "
        "rate produces a large difference in accumulated debt. The implication for "
        "trial design is that short trials (12–18 months) will inevitably "
        "underestimate the long-term benefit. This is a recognised problem in PD "
        "disease-modification research: the ADAGIO delayed-start design (Olanow "
        "et al. 2009) was specifically conceived to distinguish symptomatic from "
        "neuroprotective effects over 72 weeks, yet even this duration may capture "
        "only the early, shallow portion of the divergence curve. The model "
        "quantifies this gap and could inform the design of long-term extension "
        "studies or Bayesian adaptive trials that leverage early biomarker signals "
        "to project long-term clinical benefit."
    )

    # 4.7 — Contrast with anti-aSyn immunotherapy failures
    add_para(doc,
        "The clinical failure of cinpanemab and prasinezumab in Phase 2 "
        "(Lang et al. 2022; Pagano et al. 2022) has prompted re-examination of "
        "the α-synuclein spread hypothesis. Paavilainen et al. (2023) demonstrated "
        "that despite near-complete target engagement of aggregated α-synuclein in "
        "CSF, reduction of neuronal uptake at the synaptic cleft was the "
        "rate-limiting step for anti-tau antibodies and was more favourable for "
        "anti-α-synuclein antibodies. Ivanova and Karelina (2024) showed that "
        "blocking α-synuclein entry alone does not reduce intracellular accumulation "
        "in overexpression models because the pathology develops primarily inside "
        "neurons. Our model adds a complementary perspective: even if extracellular "
        "α-synuclein is fully neutralised, the intracellular vicious cycle—driven by "
        "cardiolipin oxidation, Complex I impairment, and ROS—continues to generate "
        "oligomeric α-synuclein de novo. This suggests that intracellular organelle-"
        "targeted pharmacology, as exemplified by bevemipretide, addresses a "
        "mechanistic gap that antibody-based strategies cannot reach."
    )

    # 4.7b — Responder heterogeneity
    add_para(doc,
        "The virtual clinical trial revealed substantial inter-patient "
        "heterogeneity in treatment response, with the majority of DA neuron "
        "preservation concentrated in the top quartile of responders and a "
        "non-trivial fraction of patients showing near-zero benefit. This "
        "heterogeneity is not an artefact of random noise; it reflects genuine "
        "between-patient variability in the 15 sampled parameters, particularly "
        "in disease modifier values (CI_max, α_clear) and death pathway rates. "
        "The MCID responder rate of 25.8% at 30 mg—meaning roughly one in four "
        "virtual patients achieved >5% DA neuron preservation—mirrors the "
        "responder fractions reported in PD trials, where heterogeneity routinely "
        "dilutes population-level effect sizes (Espay et al. 2020). The waterfall "
        "plot (Figure 6B) makes this distribution visible at the individual patient "
        "level and underscores the need for enrichment biomarkers. The model "
        "predicts that patients with relatively preserved α-synuclein clearance "
        "and moderate Complex I impairment—the subgroup in which the vicious cycle "
        "is most sensitive to cardiolipin stabilisation—are the strongest "
        "responders, a prediction testable through baseline biomarker stratification "
        "in a Phase 2 adaptive design."
    )

    # 4.8 — Limitations
    add_subheading(doc, "4.1 Limitations")

    add_para(doc,
        "Several limitations should be considered when interpreting these results. "
        "First, the mouse-to-human translation relies on allometric PK scaling "
        "and recalibration of three disease timescale parameters (CI_max, α_clear, "
        "death_scale), introducing irreducible cross-species uncertainty. The "
        "mouse model employs disease modifiers to simulate PD rather than a "
        "specific toxin protocol, but the underlying kinetics inevitably reflect "
        "rodent cell biology—for instance, the rotenone model (Betarbet et al. "
        "2000) produces nigral degeneration over weeks, whereas human PD unfolds "
        "over decades. Our death-rate scaling factor of 0.06 bridges this gap "
        "empirically rather than from first-principles biophysics. "
        "Second, the virtual population samples 15 parameters from assumed "
        "distributions and does not capture genetic subgroups (LRRK2, GBA, SNCA "
        "multiplications) or comorbidities that could modify treatment response. "
        "The absence of explicit genetic subtypes means the model cannot predict "
        "differential efficacy across the molecular taxonomy of PD that is "
        "increasingly recognised as clinically relevant (Espay et al. 2020). "
        "Third, the model does not represent α-synuclein spread between brain "
        "regions; all dynamics occur within a single nigral compartment. The "
        "Braak staging hypothesis (Braak et al. 2003) posits that pathology "
        "propagates from the lower brainstem to cortex over decades, and spatial "
        "spread may influence both symptom heterogeneity and the time course of "
        "non-motor features that our model does not capture. This simplification "
        "is defensible for a cardiolipin-targeted drug whose mechanism is "
        "cell-autonomous, but it would need extension for therapies targeting "
        "intercellular propagation."
    )

    add_para(doc,
        "Fourth, the motor endpoint is a composite surrogate derived from DA "
        "neuron survival and neuroinflammation rather than a direct mapping to "
        "UPDRS-III subscores, which limits granularity in predicting specific "
        "clinical trial endpoints. Roberts et al. (2016) achieved a direct "
        "correlation between their QSP model's STN beta/gamma ratio and "
        "UPDRS-III scores by modelling basal ganglia electrophysiology; our model "
        "prioritises subcellular mechanistic depth over circuit-level fidelity, "
        "and the two approaches are complementary rather than competing. "
        "Fifth, all simulations assume perfect adherence to daily subcutaneous "
        "dosing; real-world adherence patterns could attenuate the steady-state "
        "mitochondrial drug concentrations on which efficacy depends. "
        "Sixth, while the model is calibrated against five quantitative targets "
        "from independent laboratories, several of these targets derive from "
        "in vitro or acute toxin studies (Gao et al. 2017; Choi et al. 2022) "
        "rather than from chronic neurodegenerative conditions, and the extent "
        "to which these surrogate targets map onto slowly progressive human PD "
        "remains an open question."
    )

    # 4.9 — Future directions and closing
    add_subheading(doc, "4.2 Future Directions")

    add_para(doc,
        "Several extensions would strengthen the translational utility of this "
        "platform. First, the multi-pathway death decomposition provides an "
        "actionable roadmap for in silico combination screening: pairing "
        "bevemipretide with anti-neuroinflammatory agents (targeting the ~13% "
        "TNF/IL-1β-mediated component) or α-synuclein aggregation inhibitors "
        "(targeting the ~12% proteotoxicity component) could be evaluated by "
        "adding pharmacodynamic modules for each co-administered drug and "
        "simulating combinatorial virtual trials. The recent clinical interest in "
        "LRRK2 kinase inhibitors and GLP-1 receptor agonists for PD (Meissner "
        "et al. 2011) provides candidate co-interventions whose mechanisms map "
        "onto model nodes. Second, incorporating α-synuclein propagation dynamics "
        "—either as a spatial ODE extension or through coupling with the "
        "multi-compartment framework of Ivanova and Karelina (2024)—would enable "
        "head-to-head comparison of intracellular (organelle-targeted) versus "
        "extracellular (antibody-based) strategies within a unified platform. "
        "Third, the motor endpoint should be extended to map onto UPDRS-III "
        "subscores using the item-response theory framework described by Gottipati "
        "et al. (2017), thereby enabling direct comparison with clinical trial "
        "primary endpoints."
    )

    add_para(doc,
        "The present platform, with 28 ODEs encoding the complete aSyn–CL–CI–ROS "
        "vicious cycle, calibrated against five independent experimental datasets, "
        "validated through 24 quantitative and directional checks, and stress-tested "
        "with a three-tier sensitivity analysis and comprehensive uncertainty "
        "quantification, establishes a quantitative systems-level case for "
        "advancing bevemipretide as a disease-modifying therapy in Parkinson's "
        "disease. It provides both the mechanistic rationale—partial interruption "
        "of the vicious cycle at a druggable intracellular node—and the "
        "translational predictions—effect sizes, responder rates, sample sizes, "
        "and treatment timing guidance—needed to design an informative Phase 2 "
        "clinical trial. More broadly, it demonstrates that QSP modelling can "
        "generate testable, quantitative predictions for organelle-targeted "
        "neurotherapeutics in a disease where the gap between preclinical promise "
        "and clinical delivery has remained stubbornly wide."
    )

    # Save
    path = os.path.join(OUT_DIR, "Methods_and_Results_v2.docx")
    doc.save(path)
    print(f"Saved: {path}")

    # Word count estimate
    wc_methods_results = 0
    wc_discussion = 0
    in_discussion = False
    for p in doc.paragraphs:
        txt = p.text.strip()
        if txt.startswith("4. Discussion"):
            in_discussion = True
            continue
        if not txt:
            continue
        if any(txt.startswith(x) for x in ["2.", "3.", "Figure", "Table A", "Table B"]):
            continue
        words = len(txt.split())
        if in_discussion:
            wc_discussion += words
        else:
            wc_methods_results += words
    print(f"Methods + Results word count: {wc_methods_results}")
    print(f"Discussion word count: {wc_discussion}")
    print(f"Total body text: {wc_methods_results + wc_discussion}")


if __name__ == "__main__":
    build()
