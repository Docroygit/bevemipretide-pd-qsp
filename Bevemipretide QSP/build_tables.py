"""
Build two DOCX files with publication-ready tables for Bevemipretide QSP manuscript.
  1. Main_Tables.docx   — Tables 1–6
  2. Supplementary_Tables.docx — Tables S1–S4

Run: conda run python build_tables.py
"""

from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import os

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "New publication plots")
os.makedirs(OUT_DIR, exist_ok=True)


# ── Helpers ──────────────────────────────────────────────────────────────────

def make_doc():
    doc = Document()
    style = doc.styles["Normal"]
    font = style.font
    font.name = "Times New Roman"
    font.size = Pt(10)
    style.paragraph_format.space_after = Pt(2)
    style.paragraph_format.space_before = Pt(0)
    style.paragraph_format.line_spacing = 1.15
    return doc


def add_table_title(doc, title):
    p = doc.add_paragraph()
    run = p.add_run(title)
    run.bold = True
    run.font.size = Pt(11)
    run.font.name = "Times New Roman"
    p.paragraph_format.space_before = Pt(12)
    p.paragraph_format.space_after = Pt(4)


def add_table_note(doc, note):
    p = doc.add_paragraph()
    run = p.add_run(note)
    run.font.size = Pt(9)
    run.font.name = "Times New Roman"
    run.italic = True
    p.paragraph_format.space_before = Pt(2)
    p.paragraph_format.space_after = Pt(12)


def build_table(doc, headers, rows, col_widths=None):
    """Create a formatted table with header row and data rows."""
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    # Header
    for j, h in enumerate(headers):
        cell = table.rows[0].cells[j]
        cell.text = ""
        p = cell.paragraphs[0]
        run = p.add_run(h)
        run.bold = True
        run.font.size = Pt(9)
        run.font.name = "Times New Roman"
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        # Light grey shading
        shading = cell._element.get_or_add_tcPr()
        shd = shading.makeelement(qn("w:shd"), {
            qn("w:fill"): "D9E2F3", qn("w:val"): "clear"
        })
        shading.append(shd)

    # Data rows
    for i, row_data in enumerate(rows):
        for j, val in enumerate(row_data):
            cell = table.rows[i + 1].cells[j]
            cell.text = ""
            p = cell.paragraphs[0]
            run = p.add_run(str(val))
            run.font.size = Pt(9)
            run.font.name = "Times New Roman"
            # Left-align first column, center others
            if j == 0:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            else:
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # Column widths
    if col_widths:
        for j, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[j].width = Cm(w)

    return table


def page_break(doc):
    doc.add_page_break()


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN MANUSCRIPT TABLES
# ═══════════════════════════════════════════════════════════════════════════════

def build_main_tables():
    doc = make_doc()

    # ── TABLE 1: Model Architecture ─────────────────────────────────────────
    add_table_title(doc,
        "Table 1. QSP model architecture: modules, state variables, and inter-module coupling.")

    headers = ["Module", "Description", "State Variables", "No. ODEs", "Key Outputs"]
    rows = [
        ["M0", "Pharmacokinetics",
         "Depot, C_plasma, C_periph, C_brain, C_mito", "5",
         "C_mito → M1"],
        ["M1", "Cardiolipin dynamics",
         "CL_n, CL_ox, CL_ext, ALCAT1, TAZ", "5",
         "CL_ratio → M3; CL_ext → M2, M5"],
        ["M2", "α-Synuclein aggregation",
         "αSyn_mono, αSyn_olig, αSyn_ext", "3",
         "αSyn_olig → M1, M3; αSyn_ext → M6"],
        ["M3", "ETC bioenergetics",
         "CI_activity, SC_integrity, Δψ_m, ATP, mPTP_open", "5",
         "CI, SC → M4; ATP → M2; mPTP → M5"],
        ["M4", "Oxidative stress",
         "mROS, SOD2_act, GPx4_act", "3",
         "mROS → M1, M2, M3"],
        ["M5", "Mitophagy / Cell death",
         "PINK1_act, mito_damage, DA_neuron", "3",
         "DA_neuron → M7"],
        ["M6", "Neuroinflammation",
         "MG_active, TNF, IL1β", "3",
         "MG_active → M2; TNF, IL1β → M5"],
        ["M7", "Motor endpoint",
         "Motor_score", "1",
         "Motor_score (clinical endpoint)"],
    ]
    build_table(doc, headers, rows, col_widths=[1.5, 3.0, 5.5, 1.5, 5.0])
    add_table_note(doc,
        "Total: 28 ODEs, 28 state variables. All state variables normalised to [0, 1]. "
        "Time unit: hours. Solver: deSolve::ode() with lsoda, atol/rtol = 10⁻¹⁰. "
        "The vicious cycle (M1 → M3 → M4 → M2 → M1) is the core disease amplification mechanism.")

    page_break(doc)

    # ── TABLE 2: Key Parameters ─────────────────────────────────────────────
    add_table_title(doc,
        "Table 2. Key model parameters with values, biological roles, and sources.")

    headers = ["Parameter", "Value", "Module", "Biological Role", "Source"]
    rows = [
        # Disease modifiers
        ["CI_max", "0.74", "Disease", "Max CI repair target (PINK1/Parkin impairment)", "Estimated"],
        ["α_clear", "0.25", "Disease", "aSyn clearance capacity (GBA/aging)", "Estimated"],
        ["k_impair", "0.35", "M2→M6", "MG-driven clearance impairment", "Estimated"],
        # M0: PK
        ["k_a", "1.5 h⁻¹", "M0", "SC absorption rate", "Estimated (IP bolus)"],
        ["k_e,plasma", "0.231 h⁻¹", "M0", "Plasma elimination (t½ = 3.0 h)", "Literature [1]"],
        ["k_mito_in", "0.50 h⁻¹", "M0", "Mitochondrial uptake", "Literature [2]"],
        ["C_ref", "4.43", "M0", "Reference mitochondrial concentration", "Estimated (5 mg/kg SS)"],
        # M1: CL
        ["k_ox", "0.020 h⁻¹", "M1", "ROS-driven CL oxidation", "Literature [3]"],
        ["k_drug", "0.100 h⁻¹", "M1", "Drug-mediated CL_ox reduction", "Estimated"],
        ["k_red", "0.070 h⁻¹", "M1", "Tafazzin-mediated CL remodelling", "Literature [3]"],
        # M2: aSyn
        ["k_agg", "0.009 h⁻¹", "M2", "Basal aSyn aggregation rate", "Literature [4]"],
        ["k_clear", "0.241 h⁻¹", "M2", "Oligomer clearance (autophagy/UPS)", "Literature [5]"],
        ["k_ROS_agg", "2.0", "M2", "ROS-enhanced aggregation sensitivity", "Estimated"],
        # M3: ETC
        ["k_CI_repair", "0.050 h⁻¹", "M3", "Complex I repair rate", "Literature [6]"],
        ["K_mPTP", "0.40", "M3", "mPTP opening threshold (Hill K)", "Estimated"],
        ["k_ATP_syn", "1.20 h⁻¹", "M3", "ATP synthesis rate", "Literature [7]"],
        # M4: ROS
        ["k_ROS_basal", "0.122 h⁻¹", "M4", "Basal mitochondrial ROS production", "Estimated"],
        ["k_SOD2", "1.0 h⁻¹", "M4", "SOD2 dismutation rate", "Literature [8]"],
        ["k_shunt", "0.125 h⁻¹", "M4", "CI/SC-dependent electron leak", "Estimated"],
        # M5: Death
        ["k_death", "7×10⁻⁴ h⁻¹", "M5", "Mitochondrial death rate", "Estimated"],
        ["k_death_inflam", "2×10⁻⁴ h⁻¹", "M5", "Neuroinflammation-mediated death", "Literature [9]"],
        ["k_death_aSyn", "1.5×10⁻⁴ h⁻¹", "M5", "Direct aSyn proteotoxicity", "Literature [10]"],
        ["K_death", "0.25", "M5", "Death threshold (Hill K for mito_damage)", "Estimated"],
        ["n_death", "4", "M5", "Hill cooperativity for death", "Estimated"],
        # M6
        ["k_act_MG", "20.0 h⁻¹", "M6", "Microglial TLR2-mediated activation", "Literature [11]"],
        # M7
        ["K_motor", "0.75", "M7", "DA neuron threshold for motor deficit", "Literature [12]"],
    ]
    build_table(doc, headers, rows, col_widths=[2.2, 2.2, 1.5, 5.5, 3.5])
    add_table_note(doc,
        "Sources: [1] SBT-272 Phase 1; [2] Szeto 2014 (PMID: 24117165); "
        "[3] Schlame 2008 (PMID: 18077827); [4] Auluck 2010 (PMID: 20500090); "
        "[5] Webb 2003 (PMID: 12719433); [6] Morais 2009 (PMID: 20049710); "
        "[7] Brand 2005 (PMID: 16246006); [8] Fukai 2011 (PMID: 21473702); "
        "[9] Hirsch & Hunot 2009 (PMID: 19296921); [10] Winner 2011 (PMID: 21325059); "
        "[11] Kim 2013 (PMID: 23463005); [12] Bernheimer 1973 (PMID: 4272516). "
        "'Estimated' = fitted to quantitative calibration targets (see Table 3).")

    page_break(doc)

    # ── TABLE 3: Calibration & Validation ───────────────────────────────────
    add_table_title(doc,
        "Table 3. Calibration and validation targets with model outputs.")

    headers = ["Target", "Endpoint", "Literature Value", "Acceptable Range",
               "Model Output", "Source", "Status"]
    rows = [
        ["CI activity", "CI_activity", "~49% of WT", "0.40–0.58",
         "0.459", "Gao et al. 2017", "PASS"],
        ["CL drop", "1 − CL_n/CL_n,healthy", "15–35%", "15–35%",
         "15.7%", "Literature composite", "PASS"],
        ["mROS elevation", "mROS/mROS_basal", "~177% of basal", "150–250%",
         "161.6%", "Choi et al. 2022", "PASS"],
        ["DA neuron loss", "1 − DA_neuron", "~36% at 5 weeks", "25–45%",
         "36.0%", "Bernheimer 1973", "PASS"],
        ["TLR2 KO aSyn", "ΔaSyn_olig (k_impair=0)", "23–45% reduction", "15–50%",
         "20.3%", "Ivanova et al. 2024", "PASS"],
        ["Motor score", "Motor_score at 5 weeks", ">10%", ">10%",
         "15.6%", "Bernheimer 1973", "PASS"],
        ["Drug rescue DA", "DA_drug > DA_disease", "Directional", "Positive",
         "+12.9%", "Mechanistic", "PASS"],
        ["Drug rescue CI", "CI_drug > CI_disease", "Directional", "Positive",
         "Correct", "Mechanistic", "PASS"],
        ["Drug rescue aSyn", "aSyn_drug < aSyn_disease", "Directional", "Negative",
         "Correct", "Mechanistic", "PASS"],
        ["Drug rescue mROS", "mROS_drug < mROS_disease", "Directional", "Negative",
         "Correct", "Mechanistic", "PASS"],
        ["Drug rescue Motor", "Motor_drug < Motor_disease", "Directional", "Negative",
         "Correct", "Mechanistic", "PASS"],
        ["Dose ordering", "Monotonic dose-response", "Directional", "Monotonic",
         "Correct", "Mechanistic", "PASS"],
        ["Boundedness", "All 28 states in [0,1]", "Physical constraint", "[0, 1]",
         "All pass", "Model property", "PASS"],
        ["Healthy baseline", "DA_neuron > 0.99 at 50d", "Physiological", ">0.99",
         "0.994", "Model property", "PASS"],
    ]
    build_table(doc, headers, rows, col_widths=[2.5, 3.0, 2.5, 2.0, 2.0, 2.5, 1.2])
    add_table_note(doc,
        "Calibration performed via 4-stage systematic sweep: (1) healthy baseline, "
        "(2) CI_max × α_clear disease modifier grid, (3) TLR2 KO k_impair sweep, "
        "(4) joint fine-grid refinement. 24/24 integrated validation tests passing.")

    page_break(doc)

    # ── TABLE 4: Allometric PK Scaling ──────────────────────────────────────
    add_table_title(doc,
        "Table 4. Allometric PK scaling from mouse to human.")

    headers = ["Parameter", "Mouse", "Human", "Scaling Method",
               "t½ Mouse", "t½ Human"]
    rows = [
        ["k_a (absorption)", "1.50 h⁻¹", "0.40 h⁻¹",
         "Estimated (SC vs IP)", "0.5 h", "1.7 h"],
        ["k_e,plasma (elimination)", "0.231 h⁻¹", "0.032 h⁻¹",
         "BW⁻⁰·²⁵ allometry", "3.0 h", "21.7 h"],
        ["k_12 (distribution)", "0.30 h⁻¹", "0.041 h⁻¹",
         "BW⁻⁰·²⁵ allometry", "—", "—"],
        ["k_21 (redistribution)", "0.20 h⁻¹", "0.027 h⁻¹",
         "BW⁻⁰·²⁵ allometry", "—", "—"],
        ["k_brain_in (BBB entry)", "0.040 h⁻¹", "0.015 h⁻¹",
         "Monkey CSF:plasma ratio", "17.3 h", "46.2 h"],
        ["k_brain_out (BBB exit)", "0.016 h⁻¹", "0.010 h⁻¹",
         "BW⁻⁰·²⁵ allometry", "—", "—"],
        ["k_mito_in (mito uptake)", "0.50 h⁻¹", "0.50 h⁻¹",
         "Unscaled (physiology-independent)", "—", "—"],
        ["k_mito_out (mito efflux)", "0.020 h⁻¹", "0.020 h⁻¹",
         "Unscaled (physiology-independent)", "—", "—"],
        ["C_ref", "4.43", "2.15",
         "Computed (120-day SS at ref dose)", "—", "—"],
        ["Reference dose", "5 mg/kg IP", "30 mg SC daily",
         "FDA K_m conversion (mouse → human)", "—", "—"],
    ]
    build_table(doc, headers, rows, col_widths=[3.5, 2.2, 2.2, 3.8, 1.8, 1.8])
    add_table_note(doc,
        "Allometric scaling: BW⁻⁰·²⁵ for elimination/distribution rate constants "
        "(mouse 25 g → human 70 kg; scale factor = 7.36). "
        "FDA K_m conversion: mouse 5 mg/kg × 3 = 15 mg/m²; human = 15 × 37/70 ≈ 28.4 mg, "
        "rounded to 30 mg SC daily. Clinical dose levels: 10, 30, 60 mg SC daily (60 mg = Phase 1 MAD ceiling).")

    page_break(doc)

    # ── TABLE 5: Clinical Trial Simulation ──────────────────────────────────
    add_table_title(doc,
        "Table 5. Virtual clinical trial simulation results at 18 months (N = 244/arm).")

    headers = ["Dose", "Endpoint", "Mean Δ\nvs placebo",
               "d_unpaired", "Cliff’s δ", "MCID\nResp. %",
               "N/arm\n(80% power)"]
    rows = [
        ["10 mg", "DA neuron", "+0.027", "0.133", "+0.074", "16.4%", "889"],
        ["30 mg", "DA neuron", "+0.044", "0.221", "+0.126", "25.8%", "322"],
        ["60 mg", "DA neuron", "+0.054", "0.273", "+0.164", "31.6%", "210"],
        ["10 mg", "Motor score", "−0.013", "0.135", "+0.107", "9.0%", "862"],
        ["30 mg", "Motor score", "−0.022", "0.221", "+0.180", "17.2%", "322"],
        ["60 mg", "Motor score", "−0.028", "0.295", "+0.226", "22.1%", "181"],
    ]
    build_table(doc, headers, rows, col_widths=[1.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
    add_table_note(doc,
        "Virtual population: N = 250/arm, 15 parameters varied (log-normal), 97.6% completion rate. "
        "MCID thresholds: DA neuron >5% preservation, Motor score >3% improvement. "
        "d_unpaired = unpaired Cohen’s d (comparable to clinical trial effect sizes). "
        "Dose-response is monotonic with diminishing returns above 30 mg. "
        "60 mg matches Phase 1 MAD ceiling dose. "
        "For context: ADAGIO (rasagiline) d ≈ 0.20; SPARK (isradipine) d ≈ 0.15.")

    page_break(doc)

    # ── TABLE 6: Sensitivity Analysis Summary ───────────────────────────────
    add_table_title(doc,
        "Table 6. Sensitivity analysis: top 10 parameters ranked by three methods.")

    headers = ["Rank", "OAT Tornado\n(DA neuron, ±50%)",
               "Morris μ*\n(DA neuron)",
               "Sobol S_T\n(DA neuron)",
               "Sobol S_T\n(Motor score)"]
    rows = [
        ["1", "k_agg (49%)", "CI_max (0.31)", "k_ROS_basal (0.50)", "k_ROS_basal (0.40)"],
        ["2", "k_ROS_basal (49%)", "k_ROS_basal (0.27)", "K_mPTP (0.33)", "k_SOD2 (0.32)"],
        ["3", "CI_max (49%)", "K_mPTP (0.24)", "CI_max (0.29)", "K_mPTP (0.30)"],
        ["4", "k_SOD2 (48%)", "k_death (0.24)", "k_SOD2 (0.28)", "CI_max (0.27)"],
        ["5", "K_mPTP (47%)", "k_SOD2 (0.23)", "k_agg (0.14)", "k_agg (0.19)"],
        ["6", "k_shunt (42%)", "k_shunt (0.22)", "k_shunt (0.12)", "k_shunt (0.16)"],
        ["7", "k_clear (36%)", "k_agg (0.20)", "α_clear (0.10)", "α_clear (0.14)"],
        ["8", "α_clear (34%)", "α_clear (0.19)", "k_CL_loss (0.06)", "k_death (0.11)"],
        ["9", "k_mPTP_open (32%)", "k_mPTP_open (0.18)", "k_CI_repair (0.05)", "K_motor (0.09)"],
        ["10", "k_CI_repair (28%)", "k_CI_repair (0.17)", "k_red (0.04)", "k_clear (0.07)"],
    ]
    build_table(doc, headers, rows, col_widths=[1.0, 3.5, 3.0, 3.5, 3.5])
    add_table_note(doc,
        "OAT: one-at-a-time ±50% perturbation (68 runs). "
        "Morris: 20 trajectories (700 runs); μ* = mean absolute elementary effect. "
        "Sobol: variance-based decomposition (N = 256, 3072 runs); S_T = total-order index. "
        "Key insight: ROS module parameters (k_ROS_basal, k_SOD2, k_shunt) dominate both endpoints. "
        "Most top parameters have σ > μ* (Morris), indicating strong nonlinearity/interactions "
        "consistent with vicious cycle amplification. "
        "Disease modifiers (CI_max, α_clear) act largely through interactions (S_T ≫ S_1).")

    # Save
    path = os.path.join(OUT_DIR, "Main_Tables.docx")
    doc.save(path)
    print(f"Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  SUPPLEMENTARY TABLES
# ═══════════════════════════════════════════════════════════════════════════════

def build_supplementary_tables():
    doc = make_doc()

    # ── TABLE S1: Complete Parameter Catalogue ──────────────────────────────
    add_table_title(doc,
        "Table S1. Complete parameter catalogue for the integrated 28-ODE model.")

    headers = ["Parameter", "Value", "Unit", "Module", "Description", "Source Type", "Reference"]

    rows = [
        # Disease modifiers
        ["CI_max", "0.74", "—", "Disease", "Maximum CI repair target", "Estimated", "Calibration (Gao 2017)"],
        ["α_clear", "0.25", "—", "Disease", "aSyn clearance capacity", "Estimated", "Calibration (Choi 2022)"],
        ["k_impair", "0.35", "h⁻¹", "M2→M6", "MG-driven clearance impairment", "Estimated", "Calibration (Ivanova 2024)"],
        # M0: PK (10 params)
        ["k_a", "1.50", "h⁻¹", "M0", "Absorption rate constant", "Estimated", "IP bolus model"],
        ["k_e,plasma", "0.231", "h⁻¹", "M0", "Plasma elimination", "Literature", "SBT-272 Phase 1"],
        ["k_12", "0.30", "h⁻¹", "M0", "Central→peripheral distribution", "Literature", "SBT-272 Phase 1"],
        ["k_21", "0.20", "h⁻¹", "M0", "Peripheral→central redistribution", "Literature", "SBT-272 Phase 1"],
        ["k_brain_in", "0.040", "h⁻¹", "M0", "BBB penetration (plasma→brain)", "Literature", "Szeto 2014"],
        ["k_brain_out", "0.016", "h⁻¹", "M0", "Brain efflux", "Literature", "Szeto 2014"],
        ["k_mito_in", "0.50", "h⁻¹", "M0", "Mitochondrial uptake", "Literature", "Szeto 2014"],
        ["k_mito_out", "0.020", "h⁻¹", "M0", "Mitochondrial efflux", "Literature", "Szeto 2014"],
        ["Dose_norm", "dose/5.0", "—", "M0", "Normalised dose input", "Derived", "—"],
        ["C_ref", "4.43", "—", "M0", "Reference C_mito (5 mg/kg SS)", "Estimated", "Steady-state computation"],
        # M1: CL (12 params)
        ["k_syn_CL", "0.007", "h⁻¹", "M1", "CL synthesis rate", "Literature", "Schlame 2008"],
        ["K_CL", "1.0", "—", "M1", "CL synthesis saturation", "Assumed", "Michaelis-Menten form"],
        ["k_ox", "0.020", "h⁻¹", "M1", "ROS-driven CL oxidation", "Literature", "Kagan 2005"],
        ["k_ALCAT", "5.0", "h⁻¹", "M1", "ALCAT1-mediated oxidation", "Estimated", "Li 2010"],
        ["k_aSyn_ox", "0.020", "h⁻¹", "M1", "aSyn-driven CL oxidation", "Literature", "Robotta 2014"],
        ["k_red", "0.070", "h⁻¹", "M1", "Tafazzin-mediated CL remodelling", "Literature", "Schlame 2008"],
        ["k_drug", "0.100", "h⁻¹", "M1", "Drug-mediated CL_ox reduction", "Estimated", "SBT-272 mechanism"],
        ["k_degrad_ox", "0.010", "h⁻¹", "M1", "Oxidised CL degradation", "Assumed", "First-order clearance"],
        ["k_ext", "0.001", "h⁻¹", "M1", "CL externalisation", "Literature", "Kagan 2014"],
        ["k_clear_ext_CL", "0.050", "h⁻¹", "M1", "Externalised CL clearance", "Assumed", "Mitophagy signal"],
        ["k_act_A", "0.50", "h⁻¹", "M1", "ALCAT1 activation", "Estimated", "Li 2010"],
        ["k_deact_A", "0.45", "h⁻¹", "M1", "ALCAT1 deactivation", "Estimated", "Li 2010"],
        ["k_act_T", "0.50", "h⁻¹", "M1", "Tafazzin activation", "Estimated", "Schlame 2008"],
        ["k_deact_T", "0.556", "h⁻¹", "M1", "Tafazzin deactivation", "Estimated", "Schlame 2008"],
        # M2: aSyn (11 params)
        ["k_syn_aSyn", "0.50", "h⁻¹", "M2", "aSyn monomer synthesis", "Literature", "Bhatt 2009"],
        ["k_turn", "0.014", "h⁻¹", "M2", "Monomer turnover", "Literature", "Bhatt 2009"],
        ["k_agg", "0.009", "h⁻¹", "M2", "Basal aggregation rate", "Literature", "Auluck 2010"],
        ["k_ROS_agg", "2.0", "—", "M2", "ROS-enhanced aggregation", "Estimated", "Souza 2000"],
        ["k_CL_loss", "5.0", "—", "M2", "CL-loss seeding sensitivity", "Estimated", "Nakamura 2008"],
        ["k_seed_CL", "1.0", "—", "M2", "CL_ext seeding coefficient", "Estimated", "Kagan 2014"],
        ["k_refold", "0.010", "h⁻¹", "M2", "Oligomer refolding", "Assumed", "Slow relative to clearance"],
        ["k_clear", "0.241", "h⁻¹", "M2", "Oligomer clearance", "Literature", "Webb 2003"],
        ["k_secrete", "0.005", "h⁻¹", "M2", "aSyn secretion", "Literature", "Lee 2005"],
        ["k_clear_ext_aSyn", "0.10", "h⁻¹", "M2", "Extracellular aSyn clearance", "Literature", "Lee 2005"],
        # M3: ETC (15 params)
        ["k_CI_repair", "0.050", "h⁻¹", "M3", "Complex I repair rate", "Literature", "Morais 2009"],
        ["k_CI_ROS", "0.020", "h⁻¹", "M3", "ROS-mediated CI damage", "Literature", "Sherer 2003"],
        ["k_CI_aSyn", "0.200", "h⁻¹", "M3", "aSyn-mediated CI inhibition", "Literature", "Devi 2008"],
        ["k_SC_form", "0.10", "h⁻¹", "M3", "Supercomplex formation", "Literature", "Lapuente-Brun 2013"],
        ["k_SC_dissoc", "0.05", "h⁻¹", "M3", "Supercomplex dissociation", "Estimated", "CL-dependent"],
        ["k_resp", "1.0", "h⁻¹", "M3", "ETC respiration rate", "Assumed", "Normalised"],
        ["k_ATPase", "0.020", "h⁻¹", "M3", "F1Fo ATPase activity", "Literature", "Brand 2005"],
        ["k_leak", "0.010", "h⁻¹", "M3", "Proton leak rate", "Literature", "Brand 2005"],
        ["k_mPTP_depol", "0.50", "h⁻¹", "M3", "Depolarisation-driven mPTP", "Estimated", "Bernardi 2015"],
        ["k_ATP_syn", "1.20", "h⁻¹", "M3", "ATP synthesis rate", "Literature", "Brand 2005"],
        ["k_ATP_consume", "0.025", "h⁻¹", "M3", "ATP consumption rate", "Assumed", "Basal demand"],
        ["k_ATP_mPTP", "2.0", "—", "M3", "mPTP-induced ATP depletion", "Estimated", "Bernardi 2015"],
        ["k_mPTP_open", "0.50", "h⁻¹", "M3", "mPTP opening rate", "Estimated", "Bernardi 2015"],
        ["k_mPTP_close", "0.50", "h⁻¹", "M3", "mPTP closing rate", "Estimated", "Bernardi 2015"],
        ["K_mPTP", "0.40", "—", "M3", "mPTP Hill threshold", "Estimated", "Calibration"],
        ["n_Hill", "3", "—", "M3", "mPTP Hill coefficient", "Assumed", "Cooperative gating"],
        # M4: ROS (8 params)
        ["k_ROS_basal", "0.122", "h⁻¹", "M4", "Basal ROS production", "Estimated", "Calibration"],
        ["k_shunt", "0.125", "h⁻¹", "M4", "CI/SC-dependent electron leak", "Estimated", "Murphy 2009"],
        ["k_SOD2", "1.0", "h⁻¹", "M4", "SOD2 dismutation", "Literature", "Fukai 2011"],
        ["k_GPx4", "0.5", "h⁻¹", "M4", "GPx4 reduction", "Literature", "Brigelius-Flohe 1999"],
        ["k_ROS_other", "0.05", "h⁻¹", "M4", "Non-enzymatic ROS clearance", "Assumed", "Background"],
        ["k_SOD2_on", "0.50", "h⁻¹", "M4", "SOD2 activation", "Estimated", "Fukai 2011"],
        ["k_SOD2_off", "1.0", "h⁻¹", "M4", "SOD2 deactivation", "Estimated", "Fukai 2011"],
        ["k_GPx4_on", "0.30", "h⁻¹", "M4", "GPx4 activation", "Estimated", "Brigelius-Flohe 1999"],
        ["k_GPx4_off", "1.5", "h⁻¹", "M4", "GPx4 deactivation", "Estimated", "Brigelius-Flohe 1999"],
        # M5: Death (12 params)
        ["k_PINK1_on", "2.0", "h⁻¹", "M5", "PINK1 activation", "Literature", "Narendra 2010"],
        ["k_PINK1_CL", "20.0", "—", "M5", "CL_ext-driven PINK1 activation", "Literature", "Chu 2013"],
        ["k_PINK1_off", "0.50", "h⁻¹", "M5", "PINK1 deactivation", "Estimated", "Narendra 2010"],
        ["k_damage_mPTP", "4.0", "—", "M5", "mPTP damage sensitivity", "Estimated", "Calibration"],
        ["k_damage_energy", "1.0", "—", "M5", "Energy depletion damage", "Assumed", "Normalised"],
        ["k_damage_dpsi", "2.0", "—", "M5", "Depolarisation damage", "Estimated", "Bernardi 2015"],
        ["k_repair", "1.0", "h⁻¹", "M5", "Mitochondrial repair rate", "Assumed", "PINK1/Parkin pathway"],
        ["k_biogen", "3.0", "h⁻¹", "M5", "Mitochondrial biogenesis", "Estimated", "PGC-1α pathway"],
        ["k_death", "7×10⁻⁴", "h⁻¹", "M5", "Mitochondrial death rate", "Estimated", "Calibration"],
        ["K_death", "0.25", "—", "M5", "Death threshold (Hill K)", "Estimated", "Calibration"],
        ["n_death", "4", "—", "M5", "Death Hill cooperativity", "Estimated", "Calibration"],
        ["k_death_inflam", "2×10⁻⁴", "h⁻¹", "M5", "Inflammation-mediated death", "Literature", "Hirsch & Hunot 2009"],
        ["K_inflam_death", "0.45", "—", "M5", "Inflammation death threshold", "Estimated", "Hill gating"],
        ["n_inflam_death", "4", "—", "M5", "Inflammation death cooperativity", "Assumed", "Threshold behaviour"],
        ["k_death_aSyn", "1.5×10⁻⁴", "h⁻¹", "M5", "Direct aSyn proteotoxicity", "Literature", "Winner 2011"],
        ["K_aSyn_death", "0.20", "—", "M5", "aSyn death threshold", "Estimated", "Hill gating"],
        ["n_aSyn_death", "4", "—", "M5", "aSyn death cooperativity", "Assumed", "Threshold behaviour"],
        # M6: Neuroinflammation (7 params)
        ["k_act_MG", "20.0", "h⁻¹", "M6", "TLR2-mediated microglial activation", "Literature", "Kim 2013"],
        ["k_auto", "0.50", "h⁻¹", "M6", "Microglial autocrine activation", "Estimated", "Perry 2010"],
        ["k_deact_MG", "1.0", "h⁻¹", "M6", "Microglial deactivation", "Estimated", "Perry 2010"],
        ["k_TNF_rel", "0.20", "h⁻¹", "M6", "TNF release rate", "Literature", "Mogi 1994"],
        ["k_TNF_deg", "0.10", "h⁻¹", "M6", "TNF degradation", "Literature", "t½ ~ 7 h"],
        ["k_IL1b_rel", "0.15", "h⁻¹", "M6", "IL-1β release rate", "Literature", "Mogi 1994"],
        ["k_IL1b_deg", "0.10", "h⁻¹", "M6", "IL-1β degradation", "Literature", "t½ ~ 7 h"],
        # M7: Motor (4 params)
        ["K_motor", "0.75", "—", "M7", "DA neuron threshold for motor deficit", "Literature", "Bernheimer 1973"],
        ["n_motor", "4", "—", "M7", "Motor Hill cooperativity", "Assumed", "Threshold behaviour"],
        ["k_inflam", "0.30", "—", "M7", "Inflammation contribution to motor", "Estimated", "Calibration"],
        ["k_motor_rate", "0.05", "h⁻¹", "M7", "Motor score approach rate", "Assumed", "Timescale matching"],
    ]
    build_table(doc, headers, rows, col_widths=[2.3, 1.6, 1.0, 1.3, 4.0, 1.8, 2.5])
    add_table_note(doc,
        "Complete catalogue of all parameters in the integrated model. "
        "Source types: Literature = derived from published data with verified DOI/PMID; "
        "Estimated = fitted to quantitative calibration targets; "
        "Assumed = mechanistic justification with literature reference. "
        "All state variables normalised to [0, 1]; rate constants in h⁻¹ unless dimensionless.")

    page_break(doc)

    # ── TABLE S2: VPop Specifications ───────────────────────────────────────
    add_table_title(doc,
        "Table S2. Virtual population parameter specifications (N = 250/arm).")

    headers = ["Parameter", "Nominal Value", "CV", "Distribution", "Category", "Bounds"]
    rows = [
        ["k_ROS_basal", "0.122", "0.25", "Log-normal", "PD mechanism", "—"],
        ["K_mPTP", "0.40", "0.25", "Log-normal", "PD mechanism", "—"],
        ["k_SOD2", "1.0", "0.25", "Log-normal", "PD mechanism", "—"],
        ["k_agg", "0.009", "0.30", "Log-normal", "PD mechanism", "—"],
        ["k_shunt", "0.125", "0.20", "Log-normal", "PD mechanism", "—"],
        ["k_clear", "0.241", "0.25", "Log-normal", "PD mechanism", "—"],
        ["k_death", "4.2×10⁻⁵", "0.30", "Log-normal", "PD mechanism", "—"],
        ["K_death", "0.25", "0.25", "Log-normal", "PD mechanism", "—"],
        ["k_death_inflam", "1.2×10⁻⁵", "0.50", "Log-normal", "PD mechanism", "—"],
        ["k_death_aSyn", "9.0×10⁻⁶", "0.50", "Log-normal", "PD mechanism", "—"],
        ["k_a", "0.40", "0.30", "Log-normal", "PK", "—"],
        ["k_e,plasma", "0.032", "0.30", "Log-normal", "PK", "—"],
        ["k_brain_in", "0.015", "0.40", "Log-normal", "PK", "—"],
        ["CI_max", "0.90", "SD=0.03", "Normal (truncated)", "Disease modifier", "[0.80, 0.96]"],
        ["α_clear", "0.35", "SD=0.05", "Normal (truncated)", "Disease modifier", "[0.20, 0.55]"],
    ]
    build_table(doc, headers, rows, col_widths=[2.5, 2.2, 1.5, 2.8, 2.5, 2.0])
    add_table_note(doc,
        "Death rate nominals are human-scaled (mouse value × death_scale_human = 0.06). "
        "Log-normal sampling ensures positivity: σ = √(ln(1+CV²)), "
        "μ = ln(nominal) − σ²/2. "
        "Disease modifiers use truncated normal to prevent physiologically implausible values. "
        "4 arms: placebo, 10 mg, 30 mg, 60 mg SC daily; 547 days; 244/250 patients completed (97.6%).")

    page_break(doc)

    # ── TABLE S3: Uncertainty Quantification ────────────────────────────────
    add_table_title(doc,
        "Table S3. Uncertainty quantification: bootstrap confidence intervals and parameter envelope.")

    # Part A: Bootstrap CIs
    p = doc.add_paragraph()
    run = p.add_run("A. Bootstrap 95% confidence intervals (B = 10,000, percentile method)")
    run.bold = True
    run.font.size = Pt(10)

    headers = ["Dose", "Endpoint", "Statistic", "Point Estimate", "95% CI Lower", "95% CI Upper"]
    rows = [
        ["10 mg", "DA neuron", "d_unpaired", "0.133", "0.100", "0.170"],
        ["10 mg", "DA neuron", "Cliff’s δ", "0.074", "0.054", "0.097"],
        ["10 mg", "DA neuron", "MCID Resp. %", "16.4%", "12.3%", "21.3%"],
        ["30 mg", "DA neuron", "d_unpaired", "0.221", "0.177", "0.273"],
        ["30 mg", "DA neuron", "Cliff’s δ", "0.126", "0.103", "0.160"],
        ["30 mg", "DA neuron", "MCID Resp. %", "25.8%", "20.5%", "31.6%"],
        ["60 mg", "DA neuron", "d_unpaired", "0.273", "0.224", "0.330"],
        ["60 mg", "DA neuron", "Cliff’s δ", "0.164", "0.134", "0.198"],
        ["60 mg", "DA neuron", "MCID Resp. %", "31.6%", "26.2%", "37.3%"],
        ["30 mg", "Motor score", "d_unpaired", "0.221", "0.174", "0.273"],
        ["30 mg", "Motor score", "Cliff’s δ", "0.180", "0.147", "0.216"],
        ["30 mg", "Motor score", "MCID Resp. %", "17.2%", "13.1%", "21.7%"],
    ]
    build_table(doc, headers, rows, col_widths=[1.5, 2.0, 2.0, 2.2, 2.2, 2.2])

    doc.add_paragraph("")

    # Part B: Parameter envelope
    p = doc.add_paragraph()
    run = p.add_run("B. Joint parameter uncertainty envelope (125-combination sweep)")
    run.bold = True
    run.font.size = Pt(10)

    headers = ["Parameter", "Sweep Range", "Nominal", "Impact on d_unpaired"]
    rows = [
        ["CI_max_human", "0.86–0.94", "0.90", "Moderate"],
        ["α_clear_human", "0.28–0.42", "0.35", "Largest (±2.6× range on DA_saved)"],
        ["death_scale_human", "0.04–0.08", "0.06", "±31% relative uncertainty"],
        ["—", "—", "—", "—"],
        ["30 mg DA_saved", "95% envelope", "0.044", "[0.015, 0.069]"],
        ["30 mg d_unpaired", "Approx. range", "0.221", "[0.062, 0.403]"],
    ]
    build_table(doc, headers, rows, col_widths=[3.0, 2.5, 2.5, 5.0])
    add_table_note(doc,
        "Bootstrap: nonparametric resampling of matched patient pairs across all 4 arms. "
        "Parameter envelope: 5 × 5 × 5 grid sweep of three translational parameters. "
        "Nominal d_unpaired (0.221) sits within the credible range. "
        "α_clear dominates uncertainty because it directly gates vicious cycle entry via aSyn clearance.")

    page_break(doc)

    # ── TABLE S4: Phased Simulation ─────────────────────────────────────────
    add_table_title(doc,
        "Table S4. Phased simulation: disease milestones and treatment delay analysis.")

    # Part A: Disease milestones
    p = doc.add_paragraph()
    run = p.add_run("A. Disease progression milestones (human timescale)")
    run.bold = True
    run.font.size = Pt(10)

    headers = ["Milestone", "DA Neuron Threshold", "Time (Years)", "Literature Reference"]
    rows = [
        ["Prodromal onset", "<70% (30% loss)", "5.6", "Postuma 2012: 5–10 yr prodromal"],
        ["Diagnosis (after 12-mo delay)", "<70% + 12 mo", "6.6", "Breen 2013: median 12 mo to Dx"],
        ["Clinical PD", "<50% (50% loss)", "10.8", "Fearnley & Lees 1991"],
        ["Bernheimer threshold", "<40% (60% loss)", "14.3", "Bernheimer 1973"],
        ["Healthy baseline (20 yr)", "90.4%", "20.0", "~5%/decade age-related"],
    ]
    build_table(doc, headers, rows, col_widths=[3.5, 3.0, 2.0, 5.0])

    doc.add_paragraph("")

    # Part B: Dose x delay matrix
    p = doc.add_paragraph()
    run = p.add_run("B. Treatment efficacy by dose and diagnostic delay (18-month treatment, 30 mg)")
    run.bold = True
    run.font.size = Pt(10)

    headers = ["Delay from\nSymptoms", "DA at\nTreatment Start",
               "DA After\n18-mo Treatment", "DA Without\nTreatment",
               "ΔDA\n(drug benefit)", "Efficacy\nLoss vs 0-mo"]
    rows = [
        ["0 months", "0.700", "0.744", "0.656", "+0.088", "Reference"],
        ["6 months", "0.678", "0.720", "0.634", "+0.086", "−2%"],
        ["12 months", "0.656", "0.700", "0.612", "+0.088", "0%"],
        ["24 months", "0.612", "0.649", "0.572", "+0.077", "−12%"],
    ]
    build_table(doc, headers, rows, col_widths=[2.2, 2.2, 2.5, 2.5, 2.2, 2.2])

    doc.add_paragraph("")

    # Part C: Multi-dose matrix
    p = doc.add_paragraph()
    run = p.add_run("C. Drug benefit (ΔDA) across doses and delays")
    run.bold = True
    run.font.size = Pt(10)

    headers = ["Delay", "10 mg", "30 mg", "60 mg"]
    rows = [
        ["0 months",  "+0.052", "+0.088", "+0.107"],
        ["6 months",  "+0.051", "+0.086", "+0.105"],
        ["12 months", "+0.050", "+0.088", "+0.106"],
        ["24 months", "+0.045", "+0.077", "+0.094"],
    ]
    build_table(doc, headers, rows, col_widths=[3.0, 3.0, 3.0, 3.0])
    add_table_note(doc,
        "Phased simulation: prodromal progression → diagnostic delay → treatment initiation. "
        "State vector extracted at each phase transition via init_override parameter. "
        "Efficacy loss with 24-month delay: ~12% relative to immediate treatment. "
        "Drug benefit is robust across delays, consistent with partial vicious cycle interruption. "
        "Values approximate; exact results depend on VPop realisation.")

    # Save
    path = os.path.join(OUT_DIR, "Supplementary_Tables.docx")
    doc.save(path)
    print(f"Saved: {path}")


# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    build_main_tables()
    build_supplementary_tables()
    print("\nDone. Both DOCX files in 'New publication plots/' folder.")
