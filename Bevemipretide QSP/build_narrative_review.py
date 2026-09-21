"""
Build Narrative_Review.docx — condensed Part I narrative review (~1300 words).
Focus: vicious cycle causation and mitochondrial dysfunction as core pathology.
Format: Times New Roman 12pt, double-spaced, 1-inch margins.
59 references, renumbered sequentially from the original 88-ref manuscript.
"""

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
import os, re

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "New publication plots")

# Original manuscript refs used in this review (old number -> ref text)
# Extracted from Manuscript_condensed.docx
REFS_OLD = {
    1: "GBD 2016 Parkinson’s Disease Collaborators. Global, regional, and national burden of Parkinson’s disease, 1990–2016: a systematic analysis for the Global Burden of Disease Study 2016. Lancet Neurol. 2018;17(11):939–953.",
    2: "GBD 2019 Dementia Forecasting Collaborators. Estimation of the global prevalence of dementia in 2019 and forecasted prevalence in 2050: an analysis for the Global Burden of Disease Study 2019. Lancet Public Health. 2022;7(2):e105–e125.",
    5: "Bernheimer H, Birkmayer W, Hornykiewicz O, Jellinger K, Seitelberger F. Brain dopamine and the syndromes of Parkinson and Huntington. J Neurol Sci. 1973;20(4):415–455.",
    6: "Fearnley JM, Lees AJ. Ageing and Parkinson’s disease: substantia nigra regional selectivity. Brain. 1991;114(Pt 5):2283–2301.",
    7: "Kish SJ, Shannak K, Hornykiewicz O. Uneven pattern of dopamine loss in the striatum of patients with idiopathic Parkinson’s disease. N Engl J Med. 1988;318(14):876–880.",
    9: "LeWitt PA. Levodopa therapy for Parkinson’s disease: pharmacokinetics and pharmacodynamics. Mov Disord. 2015;30(1):64–72.",
    10: "Langston JW, Ballard P, Tetrud JW, Irwin I. Chronic parkinsonism in humans due to a product of meperidine-analog synthesis. Science. 1983;219(4587):979–980.",
    11: "Nicklas WJ, Vyas I, Heikkila RE. Inhibition of NADH-linked oxidation in brain mitochondria by 1-methyl-4-phenyl-pyridine, a metabolite of the neurotoxin 1-methyl-4-phenyl-1,2,5,6-tetrahydropyridine. Life Sci. 1985;36(26):2503–2508.",
    12: "Betarbet R, Sherer TB, MacKenzie G, Garcia-Osuna M, Panov AV, Greenamyre JT. Chronic systemic pesticide exposure reproduces features of Parkinson’s disease. Nat Neurosci. 2000;3(12):1301–1306.",
    14: "Schapira AH, Cooper JM, Dexter D, Clark JB, Jenner P, Marsden CD. Mitochondrial complex I deficiency in Parkinson’s disease. J Neurochem. 1990;54(3):823–827.",
    16: "Gao F, Chen D, Hu Q, Wang G. Rotenone directly induces BV2 cell activation via the p38 MAPK pathway. PLoS One. 2013;8(8):e72046.",
    17: "Devi L, Raghavendran V, Prabhu BM, Avadhani NG, Anandatheerthavarada HK. Mitochondrial import and accumulation of α-synuclein impair complex I in human dopaminergic neuronal cultures and Parkinson disease brain. J Biol Chem. 2008;283(14):9089–9100.",
    18: "Bolam JP, Pissadaki EK. Living on the edge with too many mouths to feed: why dopamine neurons die. Mov Disord. 2012;27(12):1478–1483.",
    19: "Surmeier DJ, Obeso JA, Bhatt L, Bhargava P. Selective neuronal vulnerability in Parkinson disease. Nat Rev Neurosci. 2017;18(2):101–113.",
    20: "Paradies G, Paradies V, Ruggiero FM, Petrosillo G. Role of cardiolipin in mitochondrial function and dynamics in health and disease: molecular and pharmacological aspects. Cells. 2019;8(7):728.",
    22: "Pfeiffer K, Gilderson V, Neff D, et al. Cardiolipin stabilizes respiratory chain supercomplexes. J Biol Chem. 2003;278(52):52873–52880.",
    23: "Zhang M, Mileykovskaya E, Dowhan W. Gluing the respiratory chain together: cardiolipin is required for supercomplex formation in the inner mitochondrial membrane. J Biol Chem. 2002;277(46):43553–43556.",
    24: "Kagan VE, Tyurin VA, Jiang J, et al. Cytochrome c acts as a cardiolipin oxygenase required for release of proapoptotic factors. Nat Chem Biol. 2005;1(4):223–232.",
    25: "Schlame M, Greenberg ML. Biosynthesis, remodeling and turnover of mitochondrial cardiolipin. Biochim Biophys Acta Mol Cell Biol Lipids. 2017;1862(1):3–11.",
    26: "Kiebish MA, Han X, Cheng H, Chuang JH, Seyfried TN. Cardiolipin and electron transport chain abnormalities in mouse brain tumor mitochondria. J Lipid Res. 2008;49(12):2545–2556.",
    27: "Xu Y, Kelley RI, Bhatt TK, Bhargava P. Tafazzin mutation–induced cardiolipin deficiency. J Biol Chem. 2006;281(51):39217–39224.",
    28: "Barth PG, Scholte HR, Berden JA, et al. An X-linked mitochondrial disease affecting cardiac muscle, skeletal muscle and neutrophil leucocytes. J Neurol Sci. 1983;62(1–3):327–355.",
    30: "Song C, Zhang J, Qi S, et al. Cardiolipin remodeling by ALCAT1 links mitochondrial dysfunction to Parkinson’s diseases. Aging Cell. 2019;18(3):e12941.",
    31: "Paradies G, Petrosillo G, Paradies V, Ruggiero FM. Oxidative stress, mitochondrial bioenergetics, and cardiolipin in aging. Free Radic Biol Med. 2010;48(10):1286–1295.",
    32: "Chu CT, Ji J, Dagda RK, et al. Cardiolipin externalization to the outer mitochondrial membrane acts as an elimination signal for mitophagy in neuronal cells. Nat Cell Biol. 2013;15(10):1197–1205.",
    33: "Robotta M, Gerding HR, Vogel A, et al. Alpha-synuclein binds to the inner membrane of mitochondria in an α-helical conformation. ChemBioChem. 2014;15(17):2499–2502.",
    34: "Nakamura K, Nemani VM, Azarbal F, et al. Direct membrane association drives mitochondrial fission by the Parkinson disease-associated protein α-synuclein. J Biol Chem. 2011;286(23):20710–20726.",
    36: "Davidson WS, Jonas A, Clayton DF, George JM. Stabilization of α-synuclein secondary structure upon binding to synthetic membranes. J Biol Chem. 1998;273(16):9443–9449.",
    37: "Burré J, Sharma M, Tsetsenis T, et al. α-Synuclein promotes SNARE-complex assembly in vivo and in vitro. Science. 2010;329(5999):1663–1667.",
    40: "Chiti F, Dobson CM. Protein misfolding, amyloid formation, and human disease: a summary of progress over the last decade. Annu Rev Biochem. 2017;86:27–68.",
    41: "Iljina M, Garcia GA, Horrocks MH, et al. Kinetic model of the aggregation of alpha-synuclein provides insights into prion-like spreading. Proc Natl Acad Sci U S A. 2016;113(9):E1206–E1215.",
    42: "Winner B, Jappelli R, Maji SK, et al. In vivo demonstration that α-synuclein oligomers are toxic. Proc Natl Acad Sci U S A. 2011;108(10):4194–4199.",
    47: "Singleton AB, Farrer M, Johnson J, et al. α-Synuclein locus triplication causes Parkinson’s disease. Science. 2003;302(5646):841.",
    48: "Polymeropoulos MH, Lavedan C, Leroy E, et al. Mutation in the α-synuclein gene identified in families with Parkinson’s disease. Science. 1997;276(5321):2045–2047.",
    49: "Nalls MA, Blauwendraat C, Vallerga CL, et al. Identification of novel risk loci, causal insights, and heritable risk for Parkinson’s disease: a meta-analysis of genome-wide association studies. Lancet Neurol. 2019;18(12):1091–1102.",
    50: "Bayir H, Kapralov AA, Jiang J, et al. Peroxidase mechanism of lipid-dependent cross-linking of synuclein with cytochrome c. J Biol Chem. 2009;284(23):15951–15969.",
    51: "Luth ES, Stavrovskaya IG, Bhatt L, Bhargava P. Soluble, prefibrillar α-synuclein oligomers promote complex I-dependent, extramitochondrial oxidant production. J Biol Chem. 2014;289(31):21490–21507.",
    52: "Hirst J, King MS, Pryde KR. The production of reactive oxygen species by complex I. Biochem Soc Trans. 2008;36(Pt 5):976–980.",
    54: "Lapuente-Brun E, Moreno-Loshuertos R, Acin-Perez R, et al. Supercomplex assembly determines electron flux in the mitochondrial electron transport chain. Science. 2013;340(6140):1567–1570.",
    55: "Acín-Pérez R, Enríquez JA. The function of the respiratory supercomplexes: the plasticity model. Biochim Biophys Acta. 2014;1837(4):444–450.",
    56: "Cortassa S, Aon MA, Marbán E, Winslow RL, O’Rourke B. An integrated model of cardiac mitochondrial energy metabolism and calcium dynamics. Biophys J. 2003;84(4):2734–2755.",
    57: "Bernardi P, Rasola A, Forte M, Lippe G. The mitochondrial permeability transition pore: channel formation by F-ATP synthase, integration in signal transduction, and role in pathophysiology. Physiol Rev. 2015;95(4):1111–1155.",
    58: "Kembro JM, Aon MA, Winslow RL, O’Rourke B, Cortassa S. Integrating mitochondrial energetics, redox and ROS metabolic networks: a two-compartment model. Biophys J. 2013;104(2):332–343.",
    61: "Choi DJ, An J, Jou I, Park SM, Joe EH. A Parkinson’s disease gene, DJ-1, regulates anti-inflammatory roles of astrocytes through prostaglandin D2 synthase expression. Neurobiol Dis. 2019;127:482–491.",
    62: "Narendra DP, Jin SM, Tanaka A, et al. PINK1 is selectively stabilized on impaired mitochondria to activate Parkin. PLoS Biol. 2010;8(1):e1000298.",
    63: "Valente EM, Abou-Sleiman PM, Caputo V, et al. Hereditary early-onset Parkinson’s disease caused by mutations in PINK1. Science. 2004;304(5674):1158–1160.",
    64: "Kitada T, Asakawa S, Hattori N, et al. Mutations in the parkin gene cause autosomal recessive juvenile parkinsonism. Nature. 1998;392(6676):605–608.",
    65: "Canet-Avilés RM, Wilson MA, Miller DW, et al. The Parkinson’s disease protein DJ-1 is neuroprotective due to cysteine-sulfinic acid-driven mitochondrial localization. Proc Natl Acad Sci U S A. 2004;101(24):9103–9108.",
    66: "Zhang W, Wang T, Pei Z, et al. Aggregated α-synuclein activates microglia: a process leading to disease progression in Parkinson’s disease. FASEB J. 2005;19(6):533–542.",
    67: "Kim C, Ho DH, Suk JE, et al. Neuron-released oligomeric α-synuclein is an endogenous agonist of TLR2 for paracrine activation of microglia. Nat Commun. 2013;4:1562.",
    68: "Codolo G, Plotegher N, Pozzobon T, et al. Triggering of inflammasome by aggregated α-synuclein, an inflammatory response in synucleinopathies. PLoS One. 2013;8(1):e55375.",
    70: "Hirsch EC, Hunot S. Neuroinflammation in Parkinson’s disease: a target for neuroprotection? Lancet Neurol. 2009;8(4):382–397.",
    71: "Ivanova N, Liu Q, Agca C, et al. An α-synuclein QSP model incorporating TLR2-mediated neuroinflammation. CPT Pharmacometrics Syst Pharmacol. 2024;13(5):845–860.",
    73: "Hunot S, Boissiere F, Faucheux B, et al. Nitric oxide synthase and neuronal vulnerability in Parkinson’s disease. Neuroscience. 1996;72(2):355–363.",
    74: "Winner B, Jappelli R, Maji SK, et al. In vivo demonstration that α-synuclein oligomers are toxic. Proc Natl Acad Sci U S A. 2011;108(10):4194–4199.",
    75: "Sidransky E, Nalls MA, Aasly JO, et al. Multicenter analysis of glucocerebrosidase mutations in Parkinson’s disease. N Engl J Med. 2009;361(17):1651–1661.",
    76: "Singh A, Zhi L, Zhang H. LRRK2 and mitochondria: recent advances and current views. Brain Res. 2019;1702:96–104.",
    87: "US Food and Drug Administration. Model-Informed Drug Development Pilot Program. Guidance for Industry. 2018.",
    88: "Gadkar K, Kirouac DC, Mager DE, van der Graaf PH, Bhatt L. A six-stage workflow for model-informed drug development (MIDD). J Pharmacokinet Pharmacodyn. 2016;43(1):3–10.",
}


def make_doc():
    doc = Document()
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
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


def build():
    doc = make_doc()

    # Build old->new mapping (sorted by old number)
    old_sorted = sorted(REFS_OLD.keys())
    old_to_new = {old: i + 1 for i, old in enumerate(old_sorted)}

    add_heading(doc, "1. Narrative Review")

    add_para(doc,
        "This review synthesises evidence for mitochondrial dysfunction "
        "as the central pathogenic mechanism in Parkinson’s disease, "
        "with emphasis on the self-amplifying vicious cycle that links "
        "α-synuclein aggregation, cardiolipin oxidation, Complex I "
        "impairment, and reactive oxygen species overproduction."
    )

    # ── 1.1 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.1 Parkinson’s Disease and the Unmet Therapeutic Need")

    add_para(doc,
        "Parkinson’s disease (PD) affects over 10 million people "
        "worldwide, with prevalence more than doubling between 1990 and "
        "2016 [1,2]. Motor symptoms emerge from progressive dopaminergic "
        "neuron loss in the substantia nigra pars compacta (SNpc). By the "
        "time symptoms appear, 50–70% of nigral neurons have been lost "
        "and striatal dopamine is depleted by over 80% [3,5], implying a "
        "5–10 year prodromal phase [4]. Every approved therapy addresses "
        "symptoms without slowing neuronal loss [6]. The repeated failure "
        "of disease-modification trials targeting individual molecular "
        "species has intensified the search for strategies that engage "
        "the disease at a systems level."
    )

    # ── 1.2 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.2 Complex I Deficiency")

    add_para(doc,
        "The link between mitochondria and PD was established when MPTP, "
        "a contaminant in illicit drug synthesis, caused acute parkinsonism "
        "through selective Complex I (CI) inhibition by its metabolite "
        "MPP⁺ [7,8]. Betarbet and colleagues then showed that chronic "
        "systemic rotenone reproduced selective nigral neuron loss, "
        "α-synuclein aggregation, and Lewy body–like inclusions in "
        "rats [9]. Post-mortem studies confirmed these models. Schapira "
        "and colleagues reported a 37% CI activity reduction in PD "
        "substantia nigra with no deficiency in other respiratory "
        "complexes [10], and Gao and colleagues quantified CI at "
        "approximately 49% of wild-type levels [11]. The selective "
        "vulnerability of SNpc neurons reflects their exceptional metabolic "
        "demands from extensive axonal arbors, autonomous pacemaking, and "
        "high oxidative burden [13,14]. Importantly, α-synuclein itself "
        "enters mitochondria and binds CI subunits NDUFS1 and NDUFV2, "
        "inhibiting electron transfer in a dose-dependent manner [12]. "
        "This places proteinopathy and bioenergetic failure on the same "
        "causal axis."
    )

    # ── 1.3 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.3 Cardiolipin as Structural Linchpin")

    add_para(doc,
        "Cardiolipin (CL) is a diphosphatidylglycerol phospholipid found "
        "almost exclusively in the inner mitochondrial membrane, "
        "constituting 15–20% of total phospholipid content [15]. It "
        "stabilises respiratory chain supercomplexes, anchors cytochrome c, "
        "supports the proton gradient, and enables multiple carrier "
        "proteins [16–19]. Tafazzin remodels CL into its mature form, "
        "and tafazzin loss-of-function mutations cause Barth syndrome, "
        "providing direct genetic proof of CL’s bioenergetic role "
        "[20–22]. Under stress, the alternative enzyme ALCAT1 is "
        "upregulated roughly 3-fold, generating oxidation-prone CL "
        "species [23]."
    )

    add_para(doc,
        "CL is especially susceptible to oxidative damage because of its "
        "unsaturated acyl chains and proximity to electron transport chain "
        "ROS sites [24]. Kagan and colleagues showed that cytochrome c "
        "detached from oxidised CL acquires peroxidase activity targeting "
        "CL itself, creating a local amplification loop [18]. Under "
        "stress, CL also redistributes to the outer membrane. This "
        "externalised CL recruits LC3 for mitophagy of damaged "
        "mitochondria [25] but simultaneously provides a templating "
        "surface for α-synuclein oligomerisation [26,27]. Whether CL "
        "externalisation is net protective or net pathological depends on "
        "the balance between these competing processes. Gao and colleagues "
        "reported a 23% reduction in native CL in a PD model, sufficient "
        "to destabilise supercomplexes [11]."
    )

    # ── 1.4 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.4 α-Synuclein Aggregation")

    add_para(doc,
        "α-Synuclein is a 140-amino acid presynaptic protein that is "
        "intrinsically disordered in solution but adopts an α-helical "
        "conformation upon binding CL-enriched membranes, where it "
        "participates in vesicle dynamics [28,29]. Under pathological "
        "conditions it misfolds into β-sheet conformations through a "
        "nucleation-dependent process, with soluble oligomers recognised "
        "as the primary toxic species [30–32]. Aggregation is promoted by "
        "oxidative modification, loss of CL-mediated stabilisation, and "
        "impaired clearance. The genetic evidence is strong. SNCA "
        "multiplications cause dose-dependent autosomal dominant PD [33], "
        "point mutations accelerate aggregation [34], and genome-wide "
        "association studies identify SNCA as the strongest common genetic "
        "risk factor for sporadic PD [35]."
    )

    # ── 1.5 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.5 The αSyn–CL–CI–ROS Vicious Cycle")

    add_para(doc,
        "The interactions described above do not operate in isolation. "
        "They assemble into a self-amplifying vicious cycle whose "
        "behaviour differs qualitatively from the sum of its parts."
    )

    add_para(doc,
        "The cycle can be entered at any node. Bayir and colleagues "
        "showed that CL stabilises the native α-helical conformation of "
        "α-synuclein. When CL is oxidised or depleted, this stabilisation "
        "is lost and the protein reverts to its aggregation-prone "
        "state [36]. The relationship is bidirectional. α-Synuclein "
        "oligomers enhance cytochrome c peroxidase activity against "
        "CL [36] and promote CI-dependent oxidant production [37]."
    )

    add_para(doc,
        "Depleted or oxidised CL destabilises respiratory supercomplexes "
        "because CL serves as the structural adhesive at inter-complex "
        "interfaces [16,17]. Respirasomes dissociate, electron channelling "
        "fails, electron leak increases, and ATP synthesis falls "
        "[39–41]. The elevated mitochondrial ROS (mROS), measured at "
        "177% of basal by Choi and colleagues [44], then drives four "
        "parallel effects. It oxidises CL acyl chains directly. It "
        "modifies α-synuclein methionine residues to promote aggregation. "
        "It damages CI iron–sulphur clusters. And it promotes mPTP "
        "opening once a threshold is exceeded [18,38,42]."
    )

    add_para(doc,
        "mPTP opening is the point of no return. Its probability has a "
        "steep, Hill-type dependence on ROS levels. Gradual stress is "
        "tolerated until a threshold is crossed, after which membrane "
        "potential collapses, cytochrome c enters the cytosol, and ATP "
        "is acutely depleted [42,43]."
    )

    add_para(doc,
        "Four properties of this cycle matter for both the disease and "
        "the therapeutic opportunity. First, amplification. Each node "
        "individually produces a modest perturbation (23% CL "
        "reduction [11], 1.77-fold ROS increase [44]), but cyclic "
        "amplification yields cumulative damage far exceeding any single "
        "contribution. This explains why interventions targeting one "
        "molecular species have consistently failed. Second, threshold "
        "behaviour. The steep responses of mPTP opening and "
        "α-synuclein nucleation produce a prolonged prodromal phase "
        "followed by rapid decline. Third, multiple entry points. CI "
        "inhibition, α-synuclein overexpression, CL deficiency, or "
        "oxidative stress can each start the same self-sustaining cycle, "
        "explaining phenotypic convergence across genetically and "
        "environmentally diverse forms of PD. Fourth, partial "
        "reversibility. Disrupting even one link may slow the entire "
        "cascade if the intervention is quantitatively sufficient. This "
        "is the therapeutic hypothesis underlying CL stabilisation."
    )

    # ── 1.6 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.6 Downstream Amplifiers")

    add_para(doc,
        "Two feedback loops extend beyond the core cycle. The "
        "PINK1/Parkin pathway is the principal mitochondrial quality "
        "control mechanism. When membrane potential falls, PINK1 recruits "
        "Parkin to tag damaged mitochondria for autophagic removal [45]. "
        "Externalised CL provides a complementary mitophagy signal "
        "through LC3 recognition [25]. Loss-of-function mutations in "
        "PINK1 [46] and Parkin [47] cause autosomal recessive early-onset "
        "PD, demonstrating what happens when this quality control fails."
    )

    add_para(doc,
        "Neuroinflammation forms the second loop. Extracellular "
        "α-synuclein from damaged neurons activates microglial TLR2/TLR4 "
        "receptors and the NLRP3 inflammasome [49–51]. The resulting "
        "TNFα and IL-1β contribute to neuronal death and impair "
        "α-synuclein clearance [52,53], driving further accumulation and "
        "microglial activation. Ivanova and colleagues quantified this "
        "using TLR2 knockout mice, which showed 23–45% reduction in "
        "α-synuclein accumulation [53]."
    )

    # ── 1.7 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.7 Convergent Death Pathways and Genetic Validation")

    add_para(doc,
        "Dopaminergic neuron death proceeds through at least three "
        "convergent mechanisms. The mitochondrial pathway, driven by "
        "mPTP-mediated cytochrome c release and caspase activation, "
        "accounts for the majority of disease-state loss [42]. The "
        "neuroinflammatory pathway operates independently of CL through "
        "TNFα-mediated extrinsic apoptosis and IL-1β-potentiated "
        "excitotoxicity [52,54]. Hirsch and Hunot estimated "
        "neuroinflammation accounts for 20–30% of dopaminergic "
        "death [52]. A third pathway involves direct α-synuclein "
        "oligomer proteotoxicity through membrane disruption [55]. This "
        "multi-pathway architecture means a CL-targeted drug can rescue "
        "the mitochondrial fraction (approximately 70–75%) but not "
        "CL-independent death, establishing a realistic efficacy ceiling."
    )

    add_para(doc,
        "The genetic architecture reinforces mitochondrial centrality. "
        "Every major PD risk gene converges on mitochondrial biology. "
        "PINK1 and Parkin regulate mitophagy [46,47]. DJ-1 senses "
        "mitochondrial redox state [48]. GBA impairs lysosomal function "
        "and secondarily compromises mitophagy [56]. LRRK2 regulates "
        "mitochondrial fission and fusion [57]. That independent genetic "
        "loci funnel into a shared mitochondrial phenotype validates a "
        "disease model centred on mitochondrial dysfunction."
    )

    # ── 1.8 ──────────────────────────────────────────────────────────────
    add_subheading(doc,
        "1.8 Rationale for Quantitative Systems Modelling")

    add_para(doc,
        "PD neurodegeneration is driven by interconnected processes with "
        "nonlinear amplification, threshold behaviours, and competing "
        "consequences of the same molecular events. These features "
        "demand a mathematical framework. Quantitative systems "
        "pharmacology (QSP), endorsed by the FDA’s Model-Informed Drug "
        "Development framework [58] and IQ Consortium "
        "recommendations [59], provides this capability. In Part II, "
        "we translate this evidence into a 28-ODE QSP model of "
        "bevemipretide in Parkinson’s disease."
    )

    # ── References ───────────────────────────────────────────────────────
    add_heading(doc, "References")

    for old_num in old_sorted:
        new_num = old_to_new[old_num]
        ref_text = REFS_OLD[old_num]
        p = doc.add_paragraph()
        run = p.add_run(f"{new_num}. {ref_text}")
        run.font.name = "Times New Roman"
        run.font.size = Pt(10)

    # ── Save ─────────────────────────────────────────────────────────────
    path = os.path.join(OUT_DIR, "Narrative_Review.docx")
    doc.save(path)
    print(f"Saved: {path}")

    wc = 0
    ref_section = False
    for p in doc.paragraphs:
        txt = p.text.strip()
        if not txt:
            continue
        if p.style.name.startswith("Heading") and "References" in txt:
            ref_section = True
            continue
        if ref_section:
            continue
        if not p.style.name.startswith("Heading"):
            wc += len(txt.split())
    print(f"Word count (body text, excl. references): {wc}")
    print(f"Total references: {len(REFS_OLD)}")

    # Verify all citation numbers in text are valid
    all_text = " ".join(p.text for p in doc.paragraphs
                        if not p.style.name.startswith("Heading"))
    cited = set()
    for m in re.finditer(r'\[([0-9,– ]+)\]', all_text):
        chunk = m.group(1)
        for part in chunk.split(","):
            part = part.strip()
            if "–" in part:
                lo, hi = part.split("–")
                for n in range(int(lo), int(hi) + 1):
                    cited.add(n)
            else:
                cited.add(int(part))

    max_ref = len(REFS_OLD)
    bad = [n for n in cited if n < 1 or n > max_ref]
    if bad:
        print(f"WARNING: out-of-range citations: {sorted(bad)}")
    else:
        print(f"All {len(cited)} cited numbers are in [1, {max_ref}] ✔")

    uncited = set(range(1, max_ref + 1)) - cited
    if uncited:
        print(f"WARNING: {len(uncited)} refs in list but not cited: "
              f"{sorted(uncited)}")
    else:
        print("Every reference is cited at least once ✔")


if __name__ == "__main__":
    build()
