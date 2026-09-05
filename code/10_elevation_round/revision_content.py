# -*- coding: utf-8 -*-
"""Content for the JTM 升华 revision. All strings plain; citations as unicode superscripts."""

TITLE = ("Chronic obstructive pulmonary disease and heart failure: a vertically integrated audit "
         "of genetic causality, global burden, and a tissue-confined inflammo-fibrotic mechanism "
         "in a cross-organ comorbidity syndrome")

ABSTRACT_BACKGROUND = ("Background Heart failure (HF) is a frequent, lethal complication of chronic obstructive "
    "pulmonary disease (COPD), but causality, burden, and mechanisms remain uncertain. We treated the COPD–HF "
    "syndrome as a prototype cross-organ comorbidity axis and triangulated genetic, population-level, and "
    "tissue evidence under an audit discipline.")

ABSTRACT_METHODS = ("Methods We performed two-sample Mendelian randomization (MR) in European, FinnGen R12, and "
    "BioBank Japan (BBJ) cohorts with meta-analysis, multivariable adjustment, mediation, and power analysis; "
    "extended replication to phenotypically enriched and HF-subtype endpoints, with winner's-curse correction "
    "and independent-exposure re-estimation; decomposed between-cohort heterogeneity by meta-regression; "
    "quantified COPD-attributable HF in 204 countries/territories (1990–2021) using GBD 2021 and Bayesian "
    "projections to 2050; and integrated bulk and disease-state single-cell transcriptomics with interorgan "
    "ligand–receptor modelling and phenome-wide safety scanning.")

ABSTRACT_RESULTS = ("Results Genetically predicted COPD increased HF risk in the discovery cohort (odds ratio [OR] "
    "1.15, 95% CI 1.08–1.22), independent of smoking, adiposity, blood pressure, and lung function; C-reactive "
    "protein mediated 4.4%. The estimate withstood winner's-curse correction (OR 1.16) and re-estimation with "
    "independent FinnGen exposure weights (OR 1.19). Replication was null in FinnGen and inverse in East Asian "
    "endpoints (BBJ 2025 HFpEF OR 0.87); ancestry, more than endpoint definition, structured the heterogeneity "
    "(meta-regression QM P = 2.4×10⁻⁴; East Asian I² = 0%). COPD-attributable HF affected 3.61 million people in "
    "2021 and was projected to reach 10.93 million by 2050, with high-SDI resurgence. Disease-tissue analyses "
    "localized a myeloid SERPINE1–integrin inflammo-fibrotic program, and an interorgan ligand–receptor map linked "
    "COPD lung cells to failing myocardium, with donor-robust cardiac receptor changes (CD36, IL1RAP, ITGB3). "
    "Circulating PAI-1 showed no causal association, with no safety signal.")

ABSTRACT_CONCLUSIONS = ("Conclusions COPD is a probable, modifiable HF driver whose genetic effect is phenotype- and "
    "ancestry-dependent, with a large, unequal, and growing burden. Tissue-confined SERPINE1–integrin signalling "
    "and its interorgan receptor landscape are translational targets requiring pharmacological validation; "
    "ancestry- and phenotype-stratified replication is now a pre-specified, falsifiable priority. The "
    "causal–burden–mechanism audit is portable to other cross-organ comorbidity syndromes.")

KEYWORDS = ("Keywords Chronic obstructive pulmonary disease; Heart failure; Mendelian randomization; Global burden "
            "of disease; Transcriptomics; SERPINE1/PAI-1; Cross-organ comorbidity; Interorgan ligand–receptor network")

INTRO_PARA3 = ("We reasoned that these three gaps are, in fact, one gap: the absence of a vertically integrated and "
    "internally audited evidentiary chain running from cause to consequence to mechanism. Beyond the COPD–HF "
    "question itself, such a chain constitutes a reusable paradigm for dissecting cross-organ comorbidity "
    "syndromes—constellations in which chronic disease of one organ generates pathology in another—whose causal "
    "architecture, population footprint, and molecular executors are rarely examined within a single framework. "
    "We therefore designed a study in which each layer constrains the others (MR defines the causal effect and its "
    "mediation, GBD analytics quantify its population footprint and future trajectory, and multiscale "
    "transcriptomics identifies and localizes the executing molecular program), executed under an audit discipline "
    "in which every reported number must trace to a deposited output file. Where the layers disagreed, between "
    "cohorts or between plasma and tissue, we report the disagreement explicitly and interrogate it with "
    "additional data, because in this study the inconsistencies turned out to be among the most informative results.")

METHODS_DATASOURCES_ADD = (" For extended replication and phenotype stratification we additionally used the BBJ 2025 "
    "HF GWAS suite (all-cause HF, 16,251 cases; HFpEF, 7,154; HFrEF, 4,254; GWAS Catalog GCST90668009–GCST90668011)⁶⁰ "
    "and FinnGen R12 phenotypically enriched endpoints (HF and coronary heart disease, 17,532 cases; HF with "
    "BMI ≥ 25, 22,965 cases; other heart diseases, 106,069 cases; other or ill-defined heart diseases, 2,192 cases). "
    "The Million Veteran Program HFpEF/HFrEF GWAS⁶¹ (dbGaP phs001672) was identified but is access-restricted and "
    "was therefore reserved as a pre-specified external validation target.")

METH_H2_WINNERS = "Winner's curse correction and heterogeneity forensics"
METH_P_WINNERS = ("Because instruments were selected at P < 5×10⁻⁸ within the discovery GWAS, the discovery effect "
    "may be inflated by winner's curse. We corrected per-SNP exposure estimates using the conditional-likelihood "
    "estimator of Zhong and Prentice⁶² and re-estimated the causal effect; as an orthogonal check, we re-computed "
    "exposure weights for the same ten variants in three independent FinnGen R12 COPD endpoints (J10_COPD, "
    "COPD_LATER, COPD_EARLY), which are free of discovery-selection bias, and repeated MR against the discovery "
    "outcome. To decompose between-cohort heterogeneity, we assembled all ten independent cohort–endpoint IVW "
    "estimates and fitted fixed- and random-effects (DerSimonian–Laird) meta-analyses overall and within ancestry "
    "and endpoint-breadth strata, an illustrative weighted meta-regression with ancestry (East Asian vs European), "
    "endpoint strictness, and HFpEF definition as moderators, and a cohort-level funnel plot with an Egger "
    "asymmetry test (interpreted illustratively given k = 10).")

METH_H2_INTERORGAN = "Interorgan ligand–receptor modelling"
METH_P_INTERORGAN = ("To render the lung–heart axis as an explicit network, we scored 31 curated ligand–receptor "
    "pairs across organs. Ligands (SERPINE1, THBS1, S100A8/A9/A12, SPP1, FN1, TGFB1, AREG, IL6, CCL2, CXCL8, IL33, "
    "PLAU) were summarized as log1p-CP10k means per sender cell type and condition in the COPD lung dataset, and "
    "receptors (PLAUR, LRP1, VTN, CD36, CD47, TLR4, AGER, ITGAV, ITGB1, ITGB3, ITGA5, CD44, TGFBR1, TGFBR2, EGFR, "
    "IL6R, IL6ST, CCR2, CXCR1, CXCR2, IL1RL1, IL1RAP) per receiver compartment and condition in the failing-heart "
    "atlas (fibroblasts, endothelium, macrophages). Edge scores (product of sender ligand and receiver receptor "
    "means) were contrasted between the COPD–failing and normal–normal states (1,428 directed edges), and ligand- "
    "and receptor-level changes were re-tested at donor level with Welch tests and Benjamini–Hochberg correction.")

METHODS_SC_ADD = (" The cis screen was subsequently extended beyond gene-expression QTLs to exon-level, "
    "transcript-usage (tx and txrev), and splicing (leafcutter) QTLs in five cardiorespiratory tissues (lung, left "
    "ventricle, atrial appendage, aorta, coronary artery) using eQTL Catalogue permuted top-variant files⁶³, "
    "testing whether any molecular layer yields a usable SERPINE1 cis instrument.")

RES_H2_WINNERS = "The discovery estimate withstands winner's-curse correction and independent re-estimation"
RES_P_WINNERS = ("Two formal checks addressed whether the discovery estimate (OR 1.15) is inflated by instrument "
    "selection within the discovery GWAS. First, conditional-likelihood correction for genome-wide selection⁶² left "
    "the estimate essentially unchanged (IVW OR 1.16, 95% CI 1.08–1.24; per-SNP shrinkage < 2% for seven of ten "
    "instruments; Fig. 8C). Second, re-estimating the same ten variants' exposure weights in three independent "
    "FinnGen COPD endpoints—free of discovery-selection bias—reproduced the effect at full strength (J10_COPD "
    "OR 1.19, 95% CI 1.10–1.28, P = 3.8×10⁻⁶; later-onset COPD OR 1.19, P = 9.2×10⁻⁶; early-onset COPD OR 1.17, "
    "P = 3.1×10⁻⁶). Winner's curse is therefore excluded as an explanation of the discovery–replication gap, "
    "localizing the discrepancy downstream of instrument strength.")

RES_H2_EXTENDED = "Extended replication across phenotype- and ancestry-diverse endpoints"
RES_P_EXTENDED = ("We widened the replication surface to seven further endpoints (Fig. 8A). Within FinnGen, "
    "phenotypically enriched composites remained null: HF with coronary heart disease OR 1.048 (95% CI "
    "0.973–1.128, P = 0.21), HF with BMI ≥ 25 OR 1.034 (0.968–1.105, P = 0.32), and the broad other-heart-disease "
    "endpoints OR 0.993 (0.953–1.035, P = 0.75) and OR 0.879 (0.728–1.061, P = 0.18; 2,192 cases). In the newly "
    "released BBJ 2025 HF suite⁶⁰, the all-cause HF estimate reproduced the legacy BBJ null (OR 0.925, "
    "0.859–0.996, P = 0.039), and subtype stratification separated cleanly by phenotype without approaching the "
    "discovery estimate: HFpEF OR 0.865 (0.776–0.964, P = 0.009; nominally inverse) and HFrEF OR 1.002 "
    "(0.862–1.165, P = 0.98). The access-restricted MVP HFpEF/HFrEF GWAS⁶¹ remains the decisive European subtype "
    "test and is reserved as a pre-specified validation target.")

RES_H2_HETEROG = "Heterogeneity is structured chiefly by ancestry"
RES_P_HETEROG = ("Pooling all ten independent cohort–endpoint estimates gave an attenuated fixed-effect summary "
    "(OR 1.013, 95% CI 0.991–1.035) with persistent heterogeneity (I² = 77%). Subgroup meta-analysis showed that "
    "heterogeneity is organized by ancestry rather than by endpoint definition: the six European estimates averaged "
    "OR 1.036 (1.011–1.062, P = 0.004; within-Europe I² = 73%, driven by the discovery cohort), whereas the four "
    "East Asian estimates were homogeneous and consistently inverse (OR 0.922, 0.878–0.968, P = 0.001; I² = 0%) "
    "(Fig. 8A). An illustrative weighted meta-regression confirmed that the ancestry, strictness, and HFpEF "
    "moderators jointly accounted for about half of the between-cohort variance (QM P = 2.4×10⁻⁴; R² analogue "
    "50%), with the East-Asian coefficient negative (−0.100 log-OR, P = 0.12); a cohort-level funnel plot was "
    "approximately symmetric (Egger P = 0.23; Fig. 8B), arguing against small-study–type bias as the source of "
    "heterogeneity.")

RES_H2_INTERORGAN = "An interorgan ligand–receptor map of the lung–heart axis"
RES_P_INTERORGAN = ("To convert the lung–heart axis from a concept into an explicit network, we scored 31 curated "
    "ligand–receptor pairs between COPD lung sender cell types and failing-heart receiver compartments (1,428 "
    "directed edges; Fig. 8D, Supporting Information Table S11). In disease, inferred signalling from lung "
    "capillary endothelium through SERPINE1 → PLAUR and SERPINE1 → VTN was attenuated (Δ log₂ edge score −2.26 "
    "and −1.71, receiver-averaged), consistent with the cell-type-specific SERPINE1 downregulation in that "
    "compartment, whereas fibroblast-originating S100A9 → TLR4 edges increased (+1.03). Receptor-side changes in "
    "the failing heart were robust to donor-level aggregation: CD36 induction in cardiac fibroblasts "
    "(P = 1.6×10⁻⁴), endothelium (P = 2.4×10⁻³), and macrophages (P = 2.4×10⁻³); macrophage IL1RAP induction "
    "(P = 1.1×10⁻³); endothelial ITGB3 induction (P = 3.7×10⁻³); and macrophage CD44 downregulation "
    "(P = 9.4×10⁻⁵) (all FDR < 0.05; Supporting Information Table S12). These donor-robust receptor changes "
    "nominate THBS1 → CD36 and IL-33 → IL1RL1/IL1RAP as the interorgan edges most likely to be operative in "
    "patients.")

DRUGGABILITY_ADD = (" The extended molecular-layer screen found no usable cis instrument either: across "
    "exon-level, transcript-usage, and splicing QTLs in five cardiorespiratory tissues, no SERPINE1 cis "
    "association survived gene-level permutation adjustment (minimum adjusted P = 0.06)⁶³, extending the absent "
    "cis-genetic anchoring from expression to isoform and splicing biology.")

DISC_PARA1 = ("This study delivers an audited, end-to-end evidentiary chain for the COPD–HF syndrome as a prototype "
    "cross-organ comorbidity axis. Four findings stand out. First, liability to COPD causally increases HF risk "
    "across a comprehensive estimator battery, independently of smoking, adiposity, blood pressure, and lung "
    "function, with 4.4% of the effect flowing through CRP; this discovery estimate withstands both "
    "winner's-curse correction and re-estimation with independent exposure weights. Second, replication is "
    "structured, not random: across ten cohort–endpoint combinations, estimates partition by ancestry (European "
    "OR 1.036; East Asian OR 0.922, I² = 0%), and East Asian HF-subtype stratification yielded a nominally inverse "
    "HFpEF estimate with a null HFrEF estimate. Third, the attributable burden is large (3.61 million cases; "
    "26-fold gradient), demographically driven, resurging in high-SDI settings, and projected to increase "
    "substantially by 2050. Fourth, the mechanistic bridge is a tissue-confined SERPINE1–integrin inflammo-fibrotic "
    "axis with an explicitly mapped interorgan receptor landscape (donor-robust cardiac CD36, IL1RAP, and ITGB3 "
    "changes), phenome-wide safe and druggable in principle.")

DISC_PARA2 = ("The replication architecture is itself a finding. The discovery estimate cannot be dismissed as "
    "winner's curse (corrected OR 1.16) or as instrument weakness (independent FinnGen exposure weights reproduce "
    "OR 1.19); the between-cohort discrepancy is therefore biological or ascertainment-based rather than a "
    "statistical artefact. Three non-exclusive explanations reconcile the estimates. Outcome ascertainment differs "
    "in kind: the discovery GWAS aggregates registry, self-report, and hospital HF, a definition enriched for "
    "chronic and possibly HFpEF-like phenotypes, whereas FinnGen's strict endpoint and BBJ's definitions capture "
    "different case mixes; direct HF-subtype testing in East Asians, however, did not rescue the effect (HFpEF "
    "OR 0.87; HFrEF OR 1.00), arguing against the hypothesis that the COPD signal is specifically an HFpEF "
    "phenomenon in that ancestry⁶⁰. Ancestry is the dominant moderator: East Asian estimates are homogeneous "
    "(I² = 0%) and inverse, whereas European estimates are positive but dispersed; instrument validity, LD "
    "structure, and the biomass- and never-smoker-enriched COPD phenotypes of East Asian cohorts are plausible "
    "contributors. These distinctions generate falsifiable predictions: (i) in European-ancestry HFpEF-stratified "
    "GWAS (e.g., the access-restricted MVP resource⁶¹), the COPD effect should exceed the corresponding HFrEF "
    "effect; (ii) multi-ancestry cohorts with harmonized HFpEF ascertainment should recover a European-type "
    "effect; and (iii) the East Asian inverse signal should concentrate in biomass-related COPD subgroups. We "
    "therefore grade the causal evidence as probable and phenotype- and ancestry-dependent rather than definitive, "
    "and we specify which cohorts will resolve it.")

DISC_PARA3 = ("In relation to prior work, our findings revise the field in three ways. Previous MR studies of COPD "
    "and cardiovascular outcomes were univariable, single-method, and unreplicated, and considered composite "
    "endpoints rather than HF specifically¹⁴; our three-cohort design with meta-analysis and power calculation—here "
    "extended to ten cohort–endpoint combinations with winner's-curse correction, independent-exposure "
    "re-estimation, and meta-regression—sets a new standard, and our MVMR is, to our knowledge, the first "
    "indication that the COPD signal survives adjustment for lung function itself. In burden research, COPD and HF "
    "have been quantified separately¹⁵,¹⁸; our attribution analysis treats their intersection as a disease entity "
    "in its own right, revealing the high-SDI resurgence and ageing-driven growth that single-disease analyses "
    "obscure. In mechanistic research, the few studies linking COPD to myocardial remodelling have been "
    "cross-sectional and imaging- or biomarker-based⁷,¹⁰; by using disease-state single-cell evidence from both "
    "organs, we move the lung–heart axis from association to cellular mechanism, map its interorgan "
    "ligand–receptor topology, and reconcile the apparent PAI-1 contradiction (pro-fibrotic in tissue yet null in "
    "circulating-protein MR) through the tissue confinement of the causal biology⁴³,⁴⁶. More broadly, the study "
    "demonstrates a portable audit paradigm for cross-organ comorbidity syndromes—causal genetics, attributable "
    "burden, and disease-tissue mechanism in one internally consistent framework—that can be reused for other "
    "organ axes.")

DISC_PARA4_ADD = (" The interorgan map adds directionality to this model: attenuation of endothelial "
    "SERPINE1 → PLAUR/VTN edges together with donor-robust induction of the cardiac THBS1 receptor CD36 and the "
    "IL-33 co-receptor IL1RAP suggests that the operative interorgan signals in patients are alarmin- and "
    "THBS1-driven rather than driven by circulating PAI-1, consistent with the plasma–tissue MR disconnect.")

DISC_LIMIT1 = ("Several limitations temper these conclusions. First, the genetic replication was inconsistent "
    "across ancestries; our causal language reflects this, and pooled estimates are sensitive to model choice and "
    "cohort composition (fixed-effect OR 1.013 across ten cohort–endpoint estimates; positive only within the "
    "European subgroup). Second, the discovery MR rests on ten instruments, adequate for robust estimation "
    "(all F > 30) but at the lower bound for some pleiotropy models⁴⁵; full CAUSE model comparison was not "
    "identifiable, MVMR conditional F-statistics were weak (1.8–2.2), and residual confounding cannot be excluded, "
    "so adjusted estimates should be read as supportive rather than definitive. Third, GBD estimates are "
    "model-based: the attribution of HF to COPD follows comparative-risk conventions, which imply a higher "
    "effective relative risk (≈2.28) than the MR estimate (1.15); we report both frameworks rather than blending "
    "them, and the 2021 rebound should be interpreted against COVID-19-era data disruption. Fourth, projections "
    "carry model uncertainty (wide credible intervals; count-space BAPC extrapolates demographic momentum; "
    "Gaussian-approximated intervals). Fifth, mediation paths for IL-6 and SERPINE1 used few instruments (1–3 SNPs "
    "on some paths), so small mediation effects cannot be ruled out. Sixth, single-cell cohorts are modest in size "
    "and cross-sectional, and bulk cohorts are small; cell-level single-cell tests inflate significance through "
    "pseudoreplication⁵⁹, and our donor-level pseudobulk re-testing confirmed only ITGAV (cardiac endothelium) "
    "among the key expression associations, so single-cell findings are presented as hypothesis-generating and the "
    "41-gene program warrants experimental perturbation; the interorgan edge scores are expression-product "
    "inferences, not physical interaction measurements. Seventh, subnational analyses for China require a "
    "dedicated GBD download not yet available. Eighth, the extended replication relies on endpoints not perfectly "
    "harmonized in ascertainment; the BBJ 2025 HFpEF/HFrEF strata are modest (7,154 and 4,254 cases), the "
    "meta-regression is illustrative (k = 10), and three pre-specified analyses—the European HFpEF/HFrEF test in "
    "MVP (dbGaP phs001672)⁶¹, an individual-level PRS→HF incidence validation in UK Biobank, and a PM2.5-exposure "
    "MR layer—are not executable with public data and are therefore reported as registered next analyses rather "
    "than results.")

DISC_NEXT = ("Genetically, the decisive next step is now precisely specified: HFpEF- versus HFrEF-stratified MR in "
    "European ancestry (MVP dbGaP phs001672)⁶¹, multi-ancestry GWAS with harmonized ascertainment, and an "
    "individual-level PRS→HF incidence validation in UK Biobank connecting population burden to personal risk "
    "stratification; a positive European HFpEF replication would convert our probable causal estimate into a "
    "definitive one, whereas confirmation of the East Asian inverse pattern would redefine the syndrome's ancestry "
    "boundaries. Mechanistically, perturbation of the SERPINE1 node (antisense, small-molecule PAI-1 inhibition, or "
    "PLAUR blockade) and of the donor-robust receptor axes (CD36, IL1RAP) in cardiac fibroblast and endothelial "
    "systems will test the causal sufficiency of the axis, and tissue- versus plasma-resolved PAI-1 "
    "quantification in patient cohorts will test the compartment model. Clinically, a pragmatic trial of "
    "natriuretic-peptide screening in older patients with COPD⁴⁹ would determine whether the burden we quantify "
    "can be intercepted earlier.")

CONCLUSIONS = ("In conclusion, COPD is a probable, modifiable, mechanistically tractable driver of heart failure, "
    "with a phenotype- and ancestry-dependent genetic effect. Its attributable burden of 3.61 million people in "
    "2021, geographically unequal and resurging precisely where cardiovascular medicine is most advanced, will "
    "increase substantially by 2050. The molecular bridge between lung and heart—a SERPINE1-centred, "
    "tissue-confined inflammo-fibrotic axis with low-grade systemic inflammation as its relay and a donor-robust "
    "cardiac receptor landscape (CD36, IL1RAP, ITGB3)—is supported by bulk and single-cell evidence in diseased "
    "human tissue, phenome-wide safe, and druggable in principle. The structured heterogeneity of the genetic "
    "effect, far from weakening the study, defines the next generation of cohorts needed and tempers translation "
    "with appropriate caution. Surveillance of older patients with COPD for heart failure, aggressive tobacco "
    "control, and tissue-targeted anti-fibrotic strategies are the actionable outputs of this work, and the "
    "causal–burden–mechanism audit on which they rest is reusable for other cross-organ comorbidity syndromes.")

RESULTS_REPLICATION_ENDING = (" We interpret this pattern as a genuine but phenotype- and ancestry-dependent causal "
    "signal and formalize its sources below.")

DATA_AVAIL = ("All source datasets are publicly available or freely accessible through their providers: IEU "
    "OpenGWAS (ebi-a-GCST90018807, ebi-a-GCST009541), FinnGen R12, BioBank Japan (including the BBJ 2025 HF suite: "
    "GWAS Catalog GCST90668009, GCST90668010, GCST90668011), GBD 2021 (IHME), GEO (GSE57148, GSE57338), CELLxGENE, "
    "GTEx v8, and deCODE. The MVP HFpEF/HFrEF GWAS is available through dbGaP (phs001672) under controlled access. "
    "All derived result tables, analysis scripts, and execution logs are provided in the Supporting Information "
    "and deposited at GitHub (github.com/wang76qiang/copd-hf-causality-burden-mechanism) and Zenodo "
    "(DOI 10.5281/zenodo.21943762).")

NEW_REFS = [
    "60. Enzan N, Miyazawa K, Koyama S, et al. Genome-wide analysis of heart failure yields insights into disease heterogeneity and enables prognostic prediction in the Japanese population. Nat Commun 2025; 16: 9680.",
    "61. Joseph J, Liu C, Hui Q, et al. Genetic architecture of heart failure with preserved versus reduced ejection fraction. Nat Commun 2022; 13: 7753.",
    "62. Zhong H, Prentice RL. Bias-reduced estimators and confidence intervals for odds ratios in genome-wide association studies. Biostatistics 2008; 9: 621–34.",
    "63. Kerimov N, Hayhurst JD, Peikova K, et al. A compendium of uniformly processed human gene expression and splicing quantitative trait loci. Nat Genet 2021; 53: 1290–9.",
]

FIG8_LEGEND = ("Figure 8. Replication architecture, heterogeneity forensics, and interorgan crosstalk. "
    "(A) Cohort–endpoint forest plot with ancestry subgroup diamonds (winner's-curse-corrected discovery shown as "
    "a non-pooled sensitivity, square marker). (B) Cohort-level funnel plot with pseudo-95% confidence limits "
    "around the fixed-effect summary (Egger test illustrative, k = 10). (C) Per-SNP winner's-curse correction of "
    "discovery exposure estimates; the IVW estimate is unchanged (OR 1.151 → 1.156) and independent "
    "FinnGen-exposure re-estimation reproduces the effect (OR 1.19). (D) Lung→heart interorgan ligand–receptor "
    "delta heatmap (receiver-averaged Δ log₂ edge score, COPD–failing vs normal–normal); only edges with "
    "non-negligible expression (max score ≥ 0.005) are shown.")
