# Statistical Analysis Plan (SAP) — Pre-registered next analyses for the COPD–HF study

Version 1.0 | Locked 2026-09-04 | To be appended to the Supporting Information as Note 5.

This SAP formally specifies three analyses that are pre-registered but not executable with publicly available data at submission. Each analysis lists hypotheses, data, instruments, estimators, and success/failure criteria, so that the falsifiable predictions in the main text carry registration force.

## Analysis R1 — European-ancestry HFpEF- versus HFrEF-stratified MR

- **Objective.** Test prediction (i): in European ancestry, the COPD causal effect on HFpEF exceeds that on HFrEF.
- **Data.** MVP HFpEF (n = 19,589 cases) and HFrEF (n = 19,495) GWAS (Joseph et al., Nat Commun 2022), via dbGaP phs001672 (controlled access; application to be filed).
- **Instruments.** The ten discovery COPD instruments (Table S1), harmonized per the pipeline rules (aligned/flipped/strand; palindromic-frequency rule |ΔEAF| ≤ 0.15); SNPs failing harmonization are reported, not imputed.
- **Estimators.** Fixed-effect IVW primary; weighted median, MR-Egger, weighted mode sensitivity; heterogeneity (Cochran Q, I²); Egger intercept.
- **Success criterion (supports the phenotype-specific hypothesis).** HFpEF IVW OR significantly > 1 AND formally larger than HFrEF (difference of log-ORs, two-sided P < 0.05, correlated-estimate adjustment via the shared-control covariance or bootstrap over SNPs).
- **Failure criterion.** HFpEF null, or HFpEF ≤ HFrEF → the fibrosis–HFpEF causal hypothesis is downgraded to "not supported in any ancestry"; the mechanism layer is re-anchored to the donor-robust integrin/receptor signals only.
- **Multiplicity.** Two endpoints, one primary contrast; no adjustment beyond the stated two-sided α = 0.05 for the contrast.

## Analysis R2 — Individual-level COPD-PRS → incident HF validation (UK Biobank)

- **Objective.** Test whether polygenic COPD risk predicts incident HF at the individual level, connecting population burden to personal risk stratification.
- **Data.** UK Biobank (requires approved application): genotyped participants without prevalent HF at baseline; incident HF ascertained from HES/death registry (ICD-10 I50); follow-up through the latest available censoring.
- **Instruments/score.** COPD PRS from the discovery GWAS weights (primary: the 10-SNP set; secondary: genome-wide PRS via LDpred2/PRS-CS using the discovery summary statistics), computed in European-ancestry participants; non-European participants analysed separately as exploratory.
- **Estimators.** Cox proportional hazards (time since baseline) of incident HF on PRS (per SD and quintiles), adjusted for age, sex, assessment centre, genotyping array, 40 PCs; sensitivity: (a) excluding prevalent COPD; (b) adjusting for smoking status/pack-years, BMI, SBP; (c) Fine–Gray with death as competing risk; (d) COPD-stratified analysis (does PRS predict HF among diagnosed COPD — the clinically actionable stratum).
- **Success criterion.** HR per SD > 1 with 95% CI excluding 1 in the primary model, and monotone quintile gradient.
- **Failure criterion.** Null or non-monotone → the individual-level risk-stratification claim is withdrawn; surveillance recommendations are re-based on the burden layer alone.
- **Multiplicity.** One primary endpoint (incident HF) and one primary score (10-SNP PRS); secondary analyses flagged exploratory.

## Analysis R3 — PM2.5-exposure MR layer

- **Objective.** Test whether ambient PM2.5 is causally upstream of COPD and HF, adding a modifiable environmental driver to the causal path (relevant to the high-SDI resurgence and East-Asian burden).
- **Data.** PM2.5 GWAS: UK Biobank ESCAPE-based residential PM2.5 (OpenGWAS ukb-b-10817, n = 423,796; requires a personal OpenGWAS token). Outcomes: discovery COPD (ebi-a-GCST90018807) and discovery HF (ebi-a-GCST009541).
- **Instruments.** P < 5×10⁻⁸ primary; if < 3 instruments pass, a pre-specified relaxed threshold (P < 1×10⁻⁵) is used and flagged; clumping r² < 0.001, 10 Mb, European 1000 Genomes reference.
- **Estimators.** IVW primary; MR-Egger and weighted median sensitivity; MVMR adjusting for smoking initiation (GCST90029014) to test independence from smoking.
- **Success criterion.** IVW P < 0.05 with direction-consistent sensitivity estimators, and persistence in the smoking-adjusted MVMR.
- **Failure criterion.** Null, or attenuation to null after smoking adjustment → the environmental-driver claim is limited to the ecological/GBD level and stated as such.
- **Multiplicity.** Two outcomes; Bonferroni α = 0.025.

## Governance

- Any deviation from this SAP (instrument set, thresholds, estimators, criteria) will be documented in a versioned amendment.
- Results, whether positive or negative, will be reported with the same audit discipline as the present manuscript (deposited scripts, logs, and result files).
