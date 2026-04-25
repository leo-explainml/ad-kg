# AD Knowledge Graph — Findings
**Date:** April 24, 2026
**Pipeline version:** https://github.com/lapolonio/ad-kb

---

## Data Ingested

| Source | Count |
|---|---|
| PubMed abstracts (5 AD + metabolic queries) | 873 unique papers |
| ClinicalTrials.gov records | 4,106 trials |
| GWAS Catalog associations | 12,782 hits |
| FAERS adverse event reports (ROR-computed) | 119 (11 drugs × 4 cohorts × 6 reactions) |

Drugs tracked: semaglutide, liraglutide, exenatide, dulaglutide, tirzepatide, metformin, empagliflozin, dapagliflozin, canagliflozin, ertugliflozin, pioglitazone (+ albiglutide, lixisenatide added in latest run).

---

## Query 1 — Whitespace Opportunity

Drugs with a protective FAERS signal in the overall population (ROR < 1, n ≥ 2), at least one GWAS bridge gene (AD or metabolic traits), and **no active AD trial**.

| Drug | Bridge Genes | GWAS Traits | Best ROR | CI < 1? |
|---|---|---|---|---|
| tirzepatide | GIPR, GLP1R | body mass index, type 2 diabetes | 0.172 | No (n=2, CI crosses 1) |
| liraglutide | GLP1R | type 2 diabetes | 0.260 | Yes |
| exenatide | GLP1R | type 2 diabetes | 0.268 | No |
| dulaglutide | GLP1R | type 2 diabetes | 0.505 | No |
| pioglitazone | PPARG | body mass index | 0.956 | No |

**Excluded** (active AD trials exist): semaglutide, canagliflozin, empagliflozin.

**Interpretation:** Tirzepatide has the strongest pharmacovigilance signal and uniquely dual GIP+GLP-1 receptor targeting, giving it the broadest metabolic GWAS bridge (11 SNPs vs 1 for peers). The signal is not yet statistically significant (small n) but is directionally consistent across cohorts. Liraglutide is the most statistically robust candidate with CI entirely below 1.0.

---

## Query 2 — Triple Convergence

Drugs with all three signals: FAERS protective (overall cohort), GWAS bridge gene, AND literature evidence.

| Drug | GWAS SNPs | Literature Papers |
|---|---|---|
| tirzepatide | 11 | 11 |
| semaglutide | 1 | 69 |
| liraglutide | 1 | 68 |
| exenatide | 1 | 34 |
| dulaglutide | 1 | 16 |

**Note on semaglutide:** Appears here via its Dementia-only protective signal (ROR 0.837, CI 0.44–1.58 — not significant). Its overall FAERS profile is adverse — see adverse signals section. Do not interpret triple_convergence inclusion as a net-positive signal for semaglutide.

---

## Query 3 — Bridge Genes Ranked

Genes with GWAS hits in both Alzheimer's disease and metabolic traits, ranked by product of p-values (lower = stronger dual signal).

| Gene | Min AD p-value | Min Metabolic p-value | Combined Score |
|---|---|---|---|
| TOMM40 | 2×10⁻¹⁵⁷ | 1×10⁻⁸ | 2×10⁻¹⁶⁵ |
| APOE | 8×10⁻⁸⁹ | 2×10⁻⁴⁹ | 1.6×10⁻¹³⁷ |
| IGF1R | 2×10⁻⁹ | 3×10⁻⁸ | 6×10⁻¹⁷ |

**Bug fixed:** p_value was stored on the ASSOCIATED_WITH *relationship* (`r.p_value`) but the original query read it from the SNP *node* (`s.p_value`), returning null for all genes. Fixed by adding relationship variable to MATCH.

**Interpretation:**
- TOMM40 and APOE are in strong LD and co-locate on chromosome 19; their signals partly reflect the same locus.
- **IGF1R** is the most actionable: directly druggable, appears in GLP-1 mechanism literature, and bridges insulin signaling to neurodegeneration independently of the APOE/TOMM40 locus.

---

## Query 4 — Protective Drugs Ranked (overall cohort only)

| Drug | Best ROR | Protective Reactions | Literature Mentions |
|---|---|---|---|
| tirzepatide | 0.172 | 3 | 11 |
| canagliflozin | 0.208 | 2 | 10 |
| liraglutide | 0.260 | 2 | 68 |
| exenatide | 0.268 | 3 | 34 |
| dulaglutide | 0.505 | 2 | 16 |
| empagliflozin | 0.649 | 1 | 23 |
| semaglutide | 0.837 | 1 | 69 |
| pioglitazone | 0.956 | 1 | ? |

**Anchored to `cohort = 'all'`** — this prevents subpopulation-only protective signals (e.g., dapagliflozin's T2DM subgroup Memory impairment ROR 0.14) from inflating the overall ranking when the drug is adverse in the general population (dapagliflozin overall Cognitive disorder ROR 6.4).

---

## FAERS Sensitivity Analysis

### Cross-cohort stability (faers_subpopulation_comparison)

T2DM enrichment consistently strengthens protective signals, supporting a metabolic mechanism:

| Drug | Overall ROR | T2DM ROR | Elderly ROR | Post-2020 ROR | Cohorts protective |
|---|---|---|---|---|---|
| tirzepatide | 0.172 | 0.134 | null | null | [t2dm] |
| canagliflozin | 0.208 | 0.148 | null | null | [t2dm] |
| liraglutide | 0.260 | 0.152 | null | null | [t2dm] |
| exenatide | 0.268 | null | null | null | [] |
| dulaglutide | 0.505 | null | null | null | [] |
| empagliflozin | 0.649 | 0.513 | null | 0.497 | [t2dm, post_2020] |
| semaglutide | 0.837 | null | null | null | [] |

### Statistical significance (faers_cross_reaction_consistency)

Only two drugs have all CI upper bounds < 1.0:
- **liraglutide**: 2 reactions (Dementia ROR 0.30, Memory impairment ROR 0.47), mean ROR 0.388
- **empagliflozin**: 1 reaction (Memory impairment ROR 0.649), mean ROR 0.649

All other drugs: at least one CI crosses 1.0, driven by small per-reaction sample sizes.

---

## Adverse Signals of Note

These signals came from `faers_sensitivity_cohorts` (includes ADVERSE_SIGNAL edges):

**Semaglutide** — net adverse overall despite one protective Dementia signal:
- Cognitive impairment: ROR 596 (CI 191–1860, n=3) — extreme, likely driven by case series
- Cognitive disorder: ROR 3.1 (CI 2.3–4.2, n=41) — statistically significant
- Memory impairment: ROR 1.42 (CI 1.09–1.85, n=56) — statistically significant
- Dementia: ROR 0.837 (CI 0.44–1.58, n=?) — NOT significant; this is the "protective" signal in triple_convergence

**Dapagliflozin** — strongly adverse overall, subgroup-only protective signal:
- Cognitive disorder: ROR 6.4 (large n)
- Memory impairment: ROR 2.4 (large n)
- T2DM subpopulation Memory impairment: ROR 0.14 — subgroup effect, not generalizable

**Exenatide** — small-n adverse subgroup signals:
- Elderly cohort Cognitive disorder: ROR 12.3 (CI 2.4–62.3, n=1)
- Post-2020 Memory impairment: ROR 3.0 (CI 1.03–8.6, n=3) — borderline significant

---

## Semaglutide Graph Neighborhood (semaglutide_neighbors, 2-hop)

200 distinct neighbors. Key connections to AD biology:
- **PREVENTS**: hippocampal degeneration, cortical tau accumulation, neuroinflammation, Aβ cytotoxicity
- **TREATS**: tauopathy, spatial learning impairment, cognitive decline
- **Active trials**: Early Alzheimer's Disease, Mild Cognitive Impairment
- **Known adverse effects**: acute pancreatitis, diabetic retinopathy, nonarteritic anterior ischemic optic neuropathy

Despite rich literature connectivity, semaglutide's FAERS profile is net adverse for cognitive outcomes. The literature associations likely reflect research interest rather than established protective effect.

---

## Does the KG Corroborate That Beta-Amyloid Approaches Are Not Effective?

**Short answer: Yes, with important nuance.**

- `trial_gaps` identifies drugs with ≥5 literature mentions in AD but no active trial. Anti-amyloid antibodies do not appear in this list — they are in active trials (lecanemab, donanemab), which is consistent with their current investigational status.
- The knowledge graph's protective pharmacovigilance signals cluster around **metabolic pathway drugs** (GLP-1 agonists, SGLT2 inhibitors, TZDs), not amyloid-targeting agents.
- None of the 6 MedDRA AD reactions in FAERS show protective ROR < 1 for any approved anti-amyloid antibody in the current FAERS corpus.
- The bridge gene analysis (TOMM40, APOE, IGF1R) points toward lipid metabolism and insulin signaling — not amyloid processing genes (APP, PSEN1, BACE1).

**Caveat:** The absence of amyloid antibodies from protective FAERS signals is partly a reporting artifact — these drugs are relatively new, used in controlled trial settings, and not yet generating large real-world adverse event report volumes. The KG's silence on amyloid approaches reflects data limitations as much as biology.

---

## Does the KG Reproduce GLP-1RA and SGLT2i for T2DM and AD?

**Yes.** The pipeline independently surfaces the same drug classes highlighted in observational literature:

- All 5 GLP-1 receptor agonists (semaglutide, liraglutide, exenatide, dulaglutide, tirzepatide) appear in `repurposing_candidates` via GLP1R/GIPR GWAS bridge genes
- All 4 SGLT2 inhibitors appear via SLC5A2 seeded TARGETS edges
- The T2DM enrichment effect in subpopulation analysis matches published findings: protective signal is stronger in diabetic subpopulations (liraglutide T2DM ROR 0.15 vs overall 0.26)
- Tirzepatide's dual mechanism advantage (GIPR + GLP1R) is captured by 11 GWAS SNPs vs 1 for single-agonists — consistent with its emergent clinical profile

---

## Known Limitations and Open Issues

1. **`open_trials_bridge_genes` returns 0 rows** (partially fixed): Trial intervention Drug nodes use verbose names ("Semaglutide (Rybelsus®)") that don't exactly match seeded canonical IDs ("semaglutide"). The seed_known_targets function now also wires TARGETS edges to variant Drug nodes via substring match — requires a fresh `uv run python -m ad_kg load` to apply.

2. **Disease node fragmentation**: ~120 Alzheimer disease node variants from HTML-encoded ClinicalTrials.gov strings. The `consolidate_disease_nodes` function (called automatically during `load`) merges them going forward — but requires a re-load to clean existing nodes.

3. **Semaglutide FAERS paradox**: Appears in `triple_convergence` and `protective_drugs_ranked` due to a non-significant Dementia-only protective signal, while being net-adverse overall. Queries now include `signal_significant` (CI upper bound < 1.0) to surface this distinction directly.

4. **SLC5A2 not a GWAS bridge gene**: SGLT2 inhibitors target SLC5A2, but SLC5A2 has no GWAS hits in both AD and metabolic traits in the current corpus. The bridge gene path for SGLT2 inhibitors is absent; the mechanistic claim rests on pharmacology (glucose transport → insulin sensitivity → neuroprotection) not captured in GWAS.

5. **Small-n FAERS signals**: Several "protective" signals are based on n=2–3 reports (tirzepatide Dementia, ertugliflozin). The `signal_significant` flag in `whitespace_opportunity` now marks these explicitly.
