# Evidence Artifact: ONCO-DRIVER-001

**Test Unit ID:** ONCO-DRIVER-001
**Algorithm:** Driver Mutation Detection (20/20 rule)
**Date Collected:** 2026-06-14

---

## Online Sources

### Vogelstein et al. (2013) — Cancer Genome Landscapes (Science)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3749880/ (PMC full text, behind CAPTCHA at fetch time) ; DOI https://doi.org/10.1126/science.1235122 (HTTP 403 on direct fetch)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, originating source of the 20/20 rule)

**How retrieved:** Web search `Vogelstein 2013 Science Cancer Genome Landscapes 20/20 rule oncogene tumor suppressor`. The PMC HTML and the science.org DOI page could not be fetched directly (CAPTCHA / 403), so the rule's exact wording was extracted from open-access secondary sources that quote the primary text verbatim (Tokheim 2020 PMC7703750; OncodriveROLE 2014; Miller 2017 Oncotarget) — see those entries. The primary-source claims below are the ones those open sources attribute directly to Vogelstein et al. 2013.

**Key Extracted Points:**

1. **20/20 rule — oncogene:** To classify a gene as an oncogene, more than 20% of the recorded mutations in the gene are at recurrent positions and are missense (activating).
2. **20/20 rule — tumor suppressor:** To classify a gene as a tumor suppressor gene (TSG), more than 20% of the recorded mutations in the gene are inactivating (truncating / loss of function).
3. **Leniency:** The rule is lenient; all well-documented cancer genes far surpass these criteria.
4. **IDH1 worked example:** Nearly all IDH1 mutations occur at the identical amino acid (codon 132, e.g. R132H); the 20/20 rule therefore classifies IDH1 as an oncogene rather than a TSG.

### Tokheim & Karchin (2020) — Somatic selection distinguishes oncogenes and tumor suppressor genes (Bioinformatics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7703750/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed; quotes the primary rule verbatim)

**How retrieved:** WebFetch of the PMC URL with a prompt asking for the exact 20/20 rule definition and inactivating mutation types.

**Key Extracted Points:**

1. **Verbatim rule:** "OGs have >20% mutations causing missense changes at recurrent positions and TSGs have >20% mutations causing inactivating changes." (20/20+ extends the original 20/20 rule with the same thresholds.)
2. **Inactivating types:** inactivating mutations are protein-truncating: nonsense and frame-shifting mutations (frameshifts).

### Schroeder et al. (2014) — OncodriveROLE classifies cancer driver genes (Bioinformatics)

**URL:** https://academic.oup.com/bioinformatics/article/30/17/i549/201062
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed reference implementation describing the rule)

**How retrieved:** WebFetch of the Oxford Academic URL asking for the rule and which mutation types count as truncating.

**Key Extracted Points:**

1. **Rule restated:** Genes with ≥20% truncating mutations are tumor suppressors (loss of function); genes with >20% of missense mutations in recurrent positions are oncogenes (activating). Introduced by Vogelstein et al. (2013).
2. **Truncating types (explicit list):** "truncating mutations include mutations causing a frameshift, a gained or lost stop codon as well as mutations in splice donor or acceptor sites."

### Miller et al. (2017) — Identification and analysis of mutational hotspots in oncogenes and tumour suppressors (Oncotarget)

**URL:** https://www.oncotarget.com/article/15514/text/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed; restates rule and hotspot definition)

**How retrieved:** WebFetch of the open-access Oncotarget text asking for the rule, hotspot/recurrent definition, the IDH1 codon, and truncating types.

**Key Extracted Points:**

1. **Rule restated:** "if 20% of all mutations observed within a gene are truncations, then the gene is likely to be a tumour suppressor. Similarly, if 20% of all missense mutations occur at a single position in the sequence, the gene is predicted to be an oncogene."
2. **Recurrent / hotspot position:** a recurrent position requires "at least two mutations of the same class" at an identical location — i.e. recurrence = the same protein position observed ≥ 2 times.
3. **IDH1:** codon 132 is the hotspot (R132H).
4. **Truncating types:** truncations encompass "nonsense mutations and frameshift insertions and deletions."

---

## Documented Corner Cases and Failure Modes

### From Tokheim & Karchin (2020) — PMC7703750

1. **Passenger truncations in oncogenes:** "random mutational processes may introduce protein-truncating mutations into OGs, which increase in frequency via genetic drift with no significant impact on tumor fitness and mislead annotations." Consequence: the 20/20 rule is a heuristic, not a statistical test; it can misclassify genes with few mutations.

### From Schroeder et al. (2014) — OncodriveROLE

1. **Low-recurrence drivers:** "the rule fails to identify drivers included in newer catalogs, mostly the lowly recurrent ones." Consequence: a gene with few mutations may satisfy neither the OG nor the TSG criterion and is left unclassified (Ambiguous), which is the correct conservative behavior for the rule.

### From Vogelstein et al. (2013) (via secondary quotes)

1. **Single-amino-acid recurrence (IDH1):** when virtually all missense mutations fall on one codon, the recurrent-missense fraction is ~1.0 ≫ 0.20, giving a clear oncogene call.

---

## Test Datasets

### Dataset: IDH1 (oncogene archetype)

**Source:** Vogelstein et al. (2013), via Miller et al. (2017) — IDH1 hotspot at codon 132.

| Parameter | Value |
|-----------|-------|
| Mutations | 10 missense, all at codon 132 (R132H) |
| Truncating fraction | 0/10 = 0.00 |
| Recurrent-missense fraction | 10/10 = 1.00 (codon 132 seen 10 times) |
| Expected classification | Oncogene (recurrent-missense 1.00 > 0.20) |

### Dataset: Tumor-suppressor archetype (dispersed truncating)

**Source:** Vogelstein et al. (2013) / OncodriveROLE (2014) — TSGs accumulate inactivating mutations dispersed along the gene.

| Parameter | Value |
|-----------|-------|
| Mutations | 5 nonsense + 2 frameshift + 1 splice = 8 truncating, all at distinct positions; 2 missense at distinct positions |
| Truncating fraction | 8/10 = 0.80 |
| Recurrent-missense fraction | 0/10 = 0.00 (no missense position seen ≥ 2 times) |
| Expected classification | TumorSuppressor (truncating 0.80 > 0.20) |

### Dataset: Mixed / boundary

**Source:** derived from the rule (threshold = 0.20 strict `>`).

| Parameter | Value |
|-----------|-------|
| 10 mutations: exactly 2 truncating, 8 missense at distinct positions | truncating fraction = 0.20 (NOT > 0.20) → not a TSG |
| 10 mutations: 3 truncating, 7 missense at distinct positions | truncating fraction = 0.30 > 0.20 → TumorSuppressor |

---

## Assumptions

1. **ASSUMPTION: Tie-break when both criteria are met** — If a gene satisfies both the oncogene (>20% recurrent missense) and the TSG (>20% truncating) criteria, the sources do not prescribe a single label. We classify by the larger of the two fractions, and report Ambiguous on an exact tie. This affects only genes that pass both thresholds; the two archetypes (IDH1, dispersed-truncating TSG) are unaffected. Justification: Vogelstein et al. note well-documented genes "far surpass" one criterion, so a dual pass is atypical; choosing the dominant signal is the least-surprising deterministic resolution.
2. **ASSUMPTION: Strict `>` 20% threshold** — Vogelstein/Tokheim/Miller state ">20%"; OncodriveROLE writes "≥20%" for TSGs. We use strict `>` 0.20 for both, matching the primary source (Vogelstein) and the verbatim Tokheim quote. A fraction of exactly 0.20 is therefore not sufficient.
3. **ASSUMPTION: ScoreDriverPotential proxy** — The checklist names CADD/SIFT/PolyPhen for ScoreDriverPotential, but those are externally trained models that cannot be retrieved and reproduced here (forbidden to fabricate). We instead return the 20/20-rule driver-signal fraction in [0,1] (max of the two criterion fractions) as a transparent, source-derived score, and document that external pathogenicity scores are caller-supplied / not implemented.

---

## Recommendations for Test Coverage

1. **MUST Test:** IDH1-style all-recurrent-missense gene classified Oncogene with recurrent-missense fraction 1.00 — Evidence: Vogelstein 2013 / Miller 2017 IDH1 codon 132.
2. **MUST Test:** Dispersed-truncating gene (nonsense+frameshift+splice) classified TumorSuppressor with truncating fraction > 0.20 — Evidence: Vogelstein 2013 / OncodriveROLE truncating types.
3. **MUST Test:** Boundary — truncating fraction exactly 0.20 is NOT a TSG (strict `>`) — Evidence: Vogelstein ">20%".
4. **MUST Test:** MatchCancerHotspots flags a mutation whose (gene, position) is in the caller-supplied hotspot set and not otherwise — Evidence: Miller 2017 recurrent-position definition.
5. **SHOULD Test:** Gene meeting neither criterion → Ambiguous (low recurrence) — Rationale: OncodriveROLE documented failure mode.
6. **SHOULD Test:** IdentifyDriverMutations returns a subset of input variants (invariant driver ⊆ somatic) — Rationale: Registry invariant.
7. **COULD Test:** Recurrence requires the same position ≥ 2 times; a singleton missense is not recurrent — Rationale: Miller hotspot definition.

---

## References

1. Vogelstein B, Papadopoulos N, Velculescu VE, Zhou S, Diaz LA Jr, Kinzler KW. 2013. Cancer Genome Landscapes. Science 339(6127):1546–1558. https://doi.org/10.1126/science.1235122 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3749880/)
2. Tokheim C, Karchin R. 2020. Somatic selection distinguishes oncogenes and tumor suppressor genes (20/20+). Bioinformatics 36(6):1712–1719. https://doi.org/10.1093/bioinformatics/btz759 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7703750/)
3. Schroeder MP, Rubio-Perez C, Tamborero D, Gonzalez-Perez A, Lopez-Bigas N. 2014. OncodriveROLE classifies cancer driver genes in loss of function and activating mode of action. Bioinformatics 30(17):i549–i555. https://doi.org/10.1093/bioinformatics/btu467 (https://academic.oup.com/bioinformatics/article/30/17/i549/201062)
4. Miller ML, Reznik E, Gauthier NP, et al. 2017. Identification and analysis of mutational hotspots in oncogenes and tumour suppressors. Oncotarget 8(20):33321–33333. https://doi.org/10.18632/oncotarget.15514 (https://www.oncotarget.com/article/15514/text/)

---

## Change History

- **2026-06-14**: Initial documentation.
