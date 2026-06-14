# Evidence Artifact: ONCO-SV-001

**Test Unit ID:** ONCO-SV-001
**Algorithm:** Somatic Complex Rearrangement Classification (Chromothripsis Inference)
**Date Collected:** 2026-06-15

---

## Online Sources

### Korbel & Campbell (2013) — "Criteria for Inference of Chromothripsis in Cancer Genomes", Cell 152:1226–1236

**URL:** https://pubmed.ncbi.nlm.nih.gov/23498933/ (PMID 23498933); full criteria discussed in the open review PMC3861665 (below)
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary paper)

**Retrieved via:** WebSearch query `"Korbel Campbell 2013 criteria hallmarks chromothripsis Cell"` → returned the Cell article and citing reviews; the criteria content was extracted by fetching the open-access review PMC3861665 (Maher & Wilson 2012, "Chromothripsis and beyond"), which enumerates the Korbel & Campbell hallmark criteria verbatim.

**Key Extracted Points:**

1. **Six hallmark criteria:** (A) clustering of breakpoints; (B) regularity of oscillating copy-number states; (C) interspersed regions of retained / lost heterozygosity; (D) prevalence of rearrangements affecting a specific haplotype; (E) randomness of DNA fragment order and fragment joins; (F) ability to walk the derivative chromosome (invariant alternation between head and tail sequences).
2. **Clustering test (criterion A):** "the null hypothesis of random breakpoints predicts that the distance between breakpoints should be distributed exponentially." (verbatim, via PMC3861665) — i.e. clustering is inferred when observed inter-breakpoint distances deviate from the exponential null toward many short distances.
3. **Oscillating copy-number states (criterion B):** chromothripsis shows "an alternating copy number profile" with the "hallmark two copy number states observed in chromothripsis" (verbatim, via PMC3861665). The CN profile oscillates between (canonically) two states rather than showing progressively rising amplification levels.
4. **Heterozygosity (criterion C):** interspersed retention and loss of heterozygosity along the affected region (verified from SNP-array / sequencing genotype information).
5. **Clustering necessary-but-not-sufficient:** clustering alone is "necessary but not sufficient" for a chromothripsis diagnosis (PMC3861665).

### Maher & Wilson / "Chromothripsis and beyond" review (PMC3861665)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3861665/
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed review citing primaries)

**Retrieved via:** WebFetch of https://pmc.ncbi.nlm.nih.gov/articles/PMC3861665/ with a prompt requesting the exact CN-switch threshold and clustering test.

**Key Extracted Points:**

1. **First-pass operational threshold (Magrangeas et al. 2011):** "a common sense approach was to adopt operational definitions that would require, say, 10, 20, or 50 oscillating copy number changes" for first-pass screening (verbatim). The lowest commonly cited first-pass cutoff is **10 oscillating copy-number changes**.
2. **Two-state hallmark:** the canonical chromothripsis signature is an alternating profile between **two** copy-number states (loss and retention), contrasting with gradual models that produce multiple amplified levels.
3. **Exponential null for breakpoints:** "The degree of clustering can be assessed by simple statistical methods because the null hypothesis of random breakpoints predicts that the distance between breakpoints should be distributed exponentially."

### Cortés-Ciriano et al. (2020) — "Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing", Nature Genetics 52:331–341

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7058534/ (PMCID PMC7058534); DOI 10.1038/s41588-019-0576-7
**Accessed:** 2026-06-15
**Authority rank:** 1 (peer-reviewed primary, PCAWG pan-cancer)

**Retrieved via:** WebSearch query `Cortes-Ciriano 2020 ... 2,658 cancer genomes copy number oscillating states`, then WebFetch of PMC7058534 requesting the exact high/low-confidence segment-count and SV-count thresholds.

**Key Extracted Points:**

1. **High-confidence call:** chromothripsis high-confidence calls "display oscillations between two states in at least seven adjacent segments" (verbatim). → **≥ 7 adjacent oscillating segments**.
2. **Low-confidence call:** low-confidence calls "involve between four and six segments" (verbatim). → **4–6 adjacent oscillating segments**.
3. **Canonical event:** "more than 60% of the CN segments in the affected region oscillated between two states" (verbatim) characterises canonical (two-state) chromothripsis.
4. **Minimum SV burden:** focal events "comprising fewer than six SVs" were excluded → a practical minimum of **6 clustered intrachromosomal SVs**.
5. **Multi-state caveat:** "a considerable fraction of the events involves multiple chromosomes as well as additional structural alterations"; canonical profiles oscillate between two CN states but real events may show more.

---

## Documented Corner Cases and Failure Modes

### From Korbel & Campbell 2013 / PMC3861665

1. **Clustering necessary but not sufficient:** a cluster of breakpoints alone does not establish chromothripsis; the oscillating-CN and randomness criteria must also hold. A profile with clustered breakpoints but progressively rising CN (>2 ascending states, gradual amplification) is NOT chromothripsis (it suggests BFB / progressive amplification).
2. **Too few segments:** fewer than the screening minimum of oscillating CN changes / adjacent oscillating segments cannot be called chromothripsis (insufficient evidence).

### From Cortés-Ciriano 2020

1. **Confidence tiers:** 4–6 oscillating segments → low-confidence; ≥7 → high-confidence. Below 4 → not called.
2. **Focal-SV minimum:** events with <6 SVs are excluded from chromothripsis calling.

---

## Test Datasets

### Dataset: Korbel & Campbell screening threshold (oscillating CN changes)

**Source:** Korbel & Campbell 2013; Magrangeas et al. 2011 (first-pass cutoff), via PMC3861665.

| Parameter | Value |
|-----------|-------|
| Minimum oscillating CN changes (first-pass screen) | 10 |
| Canonical number of CN states | 2 |
| Breakpoint clustering null | exponential inter-breakpoint distances |

### Dataset: Cortés-Ciriano 2020 segment / SV thresholds

**Source:** Cortés-Ciriano et al. 2020, Nat Genet 52:331–341 (PMC7058534).

| Parameter | Value |
|-----------|-------|
| High-confidence adjacent oscillating segments | ≥ 7 |
| Low-confidence adjacent oscillating segments | 4–6 |
| Fraction of segments oscillating between two states (canonical) | > 60% |
| Minimum clustered intrachromosomal SVs | 6 |

### Dataset: Hand-constructed two-state oscillating CN profile (worked example)

**Source:** Derived from criterion B (two-state oscillation) + the ≥10 oscillation screen (Korbel & Campbell 2013 / Magrangeas et al. 2011).

| Profile (per-segment integer CN) | Oscillation count | States | Expected |
|-----------|-------|--------|----------|
| 2,1,2,1,2,1,2,1,2,1,2 (11 segments) | 10 | {1,2} (2 states) | Chromothripsis (10 ≥ 10, 2 states) |
| 2,1,2,1,2,1 (6 segments) | 5 | {1,2} | NotComplex (5 < 10) |
| 2,3,4,5,6,7 (monotone rising) | 0 (no down-up reversals) | 6 states | NotComplex (progressive amplification, not 2-state oscillation) |

---

## Assumptions

1. **ASSUMPTION: Oscillation = adjacent-segment state reversal count.** Korbel & Campbell / Magrangeas count "oscillating copy number changes". We operationalise an *oscillation* as a segment whose CN state differs from its predecessor (a state transition), and we additionally require the profile to alternate between a bounded number of states (≤ 3, canonically 2) per criterion B. The exact bracketing used by each tool varies; the count of state-transitions is the directly source-supported quantity. Justification: the screen is explicitly a "first-pass" count of CN changes; transition counting is the minimal faithful realisation.
2. **ASSUMPTION: Clustering test summarised as fraction of short inter-breakpoint gaps vs exponential null.** Korbel & Campbell only state that random breakpoints give exponentially distributed inter-breakpoint distances; they do not fix one goodness-of-fit statistic. We expose the inter-breakpoint distances and the mean, and flag clustering when the coefficient of variation exceeds 1 (the exponential distribution has CV = 1; over-dispersion toward short gaps with a few large gaps gives CV > 1). This is a transparent, source-anchored summary, not a clinical caller.

---

## Recommendations for Test Coverage

1. **MUST Test:** Two-state oscillating profile with exactly 10 transitions → classified Chromothripsis. — Evidence: Korbel & Campbell 2013 criterion B + ≥10 first-pass screen (PMC3861665).
2. **MUST Test:** Profile with 5 transitions → NotComplex (below 10 screen). — Evidence: PMC3861665.
3. **MUST Test:** Monotone rising CN (no reversals, >3 states) → NotComplex even if many segments. — Evidence: two-state hallmark, criterion B.
4. **MUST Test:** Oscillation counting matches per-segment state-transition count on hand inputs. — Evidence: Magrangeas / Korbel & Campbell.
5. **MUST Test:** SV burden below 6 → not eligible for chromothripsis (focal exclusion). — Evidence: Cortés-Ciriano 2020.
6. **SHOULD Test:** Confidence tier ≥7 oscillating segments → high vs 4–6 → low. — Rationale: Cortés-Ciriano 2020 tiering.
7. **SHOULD Test:** Breakpoint clustering CV>1 flagged, regular spacing (CV≈0) not flagged. — Rationale: exponential null (Korbel & Campbell).
8. **COULD Test:** Empty / null / single-segment inputs handled with explicit validation. — Rationale: API robustness.

---

## References

1. Korbel JO, Campbell PJ. (2013). Criteria for Inference of Chromothripsis in Cancer Genomes. Cell 152(6):1226–1236. https://doi.org/10.1016/j.cell.2013.02.023 (PMID 23498933)
2. Cortés-Ciriano I, Lee JJ-K, Xi R, et al. (2020). Comprehensive analysis of chromothripsis in 2,658 human cancers using whole-genome sequencing. Nature Genetics 52:331–341. https://doi.org/10.1038/s41588-019-0576-7 (PMCID PMC7058534)
3. Magrangeas F, Avet-Loiseau H, Munshi NC, Minvielle S. (2011). Chromothripsis identifies a rare and aggressive entity among newly diagnosed multiple myeloma patients. Blood 118(3):675–678. https://doi.org/10.1182/blood-2011-03-344069 (first-pass ≥10 oscillating-CN-change screen, cited via PMC3861665)
4. Maher CA, Wilson RK. (2012). Chromothripsis and human disease: piecing together the shattering process / "Chromothripsis and beyond". (open-access review enumerating the Korbel & Campbell criteria). https://pmc.ncbi.nlm.nih.gov/articles/PMC3861665/

---

## Change History

- **2026-06-15**: Initial documentation (ONCO-SV-001).
