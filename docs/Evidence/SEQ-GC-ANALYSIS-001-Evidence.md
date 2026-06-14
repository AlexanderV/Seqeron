# Evidence Artifact: SEQ-GC-ANALYSIS-001

**Test Unit ID:** SEQ-GC-ANALYSIS-001
**Algorithm:** Comprehensive GC Analysis (GC content, GC skew, AT skew, windowed profiles, compositional variance)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia: GC-content (citing Madigan & Martinko, Brock Biology of Microorganisms)

**URL:** https://en.wikipedia.org/wiki/GC-content
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4 (Wikipedia article; the numeric formula cites a primary textbook, ref [7] Madigan & Martinko 2003)

**Key Extracted Points:**

1. **GC-content definition:** "GC-content … is the percentage of nitrogenous bases in a DNA or RNA molecule that are either guanine (G) or cytosine (C)."
2. **GC% formula (verbatim notation):** `GC% = (G + C) / (A + T + G + C) × 100%`. The denominator is the total number of bases (all four bases A, T, G, C). Cited to Madigan MT, Martinko JM (2003) *Brock Biology of Microorganisms*, 10th ed., Pearson-Prentice Hall.

### Wikipedia: GC skew (citing Lobry 1996 and Grigoriev 1998)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4 (Wikipedia article citing primaries Lobry 1996 and Grigoriev 1998)

**Key Extracted Points:**

1. **GC skew formula (verbatim):** `GC skew = (G - C)/(G + C)`.
2. **Value range:** the nucleotide composition skew spectrum ranges from **−1 to +1**; +1 corresponds to C = 0, −1 corresponds to G = 0. Positive GC skew indicates guanine richness over cytosine; negative skew indicates cytosine richness.
3. **Primary attribution:** Lobry (1996) first reported compositional asymmetry in *E. coli*, *Bacillus subtilis*, *Haemophilus influenzae*; Grigoriev (1998) for the cumulative-skew diagram.
4. **Use:** indicator of leading/lagging strand, replication origin and terminus; a switch in GC-skew sign occurs at these boundaries.

### Biopython Bio.SeqUtils (reference implementation, v1.84)

**URL:** https://biopython.org/docs/1.84/api/Bio.SeqUtils.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation in an established library)

**Key Extracted Points:**

1. **GC_skew formula (verbatim doc):** "Calculate GC skew (G-C)/(G+C) for multiple windows along the sequence." Matches Lobry (G−C)/(G+C).
2. **GC_skew zero-division handling:** "Returns 0 for windows without any G/C by handling zero division errors." → when G + C = 0 the skew is defined as 0.
3. **GC_skew ambiguous bases:** does NOT look at ambiguous nucleotides — only G and C contribute to the skew; other symbols are ignored in the numerator/denominator.
4. **gc_fraction return:** "Calculate G+C percentage in seq (float between 0 and 1)." Biopython returns a *fraction* in [0,1]; the GC-content textbook convention (above) is the same quantity expressed as a percentage (×100).

### Cuemath: Population Variance (statistics reference, worked example)

**URL:** https://www.cuemath.com/data/population-variance/
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4 (math reference with explicit formula + worked numeric example)

**Key Extracted Points:**

1. **Population variance formula (verbatim):** `σ² = Σ(xᵢ - μ)² / n` — the sum of squared deviations from the mean divided by the total number of observations N (population variance, NOT Bessel-corrected n−1).
2. **Worked example:** data {12, 13, 12, 14, 19}: mean μ = 14; squared deviations 4 + 1 + 4 + 0 + 25 = 34; population variance = 34 / 5 = **6.8**. This anchors the variance derivation independently of the implementation.

---

## Documented Corner Cases and Failure Modes

### From Biopython Bio.SeqUtils (GC_skew)

1. **Window with no G/C:** denominator G + C = 0 → skew returns 0 (handled zero-division), not NaN/exception.
2. **Ambiguous / non-ACGT symbols:** ignored by the skew; only G and C counted.

### From Wikipedia GC skew / GC-content

1. **Pure-G window:** C = 0 → GC skew = +1 (upper bound). Pure-C window: G = 0 → GC skew = −1 (lower bound).
2. **GC content of a sequence with no G/C:** numerator 0 → GC% = 0.

### From the implementation contract (windowing)

1. **Sequence shorter than the window:** no full window fits → windowed lists are empty → window-derived variances are 0 (no data); overall scalar metrics are still defined over the whole sequence.

---

## Test Datasets

### Dataset: Hand-derived worked example "GGGCCAT" (window = whole sequence)

**Source:** derivation from GC-content (Wikipedia/Brock) and GC-skew (Lobry/Wikipedia) formulas.

Sequence `GGGCCAT`: G = 3, C = 2, A = 1, T = 1, length = 7.

| Parameter | Value | Derivation |
|-----------|-------|------------|
| OverallGcContent (%) | 5/7 × 100 = 71.42857142857143 | (G+C)/(A+T+G+C)×100 = (3+2)/7×100 |
| OverallGcSkew | (3−2)/(3+2) = 0.2 | (G−C)/(G+C) |
| OverallAtSkew | (1−1)/(1+1) = 0.0 | (A−T)/(A+T) |

### Dataset: Population variance anchor {12, 13, 12, 14, 19}

**Source:** Cuemath Population Variance worked example.

| Parameter | Value |
|-----------|-------|
| Mean μ | 14 |
| Σ(xᵢ−μ)² | 34 |
| Population variance (÷N=5) | 6.8 |

### Dataset: Windowed-variance example sequence

**Source:** derivation; window = 2, step = 2 over `GGCC` gives windows `GG` (skew +1, GC% 100) and `CC` (skew −1, GC% 100).

| Window | GC% | GC skew |
|--------|-----|---------|
| GG (0..1) | 100 | +1 |
| CC (2..3) | 100 | −1 |
| GcSkewVariance | mean = 0; ((1−0)²+(−1−0)²)/2 = 1.0 |
| GcContentVariance | mean = 100; ((100−100)²+(100−100)²)/2 = 0.0 |

---

## Assumptions

1. **ASSUMPTION: GC content reported as a percentage (×100) rather than Biopython's [0,1] fraction.** The two differ only by a factor of 100 and both are documented conventions (Brock textbook uses ×100%; Biopython `gc_fraction` uses [0,1]). The repository convention (matching the existing `GcAnalysisResult` and sibling `CalculateGcContent` which already ×100) is percentage. Not correctness-affecting at the formula level — a labeling/units choice — but tests pin the exact percentage value so the convention is locked.
2. **ASSUMPTION: window-set "variability" is the population variance (÷N) of the per-window metric.** The checklist names "variability" without specifying the estimator. The implementation uses population variance Σ(x−μ)²/N (consistent across the dataset, not a sample of a larger one), matching the Cuemath population-variance definition. Sample variance (÷N−1) would change the value; population variance is the documented choice and is the natural estimator when the windows ARE the entire population of windows.

---

## Recommendations for Test Coverage

1. **MUST Test:** OverallGcContent equals (G+C)/total×100 on a known sequence — Evidence: GC-content (Brock/Wikipedia).
2. **MUST Test:** OverallGcSkew equals (G−C)/(G+C) on a known sequence — Evidence: GC skew (Lobry/Wikipedia/Biopython).
3. **MUST Test:** OverallAtSkew equals (A−T)/(A+T) — Evidence: AT skew (Charneski/Lobry, sibling SEQ-ATSKEW-001).
4. **MUST Test:** GcSkewVariance / GcContentVariance equal the population variance Σ(x−μ)²/N of the windowed values — Evidence: Cuemath population variance worked example (6.8 anchor) + windowed-variance dataset.
5. **MUST Test:** windowed lists contain the expected number of windows with correct WindowStart/WindowEnd/Position — Evidence: sliding-window definition (Biopython GC_skew "multiple windows").
6. **SHOULD Test:** sequence shorter than window → empty windowed lists, variances 0, overall scalars still computed — Rationale: documented windowing corner case.
7. **SHOULD Test:** null DnaSequence → ArgumentNullException; null/empty string → zero/empty result — Rationale: failure-mode parity with sibling methods.
8. **COULD Test:** string and DnaSequence overloads produce identical results (delegation) — Rationale: API-shape parity.

---

## References

1. Lobry JR. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution* 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626 (primary, cited via Wikipedia GC skew — formula confirmed in retrieved Wikipedia/Biopython text)
2. Grigoriev A. (1998). Analyzing genomes with cumulative skew diagrams. *Nucleic Acids Research* 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286 (primary, cited via Wikipedia GC skew)
3. Madigan MT, Martinko JM. (2003). *Brock Biology of Microorganisms*, 10th ed. Pearson-Prentice Hall. (GC% formula, cited via Wikipedia GC-content: https://en.wikipedia.org/wiki/GC-content)
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14)
5. Wikipedia. GC-content. https://en.wikipedia.org/wiki/GC-content (accessed 2026-06-14)
6. Biopython contributors. Bio.SeqUtils package (v1.84): `GC_skew`, `gc_fraction`. https://biopython.org/docs/1.84/api/Bio.SeqUtils.html (accessed 2026-06-14)
7. Population Variance — Definition, Formula, Examples. Cuemath. https://www.cuemath.com/data/population-variance/ (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
