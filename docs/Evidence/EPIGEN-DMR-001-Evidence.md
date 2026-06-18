# Evidence Artifact: EPIGEN-DMR-001

**Test Unit ID:** EPIGEN-DMR-001
**Algorithm:** Differentially Methylated Region (DMR) detection (tiling-window + Fisher's exact test, methylKit model)
**Date Collected:** 2026-06-13

---

## Online Sources

### methylKit paper (Akalin et al. 2012, Genome Biology) — PMC full text

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3491415/
**Retrieved by:** WebFetch of the PMC article (PMC3491415, the open-access copy of Genome Biology 2012, 13:R87), after a WebSearch for "methylKit differentially methylated region DMR definition Fisher exact test logistic regression methylation difference cutoff Akalin 2012".
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation, peer-reviewed) / 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **DMC/DMR definition (verbatim):** "a base/region would subsequently be classified as a differentially methylated cytosine (DMC) or region (DMR)" when the null hypothesis is rejected. A DMR is a genomic region of adjacent CpG sites that are differentially methylated.
2. **Per-site methylation level (verbatim):** methylation is computed as "the ratio of C/(C+T) at each base" — i.e. methylated cytosine reads over total covered reads.
3. **Default difference + q-value cutoff (verbatim):** "By default, it will extract bases/regions with a q-value <0.01 and %methylation difference >25%."
4. **Hyper vs hypo (verbatim):** "Users can specify if they want hyper-methylated bases/regions (bases/regions with higher methylation compared to control samples) or hypo-methylated bases/regions (bases/regions with lower methylation compared to control samples)." Hyper = higher methylation than control; hypo = lower.

### methylKit `tileMethylCounts` man page (al2na/methylKit, GitHub master)

**URL:** https://github.com/al2na/methylKit/blob/master/man/tileMethylCounts-methods.Rd
**Retrieved by:** WebFetch of the raw Rd man page on the al2na/methylKit GitHub master branch, after a WebSearch for "methylKit tileMethylCounts win.size 1000 step.size 1000 default".
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Default tiling parameters (verbatim):** `win.size=1000`, `step.size=1000`, `cov.bases=0`, `mc.cores=1`.
2. **Function purpose (verbatim):** "The function summarizes methylated/unmethylated base counts over tilling windows accross genome."

### methylKit `getMethylDiff` / `get.methylDiff` source (diffMeth.R)

**URL:** https://rdrr.io/github/vd4mmind/methylkit/src/R/diffMeth.R
**Retrieved by:** WebFetch of the diffMeth.R source mirror, after the same tileMethylCounts WebSearch surfaced the diffMeth.R link.
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation source)

**Key Extracted Points:**

1. **Signature (verbatim):** `get.methylDiff(.Object, difference=25, qvalue=0.01, type="all")`.
2. **Hyper selection (verbatim R):** `.Object$qvalue<qvalue & (.Object$meth.diff) > difference`.
3. **Hypo selection (verbatim R):** `.Object$qvalue<qvalue & (.Object$meth.diff) < -1*difference`.
4. **Threshold is strict greater-than:** a region is hyper only if `meth.diff > 25`, hypo only if `meth.diff < -25` (default), AND `qvalue < 0.01`. `difference` is in percentage units (25 = 25 percentage points).

### methylKit `calculateDiffMeth` man page (al2na/methylKit, GitHub master)

**URL:** https://github.com/al2na/methylKit/blob/master/man/calculateDiffMeth-methods.Rd
**Retrieved by:** WebFetch of the raw Rd man page, after a WebSearch for 'methylKit calculateDiffMeth "Fisher's exact test" ... one sample per group ... no replicates'.
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Test selection (verbatim):** "If there is one sample in each group, e.g. after applying the pooling samples, the Fisher's exact test will be applied for differential methylation." Corroborated by the user guide / search: no replicates (one sample per group) → Fisher's exact test; multiple samples per group → logistic regression.

### Fisher's exact test — Wikipedia (citing Fisher 1922/1935 primaries)

**URL:** https://en.wikipedia.org/wiki/Fisher's_exact_test
**Retrieved by:** WebFetch of the Wikipedia article, after a WebSearch for "Fisher's exact test two-sided p-value 2x2 contingency table hypergeometric distribution formula".
**Accessed:** 2026-06-13
**Authority rank:** 4 (Wikipedia citing primary literature) — the hypergeometric formula is the primary mathematical definition.

**Key Extracted Points:**

1. **Hypergeometric probability of a single 2×2 table (verbatim):** `p = (a+b choose a)(c+d choose c)/(n choose a+c) = (a+b)! (c+d)! (a+c)! (b+d)! / (a! b! c! d! n!)`, for cells a,b,c,d with row sums (a+b),(c+d), column sums (a+c),(b+d), total n.
2. **Two-sided p-value rule (verbatim paraphrase):** "p-value by summing the probabilities for all tables with probabilities less than or equal to that of the observed table." Sum over all tables (same fixed margins) whose hypergeometric probability ≤ the observed table's probability.
3. **Worked example (verbatim cells + value):** table Studying = {Men 1, Women 9} (row total 10), Not-studying = {Men 11, Women 3} (row total 14), column totals 12/12, n=24. Single-table probability **p ≈ 0.001346076**. The most-extreme same-direction table has p ≈ 0.000033652; the one-tailed sum ≈ 0.001379728.

---

## Documented Corner Cases and Failure Modes

### From methylKit (tileMethylCounts / calculateDiffMeth)

1. **Empty input:** With no covered positions there are no tiles and no DMRs — the function yields an empty result set.
2. **Sparse tile (low coverage):** `cov.bases` filters tiles with too few covered bases; default `cov.bases=0` keeps all tiles. A tile must contain at least one informative cytosine to be tested.
3. **One sample per group:** Fisher's exact test is applied (the implemented regime here). With replicates, logistic regression is used instead (out of scope).

### From Fisher's exact test (Wikipedia / primary)

1. **Zero-coverage group:** If one group has no reads in a window (row total 0), the 2×2 table is degenerate; the only possible table given the margins is the observed one, so the exact two-sided p-value is 1.0 (no evidence of difference).
2. **Identical proportions / zero marginal:** When a row or column total is 0, or both groups have identical methylated/unmethylated counts producing the only feasible table, p = 1.0.

---

## Test Datasets

### Dataset: Fisher's exact test worked example (Wikipedia, "studying by gender")

**Source:** Fisher's exact test, Wikipedia (citing Fisher 1922/1935). https://en.wikipedia.org/wiki/Fisher's_exact_test

| Parameter | Value |
|-----------|-------|
| a (Studying, Men) | 1 |
| b (Studying, Women) | 9 |
| c (Not-studying, Men) | 11 |
| d (Not-studying, Women) | 3 |
| n | 24 |
| Single-table hypergeometric probability p | ≈ 0.001346076 |

This dataset validates the hypergeometric single-table probability used inside the Fisher's exact p-value computation.

### Dataset: DMR window — hyper-methylated, derived from methylKit model

**Source:** methylKit semantics (C/(C+T) per base; meth.diff = group2% − group1%; Fisher exact 2×2 of pooled methylated/unmethylated counts).

| Parameter | Value |
|-----------|-------|
| Group1 (control) per site level / coverage | 0.0 / 20 (all sites) |
| Group2 (treatment) per site level / coverage | 1.0 / 20 (all sites) |
| Pooled 2×2 table (3 sites) | methylated: g1=0, g2=60; unmethylated: g1=60, g2=0 |
| meth.diff (group2% − group1%) | +100 (percentage points) |
| Classification | Hyper-methylated (meth.diff > 25, q < 0.01) |
| Two-sided Fisher exact p | ≈ 0 (the only equally/more extreme table is the observed; complete separation) |

---

## Assumptions

1. **ASSUMPTION: Per-window pooling of single-sample sites into one 2×2 table.** methylKit's Fisher path operates on per-base coverage counts; this unit pools the covered cytosines inside a tiling window into one 2×2 table (sum of methylated reads vs sum of unmethylated reads, each group) and applies Fisher's exact test to the window. This mirrors `tileMethylCounts` (sum counts over the window) followed by Fisher's exact on the tiled counts; it is the documented methylKit tile→test pipeline, so it is treated as evidence-backed rather than a free assumption.
2. **ASSUMPTION: numC and numT are derived from `MethylationLevel × Coverage` and `(1 − MethylationLevel) × Coverage`, rounded to the nearest integer.** The repository `MethylationSite` stores a fractional level and an integer coverage rather than raw C/T counts; reconstructing integer counts is required to build the 2×2 table. methylKit stores numCs/numTs directly; the rounding is a representation detail, not a change to the C/(C+T) definition.

---

## Recommendations for Test Coverage

1. **MUST Test:** Hyper-methylated window (group2 fully methylated, group1 fully unmethylated) is reported with positive meth.diff and "Hypermethylated" annotation — Evidence: methylKit getMethylDiff hyper rule (`meth.diff > difference`) and hyper="higher methylation than control".
2. **MUST Test:** Hypo-methylated window (group2 unmethylated, group1 methylated) is reported with negative meth.diff and "Hypomethylated" annotation — Evidence: methylKit getMethylDiff hypo rule (`meth.diff < -difference`).
3. **MUST Test:** A window whose absolute meth.diff is at/under the cutoff is NOT reported (strict `>` threshold) — Evidence: methylKit `meth.diff > difference` (strict).
4. **MUST Test:** Tiling boundary — sites farther apart than win.size fall into separate windows — Evidence: tileMethylCounts win.size=1000, step.size=1000.
5. **MUST Test:** Fisher's exact single-table hypergeometric probability matches the Wikipedia worked example (a=1,b=9,c=11,d=3 → ≈0.001346076) — Evidence: Wikipedia hypergeometric formula + worked value.
6. **MUST Test:** Empty input → no DMRs — Evidence: methylKit empty-tile behavior.
7. **SHOULD Test:** Window with fewer covered sites than `minCpGCount` is not reported — Rationale: a DMR is a region of *adjacent* CpG sites; methylKit's `cov.bases` / region definition requires multiple sites.
8. **SHOULD Test:** Degenerate group (zero coverage in one group within a window) → p-value 1.0, not reported — Rationale: Fisher exact on a degenerate margin.
9. **COULD Test:** Determinism — same inputs yield identical DMR list and ordering by start position — Rationale: deterministic tiling.

---

## References

1. Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE. (2012). methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. Genome Biology 13:R87. https://doi.org/10.1186/gb-2012-13-10-r87 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC3491415/)
2. methylKit reference manual — `tileMethylCounts`, `calculateDiffMeth`, `getMethylDiff` methods. al2na/methylKit, GitHub. https://github.com/al2na/methylKit/blob/master/man/tileMethylCounts-methods.Rd ; https://github.com/al2na/methylKit/blob/master/man/calculateDiffMeth-methods.Rd ; diffMeth.R: https://rdrr.io/github/vd4mmind/methylkit/src/R/diffMeth.R
3. Fisher's exact test (hypergeometric probability of a 2×2 table, two-sided p-value, worked example). Wikipedia, citing Fisher RA (1922, 1935). https://en.wikipedia.org/wiki/Fisher's_exact_test

---

## Change History

- **2026-06-13**: Initial documentation.
