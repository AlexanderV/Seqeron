# Evidence Artifact: ASSEMBLY-STATS-001

**Test Unit ID:** ASSEMBLY-STATS-001
**Algorithm:** Assembly Statistics (N50 / L50 / Nx / Lx / auN, gap detection, contiguity summary)
**Date Collected:** 2026-06-13

---

## Online Sources

### Miller, Koren & Sutton (2010) — "Assembly algorithms for next-generation sequencing data"

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2874646/
**Accessed:** 2026-06-13 (fetched via WebFetch of the PMC full-text page; redirect from ncbi.nlm.nih.gov → pmc.ncbi.nlm.nih.gov followed)
**Authority rank:** 1 (peer-reviewed review, *Genomics* 95(6):315-327)

**Key Extracted Points:**

1. **N50 definition (verbatim, §1.2 "What is an Assembly?"):** "The contig N50 is the length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly." → the cumulative threshold is **at least 50%** (inclusive), contigs taken largest-first.

### Wikipedia — "N50, L50, and related statistics" (cites Miller 2010 and Earl 2011 primaries)

**URL:** https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics (also fetched the `action=raw` wikitext at https://en.wikipedia.org/w/index.php?title=N50,_L50,_and_related_statistics&action=raw)
**Accessed:** 2026-06-13 (fetched via WebFetch of both the rendered article and the raw wikitext)
**Authority rank:** 4 (Wikipedia; used for its worked example and as a cross-check of the Miller primary)

**Key Extracted Points:**

1. **N50 (verbatim):** "The N50 is defined as the sequence length of the shortest contig at 50% of the total assembly length." Restated as a weighted median: "50% of the entire assembly is contained in contigs or scaffolds equal to or larger than this value."
2. **L50 (verbatim):** "The L50 is defined as count of smallest number of contigs whose length sum makes up half of genome size." L50 is a **count**, N50 is a **length**.
3. **N90 (verbatim):** "it is the length for which the collection of all contigs of that length or longer contains at least 90% of the sum of the lengths of all contigs." N90 ≤ N50.
4. **Worked example — Assembly A (verbatim):** "Assembly A contains six contigs of lengths 80 kbp, 70 kbp, 50 kbp, 40 kbp, 30 kbp, and 20 kbp. The sum size of assembly A is 290 kbp, the N50 contig length is 70 kbp because 80 + 70 is greater than 50% of 290, and the L50 contig count is 2 contigs." → **N50 = 70, L50 = 2** for lengths {80,70,50,40,30,20}, total 290.
5. **Worked example — Assembly B:** adding contigs of 10 and 5 kbp (sum = 305 kbp) gives **N50 = 50 kbp** (80 + 70 + 50 = 200 > 152.5) and **L50 = 3**.

### QUAST reference implementation — `quast_libs/N50.py` (ablab/quast, master)

**URL:** https://raw.githubusercontent.com/ablab/quast/master/quast_libs/N50.py
**Accessed:** 2026-06-13 (downloaded the raw source file with curl; full function bodies read)
**Authority rank:** 3 (reference implementation in an established assembly-evaluation tool)

**Key Extracted Points:**

1. **Cumulative comparison (verbatim core of `NG50_and_LG50`):**
   ```python
   s = reference_length
   limit = reference_length * (100.0 - percentage) / 100.0
   lg50 = 0
   for l in numlist:
       s -= l
       lg50 += 1
       if s <= limit:
           ng50 = l
           return ng50, lg50
   ```
   For N50 (`percentage=50`, `reference_length = sum(numlist)`), `limit = total/2` and the stop test `s <= limit` is `total − cumulative ≤ total/2`, i.e. **cumulative ≥ total/2** — the inclusive "at least 50%" convention. `N50()` = `NG50(numlist, sum(numlist))`; `L50()` = `LG50(...)` returns the iteration count.
2. **auN / au_metric (verbatim):** `return float(sum([n ** 2 for n in numlist])) / denum` where `denum = sum(numlist)` (or a supplied reference length). → **auN = Σᵢ Lᵢ² / Σⱼ Lⱼ**. Empty list returns `None`.

### Heng Li (2020) — "auN: a new metric to measure assembly contiguity"

**URL:** https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
**Accessed:** 2026-06-13 (fetched via WebFetch)
**Authority rank:** 3 (originating definition of the auN metric by the samtools/minimap2 author)

**Key Extracted Points:**

1. **Nx (verbatim):** "Contigs no shorter than Nx covers x% of the assembly," x in 0..100.
2. **auN formula (verbatim):** "auN = ∑_i L_i · (L_i / ∑_j L_j) = ∑_i L_i² / ∑_j L_j." auN is the area under the Nx curve; it is more stable than the discrete N50.

---

## Documented Corner Cases and Failure Modes

### From QUAST `N50.py`

1. **Empty input:** `au_metric` asserts `len(numlist) > 0` and returns `None`; `NG50_and_LG50` returns `(None, None)` when the loop never reaches the limit / list is empty. → Empty assembly has no defined N50/L50/auN; the repository returns 0 for these instead of throwing (documented as a deliberate, non-correctness-affecting choice in the algorithm doc).

### From Miller et al. (2010) / Wikipedia

1. **Inclusive threshold boundary:** N50 is the *first* contig (largest-first) at which cumulative length reaches **at least** 50% of the total. The boundary is inclusive (≥), confirmed independently by Miller's "at least 50%" and QUAST's `s <= limit`.
2. **N90 ≤ N50:** raising the threshold cannot increase Nx; L90 ≥ L50.

---

## Test Datasets

### Dataset: Wikipedia Assembly A (canonical N50/L50 worked example)

**Source:** Wikipedia "N50, L50, and related statistics", worked example (cites Miller 2010); cross-checked against QUAST `N50()`/`L50()`.

| Parameter | Value |
|-----------|-------|
| Contig lengths | {80, 70, 50, 40, 30, 20} |
| Total length | 290 |
| N50 | 70 |
| L50 | 2 |
| N90 | 30 (cumulative 80+70+50+40+30 = 270 ≥ 261 = 90% of 290) |
| L90 | 5 |
| auN | (80²+70²+50²+40²+30²+20²)/290 = (6400+4900+2500+1600+900+400)/290 = 16700/290 ≈ 57.5862 |

### Dataset: Wikipedia Assembly B (threshold shift example)

**Source:** Wikipedia worked example (same article).

| Parameter | Value |
|-----------|-------|
| Contig lengths | {80, 70, 50, 40, 30, 20, 10, 5} |
| Total length | 305 |
| N50 | 50 |
| L50 | 3 |

### Dataset: QUAST `au_metric` formula check

**Source:** QUAST `N50.py` `au_metric`; lh3 (2020).

| Parameter | Value |
|-----------|-------|
| Lengths | {100, 80, 60, 40, 20} |
| Total | 300 |
| auN = Σl²/Σl | (10000+6400+3600+1600+400)/300 = 22000/300 = 73.333… |
| N50 (100+80=180 ≥ 150) | 80 |
| L50 | 2 |

---

## Assumptions

1. **ASSUMPTION: Empty-input return value (non-correctness-affecting).** Authoritative sources (QUAST) return `None` for an empty contig set; the repository returns the all-zero `AssemblyStatistics` and `Nx=Lx=0` / `auN=0`. This is an API-shape choice (no valid N50 exists for an empty assembly); it does not change any defined N50/L50/auN value and is documented in the algorithm doc §6.1.
2. **ASSUMPTION: Median convention in `CalculateStatistics.MedianLength`.** The N50 sources do not define an assembly "median contig length"; the repository reports the upper median (`lengths[count/2]` over the descending-sorted list). This field is auxiliary and outside the cited N50/L50/Nx/auN contract; the canonical statistics (N50, L50, N90, L90, largest, smallest, totals, GC, gaps) are source-backed. Tested as implemented behavior, flagged as not source-derived.

---

## Recommendations for Test Coverage

1. **MUST Test:** N50=70, L50=2 on Assembly A {80,70,50,40,30,20} — Evidence: Wikipedia worked example + Miller 2010 §1.2 + QUAST.
2. **MUST Test:** N50=50, L50=3 on Assembly B {80,70,50,40,30,20,10,5} — Evidence: Wikipedia worked example.
3. **MUST Test:** inclusive boundary (cumulative reaching exactly 50% triggers Nx) — Evidence: Miller "at least 50%" / QUAST `s <= limit`.
4. **MUST Test:** auN = Σl²/Σl exact (73.333… on {100,80,60,40,20}) — Evidence: lh3 + QUAST `au_metric`.
5. **MUST Test:** N90 ≤ N50 and L90 ≥ L50 (monotonicity) — Evidence: Wikipedia N90 definition.
6. **MUST Test:** FindGaps detects N-runs with exact 0-based inclusive [Start,End] and length — Evidence: implementation contract (gap = maximal run of N/n).
7. **SHOULD Test:** CalculateStatistics aggregates N50/L50/N90/L90/largest/smallest/total/GC/gaps on a multi-contig assembly — Rationale: canonical summary method.
8. **SHOULD Test:** FindGaps minGapLength filter, leading/trailing gaps — Rationale: documented boundary handling.
9. **COULD Test:** empty input returns zeros without throwing — Rationale: documented edge case.

---

## References

1. Miller JR, Koren S, Sutton G (2010). Assembly algorithms for next-generation sequencing data. *Genomics* 95(6):315-327. https://pmc.ncbi.nlm.nih.gov/articles/PMC2874646/ (DOI: 10.1016/j.ygeno.2010.03.001)
2. Wikipedia contributors. N50, L50, and related statistics. https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics (accessed 2026-06-13)
3. Gurevich A, Saveliev V, Vyahhi N, Tesler G (2013). QUAST: quality assessment tool for genome assemblies — `quast_libs/N50.py`. https://raw.githubusercontent.com/ablab/quast/master/quast_libs/N50.py (accessed 2026-06-13; DOI of tool: 10.1093/bioinformatics/btt086)
4. Li H (2020). auN: a new metric to measure assembly contiguity. https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation (ASSEMBLY-STATS-001).
