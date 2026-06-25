# Validation Report: CHROM-SYNT-001 — Synteny Analysis (collinear blocks & rearrangement detection)

- **Validated:** 2026-06-24   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.FindSyntenyBlocks(orthologPairs, minGenes=3, maxGap=10)`, `ChromosomeAnalyzer.DetectRearrangements(syntenyBlocks)`
- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (synteny region lines 638–823)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Synteny_Tests.cs` (19 tests)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

This is an independent re-validation (fresh context). The source has **not changed** since the prior
archived report (`git diff cb113ce..HEAD` over the synteny region is empty; last code touch
`f0568862` only strengthened tests). Conclusions below were re-derived, not copied.

---

## Stage A — Description

### Sources opened & what they confirm

- **MCScanX (Wang et al. 2012, Nucleic Acids Res. 40(7):e49)** — re-read via PMC3326336/Oxford
  Academic. A **collinear (synteny) block** is a chain of pair-wise collinear *anchor genes*
  preserving gene order between two genomes. Two parameters govern block calling:
  **MATCH_SIZE** (min number of anchor genes to call synteny, default 5) and **MAX_GAPS**
  (max gaps allowed, default 20). This directly corroborates the spec's two parameters —
  `minGenes` (min anchors) and `maxGap` (gap tolerance). Inverted (anti-diagonal /
  reverse-collinear) blocks are standard MCScanX dot-plot output and denote inversions.
- **Wikipedia "Synteny"** — modern synteny ≈ *collinearity*: "conservation of blocks of order
  within two sets of chromosomes"; syntenic blocks preserve homologous gene order from a common
  ancestor. (The article is **silent on strand/inversion** — see Note 1.)
- **Wikipedia "Chromosomal rearrangement"** — confirms the four structural types used by
  `DetectRearrangements`: **Inversion** (reversed segment), **Translocation** (segment moved
  between non-homologous chromosomes), **Deletion** (segment removed), **Duplication** (segment
  copied).

### Block definition (validated)

A **synteny block** is a maximal run of homologous marker pairs that are **collinear** —
consecutive markers preserve consistent ordering on both genomes. **Forward ('+')** when target
positions increase with reference positions; **inverted ('-')** when target positions decrease.
A gap exceeding the threshold, or an order reversal, breaks the run; runs with `< minGenes`
anchors are dropped. `SequenceIdentity = NaN` because identity is not computable from
coordinate-only input (MCScanX derives it from BLAST alignments) — sound.

### Edge-case semantics (all defined & sourced)

empty → empty; below `minGenes` → empty; fully collinear → 1 forward block; fully inverted →
1 inverted block; gap > `maxGap` → split; multiple chromosome pairs → separate blocks; single
marker → empty.

### Independent cross-check (hand-computed)

Reference order genome1 = 1,2,3,4; target order = "genome2":

| genome2 | Expected | Reasoning |
|---------|----------|-----------|
| 1,2,3,4 | 1 forward block of 4 ('+') | target strictly increasing → collinear forward |
| 4,3,2,1 | 1 inverted block of 4 ('-') | target strictly decreasing → reverse-collinear |
| 1,2,5,3,4 | 1 forward block of 3 {1,2,5}; trailing {3,4} dropped (< minGenes=3) | the 5→3 step reverses direction → break |

These match spec M1/M2 and the gap/break semantics.

### Stage A findings / divergences (notes)

1. **Strand/inverted blocks are not in the Wikipedia "Synteny" article** (the cited label for
   M2). The forward/inverted ('+'/'-') classification is nonetheless standard and authoritative
   via **MCScanX dot-plots and SyRI inversion calls**. Definition is correct; only the *citation
   label* for M2 is weak. Documentation-only, minor.
2. The spec's collinearity rule is the simpler "consecutive consistent ordering," not MCScanX's
   full scoring DP/LIS chaining. Both are valid block definitions; the spec's is a legitimate
   documented simplification.

**Stage A verdict: PASS-WITH-NOTES** — biology/definitions correct; inverted-block concept is
properly attributable to MCScanX/SyRI rather than the (strand-silent) Wikipedia Synteny page.

---

## Stage B — Implementation

### Code path reviewed

- `FindSyntenyBlocks` (lines 643–716): groups pairs by `(Chr1, Chr2)`; sorts by `Start1`; walks
  consecutive pairs. `currentForward = curr.Start2 > prev.End2`; the first step (i==1) fixes the
  block direction `isForward`; a pair stays collinear iff `currentForward == isForward` **and**
  `gap1 ≤ maxGap·1e6` **and** `gap2 ≤ maxGap·1e6` (`gap2 = Math.Abs(...)`). On break or at the
  last index, emits a block when `geneCount ≥ minGenes`; strand = `isForward ? '+' : '-'`;
  Species2 coords use `Min(first,last)` / `Max(first,last)` to stay valid under inversion;
  `SequenceIdentity = NaN`.
- `DetectRearrangements` (lines 721–823): sorts blocks by `(Chr1, Start1)`; for adjacent
  same-Chr1 pairs — different Chr2 → **Translocation**; same Chr2 + different strand →
  **Inversion**; same Chr2 + same strand + asymmetric gap (`gap1 > 0 && gap2 ≥ 0 && gap1 > 2·gap2`)
  → **Deletion** (Size = `gap1 − max(0,gap2)`). A separate O(n²) pass flags overlapping Species1
  regions mapping to different Species2 targets → **Duplication**.

### Formula realised correctly? (traced vs code)

- **M1 forward (4 genes)** → 1 block, '+', GeneCount 4, Species1 1000–8000, Species2 1000–8000. ✔
- **M2 reverse (4 genes)** → 1 block, '-', GeneCount 4; Species2 `Min(8000,2000)=2000` /
  `Max(9000,3000)=9000`. ✔ (re-traced step by step this session)
- **M5 exactly minGenes (3)** → 1 block, GeneCount 3, '+', 1000–6000. ✔
- **M6 two chr pairs** → 2 blocks, 3 genes each. ✔
- **M16 gap split** (maxGap=2 → 2 MB; gene3→gene4 gap = 2,994,000 bp > 2 MB) → 2 blocks of 3
  (1000–6000 and 3,000,000–8,000,000). ✔
- **M9 Inversion** → Pos1=50000, Pos2=60000, Size=10000. ✔
- **M10 Translocation** → Chr2="chrB", Pos1=50000, Pos2=1000, Size=null. ✔
- **M14 Deletion** → gap1=100000, gap2=5000, 100000 > 2·5000 → Size = 95000. ✔
- **M15 Duplication** → overlap 20000–50000, different Species2 → Size = 30000. ✔
- Edge: empty→empty; single marker→empty; single block→empty; collinear genome→empty. ✔

### Cross-verification table (recomputed vs actual code via test run)

`--filter FullyQualifiedName~ChromosomeAnalyzer_Synteny_Tests` → **19 passed, 0 failed**.
Every asserted value (counts, coords, strand, sizes) reproduced.

### Variant/delegate consistency

No `*Fast`/delegate variants; two independent canonical methods. Consistent.

### Numerical robustness

`maxGap·1e6` evaluated as `int` — fits in `int.MaxValue` (~2.1e9) for `maxGap` up to ~2147;
default 10 safe. `gap2` uses `Math.Abs` to guard inverted runs. No div-by-zero.

### Test quality audit

19 tests, exact hand-calculated assertions (counts, coordinates, strand, sizes), deterministic,
covering forward/inverted/gap-split/min-genes/multi-chromosome and all four rearrangement types
plus invariants (strand ∈ {+,−}, Start ≤ End, NaN identity, valid type strings, Position1 set).
Strong; not tautological. (Spec §4 says 20 tests; the file holds 19 — a documentation off-by-one,
not a coverage gap.)

### Stage B findings / defects (notes, not blocking)

1. **Greedy pairwise chaining, not DP/LIS** — markers that stay *locally* monotonic but are
   globally interleaved are absorbed into a block rather than broken out. Faithful to the spec's
   "consecutive consistent ordering" definition; no test asserts the contrary. Not a defect vs.
   the validated spec.
2. **`DetectRearrangements` is adjacency-only with heuristic thresholds** (e.g. Deletion uses
   `gap1 > 2·gap2`). Reasonable coordinate-only heuristics sourced to the structural definitions;
   internally consistent and tested. Not the full SyRI classifier.

**Stage B verdict: PASS** — code faithfully realises the validated block/rearrangement
definitions; every worked example and cross-check value reproduced exactly. Divergences are
documented simplifications, not defects.

---

## Verdict & follow-ups

- **State: CLEAN** — no defect found; implementation matches the Stage-A-validated description.
  No code changes required.
- Tests: synteny filter → **19 passed, 0 failed**. (Full suite not re-run: no code touched.)
- Follow-ups (documentation-only, non-blocking):
  - Re-attribute the inverted-block (M2) and gap-split (M16) tests to **MCScanX/SyRI** rather
    than the strand-silent Wikipedia "Synteny" page.
  - Reconcile the spec/Evidence "20 tests" with the actual **19** in the test file.
