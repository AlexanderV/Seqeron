# Validation Report: SV-BREAKPOINT-001 — Breakpoint Detection from Split (soft-clipped) Reads

- **Validated:** 2026-06-15   **Area:** StructuralVar
- **Canonical method(s):** `StructuralVariantAnalyzer.FindBreakpoints(splitReads, clusterTolerance, minSupport)`, `StructuralVariantAnalyzer.RefineBreakpoint(chromosome, regionStart, regionEnd, splitReads)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened this session (retrieved, not just cited)

| Source | URL | What it confirmed (verbatim) |
|--------|-----|------------------------------|
| ClipCrop (Suzuki et al. 2011) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3287472/ | "The marginal point between a clipped sequence and matched sequence is denoted as a breakpoint." · breakpoints are "sorted and clustered within 5-base differences." · "When the left side of the breakpoint is clipped, it is denoted as an L-breakpoint, and R-breakpoint in the opposite case." · no explicit minimum read count (uses a reliability score). |
| SoftSearch (Hart et al. 2013) | https://pmc.ncbi.nlm.nih.gov/articles/PMC3865185/ | "A putative breakpoint is defined when there is at least x soft-clipped reads beginning at position y." · default 5 ("SoftSearch_2 is default parameters" = 5; "SoftSearch_1 … decreased the minimum number of softclipped reads from 5 to 2"). · "Soft-clipped reads supporting a break point event are combined if the left/right orientation is in the same direction." · "Soft-clipped reads with more than 5 unmapped bases are passed through for further analysis." |
| Tattini, D'Aurizio & Magi (2015) | https://pmc.ncbi.nlm.nih.gov/articles/PMC4479793/ | "one end is anchored to the reference genome and the other end maps imprecisely owing to the presence of an underlying structural variant or indel breakpoint" · "SR methods allow for the detection of SVs with single base-pair resolution." · "Splitread searches for clusters of split reads…" |
| SAM/BAM Spec (samtools/hts-specs, SAMv1.tex) | https://raw.githubusercontent.com/samtools/hts-specs/master/SAMv1.tex | POS = "1-based leftmost mapping POSition of the first CIGAR operation that consumes a reference base." · consume table M(q+r), I(q), D(r), N(r), S(q only), H(neither). · "S may only have H operations between them and the ends of the CIGAR string." · "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ." |

All four sources match the repo's TestSpec §1.2 and Evidence doc **verbatim**; no citation-label drift. Authority ranks (1 peer-reviewed; 2 spec) are correct.

### Formula / definition check
- Per-read breakpoint = aligned/clipped junction (ClipCrop). Implemented as `SplitRead.SupplementaryPosition`. ✅
- Cluster membership = same chromosome AND adjacent sorted junctions within `clusterTolerance` (ClipCrop "within 5-base differences"; SAM POS per-contig). ✅
- Support = cluster size; report only if ≥ `minSupport` (SoftSearch "at least x … beginning at position y"). ✅
- Defaults: tolerance 5 (ClipCrop), minSupport 2 (SoftSearch's documented lowered minimum, matching the sibling `ClusterSplitReads`/BreakDancer -r=2 repo convention; Open-Question §7 — both source-backed). ✅

### Edge-case semantics
Empty→empty; null→ArgumentNullException; below-support→none; gap>tol→split; cross-chromosome→never merged; tolerance boundary inclusive (gap == tol clusters). All defined and sourced.

### Independent cross-check (hand-computed, sourced rules)
- M1 {5000,5002,5004}, tol 5, minSup 2 → 1 BP, support 3, mean (5000+5002+5004)/3 = **5002**. ✅
- M4 {5000,5002}/{9000,9002} → 2 BPs, means **5001 / 9001**, support 2 each. ✅
- S1 gap 5 == tol → cluster; S2 gap 6 > tol → split. ✅
- M8 {5000,5000,5004} → mode **5000**. M8b tie {5000,5004} → round((5000+5004)/2)=**5002**. ✅

### Findings / divergences
- **ASSUMPTION A1 (documented):** cluster summary coordinate = rounded mean of member junctions. Sources fix the per-read junction and the tolerance window but not the summary statistic. Sub-tolerance only; does not affect membership or support. Acceptable.
- **INV-03 is overstated** for the chosen clustering: single-linkage on adjacent gaps can chain members spanning > tolerance, and the reported mean can then sit > tolerance from an extreme member. This does not contradict ClipCrop (whose "within 5-base differences" is between adjacent sorted positions), so the *clustering* is correct; only the invariant's wording is too strong. Logged BY-DESIGN in FINDINGS_REGISTER §D.

Stage A is biologically/mathematically sound. **PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs`
- `FindBreakpoints` L577-585 (null-guard + iterator split, so ArgumentNullException is eager — verified by test M7).
- `FindBreakpointsIterator` L587-629: sort by (chromosome ordinal, junction); single-linkage on adjacent gap ≤ tolerance with chromosome equality; emit clusters with size ≥ minSupport.
- `CreateBreakpoint` L631-645: position = rounded mean of member junctions (A1).
- `RefineBreakpoint` L665-700: filters by chromosome AND inclusive [regionStart, regionEnd]; returns null if none; else mode, tie-broken by rounded mean of modal coordinates.

### Formula realised correctly?
Yes. Clustering, support gating, per-contig separation, mean coordinate, and the mode/tie consensus all match the validated description. The `Position1==Position2` symmetric `+/-` strand fields are a placeholder (junction-only model carries no second locus); not relied on by the spec.

### Cross-verification table recomputed vs code (full suite run)
| Case | Expected (sourced/hand) | Code result | Match |
|------|------------------------|-------------|-------|
| M1 | 1 BP @5002, support 3 | 1 @5002, 3 | ✅ |
| M2 | none | empty | ✅ |
| M3 | none (gap 100>5) | empty | ✅ |
| M4 | 2 BP @5001/@9001, support 2 | as expected | ✅ |
| M5 | 2 BP (chr1, chr2) | as expected | ✅ |
| M6/M7 | empty / throws | as expected | ✅ |
| S1 (gap 5) | 1 BP support 2 | ✅ | ✅ |
| S2 (gap 6) | empty | ✅ | ✅ |
| chaining {5000,5005,5010} | 1 BP @5005 support 3 | ✅ | ✅ |
| M8 mode | 5000 | ✅ | ✅ |
| M8b tie {5000,5004} | 5002 | ✅ | ✅ |
| M8c other-chr excluded | 5000 | ✅ | ✅ |
| S3 no reads | null | ✅ | ✅ |
| S3b inclusive bounds {4990,5010} | 5000 | ✅ | ✅ |
| C1 null | throws | ✅ | ✅ |

### Variant/delegate consistency
`ClusterSplitReads` (sibling, clusters on `PrimaryPosition` with distance 10) is a distinct method not under this unit; `FindBreakpoints` correctly clusters on the junction (`SupplementaryPosition`). No delegate divergence.

### Test-quality gate
- **Sourced, not code-echoes:** every expected value traces to ClipCrop/SoftSearch/SAM or a hand computation; a wrong mean/mode/tolerance would fail.
- **No green-washing:** exact `Is.EqualTo` on positions/support and `Is.Empty`/`Is.Null`/`Throws`; no weakened assertions, no widened tolerances, no skips.
- **Coverage gaps found and fixed:** the original 12 tests omitted (a) the RefineBreakpoint mode-TIE branch, (b) RefineBreakpoint cross-chromosome exclusion, (c) inclusive region-bound semantics, (d) single-linkage chaining. Added 4 tests (`AdjacentGapsWithinTolerance_ChainIntoOneBreakpoint`, `ModeTie_ReturnsRoundedMeanOfModalJunctions`, `OtherChromosomeJunctions_Excluded`, `JunctionsOnRegionBounds_AreInclusive`) → **16 tests, all green**.
- **Honest green:** FULL unfiltered suite **6489 passed, 0 failed** (1 `[Explicit]` benchmark skipped). Changed test file builds **warning-free** (the 4 NUnit2007 warnings are pre-existing in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

**Gate result: PASS** (after adding 4 coverage tests).

### Findings / defects
No code defect. One invariant-wording note (INV-03, BY-DESIGN) and the documented ASSUMPTION A1.

## Verdict & follow-ups
- **Stage A: PASS** — description matches all four authoritative sources verbatim and hand checks.
- **Stage B: PASS-WITH-NOTES** — implementation faithful; tests strengthened from 12→16 to cover the tie-break, cross-chromosome refine, inclusive bounds, and single-linkage chaining branches.
- **End-state: CLEAN** — no defect; coverage gaps fully fixed in-session; build 0 errors, suite 0 failures.
- Note carried to FINDINGS_REGISTER §D: INV-03 wording overstates the single-linkage+mean guarantee (clustering itself is correct per ClipCrop adjacent-difference semantics).
