# Test Specification: ASSEMBLY-COVER-001

**Test Unit ID:** ASSEMBLY-COVER-001
**Area:** Assembly
**Algorithm:** Coverage (Depth) Calculation — per-base sequencing depth over a reference
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Illumina — Sequencing Coverage for NGS Experiments | 2 | https://sapac.illumina.com/science/technology/next-generation-sequencing/plan-experiments/coverage.html | 2026-06-13 |
| 2 | Cook, D.E. — Calculate Depth and Breadth of Coverage From a bam File | 3 | https://www.danielecook.com/calculate-depth-and-breadth-of-coverage-from-a-bam-file/ | 2026-06-13 |
| 3 | Metagenomics Wiki — SAMtools: get breadth of coverage | 3 | https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage | 2026-06-13 |
| 4 | Daley et al. (2020) — Predicting Bases for Sufficient Coverage (PMC7398442) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7398442/ | 2026-06-13 |
| 5 | Lander & Waterman (1988) Genomics 3:231-239 | 1 | https://doi.org/10.1016/0888-7543(88)90007-9 | 2026-06-13 (restatement) |
| 6 | Bowtie 2 MANUAL (`--end-to-end`, `--nofw/--norc`, default reporting) | 2 | https://raw.githubusercontent.com/BenLangmead/bowtie2/master/MANUAL.markdown | 2026-11-25 |
| 7 | Bowtie 2 `unique.h` (BowtieMapq2) | 2 | https://raw.githubusercontent.com/BenLangmead/bowtie2/master/unique.h | 2026-11-25 |
| 8 | samtools-depth(1) | 2 | https://raw.githubusercontent.com/samtools/samtools/develop/doc/samtools-depth.1 | 2026-11-25 |
| 9 | pysam 0.24.1 `samtools depth -a [-Q 2]` (numeric oracle) | 2 | Evidence §2026-11-25 | 2026-11-25 |

### 1.2 Key Evidence Points

1. Per-base depth at a reference position = the number of reads covering / mapping to that position — Metagenomics Wiki [3]; Cook [2].
2. Average depth (coverage) = sum of per-base depths / genome size; also expressible as C = LN/G — Cook [2]; Illumina [1].
3. Breadth of coverage = fraction of reference bases covered by at least one read — Metagenomics Wiki [3]; Cook [2].
4. Under Lander-Waterman uniform-Poisson placement, P(base uncovered) = e^−c and breadth = 1 − e^−c — Lander & Waterman [5]; PMC7398442 [4].
5. Default aligner placement is end-to-end (no soft clipping), on both strands, one reported alignment per read; equal-best multi-locus reads get MAPQ ≤ 1 — [6][7]; samtools depth counts each primary record over its aligned span, `-Q` filters MAPQ — [8].

### 1.3 Documented Corner Cases

- Position with zero overlapping reads → depth 0 (does not count toward breadth) [2][3].
- Read extending past reference end → cannot align end-to-end (Bowtie 2 default, no soft clipping) [6]; contributes 0 (M4).
- No reads aligned → all per-base depths 0; average depth 0; breadth 0 [2][3].

### 1.4 Known Failure Modes / Pitfalls

1. Treating coverage as a single scalar rather than a per-position array hides uneven coverage — Metagenomics Wiki [3].
2. Lander-Waterman uniform-Poisson expectation is only a modeling caveat; the per-base depth array itself is exact regardless of placement uniformity — PMC7398442 [4].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateCoverage(string reference, IReadOnlyList<string> reads, int minOverlap = 20)` | SequenceAssembler | Canonical | Returns `int[]` of length `reference.Length`; element i = number of placed reads spanning position i. Both strands (Bowtie 2/BWA default; = 5-arg overload with `includeReverseComplement: true`); forward-only via `includeReverseComplement: false`. |
| `CalculateCoverage(reference, reads, minOverlap, bool includeReverseComplement, bool excludeMultiMapped = false)` | SequenceAssembler | Canonical | Both strands / multi-mapping filter. Placement delegates to `ApproximateMatcher.FindBestMatch` / `FindWithMismatches`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Returned array length = `reference.Length`. | Yes | Per-base depth defined per reference position [3]. |
| INV-2 | Every element ≥ 0. | Yes | Depth is a count of reads [2][3]. |
| INV-3 | Σ(depth) = Σ over placed reads of (overlap length with reference). | Yes | Σ per-base depths = total bases mapped [2]. |
| INV-4 | A read placed at pos p (length L) increments exactly positions [p, p+L) ⊆ [0, refLen). | Yes | Per-base depth = reads spanning the position [3][8]; end-to-end placement [6]. |
| INV-5 | A read that fails to place (best match < minOverlap) contributes 0 everywhere. | Yes | Unaligned reads add no depth [2][3]. |
| INV-6 | Breadth = #{i : depth[i] ≥ 1} / refLen ∈ [0,1]; average = Σ(depth)/refLen. | Yes | Breadth and average-depth definitions [2][3]. |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Worked dataset depth array | ref `ACGTTGCAAT`, reads `ACGTT`/`TTGCA`/`GCAAT` (unique, place at 0/3/5), minOverlap 5 | `[1,1,1,2,2,2,2,2,1,1]` | Per-base depth = reads spanning position [3]; Evidence dataset |
| M2 | Array length = reference length | Any reference, any reads | `result.Length == reference.Length` | INV-1 [3] |
| M3 | Single read overlap interval | ref len 10, one read `AAAAA` at pos 2, minOverlap 5 | depth 1 on [2,7), 0 elsewhere | INV-4 [3] |
| M4 | Read longer than reference | ref len 4 `AAAA`, one read `AAAAAA` (len 6), minOverlap 4 | all-zero array; read cannot align end-to-end | INV-4/INV-5 boundary [3][6] |
| M5 | Unmatched read contributes 0 | ref `AAAAAAAAAA`, read `GGGGG`, minOverlap 5 | all-zero array | INV-5 [2][3] |
| M6 | Empty reads list | ref `AAAAAAAAAA`, no reads | all-zero array length 10 | INV-5 / no aligned reads → 0 depth [2][3] |
| M7 | Average depth from array | M1 dataset | Σ=15, average = 15/10 = 1.5 | Average = Σ depth / genome size [2] |
| M8 | Breadth from array | M1 dataset | #covered = 10, breadth = 10/10 = 1.0 | Breadth = covered bases / refLen [3] |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Partial coverage breadth | ref len 10, one read `AAAAA` at pos 0 | depth `[1,1,1,1,1,0,0,0,0,0]`; breadth 0.5; average 0.5 | Uneven coverage; breadth < 1 [2][3] |
| S2 | All elements non-negative | M1 dataset | every element ≥ 0 | INV-2 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Case-insensitive placement | ref `AAAAAAAAAA`, read `aaaaa` minOverlap 5 | maps; depth 1 over its interval | Implementation upper-cases both sides |

### 4.4 Placement tests (2026-11-25, samtools depth oracle [8][9])

| ID | Test Case | Expected | Evidence | Test |
|----|-----------|----------|----------|------|
| P1 | Reverse-strand read `TTGACC` on `ACGGTCAATGCC`, minOverlap 6 | 3-arg default and `includeReverseComplement: true` → `[0,0,1,1,1,1,1,1,0,0,0,0]`; explicit `includeReverseComplement: false` → all 0; + `GGTCAA` lower/upper → 2× (pysam `samtools depth -a` oracle) | [6][9] | `CalculateCoverage_ReverseStrandRead_PlacedOnlyWhenBothStrandsSearched`; MCP `CalculateCoverage_StrandAndMultiMapOptions_DelegateToCore` (default both strands; `excludeMultiMapped` on both strands → all 0, forward-only → `[0,0,1,1,1,1,1,0,0,0]`) |
| P2 | RC strictly better (`TTTCC` on `GGGAAACCC`) | `[0,1,1,1,1,1,0,0,0]` | [6][9] | `CalculateCoverage_ReverseComplementStrictlyBetter_IsChosen` |
| P3 | Multi-mapped reads counted once / excluded | see Evidence table | [6][7][8][9] | `CalculateCoverage_MultiMappedRead_CountedOnceOrExcluded` |
| P4 | Reverse palindrome `ACGT` both strands same span = 1 locus | `[0,0,1,1,1,1,0,0]` | [7][9] | `CalculateCoverage_ReversePalindrome_BothStrandsSameSpan_IsUnique` |
| P5 | Mismatch tolerance = length − Hamming ≥ minOverlap | `AGCTCT`: 5 → [4,10); 6 → 0 | [9] | `CalculateCoverage_MismatchedRead_PlacedWhenMatchesReachMinOverlap` |
| P6 | minOverlap < 1 throws (0, −3) | `ArgumentOutOfRangeException` | zero matches ≠ alignment | `CalculateCoverage_MinOverlapBelowOne_Throws` |
| P7 | Null read element | `ArgumentException` | validation | `CalculateCoverage_NullReadElement_ThrowsArgumentException` |
| P8 | 2000 random cases, both flags | equals brute-force oracle (port of the samtools-checked Python oracle) | [9] | Fuzzing `CalculateCoverage_BothStrandsAndMultiMapping_MatchBruteForceOracle` |
| P9 | MCP delegation of options + validation | as P1/P3/P6 | — | Mcp `CalculateCoverage_StrandAndMultiMapOptions_DelegateToCore` |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `CalculateCoverage` and `SequenceAssembler_*Coverage*`. No existing test file for this method. Sibling assembly units have `SequenceAssembler_MergeContigs_Tests.cs` and `SequenceAssembler_Scaffold_Tests.cs`.
- Implementation present: `SequenceAssembler.CalculateCoverage` in `src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit; no prior test |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| M8 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_CalculateCoverage_Tests.cs` — all cases for this unit.
- **Remove:** none (no prior tests for this method).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_CalculateCoverage_Tests.cs` | Canonical unit tests | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented worked-dataset depth array test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented array-length invariant test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented single-read interval test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented read-longer-than-reference test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented unmatched-read test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented empty-reads test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented average-depth test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented breadth test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented partial-coverage test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented non-negativity test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented case-insensitive test | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `CalculateCoverage_WorkedDataset_ReturnsExpectedDepthArray` |
| M2 | ✅ Covered | `CalculateCoverage_AnyInput_ArrayLengthEqualsReferenceLength` |
| M3 | ✅ Covered | `CalculateCoverage_SingleRead_IncrementsExactInterval` |
| M4 | ✅ Covered | `CalculateCoverage_ReadLongerThanReference_ContributesZero` |
| M5 | ✅ Covered | `CalculateCoverage_UnmatchedRead_ContributesZero` |
| M6 | ✅ Covered | `CalculateCoverage_EmptyReads_ReturnsAllZeroArray` |
| M7 | ✅ Covered | `CalculateCoverage_WorkedDataset_AverageDepthIsOnePointFive` |
| M8 | ✅ Covered | `CalculateCoverage_WorkedDataset_BreadthIsOne` |
| S1 | ✅ Covered | `CalculateCoverage_PartialCoverage_BreadthAndAverageAreHalf` |
| S2 | ✅ Covered | `CalculateCoverage_WorkedDataset_AllDepthsNonNegative` |
| C1 | ✅ Covered | `CalculateCoverage_LowercaseRead_MapsCaseInsensitively` |

**Total in-scope cases:** 11 — **✅ count:** 11

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Read placement is ungapped end-to-end (Bowtie 2 `--end-to-end` without gaps; mismatch cost 1; `--score-min` replaced by `minOverlap` matches); tie → forward strand, leftmost (Bowtie 2: pseudo-random). Cross-checked against samtools depth on 2014 cases. | Read placement in all cases |

---

## 7. Open Questions / Decisions

1. Sources define depth given an alignment but not the aligner itself. Decision: test the depth-counting arithmetic (source-defined) using exact-match placements where the repository's best-match scan is unambiguous. Average depth and breadth are derived in tests from the returned per-base array per [2][3].
