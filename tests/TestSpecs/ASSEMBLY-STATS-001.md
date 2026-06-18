# Test Specification: ASSEMBLY-STATS-001

**Test Unit ID:** ASSEMBLY-STATS-001
**Area:** Assembly
**Algorithm:** Assembly Statistics (N50 / L50 / Nx / Lx / auN, gap detection, contiguity summary)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Miller, Koren & Sutton (2010), *Genomics* 95(6):315-327, §1.2 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2874646/ | 2026-06-13 |
| 2 | Wikipedia, "N50, L50, and related statistics" (worked example) | 4 | https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics | 2026-06-13 |
| 3 | QUAST `quast_libs/N50.py` (ablab/quast) | 3 | https://raw.githubusercontent.com/ablab/quast/master/quast_libs/N50.py | 2026-06-13 |
| 4 | Li H (2020), "auN: a new metric to measure assembly contiguity" | 3 | https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity | 2026-06-13 |

### 1.2 Key Evidence Points

1. N50 = length of the smallest contig in the smallest set of largest contigs whose combined length is **at least 50%** of the assembly — Miller 2010 §1.2 (source 1).
2. L50 = the **count** of those contigs; N50 is a **length** — Wikipedia (source 2).
3. Cumulative threshold is inclusive (≥); QUAST stops when `total − cumulative ≤ total/2`, i.e. cumulative ≥ total/2 (source 3).
4. auN = Σᵢ Lᵢ² / Σⱼ Lⱼ — lh3 (source 4) and QUAST `au_metric` (source 3).
5. N90 ≤ N50, L90 ≥ L50 (raising the threshold cannot increase Nx) — Wikipedia (source 2).
6. Worked example A {80,70,50,40,30,20}, total 290 → N50=70, L50=2; B {…,10,5}, total 305 → N50=50, L50=3 — Wikipedia (source 2).

### 1.3 Documented Corner Cases

- Inclusive boundary: the first largest-first contig at which cumulative ≥ 50% of total defines N50 (sources 1, 3).
- Empty input: QUAST returns `None`; repository returns all-zero statistics / `Nx=Lx=0` / `auN=0` (source 3; documented as non-correctness-affecting).
- Monotonicity: N90 ≤ N50 (source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Floating-point rounding of the cutoff (`(long)(total*0.5)`) can shift the boundary on odd totals; the integer-exact test `cumulative*100 ≥ total*threshold` avoids it — derived from QUAST's exact `s <= limit` (source 3).
2. Confusing N50 (length) with L50 (count) — Wikipedia (source 2).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateStatistics(sequences)` | GenomeAssemblyAnalyzer | Canonical | Aggregates N50/L50/N90/L90/largest/smallest/totals/GC/gaps |
| `CalculateNx(sortedLengths, totalLength, threshold)` | GenomeAssemblyAnalyzer | Canonical | Core Nx/Lx; inclusive ≥ threshold% |
| `CalculateNx(lengths, threshold)` | GenomeAssemblyAnalyzer | Delegate | Sorts + totals, delegates to core overload |
| `CalculateN50(lengths)` | GenomeAssemblyAnalyzer | Delegate | `CalculateNx(lengths, 50).Nx` |
| `CalculateAuN(lengths)` | GenomeAssemblyAnalyzer | Canonical | Σl²/Σl |
| `FindGaps(sequences, minGapLength)` | GenomeAssemblyAnalyzer | Canonical | Maximal N-runs, 0-based inclusive [Start,End] |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | N50 is the length of the shortest contig in the largest-first prefix whose cumulative length is ≥ 50% of total | Yes | Source 1 §1.2; Source 3 |
| INV-2 | L50 = number of contigs in that prefix; N50 is a length, L50 a count | Yes | Source 2 |
| INV-3 | N90 ≤ N50 and L90 ≥ L50 (Nx non-increasing in x; Lx non-decreasing) | Yes | Source 2 |
| INV-4 | auN = Σl²/Σl | Yes | Source 4; Source 3 |
| INV-5 | Threshold boundary is inclusive: cumulative reaching exactly x% of total selects that contig as Nx | Yes | Source 1 "at least 50%"; Source 3 `s <= limit` |
| INV-6 | FindGaps yields one entry per maximal run of N/n with 0-based inclusive [Start,End] and Length = End−Start+1 | Yes | Implementation contract |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | N50/L50 Assembly A | CalculateNx({80,70,50,40,30,20}, total=290, 50) | Nx=70, Lx=2 | Source 2 worked ex.; Source 1 |
| M2 | N50/L50 Assembly B | CalculateNx({80,70,50,40,30,20,10,5}, 50) | Nx=50, Lx=3 | Source 2 worked ex. |
| M3 | N90/L90 Assembly A | CalculateNx(A, 90) | Nx=30, Lx=5 | Source 2 (cum 270 ≥ 261) |
| M4 | Monotonicity | N90 ≤ N50 and L90 ≥ L50 on Assembly A | 30 ≤ 70 and 5 ≥ 2 | Source 2 (INV-3) |
| M5 | Inclusive boundary | Lengths {50,50} total 100; CalculateNx(50): cumulative 50 = exactly 50% | Nx=50, Lx=1 | Source 1/3 (INV-5) |
| M6 | auN exact | CalculateAuN({100,80,60,40,20}) | 22000/300 = 73.3333… | Source 4; Source 3 |
| M7 | CalculateN50 delegate | CalculateN50({80,70,50,40,30,20}) (unsorted input shuffled) | 70 | Source 2; M1 |
| M8 | CalculateNx 2-arg delegate | CalculateNx({20,80,50,30,70,40}, 50) on shuffled input | Nx=70, Lx=2 | Source 2; M1 |
| M9 | Statistics aggregation | CalculateStatistics of A as contigs (G/C filled to set lengths) | N50=70, L50=2, N90=30, L90=5, Largest=80, Smallest=20, TotalLength=290 | Source 2 |
| M10 | FindGaps single gap | "ACGTNNNNACGT" | one gap Start=4, End=7, Length=4 | INV-6 |
| M11 | FindGaps leading gap | "NNNNACGT" | one gap Start=0, End=3, Length=4 | INV-6 |
| M12 | FindGaps trailing gap | "ACGTNNNN" | one gap Start=4, End=7, Length=4 | INV-6 |
| M13 | FindGaps minGapLength filter | "ACGTNNACGTNNNNNNACGT", minGapLength=5 | one gap Length=6, Start=10, End=15 | INV-6 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Single contig | CalculateNx({100}, 50) | Nx=100, Lx=1 | One contig is the whole assembly |
| S2 | FindGaps no gaps | "ACGTACGT" | empty | No N runs |
| S3 | FindGaps multiple gaps | "ACGTNNNACGTNNNNNNACGT" | two gaps | Separated N runs |
| S4 | Statistics GC content | all-GC contig | GcContent=1.0 | GC fraction over non-N bases |
| S5 | Statistics gap stats | "ACGTNNNNACGT" | TotalGaps=1, TotalGapLength=4, GapPercentage=100·4/12 | Gap aggregation |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty statistics | CalculateStatistics(∅) | all-zero record | Documented edge (ASSUMPTION 1) |
| C2 | Empty Nx | CalculateNx(∅, 0, 50) | Nx=0, Lx=0 | Documented edge |
| C3 | Empty auN | CalculateAuN(∅) | 0 | Documented edge |
| C4 | All-N contig | CalculateStatistics("N"×100) | TotalGapLength=100, GapPercentage=100 | Edge |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzerTests.cs` — pre-template fixture covering the whole `GenomeAssemblyAnalyzer` class. In-scope regions for this unit: "Basic Statistics", "Nx Statistics", "Gap Analysis", plus the gap/AllNs edge tests. Out-of-scope regions (scaffolds, completeness, repeats, comparison, quality, utility) belong to other (future) units and are retained.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 N50/L50 A | ⚠ Weak | Old `CalculateNx_N50` uses different lengths and no assertion messages |
| M2 N50/L50 B | ❌ Missing | Not covered |
| M3 N90/L90 A | ⚠ Weak | Old `CalculateNx_N90` uses `GreaterThanOrEqualTo`, no exact Nx |
| M4 Monotonicity | ❌ Missing | Not covered |
| M5 Inclusive boundary | ❌ Missing | Not covered |
| M6 auN exact | ⚠ Weak | Old `CalculateAuN` uses `.Within(1)` (permissive) |
| M7 CalculateN50 delegate | ❌ Missing | Method newly added |
| M8 CalculateNx 2-arg | ❌ Missing | Overload newly added |
| M9 Statistics aggregation | ⚠ Weak | Old stats tests check only Total/Largest/Smallest, not N50/L50/N90 |
| M10 FindGaps single | ⚠ Weak | Old test lacks messages, no TestSpec ID |
| M11 FindGaps leading | ⚠ Weak | Old `FindGaps_StartingWithGap` checks only Start |
| M12 FindGaps trailing | ⚠ Weak | Old `FindGaps_EndingWithGap` checks only End |
| M13 FindGaps minLength | ⚠ Weak | Old test lacks Start/End assertions |
| S1 Single contig | ❌ Missing | Not covered |
| S2 FindGaps no gaps | ⚠ Weak | Old `FindGaps_NoGaps`, no message |
| S3 FindGaps multiple | ⚠ Weak | Old test checks count only |
| S4 Statistics GC | ⚠ Weak | Old `CalculateStatistics_GcContent`, `.Within(0.01)` |
| S5 Statistics gap stats | ⚠ Weak | Old `CalculateStatistics_WithGaps` uses `GreaterThan` |
| C1 Empty statistics | ⚠ Weak | Old `CalculateStatistics_EmptyInput`, partial |
| C2 Empty Nx | ⚠ Weak | Old `CalculateNx_EmptyInput`, no message |
| C3 Empty auN | ⚠ Weak | Old `CalculateAuN_EmptyInput`, no message |
| C4 All-N contig | ⚠ Weak | Old `CalculateStatistics_AllNs`, no message |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs` — all evidence-based MUST/SHOULD/COULD cases for N50/L50/Nx/Lx/auN/Statistics/FindGaps.
- **Remove:** the in-scope Weak/Duplicate tests from `GenomeAssemblyAnalyzerTests.cs` (the "Basic Statistics Tests", "Nx Statistics Tests", "Gap Analysis Tests" gap cases, and the `CalculateStatistics_AllNs` / `FindGaps_*Gap` edge cases). Retain the out-of-scope regions (scaffolds, completeness, repeats, comparison, quality, utility) in that file unchanged.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs` | Canonical ASSEMBLY-STATS-001 | 22 |
| `GenomeAssemblyAnalyzerTests.cs` | Out-of-scope methods (other units) | retained, in-scope tests removed |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewritten with exact values + message | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewritten exact | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented (property/monotonicity) | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented (boundary) | ✅ Done |
| 6 | M6 | ⚠ Weak | Rewritten with `.Within(1e-10)` | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ⚠ Weak | Rewritten exact aggregation | ✅ Done |
| 10 | M10 | ⚠ Weak | Rewritten exact + message | ✅ Done |
| 11 | M11 | ⚠ Weak | Rewritten exact | ✅ Done |
| 12 | M12 | ⚠ Weak | Rewritten exact | ✅ Done |
| 13 | M13 | ⚠ Weak | Rewritten exact Start/End | ✅ Done |
| 14 | S1 | ❌ Missing | Implemented | ✅ Done |
| 15 | S2 | ⚠ Weak | Rewritten | ✅ Done |
| 16 | S3 | ⚠ Weak | Rewritten exact | ✅ Done |
| 17 | S4 | ⚠ Weak | Rewritten exact | ✅ Done |
| 18 | S5 | ⚠ Weak | Rewritten exact | ✅ Done |
| 19 | C1 | ⚠ Weak | Rewritten | ✅ Done |
| 20 | C2 | ⚠ Weak | Rewritten | ✅ Done |
| 21 | C3 | ⚠ Weak | Rewritten | ✅ Done |
| 22 | C4 | ⚠ Weak | Rewritten | ✅ Done |

**Total items:** 22
**✅ Done:** 22 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact Nx=70, Lx=2 |
| M2 | ✅ | Exact Nx=50, Lx=3 |
| M3 | ✅ | Exact Nx=30, Lx=5 |
| M4 | ✅ | Monotonicity asserted |
| M5 | ✅ | Inclusive boundary asserted |
| M6 | ✅ | auN 73.3333… within 1e-10 |
| M7 | ✅ | CalculateN50 on shuffled input |
| M8 | ✅ | CalculateNx 2-arg on shuffled input |
| M9 | ✅ | Full statistics aggregation |
| M10 | ✅ | FindGaps exact |
| M11 | ✅ | Leading gap exact |
| M12 | ✅ | Trailing gap exact |
| M13 | ✅ | minGapLength filter exact |
| S1 | ✅ | Single contig |
| S2 | ✅ | No gaps empty |
| S3 | ✅ | Multiple gaps exact |
| S4 | ✅ | GC content exact |
| S5 | ✅ | Gap stats exact |
| C1 | ✅ | Empty statistics |
| C2 | ✅ | Empty Nx |
| C3 | ✅ | Empty auN |
| C4 | ✅ | All-N contig |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Empty input returns all-zero statistics / Nx=Lx=0 / auN=0 instead of `None`/throw (non-correctness-affecting API choice) | C1, C2, C3 |
| 2 | `AssemblyStatistics.MedianLength` uses the upper-median convention; not part of the cited N50/L50/Nx/auN contract | (not asserted as a canonical value) |

---

## 7. Open Questions / Decisions

1. **Project placement:** `GenomeAssemblyAnalyzer` lives in `Seqeron.Genomics.Chromosome` (not `Seqeron.Genomics.Alignment` as the Area→Project table suggests). Per the prompt, code is placed in the same project/class the sibling Registry entry uses; `Seqeron.Genomics.Tests` already references the Chromosome project. Decision: keep it in Chromosome.
2. **By-area `CalculateNx(contigs, threshold)` / `CalculateN50(contigs)` signatures:** the checklist names 1-argument-style entry points; these are added as overloads/delegates over the existing 3-arg core. No conflict with evidence.
