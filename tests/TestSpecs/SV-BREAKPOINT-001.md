# Test Specification: SV-BREAKPOINT-001

**Test Unit ID:** SV-BREAKPOINT-001
**Area:** StructuralVar
**Algorithm:** Breakpoint Detection from Split (soft-clipped) Reads
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | SAM/BAM Format Specification (samtools/hts-specs), CIGAR semantics | 2 | https://samtools.github.io/hts-specs/SAMv1.pdf | 2026-06-13 |
| 2 | Tattini, D'Aurizio & Magi (2015), Front Bioeng Biotechnol 3:92 | 1 | https://doi.org/10.3389/fbioe.2015.00092 | 2026-06-13 |
| 3 | Suzuki et al. (2011) ClipCrop, BMC Bioinformatics 12(S14):S7 | 1 | https://doi.org/10.1186/1471-2105-12-S14-S7 | 2026-06-13 |
| 4 | Hart et al. (2013) SoftSearch, PLoS ONE 8(12):e83356 | 1 | https://doi.org/10.1371/journal.pone.0083356 | 2026-06-13 |

### 1.2 Key Evidence Points

1. A split read has one anchored end and one end mapping imprecisely at the breakpoint — Tattini et al. 2015.
2. The breakpoint is the "marginal point between a clipped sequence and matched sequence" (the aligned/clipped junction), giving single-base resolution — ClipCrop (Suzuki et al. 2011); Tattini et al. 2015.
3. A breakpoint is called when "at least x soft-clipped reads beginning at position y" share a clip position; support = number of clipped reads at that position — SoftSearch (Hart et al. 2013).
4. Clipped reads are combined into one breakpoint only when their clip side (left/right orientation) is the same and their positions agree within a small tolerance ("clustered within 5-base differences") — SoftSearch; ClipCrop.
5. Soft-clipped (S) bases consume the read but NOT the reference; POS is the 1-based leftmost reference position of the first reference-consuming operation — SAM spec.
6. A breakpoint position is chromosome-local (SAM POS is per-contig); reads on different chromosomes belong to different breakpoints — SAM spec; ClipCrop (sort/cluster by position within a contig).

### 1.3 Documented Corner Cases

- Clip positions with fewer than the minimum supporting clipped reads are not reported (SoftSearch default x = 5, configurable to 2).
- Reads with ≤ 5 clipped bases are too short to be a reliable breakpoint signal (SoftSearch).
- Mapping imprecision spreads a true breakpoint across nearby positions; clustering within a tolerance merges them (ClipCrop, 5 b).

### 1.4 Known Failure Modes / Pitfalls

1. Merging two distinct clip stacks separated by more than the cluster tolerance into one breakpoint — they are separate events (ClipCrop clusters only within tolerance).
2. Clustering reads across chromosomes — a breakpoint coordinate is per-chromosome (SAM POS).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindBreakpoints(splitReads, clusterTolerance, minSupport)` | StructuralVariantAnalyzer | Canonical | Clusters split-read junctions per chromosome within tolerance and reports breakpoints meeting minimum support. |
| `RefineBreakpoint(chromosome, regionStart, regionEnd, splitReads)` | StructuralVariantAnalyzer | Refinement | Returns the consensus junction coordinate (mode, then mean tie-break) of split reads whose junction falls inside the region. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Every reported breakpoint has SupportingReads ≥ minSupport. | Yes | SoftSearch ("at least x soft-clipped reads") |
| INV-02 | Two split reads cluster into one breakpoint only if same chromosome AND junction positions within clusterTolerance. | Yes | ClipCrop ("clustered within 5-base differences"); SAM POS per-contig |
| INV-03 | A reported breakpoint's position lies within [min, max] junction of its member reads (so within clusterTolerance of every member). | Yes | ClipCrop (junction = breakpoint; tolerance window) |
| INV-04 | SupportingReads of a breakpoint equals the number of split reads in its cluster. | Yes | SoftSearch (support = clipped reads at the position) |
| INV-05 | The number of reported breakpoints never exceeds the number of input split reads. | Yes | Each read joins at most one cluster (partition) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Cluster of agreeing junctions | 3 split reads, same chr, junctions {5000, 5002, 5004} within tolerance 5, minSupport 2 | One breakpoint at chr position 5002 (mean), SupportingReads = 3 | ClipCrop (junction=breakpoint, cluster within 5 b); SoftSearch (support count) |
| M2 | Below minimum support | 1 isolated split read at junction 5000, minSupport 2 | No breakpoint reported | SoftSearch ("at least x soft-clipped reads"); default min support |
| M3 | Gap exceeds tolerance | 2 reads at junctions 5000 and 5100 (gap 100 > tolerance 5), minSupport 2 | No breakpoint (two singleton groups, each below support) | ClipCrop (cluster only within tolerance) |
| M4 | Two distinct breakpoints | 2 reads at ~5000 and 2 reads at ~9000, minSupport 2 | Two breakpoints (≈5000 and ≈9000), each SupportingReads = 2 | ClipCrop (sort/cluster by position); SoftSearch (support) |
| M5 | Different chromosomes not merged | 2 reads at junction 5000 on chr1, 2 reads at 5000 on chr2, minSupport 2 | Two breakpoints, one per chromosome | SAM POS per-contig; ClipCrop per-contig clustering |
| M6 | Empty input | No split reads | Empty result | Defined trivial behavior |
| M7 | Null input | `splitReads = null` | Throws ArgumentNullException | Input-validation contract (sibling methods) |
| M8 | RefineBreakpoint consensus | Region [4990,5010] over reads with junctions {5000,5000,5004} | Returns 5000 (mode/consensus junction) | ClipCrop (breakpoint = junction); SoftSearch (reads at a position) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Tolerance boundary (inclusive) | 2 reads exactly `clusterTolerance` apart (5000 and 5005, tol 5), minSupport 2 | One breakpoint, SupportingReads = 2 | Adjacency gap = tolerance is within window |
| S2 | Tolerance boundary (exclusive) | 2 reads `clusterTolerance + 1` apart (5000 and 5006, tol 5), minSupport 2 | No breakpoint (two singletons) | Gap > tolerance splits the cluster |
| S3 | RefineBreakpoint no reads in region | Region with no member junctions | Returns null (no consensus) | Refinement undefined without support |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | RefineBreakpoint null input | `splitReads = null` | Throws ArgumentNullException | Input-validation contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `FindBreakpoints` / `RefineBreakpoint`. Existing SV tests: `StructuralVariantAnalyzerTests.cs` (legacy, broad) and `StructuralVariantAnalyzer_DetectSVs_Tests.cs` (SV-DETECT-001). Neither tests `FindBreakpoints` or `RefineBreakpoint` (these methods did not exist before this unit). No prior canonical file for SV-BREAKPOINT-001.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New method, no prior test |
| M2 | ❌ Missing | New method, no prior test |
| M3 | ❌ Missing | New method, no prior test |
| M4 | ❌ Missing | New method, no prior test |
| M5 | ❌ Missing | New method, no prior test |
| M6 | ❌ Missing | New method, no prior test |
| M7 | ❌ Missing | New method, no prior test |
| M8 | ❌ Missing | New method, no prior test |
| S1 | ❌ Missing | New method, no prior test |
| S2 | ❌ Missing | New method, no prior test |
| S3 | ❌ Missing | New method, no prior test |
| C1 | ❌ Missing | New method, no prior test |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_FindBreakpoints_Tests.cs` — all SV-BREAKPOINT-001 cases.
- **Remove:** none. Legacy `StructuralVariantAnalyzerTests.cs` does not cover these methods; left untouched (out of scope).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `StructuralVariantAnalyzer_FindBreakpoints_Tests.cs` | Canonical for SV-BREAKPOINT-001 | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented cluster-of-agreeing-junctions test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented below-min-support test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented gap-exceeds-tolerance test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented two-distinct-breakpoints test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented different-chromosomes test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented empty-input test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented null-input test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented RefineBreakpoint consensus test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented tolerance-boundary-inclusive test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented tolerance-boundary-exclusive test | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented RefineBreakpoint-no-reads test | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented RefineBreakpoint null test | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindBreakpoints_AgreeingJunctions_ReturnsOneBreakpointWithSupport` |
| M2 | ✅ Covered | `FindBreakpoints_BelowMinSupport_ReturnsEmpty` |
| M3 | ✅ Covered | `FindBreakpoints_GapExceedsTolerance_ReturnsEmpty` |
| M4 | ✅ Covered | `FindBreakpoints_TwoSeparateClusters_ReturnsTwoBreakpoints` |
| M5 | ✅ Covered | `FindBreakpoints_DifferentChromosomes_NotMerged` |
| M6 | ✅ Covered | `FindBreakpoints_EmptyInput_ReturnsEmpty` |
| M7 | ✅ Covered | `FindBreakpoints_NullInput_Throws` |
| M8 | ✅ Covered | `RefineBreakpoint_ConsensusJunction_ReturnsMode` |
| S1 | ✅ Covered | `FindBreakpoints_JunctionsExactlyToleranceApart_Cluster` |
| S2 | ✅ Covered | `FindBreakpoints_JunctionsBeyondTolerance_DoNotCluster` |
| S3 | ✅ Covered | `RefineBreakpoint_NoReadsInRegion_ReturnsNull` |
| C1 | ✅ Covered | `RefineBreakpoint_NullInput_Throws` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Cluster summary coordinate is the rounded mean of member junctions (sources fix per-read junction + tolerance window but not the summary statistic). Sub-tolerance only; does not affect membership or support. | M1 (reported position), INV-03 |

---

## 7. Open Questions / Decisions

1. Decision: minimum-support default is set to 2 (SoftSearch's documented configurable minimum) rather than its default 5, to match the sibling `ClusterSplitReads`/`DetectSVs` convention in the same class (BreakDancer -r = 2). Both are source-backed; the lower value is the established repository default and is overridable by the caller.
