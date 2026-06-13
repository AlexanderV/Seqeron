# Test Specification: ASSEMBLY-OLC-001

**Test Unit ID:** ASSEMBLY-OLC-001
**Area:** Assembly
**Algorithm:** Overlap-Layout-Consensus (overlap detection + OLC assembly)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Compeau, Pevzner & Tesler (2011), "How to apply de Bruijn graphs to genome assembly", Nat Biotechnol 29:987–991 | 1 | https://doi.org/10.1038/nbt.2023 | 2026-06-13 |
| 2 | Langmead B., "Overlap Layout Consensus assembly" (JHU lecture notes) | 3 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf | 2026-06-13 |
| 3 | Langmead B., "Assembly & Shortest Common Superstring" (JHU lecture notes) | 3 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/16_assembly_scs_v2.pdf | 2026-06-13 |

### 1.2 Key Evidence Points

1. Overlap graph: one node per read; directed edge A→B when a suffix of A equals a prefix of B with overlap length above a minimum threshold — Source 1 (overlap graph); Source 2 p.5; Source 3 p.16, p.23.
2. Report only the longest suffix-prefix match for an ordered pair — Source 2 p.10.
3. Edge weight = overlap length; for `GTACGTACGAT` 6-mers with minOverlap 4 the graph has 12 directed edges with weights 4 and 5 — Source 3 p.24–25 (re-derived in §1 Evidence dataset).
4. Layout = find a Hamiltonian path; exact layout is NP-complete → heuristic layout in practice — Source 1.
5. Layout uses transitive-edge reduction then emits contigs along non-branching stretches; consensus = per-column majority vote — Source 2 p.21–25, p.28.
6. Greedy maximal-overlap merging is a heuristic and not necessarily optimal — Source 3 p.45, p.57.
7. All-pairs overlap detection is O(N²) (O(d²n²)); suffix-tree overlap is O(N + a) — Source 2 p.10, p.16.

### 1.3 Documented Corner Cases

- Reads with no above-threshold overlap → no edges → each read is its own singleton contig (Source 1 overlap-graph definition; Source 2 layout).
- Unresolvable repeats split the layout into multiple contigs (Source 2 p.25).
- Report only the longest overlap per ordered pair (Source 2 p.10).
- Empty read set: no source specifies behavior → treated as the trivial empty→empty identity case (see §6 ASSUMPTION-2).

### 1.4 Known Failure Modes / Pitfalls

1. Greedy layout can produce a longer-than-minimal superstring (suboptimal) — Source 3 p.57.
2. Repeats longer than the read/k-mer length collapse and foil reconstruction — Source 3 p.58–62.
3. Sequencing errors create spurious dead-end branches that must be pruned — Source 2 p.26.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AssembleOLC(reads, parameters)` | SequenceAssembler | Canonical | Full OLC: overlap → greedy layout → chain consensus → stats |
| `FindAllOverlaps(reads, minOverlap, minIdentity)` | SequenceAssembler | Canonical | Builds the directed overlap graph (edges with overlap length) |
| `FindOverlap(s1, s2, minOverlap, minIdentity)` | SequenceAssembler | Internal | Longest suffix-prefix overlap for one ordered pair; tested via FindAllOverlaps + 1 direct case |
| `FindAllOverlaps(reads, minOverlap, minIdentity, CancellationToken, IProgress)` | SequenceAssembler | Delegate | Cancellable variant; smoke test proving same results + cancellation |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `FindAllOverlaps` never reports a self-overlap: `ReadIndex1 != ReadIndex2` for every edge | Yes | Source 1 (nodes are distinct reads); Source 2 p.5 |
| INV-2 | Every reported `OverlapLength` ≥ `minOverlap` and ≤ min(len(read1), len(read2)) | Yes | Source 2 p.5 (`l`), p.10 (longest match) |
| INV-3 | The reported overlap for an ordered pair is the **longest** suffix-prefix match (no longer valid overlap exists) | Yes | Source 2 p.10 |
| INV-4 | For an unambiguous tiling, the assembled contig is a superstring containing every input read as a substring; its length ≤ Σ read lengths and ≥ the longest single read | Yes | Source 1 (path merges reads in order); Source 3 p.26 (superstring) |
| INV-5 | Reads with no above-threshold overlap each appear as their own contig (count of contigs = count of reads when the overlap graph is edgeless) | Yes | Source 1 (no edge below threshold); Source 2 layout |
| INV-6 | `AssemblyResult.TotalLength` = Σ over emitted contigs of contig length; `LongestContig` = max contig length | Yes | Definitional (assembly statistics) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FindAllOverlaps_6mers_ExactGraph | `FindAllOverlaps` on the 6 distinct 6-mers of `GTACGTACGAT`, minOverlap 4, identity 1.0 | Exactly the 12 directed edges (A→B,len) from the Evidence table; 12 edges total | Source 3 p.24–25 (re-derived) |
| M2 | FindAllOverlaps_NoSelfOverlap | Overlap graph contains no self edges | Every edge has `ReadIndex1 != ReadIndex2` (INV-1) | Source 1; Source 2 p.5 |
| M3 | FindOverlap_LongestSuffixPrefix | `FindOverlap("CTCTAGGCC","TAGGCCCTC",3,1.0)` | length 6, pos1 = 3, pos2 = 0 (longest suffix-prefix `TAGGCC`) | Source 2 p.5, p.10 |
| M4 | AssembleOLC_UnambiguousChain_SingleContig | `AssembleOLC(["AAAAACCCCC","CCCCCGGGGG","GGGGGTTTTT"], minOverlap 5, identity 1.0, minLen 10)` | One contig `AAAAACCCCCGGGGGTTTTT` (length 20); TotalReads 3 | Source 2 p.5 + consensus chain merge; INV-4 |
| M5 | AssembleOLC_NoOverlaps_Singletons | `AssembleOLC(["AAAAAAAAAA","CCCCCCCCCC","GGGGGGGGGG"], minOverlap 5, identity 1.0, minLen 5)` | 3 contigs, each equal to an input read; TotalLength 30 | Source 1 (edgeless graph); INV-5 |
| M6 | AssembleOLC_EmptyReads_EmptyResult | `AssembleOLC([])` | Empty contigs; TotalLength 0; TotalReads 0 | §6 ASSUMPTION-2 (trivial identity) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindOverlap_IdentityThreshold | overlap with 1 mismatch in 8 (7/8 = 0.875) | accepted at minIdentity 0.85, rejected at 0.95 | Source 2 p.11–15 (approximate overlap) |
| S2 | FindOverlap_MinOverlapBoundary | overlap exactly = minOverlap accepted; one below rejected | accept at L = minOverlap, null below | Source 2 p.5 (`l`) |
| S3 | FindAllOverlaps_BelowThreshold_NoEdges | reads sharing only a 3-base overlap, minOverlap 4 | empty overlap list | INV-2 |
| S4 | AssembleOLC_ContigLengthBounds | contig length ≤ Σ read lengths and ≥ longest read, for the unambiguous chain | bounds hold (property) | INV-4 |
| S5 | FindAllOverlaps_Cancellable_SameResult | cancellable overload with `CancellationToken.None` yields same edge set as the basic overload | edge sets equal | Delegate smoke |
| S6 | FindAllOverlaps_Cancellable_Throws | cancellable overload with an already-cancelled token | throws `OperationCanceledException` | Delegate smoke |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | AssembleOLC_RepeatLimitation | reads with an internal repeat are not collapsed into a contig shorter than the longest read | result respects INV-4 lower bound | Source 3 p.58–62 (documents heuristic limitation, not exactness) |
| C2 | FindAllOverlaps_CaseInsensitive | lowercase reads produce the same overlaps as uppercase | identical edge set | `CalculateIdentity` is case-insensitive |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssemblerTests.cs` — a broad fixture covering many `SequenceAssembler` methods (FindOverlap, FindAllOverlaps, CalculateIdentity, AssembleOLC, AssembleDeBruijn, MergeContigs, CalculateStats, Scaffold, CalculateCoverage, ComputeConsensus, QualityTrimReads, ErrorCorrectReads). The OLC-relevant tests there are permissive.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 FindAllOverlaps exact 12-edge graph | ❌ Missing | existing `FindAllOverlaps_ReturnsAllValidOverlaps` only asserts `Count > 0` (Weak) |
| M2 FindAllOverlaps no self-overlap | ⚠ Weak | existing `FindAllOverlaps_NoSelfOverlaps` checks the predicate but on a 2-read input with no evidence anchor |
| M3 FindOverlap longest suffix-prefix exact | ⚠ Weak | existing `FindOverlap_PerfectOverlap_ReturnsCorrect` checks length but not pos1/pos2; not the evidence string |
| M4 AssembleOLC unambiguous chain → exact single contig | ❌ Missing | existing `AssembleOLC_PerfectOverlappingReads_MergesIntoOne` asserts only `TotalLength >= 15` (Weak) |
| M5 AssembleOLC no-overlap singletons | ⚠ Weak | existing `AssembleOLC_NoOverlaps_ReturnsSingletonContigs` checks count 3 but not contig contents/total length |
| M6 AssembleOLC empty → empty | ✅ Covered | existing `AssembleOLC_EmptyReads_ReturnsEmptyResult` adequate (will be re-expressed in canonical file) |
| S1 identity threshold | ⚠ Weak | existing `FindOverlap_WithMismatches_RespectsIdentityThreshold` ok but uses non-evidence strings; will be rewritten in canonical file |
| S2 minOverlap boundary | ❌ Missing | none |
| S3 below-threshold → no edges | ❌ Missing | none |
| S4 contig length-bound property | ❌ Missing | none |
| S5 cancellable same result | ❌ Missing | none |
| S6 cancellable throws | ❌ Missing | none |
| C1 repeat limitation | ❌ Missing | none |
| C2 case-insensitive overlaps | ❌ Missing | none |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleOLC_Tests.cs` — all evidence-based OLC tests (M1–M6, S1–S6, C1–C2) for `AssembleOLC`, `FindAllOverlaps`, and the internal `FindOverlap` used by them.
- **Remove:** the OLC-scoped tests in `SequenceAssemblerTests.cs` that this unit supersedes (the `FindOverlap`, `FindAllOverlaps`, and `AssembleOLC` regions), to avoid duplicate/weak coverage. Non-OLC tests in that file (CalculateIdentity, AssembleDeBruijn, MergeContigs, CalculateStats, Scaffold, CalculateCoverage, ComputeConsensus, QualityTrimReads, ErrorCorrectReads) belong to other units and are left untouched.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_AssembleOLC_Tests.cs` | Canonical ASSEMBLY-OLC-001 fixture | 14 |
| `SequenceAssemblerTests.cs` | Other-unit tests (OLC regions removed) | reduced (non-OLC only) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact 12-edge assertion | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote on evidence 6-mer input | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote with evidence string + pos1/pos2 | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented exact single-contig assertion | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote with contig contents + total length | ✅ Done |
| 6 | M6 | ✅ Covered | Re-expressed in canonical file | ✅ Done |
| 7 | S1 | ⚠ Weak | Rewrote in canonical file | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented boundary test | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented below-threshold test | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented property test (INV-4) | ✅ Done |
| 11 | S5 | ❌ Missing | Implemented cancellable-equality smoke | ✅ Done |
| 12 | S6 | ❌ Missing | Implemented cancellation-throws smoke | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented repeat-limitation property | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented case-insensitive test | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact 12-edge overlap graph asserted |
| M2 | ✅ Covered | no self-overlap on evidence input |
| M3 | ✅ Covered | longest suffix-prefix length + positions |
| M4 | ✅ Covered | exact single contig `AAAAACCCCCGGGGGTTTTT` |
| M5 | ✅ Covered | 3 singleton contigs + total length 30 |
| M6 | ✅ Covered | empty → empty |
| S1 | ✅ Covered | identity threshold 0.85 vs 0.95 |
| S2 | ✅ Covered | minOverlap boundary accept/reject |
| S3 | ✅ Covered | below-threshold → no edges |
| S4 | ✅ Covered | contig-length bounds property |
| S5 | ✅ Covered | cancellable overload equality |
| S6 | ✅ Covered | cancellable throws on cancelled token |
| C1 | ✅ Covered | repeat limitation lower-bound property |
| C2 | ✅ Covered | case-insensitive overlaps |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Canonical numeric cases use exact-match reads (identity 1.0); approximate-overlap path is exercised separately via the identity threshold | M1, M3, M4, S1 |
| 2 | Empty read set → empty `AssemblyResult` (no source specifies; trivial identity) | M6 |

---

## 7. Open Questions / Decisions

1. The greedy layout (`BuildContigsFromOverlaps`) is a heuristic; exact OLC layout is NP-complete (Source 1). Tests assert the **unambiguous-chain** contract (M4, INV-4) and the documented repeat limitation (C1) rather than exact reconstruction of repeat-containing genomes — consistent with the sources. No correctness-affecting constants exist in the OLC path (no scoring tables/penalties), so no values require source-backing beyond the threshold parameters, which are caller-supplied.
