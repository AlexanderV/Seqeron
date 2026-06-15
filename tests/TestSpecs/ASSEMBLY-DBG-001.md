# Test Specification: ASSEMBLY-DBG-001

**Test Unit ID:** ASSEMBLY-DBG-001
**Area:** Assembly
**Algorithm:** De Bruijn graph assembly — graph construction (`BuildDeBruijnGraph`) and Eulerian-walk reconstruction (`AssembleDeBruijn`)
**Status:** ☑ Validated (ASSEMBLY-DBG-001, 2026-06-15)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Langmead B. "De Bruijn Graph assembly" (JHU lecture notes) | 3 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_dbg.pdf | 2026-06-13 |
| 2 | Jones NC, Pevzner PA. 2004. *An Introduction to Bioinformatics Algorithms*, MIT Press (Theorems 8.1, 8.2; §8.8-8.9) | 1 | https://eclass.uoa.gr/modules/document/file.php/NURS565/BioinformaticsAlgsBook.pdf | 2026-06-13 |
| 3 | Compeau PEC, Pevzner PA, Tesler G. 2011. How to apply de Bruijn graphs to genome assembly. *Nat Biotechnol* 29:987-991 | 1 | https://doi.org/10.1038/nbt.2023 | 2026-06-13 |

### 1.2 Key Evidence Points

1. Nodes are (k-1)-mers; each input k-mer is a directed edge from its prefix (left k-1-mer) to its suffix (right k-1-mer); repeated k-mers create multiedges — Langmead DBG p.5-9.
2. Construction iterates `i in [0, len(read)-(k-1))`, `kmer = read[i:i+k]`, edge `read[i:i+k-1] → read[i+1:i+k]`; reads shorter than k contribute no edges — Langmead DBG p.16.
3. Assembly = an Eulerian walk that visits each edge exactly once; a connected graph is Eulerian (has an Eulerian path) iff it has at most two semibalanced nodes, all others balanced — Langmead DBG p.10; Jones & Pevzner Theorems 8.1, 8.2.
4. With perfect sequencing the construction is Eulerian: left-end node has one extra out-edge, right-end node one extra in-edge, all others balanced — Langmead DBG p.15.
5. Reconstruction spells the walk: `path[0] + concat(node[-1] for node in path[1:])` — Langmead DBG p.18-19.
6. Worked reconstructions (exact): `AAABBBA` k=3 → `AAABBBA`; `a_long_long_long_time` k=5 → itself; `to_every…season` k=4 → itself — Langmead DBG p.11, p.18, p.19/22.

### 1.3 Documented Corner Cases

- Repeats ≥ k-1 create shared nodes joining cycles ⇒ multiple Eulerian walks, only one is the true genome (Langmead p.21); a repeat unresolvable at small k may be resolvable at larger k (p.22).
- Coverage gaps disconnect the graph ⇒ multiple contigs (p.24-25); coverage excess / errors make the graph non-Eulerian (p.26-27).
- A read shorter than k yields no k-mers (p.16).

### 1.4 Known Failure Modes / Pitfalls

1. Non-unique Eulerian walk when a (k-1)-mer repeats — reconstruction may not equal the source genome (Langmead p.21-22).
2. Disconnected / unbalanced graph from real (errored, uneven-coverage) data has no single Eulerian walk — Langmead p.24-27; Jones & Pevzner Theorem 8.2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `BuildDeBruijnGraph(reads, k)` | SequenceAssembler | Canonical | (k-1)-mer multigraph adjacency; node/edge correctness |
| `AssembleDeBruijn(reads, parameters)` | SequenceAssembler | Canonical | Eulerian-walk reconstruction into contigs + stats |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Every emitted graph edge goes from a (k-1)-mer prefix to the (k-1)-mer suffix of some input k-mer | Yes | Langmead DBG p.6-7 |
| INV-02 | The total number of directed edges equals the total number of k-mers chopped from all reads (multigraph; repeats counted) | Yes | Langmead DBG p.7-8 |
| INV-03 | A read of length < k contributes zero edges | Yes | Langmead DBG p.16 |
| INV-04 | For a connected graph with a unique Eulerian walk, the reconstruction equals `path[0] + concat(last char of each subsequent node)` and is the source genome | Yes | Langmead DBG p.18; J&P Thm 8.2 |
| INV-05 | Every input k-mer's spelled string appears as a substring of some emitted contig | Yes | Langmead DBG p.7, p.18 |
| INV-06 | Result statistics are consistent: `TotalLength` = sum of contig lengths, `LongestContig` = max contig length | Yes | definition (mirrors `AssembleOLC`) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | BuildDeBruijnGraph AAABBBA | k=3 graph of `AAABBBA` | Nodes `{AA,AB,BB,BA}`; edges `AA→AA, AA→AB, AB→BB, BB→BB, BB→BA` (5 edges) | Langmead DBG p.5-11 |
| M2 | BuildDeBruijnGraph multiedge | k=3 graph of `AAABBBBA` | Two `BB→BB` edges (multiedge); 6 edges total | Langmead DBG p.8 |
| M3 | BuildDeBruijnGraph edge count | edges == #k-mers (INV-02) for `AAABBBA` | 5 edges == (7-3+1)=5 k-mers | Langmead DBG p.7 |
| M4 | AssembleDeBruijn AAABBBA | k=3, single read `AAABBBA` | one contig `AAABBBA` | Langmead DBG p.11 |
| M5 | AssembleDeBruijn a_long… | k=5, `a_long_long_long_time` | one contig `a_long_long_long_time` | Langmead DBG p.18 |
| M6 | AssembleDeBruijn to_every k=4 | k=4, `to_every_thing_turn_turn_turn_there_is_a_season` | one contig == input | Langmead DBG p.19, p.22 |
| M7 | AssembleDeBruijn DNA unique | k=4, `ATGGCGTGCA` | one contig `ATGGCGTGCA` | construction+spelling rule, p.18 |
| M8 | Empty reads | `AssembleDeBruijn([])` | empty `AssemblyResult` (0 contigs/reads/length) | ASSUMPTION-2 |
| M9 | Null reads | `AssembleDeBruijn(null)` | empty result, no exception | contract (mirrors OLC) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Read shorter than k | `BuildDeBruijnGraph(["AC"], 5)` | empty graph (INV-03) | Langmead p.16; ASSUMPTION-3 |
| S2 | Empty graph build | `BuildDeBruijnGraph([], 3)` | empty graph | trivial identity |
| S3 | INV-05 substring | each input k-mer of `ATGGCGTGCA` (k=4) is a substring of the reconstructed contig | all true | Langmead p.7, p.18 |
| S4 | Stats consistency | `AssembleDeBruijn` on `a_long…` (k=5): TotalLength == contig length, LongestContig == that length | INV-06 holds | definition |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Repeat → branch node | `to_every…` graph at k=3 | some node has out-degree ≥ 2 (branching ⇒ multiple Eulerian walks) | Langmead DBG p.21-22 |
| C2 | Property: edges == #kmers | random-free deterministic multi-read set; |edges| equals total chopped k-mers (INV-02) | equality holds | Langmead p.7 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file targets `AssembleDeBruijn` or `BuildDeBruijnGraph`. The sibling `SequenceAssembler_AssembleOLC_Tests.cs` covers only the OLC methods (ASSEMBLY-OLC-001). The DBG methods existed in `SequenceAssembler.cs` (`AssembleDeBruijn`, private `BuildDeBruijnGraph`/`FindContigs`/`TraceContig`) but were untested.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new |
| M2 | ❌ Missing | new |
| M3 | ❌ Missing | new |
| M4 | ❌ Missing | new |
| M5 | ❌ Missing | new |
| M6 | ❌ Missing | new |
| M7 | ❌ Missing | new |
| M8 | ❌ Missing | new |
| M9 | ❌ Missing | new |
| S1 | ❌ Missing | new |
| S2 | ❌ Missing | new |
| S3 | ❌ Missing | new |
| S4 | ❌ Missing | new |
| C1 | ❌ Missing | new |
| C2 | ❌ Missing | new |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleDeBruijn_Tests.cs` — all M/S/C cases above.
- **Remove:** none (no prior DBG tests exist).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_AssembleDeBruijn_Tests.cs` | canonical DBG fixture | 15 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented graph node/edge assertion | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented multiedge assertion | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented edge-count == #kmers | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented AAABBBA reconstruction | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented a_long… reconstruction | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented to_every… k=4 reconstruction | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented ATGGCGTGCA reconstruction | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented empty-reads result | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented null-reads result | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented read<k empty graph | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented empty-reads empty graph | ✅ Done |
| 12 | S3 | ❌ Missing | Implemented INV-05 substring property | ✅ Done |
| 13 | S4 | ❌ Missing | Implemented stats consistency | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented repeat→branch-node structural assertion | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented edges==#kmers property | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `BuildDeBruijnGraph_AaabbbaK3_ProducesExpectedNodesAndEdges` |
| M2 | ✅ Covered | `BuildDeBruijnGraph_RepeatedKmer_ProducesMultiedge` |
| M3 | ✅ Covered | `BuildDeBruijnGraph_EdgeCount_EqualsNumberOfKmers` |
| M4 | ✅ Covered | `AssembleDeBruijn_Aaabbba_ReconstructsSingleContig` |
| M5 | ✅ Covered | `AssembleDeBruijn_ALongRepeat_ReconstructsInput` |
| M6 | ✅ Covered | `AssembleDeBruijn_ToEveryK4_ReconstructsInput` |
| M7 | ✅ Covered | `AssembleDeBruijn_DnaUniqueWalk_ReconstructsInput` |
| M8 | ✅ Covered | `AssembleDeBruijn_EmptyReads_ReturnsEmptyResult` |
| M9 | ✅ Covered | `AssembleDeBruijn_NullReads_ReturnsEmptyResult` |
| S1 | ✅ Covered | `BuildDeBruijnGraph_ReadShorterThanK_ProducesEmptyGraph` |
| S2 | ✅ Covered | `BuildDeBruijnGraph_EmptyReads_ProducesEmptyGraph` |
| S3 | ✅ Covered | `AssembleDeBruijn_DnaUniqueWalk_EveryKmerIsSubstringOfContig` |
| S4 | ✅ Covered | `AssembleDeBruijn_ALongRepeat_StatisticsAreConsistent` |
| C1 | ✅ Covered | `BuildDeBruijnGraph_ToEveryK3_HasBranchNode` |
| C2 | ✅ Covered | `BuildDeBruijnGraph_MultiRead_EdgeCountEqualsKmerCount` |

### 5.7 Validation Addendum (ASSEMBLY-DBG-001, 2026-06-15)

Validation added three tests for documented Stage-A branches that the original fixture
(17 tests) left untested. Values were re-derived from the fetched Langmead PDF + an
independent Hierholzer reference, not from the implementation's output:

| Test | Branch covered | Expected value | Source |
|------|----------------|----------------|--------|
| `AssembleDeBruijn_DisconnectedGraph_OneContigPerComponent` | weakly-connected components → one contig each | `{ATGGCGTGCA, GATTACAGGTC}` (reads share no 3-mer; each unique-walk) | Langmead DBG p.24; algorithm step 3-4 |
| `AssembleDeBruijn_MinContigLength_FiltersShortContigs` | `MinContigLength` discard filter | only `GATTACAGGTC` survives at MinContigLength=11 | contract §3.1 |
| `BuildDeBruijnGraph_NullReads_ProducesEmptyGraph` | null-reads guard | empty graph, no exception | contract §3.3 |

Fixture count after validation: **20** (17 prior + 3). Full unfiltered suite: 6497 passed, 0 failed.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Walk-selection among multiple Eulerian walks is unspecified; exact reconstruction asserted only on unique-walk inputs | M4-M7, C1 (C1 asserts the branching structure, not a specific walk) |
| 2 | Empty/null reads → empty result | M8, M9 |
| 3 | Read shorter than k contributes no k-mers | S1 |

---

## 7. Open Questions / Decisions

1. **Decision:** `BuildDeBruijnGraph` is promoted from private to public (per the unit's Methods table listing it as a tested method) returning the (k-1)-mer adjacency multigraph.
2. **Decision:** the prior `AssembleDeBruijn` (greedy non-branching-stretch tracer via `FindContigs`/`TraceContig`) did not compute an Eulerian walk and could not reconstruct the published examples (e.g. `a_long_long_long_time`). It is replaced by an Eulerian-path (Hierholzer) reconstruction conforming to Langmead p.18 / J&P Theorem 8.2. The checklist text ("Find Eulerian paths / contigs") is consistent with this correction.
