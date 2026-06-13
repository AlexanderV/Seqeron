# Test Specification: COMPGEN-REARR-001

**Test Unit ID:** COMPGEN-REARR-001
**Area:** Comparative
**Algorithm:** Genome Rearrangement Detection by Breakpoints (signed gene-order comparison)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Hunter College CompBio Lecture 16 — Genome rearrangements, sorting by reversals | 1 | https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf | 2026-06-13 |
| 2 | Tannier et al. — breakpoint distance (via PMC "On the Complexity of Rearrangement Problems under the Breakpoint Distance") | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3887456/ | 2026-06-13 |
| 3 | Bafna & Pevzner (1998) — Sorting by Transpositions, SIAM J. Discrete Math. 11(2):224–240 | 1 | https://doi.org/10.1137/S089548019528280X | 2026-06-13 |

### 1.2 Key Evidence Points

1. A breakpoint is a position adjacent in one genome but not the other — Tannier (PMC3887456); Hunter Lecture 16.
2. Breakpoint criterion on extended signed permutation: pair `(x,y)` consecutive in α is a breakpoint iff **neither `(x,y)` nor `(−y,−x)`** is an adjacency of the identity β — Hunter Lecture 16.
3. Extend the permutation with `π₀ = 0` and `π_{n+1} = n+1` before counting — Bafna & Pevzner (1998); Hunter Lecture 16.
4. Worked example: `α=(−2,−3,+1,+6,−5,−4)` has exactly **6 breakpoints**; `(−5,−4)` is NOT a breakpoint because `(4,5) ∈ β` — Hunter Lecture 16.
5. Identity / collinear order has `b(β) = 0` — Hunter Lecture 16.
6. A reversal (inversion) `α[i,j]` reverses a segment **and negates its signs**; lower bound `d ≥ b(α)/2` — Hunter Lecture 16.
7. A transposition moves a block to a new location **preserving orientation** (no sign change) — Bafna & Pevzner (1998).
8. Operation classes: inversions, transpositions, deletions, insertions, duplications — Bafna & Pevzner (1998).
9. Breakpoint distance formula `d(π₁,π₂) = n − sim(π₁,π₂)` (sim = common adjacencies) — Tannier (PMC3887456).

### 1.3 Documented Corner Cases

- Identity / collinear order → 0 breakpoints (Hunter).
- Sign-consecutive descending pair (e.g. `(−5,−4)`) → not a breakpoint (Hunter, criterion tests `(−y,−x)`).
- Extended endpoints `(0,π₁)` and `(πₙ,n+1)` may themselves be breakpoints (Hunter worked example).
- No common adjacencies → `d_BP = n` (Tannier).

### 1.4 Known Failure Modes / Pitfalls

1. Forgetting the `(−y,−x)` half of the criterion → over-counting descending runs as breakpoints — Hunter Lecture 16.
2. Forgetting to extend with `0` and `n+1` → missing boundary breakpoints — Bafna & Pevzner (1998).
3. Treating a block relocation (transposition) as an inversion when no sign flip occurred — Bafna & Pevzner (1998).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectRearrangements(genome1Genes, genome2Genes, orthologMap)` | ComparativeGenomics | **Canonical** | Returns one `RearrangementEvent` per breakpoint in the signed gene-order comparison. |
| `ClassifyRearrangement(RearrangementEvent)` | ComparativeGenomics | **Canonical** | Maps a detected breakpoint's local signature to `RearrangementType` (Inversion vs Transposition). |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Identical signed gene order ⇒ 0 breakpoint events (`b(β)=0`). | Yes | Hunter Lecture 16 |
| INV-2 | Breakpoint pair `(x,y)` excludes both `(x,y)` and `(−y,−x)` matching identity adjacency. | Yes | Hunter Lecture 16 |
| INV-3 | The permutation is extended with `0` and `n+1`; boundary pairs participate. | Yes | Bafna & Pevzner (1998) |
| INV-4 | Number of breakpoint events lies in `[0, n+1]` (n+1 internal pairs of the extended permutation of n markers). | Yes | Hunter Lecture 16 (extension) |
| INV-5 | `ClassifyRearrangement` returns Inversion iff the disrupted pair shows a local sign reversal; Transposition for an orientation-preserving relocation. | Yes | Hunter (reversal negates signs); Bafna & Pevzner (transposition preserves orientation) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Hunter worked example | α=(−2,−3,+1,+6,−5,−4) vs identity | exactly 6 breakpoint events | Hunter Lecture 16 lines 287–289 |
| M2 | Identity order | collinear, identical strands | 0 events | Hunter "b_β(β)=0" |
| M3 | Single reversed block | α=(+1,−4,−3,−2,+5) | exactly 2 breakpoint events | reversal def (Hunter) + criterion |
| M4 | Sign-consecutive descending pair | adjacency `(−5,−4)` present | that pair is NOT counted as a breakpoint | Hunter "(−5,−4) is not a breakpoint" |
| M5 | Classify inversion | event from a sign-flipped local reversal | `RearrangementType.Inversion` | Hunter reversal negates signs |
| M6 | Classify transposition | event from same-strand block relocation | `RearrangementType.Transposition` | Bafna & Pevzner transposition |
| M7 | Fewer than 2 mappable orthologs | one or zero anchors | 0 events (no internal adjacency) | breakpoints over consecutive pairs |
| M8 | Null genome1 (Detect) | first arg null | `ArgumentNullException` | contract / sibling convention |
| M9 | Null orthologMap (Detect) | map arg null | `ArgumentNullException` | contract |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Dangling ortholog | target gene absent in genome2 | anchor skipped, no crash; remaining pairs evaluated | robustness, mirrors synteny unit |
| S2 | Breakpoint-distance consistency | small case cross-checked vs `n − common adjacencies` | event count equals `d_BP` | Tannier formula |
| S3 | Empty genomes (Detect) | both lists empty | 0 events, no exception | contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property: event count bound | random-but-fixed permutations | count ∈ [0, n+1] | INV-4 value bound |
| C2 | Property: identity idempotence | identical inputs of varying size | always 0 events | INV-1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Legacy file `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomicsTests.cs` contains a `DetectRearrangements` region (3 tests) with permissive assertions (`Any(...)`, `GreaterThanOrEqualTo(0)`) against the prior ad-hoc heuristic, plus `RearrangementEvent`/`RearrangementType` record/enum smoke tests.
- No `ComparativeGenomics_DetectRearrangements_Tests.cs` (the canonical per-unit file) exists.
- No `ClassifyRearrangement` method or tests exist.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (6-breakpoint example) | ❌ Missing | no evidence-based count test exists |
| M2 (identity → 0) | ⚠ Weak | legacy `CollinearGenomes_ReturnsEmpty` checks only inversion-count==0, not total |
| M3 (single reversal → 2) | ❌ Missing | |
| M4 (sign-consecutive not a breakpoint) | ❌ Missing | |
| M5 (classify inversion) | ❌ Missing | method does not exist |
| M6 (classify transposition) | ❌ Missing | method does not exist |
| M7 (<2 orthologs → 0) | ❌ Missing | |
| M8 (null genome1) | ❌ Missing | |
| M9 (null map) | ❌ Missing | |
| S1 (dangling ortholog) | ⚠ Weak | legacy `MissingGene_DetectsDeletion` uses `GreaterThanOrEqualTo(0)` (vacuous) |
| S2 (d_BP consistency) | ❌ Missing | |
| S3 (empty genomes) | ❌ Missing | |
| C1 (count bound property) | ❌ Missing | |
| C2 (identity idempotence) | ❌ Missing | |
| legacy `InvertedRegion_DetectsInversion` | 🔁 Duplicate | superseded by M1/M3/M5 against corrected behavior |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_DetectRearrangements_Tests.cs` — all M/S/C cases for `DetectRearrangements` and `ClassifyRearrangement`.
- **Remove:** the `#region DetectRearrangements Tests` (3 tests) from `ComparativeGenomicsTests.cs` — they assert the prior ad-hoc heuristic with permissive checks and conflict with the corrected breakpoint behavior. The `RearrangementEvent`/`RearrangementType` record/enum smoke tests remain (they test the data types, not the algorithm).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ComparativeGenomics_DetectRearrangements_Tests.cs` | Canonical unit tests | 14 |
| `ComparativeGenomicsTests.cs` | Legacy data-type smoke tests only (DetectRearrangements region removed) | (unchanged minus 3) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented exact 6-breakpoint test | ✅ Done |
| 2 | M2 | ⚠ Weak | rewritten: asserts total events == 0 | ✅ Done |
| 3 | M3 | ❌ Missing | implemented exact 2-breakpoint test | ✅ Done |
| 4 | M4 | ❌ Missing | implemented sign-consecutive exclusion test | ✅ Done |
| 5 | M5 | ❌ Missing | implemented classify-inversion test | ✅ Done |
| 6 | M6 | ❌ Missing | implemented classify-transposition test | ✅ Done |
| 7 | M7 | ❌ Missing | implemented <2-ortholog test | ✅ Done |
| 8 | M8 | ❌ Missing | implemented null-genome1 test | ✅ Done |
| 9 | M9 | ❌ Missing | implemented null-map test | ✅ Done |
| 10 | S1 | ⚠ Weak | rewritten: exact remaining-event check | ✅ Done |
| 11 | S2 | ❌ Missing | implemented d_BP consistency test | ✅ Done |
| 12 | S3 | ❌ Missing | implemented empty-genomes test | ✅ Done |
| 13 | C1 | ❌ Missing | implemented count-bound property test | ✅ Done |
| 14 | C2 | ❌ Missing | implemented identity-idempotence property test | ✅ Done |
| 15 | legacy InvertedRegion / CollinearGenomes / MissingGene | 🔁 Duplicate / ⚠ Weak | removed legacy DetectRearrangements region | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 6 breakpoints |
| M2 | ✅ | 0 events |
| M3 | ✅ | exact 2 breakpoints |
| M4 | ✅ | sign-consecutive pair excluded |
| M5 | ✅ | Inversion |
| M6 | ✅ | Transposition |
| M7 | ✅ | 0 events for <2 anchors |
| M8 | ✅ | ArgumentNullException |
| M9 | ✅ | ArgumentNullException |
| S1 | ✅ | dangling ortholog skipped |
| S2 | ✅ | event count == d_BP |
| S3 | ✅ | empty → 0, no throw |
| C1 | ✅ | count ∈ [0, n+1] |
| C2 | ✅ | identity → 0 across sizes |
| legacy DetectRearrangements region | ✅ | removed |

All in-scope cases ✅. Count of ✅ (14 planned cases) equals total in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Anchors supplied as ordered ortholog mapping (anchor generation delegated). | §2 method inputs |
| 2 | `Gene.Strand` `'+' / '-'` encodes the permutation sign. | M3, M4, M5 |
| 3 | Translocation/Deletion/Insertion/Duplication not classified by these two methods (no single-permutation rule). | §4.3 scope, M5/M6 |

---

## 7. Open Questions / Decisions

1. **Decision:** `DetectRearrangements` is reframed from the prior ad-hoc heuristic (which reported Inversion/Deletion/Insertion from gaps without a source basis) to **breakpoint detection**, the formally defined quantity. The checklist abbreviation "DetectRearrangements(blocks)" is satisfied by detecting rearrangement boundaries (breakpoints) over the ordered ortholog markers. Conflict with the prior checklist text noted; behavior follows the external evidence (Tannier, Hunter, Bafna–Pevzner).
2. **Decision:** classification limited to Inversion and Transposition (the two operations definable from a signed in-order permutation); other enum values out of scope per Assumption 3.
