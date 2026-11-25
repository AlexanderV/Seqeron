# Test Specification: ASSEMBLY-CONSENSUS-001

**Test Unit ID:** ASSEMBLY-CONSENSUS-001
**Area:** Assembly
**Algorithm:** Consensus Computation (column-wise majority/threshold consensus)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25 (B10 review)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Biopython `dumb_consensus` (v1.79) | 3 | https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py | 2026-06-13 |
| 2 | EMBOSS `cons` documentation | 3 | https://emboss.sourceforge.net/apps/cvs/emboss/apps/cons.html | 2026-06-13 |
| 3 | Wikipedia "Consensus sequence" (defn + IUPAC) | 4 | https://en.wikipedia.org/wiki/Consensus_sequence | 2026-06-13 |

### 1.2 Key Evidence Points

1. Consensus = most frequent residue per alignment column — Wikipedia (def'n).
2. Decision rule: emit residue only if exactly one residue holds the max count AND `max_size/num_atoms >= threshold`; else emit ambiguous — Biopython `dumb_consensus`.
3. Gap characters `-` and `.` are skipped from the per-column tally; `num_atoms` counts only non-gap residues — Biopython.
4. Ties (≥2 residues at max count) → ambiguous symbol, not an arbitrary pick — Biopython (`len(max_atoms)==1` guard).
5. Consensus length = full alignment length (longest read); shorter reads contribute nothing past their end — Biopython (`con_len = get_alignment_length`).
6. All-gap / empty column (`num_atoms==0`) → ambiguous, no division by zero (short-circuit) — Biopython.
7. Below-plurality columns have no consensus residue — EMBOSS `cons`.

### 1.3 Documented Corner Cases

Gaps skipped; residue tie → ambiguous; sub-threshold majority → ambiguous; all-gap column → ambiguous; ragged reads span longest read. (All from Biopython `dumb_consensus`; below-plurality from EMBOSS `cons`.)

### 1.4 Known Failure Modes / Pitfalls

1. Naive `MaxBy` tie-break would silently pick an arbitrary residue, violating the source rule (tie → ambiguous) — Biopython point 4.
2. Dividing by `num_atoms` without an all-gap guard risks division by zero — avoided by the `len(max_atoms)==1` short-circuit — Biopython point 6.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ComputeConsensus(IReadOnlyList<string> alignedReads, double threshold = 0.5, char ambiguous = 'N')` | SequenceAssembler | Wrapper | Case-insensitive, `require_multiple=False`; used by AssembleOLC. |
| `ComputeConsensus(IReadOnlyList<string>, double threshold, char ambiguous, bool requireMultiple, bool ignoreCase = true)` | SequenceAssembler | Canonical | Full Biopython `dumb_consensus`; `ignoreCase: false` = Biopython 1.85 exactly. |
| MCP `compute_consensus` | AlignmentTools | Wrapper | Delegates; optional threshold/ambiguous/requireMultiple/ignoreCase. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Output length equals the longest input read's length (alignment length). | Yes | Biopython `con_len = get_alignment_length` |
| INV-02 | A column whose single most-common residue has frequency (among non-gap residues) `>= threshold` emits that residue; otherwise ambiguous. | Yes | Biopython decision rule |
| INV-03 | A column with ≥2 residues tied for the maximum count emits the ambiguous symbol. | Yes | Biopython `len(max_atoms)==1` guard |
| INV-04 | Gap characters `-` and `.` never appear in the output and never count toward `num_atoms`. | Yes | Biopython tally |
| INV-05 | Empty read list → empty string. | Yes | Trivial (no columns) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Unanimous columns | All reads identical (`ACGT`×3) | `ACGT` | Biopython: freq 1.0 ≥ threshold |
| M2 | Majority above threshold | Col0 = A,A,A,T (3/4=0.75 ≥ 0.7) at threshold 0.7 | `A...` | Biopython threshold formula |
| M3 | Sub-threshold majority → ambiguous | Col = A,A,T (2/3≈0.667 < 0.7) at threshold 0.7 | `N` at that column | Biopython point 5 |
| M4 | Tie → ambiguous | Col = A,G (1/1 tie) | `N` at that column | Biopython point 4 |
| M5 | Gaps skipped | `A-GT`,`ACGT`,`ACGT` | `ACGT` | Biopython tally (gap skip) |
| M6 | Ragged reads span longest | `ACGT`,`ACG` (lengths 4,3) | length 4; col3 from the only contributing read | Biopython INV-01 |
| M7 | All-gap column → ambiguous | `A-T`,`A-T` (col1 all gaps) | `ANT` | Biopython point 6 |
| M8 | Empty read list | `[]` | `""` | Trivial INV-05 |
| M9 | Threshold reproduces Biopython 0.7 | Col = A,A,A,A,T (4/5=0.8) vs A,A,T (2/3) at threshold 0.7 | first→A, second→N | Biopython default 0.7 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null input throws | `null` read list | `ArgumentNullException` | Repository contract (sibling methods) |
| S2 | Lowercase normalized | `acgt`,`ACGT` | `ACGT` | `char.ToUpperInvariant` per existing code |

### 4.2b B10 review additions (2026-11-25)

| ID | Test Case | Expected | Evidence / Test |
|----|-----------|----------|-----------------|
| B1 | 11 Biopython 1.85 hand cases (require_multiple, case, N residue, thresholds 0/1.5/−1, ragged) | Biopython output | `ComputeConsensus_MatchesBiopython185DumbConsensus` |
| B2 | `A-`,`AC` with requireMultiple | `AN` (default `AC`) | `ComputeConsensus_RequireMultiple_SingleResidueColumnIsAmbiguous` |
| B3 | Case-sensitive mode | `a`,`A` → `N`; `aa`,`Aa`,`AA` → `Aa`; default `AA` | `ComputeConsensus_IgnoreCaseFalse_ComparesCaseSensitively` |
| B4 | `N` is a residue | `N`,`N`,`A` (amb X) → `N` | `ComputeConsensus_NResidue_IsCountedNotSkipped` |
| B5 | Null read element | `ArgumentException` | `ComputeConsensus_NullReadElement_ThrowsArgumentException` |
| B6 | NaN threshold | `ArgumentOutOfRangeException` (both overloads) | `ComputeConsensus_NaNThreshold_Throws` |
| B7 | MCP options + description | Biopython values; `ambiguous:"XY"`/NaN rejected | `ComputeConsensus_BiopythonOptions_DelegateToCanonicalConsensus` (Mcp.Alignment) |

Oracle: Biopython 1.85 3012/3012 (Evidence §B10 review).

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Custom ambiguous symbol | tie column with `ambiguous: 'X'` | `X` at that column | Biopython default symbol |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssemblerTests.cs` — `#region ComputeConsensus Tests` (4 tests: `ComputeConsensus_IdenticalReads_ReturnsRead`, `ComputeConsensus_MajorityVote`, `ComputeConsensus_IgnoresGaps`, `ComputeConsensus_EmptyReads_ReturnsEmpty`).
- No canonical `SequenceAssembler_ComputeConsensus_Tests.cs` existed.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (unanimous) | ⚠ Weak | Existing `_IdenticalReads` has no assertion message; acceptable case but rewrite for consistency. |
| M2 (majority ≥ threshold) | ⚠ Weak | Existing `_MajorityVote` checks only `consensus[0]`, no message, no threshold notion. |
| M3 (sub-threshold) | ❌ Missing | No threshold logic existed in tests. |
| M4 (tie → ambiguous) | ❌ Missing | Old code used `MaxBy` (arbitrary) — never tested. |
| M5 (gaps skipped) | ⚠ Weak | Existing `_IgnoresGaps` has no message; only tests `-` not `.`. |
| M6 (ragged) | ❌ Missing | Not tested. |
| M7 (all-gap col) | ❌ Missing | Not tested. |
| M8 (empty list) | ✅ Covered | `_EmptyReads_ReturnsEmpty` (no message). Rewrite for message. |
| M9 (threshold param) | ❌ Missing | Parameter did not exist. |
| S1 (null) | ❌ Missing | Not tested. |
| S2 (lowercase) | ❌ Missing | Not tested. |
| C1 (custom ambiguous) | ❌ Missing | Parameter did not exist. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_ComputeConsensus_Tests.cs` — all cases M1–M9, S1–S2, C1.
- **Remove:** the `#region ComputeConsensus Tests` block (4 weak tests) from `SequenceAssemblerTests.cs` — superseded by the canonical file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_ComputeConsensus_Tests.cs` | Canonical | 12 |
| `SequenceAssemblerTests.cs` (ComputeConsensus region) | Removed | 0 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote with message | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote, full-string assert + message | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote (covers `-` and `.`) | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | M8 | ✅ Covered | Rewrote with message | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | S1 | ❌ Missing | Implemented | ✅ Done |
| 11 | S2 | ❌ Missing | Implemented | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | `ComputeConsensus_UnanimousColumns_ReturnsThatResidue` |
| M2 | ✅ | `ComputeConsensus_MajorityAboveThreshold_ReturnsMajorityResidue` |
| M3 | ✅ | `ComputeConsensus_SubThresholdMajority_ReturnsAmbiguous` |
| M4 | ✅ | `ComputeConsensus_TiedResidues_ReturnsAmbiguous` |
| M5 | ✅ | `ComputeConsensus_GapCharacters_AreSkipped` |
| M6 | ✅ | `ComputeConsensus_RaggedReads_SpansLongestRead` |
| M7 | ✅ | `ComputeConsensus_AllGapColumn_ReturnsAmbiguous` |
| M8 | ✅ | `ComputeConsensus_EmptyReadList_ReturnsEmpty` |
| M9 | ✅ | `ComputeConsensus_Threshold070_ReproducesBiopython` |
| S1 | ✅ | `ComputeConsensus_NullInput_Throws` |
| S2 | ✅ | `ComputeConsensus_LowercaseReads_AreNormalized` |
| C1 | ✅ | `ComputeConsensus_CustomAmbiguousSymbol_IsEmitted` |

In-scope cases: 12. ✅ count: 12.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Ambiguous default `N` (DNA IUPAC) vs Biopython `X` — presentation only, configurable | Default `ambiguous` parameter |
| 2 | Default threshold `0.5` (majority) vs Biopython documented `0.7` — configurable; 0.7 reproducible | Default `threshold` parameter; M9 |

Both assumptions are parameter defaults (presentation/configuration), not the decision rule, which is fully source-backed. Neither default changes the rule; all source-defined values are reachable via parameters.

---

## 7. Open Questions / Decisions

1. (B10) `MotifFinder.GenerateDumbConsensus` (Analysis, B05) re-implements the same rule; cross-batch request to delegate to the 5-arg canonical overload with `ignoreCase: false` after its equal-length check.
2. The decision rule (strict `>= threshold`, tie→ambiguous, gap-skipping, all-gap guard, length=longest read) is taken verbatim from Biopython `dumb_consensus`. Only the two defaults are documented assumptions, with every source value reachable via parameters.
