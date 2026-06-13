# Test Specification: ASSEMBLY-CORRECT-001

**Test Unit ID:** ASSEMBLY-CORRECT-001
**Area:** Assembly
**Algorithm:** K-mer spectrum (two-sided) read error correction
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Liu, Schmidt & Maskell (2013), Musket, *Bioinformatics* 29(3):308-315 | 1 / 3 | https://doi.org/10.1093/bioinformatics/bts690 | 2026-06-13 |
| 2 | Kelley, Schatz & Salzberg (2010), Quake, *Genome Biology* 11:R116 | 1 | https://doi.org/10.1186/gb-2010-11-11-r116 | 2026-06-13 |
| 3 | Song & Florea (2018), Mining statistically-solid k-mers, PMC6311904 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6311904/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. A k-mer is *trusted/solid* when its multiplicity ≥ a coverage cut-off; *untrusted/weak* otherwise — Musket [1]; Quake [2]; [3].
2. A base is *trusted* if it is covered by any trusted k-mer; trusted bases are not modified — Musket [1]; Quake [2].
3. Two-sided correction finds a *unique* alternative base making all k-mers covering position *i* trusted; corrections are single-base substitutions (≤1 error per k-mer) — Musket [1]; Quake [2].
4. If more than one alternative makes the covering k-mers trusted, the base is left unchanged (ambiguity) — Musket [1]; [3].

### 1.3 Documented Corner Cases

- Ambiguous position (>1 valid alternative) → base unchanged [1].
- No valid correcting base → base unchanged / read uncorrected [2].
- ≤1-error-per-k-mer assumption limits multi-error k-mers [1].

### 1.4 Known Failure Modes / Pitfalls

1. A solid k-mer may still contain an error, and a weak k-mer may be error-free — the cut-off cannot perfectly separate them [2][3].
2. Multiple substitution errors inside a single k-mer may not be correctable under the two-sided rule [1].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ErrorCorrectReads(reads, kmerSize, minKmerFrequency)` | SequenceAssembler | Canonical | Two-sided k-spectrum substitution correction |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Output read count equals input read count | Yes | Substitution-only model [1][2] |
| INV-2 | Each output read has the same length as its input (no indels) | Yes | Single-base substitution edits [1][2] |
| INV-3 | A position covered by any trusted k-mer is never modified | Yes | Trusted-base rule [1][2] |
| INV-4 | A base is changed only to a base that makes all covering k-mers trusted, and only when that base is unique | Yes | Two-sided unique-alternative rule [1] |
| INV-5 | Determinism: same inputs → same output | Yes | Fixed A,C,G,T candidate order; spectrum built once |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single-substitution correction | k=3, cut-off=2; 3×`ACGTACGT` + 1×`ACGTTCGT` | Error read → `ACGTACGT`; all four outputs `ACGTACGT` | [1] two-sided unique alternative; derived spectrum |
| M2 | Trusted base unchanged | All reads identical (`ACGTACGT` ×3), k=3, cut-off=2 | Every output equals input (no modification) | [1][2] trusted-base rule |
| M3 | Ambiguous position unchanged | k=1, cut-off=2; reads `A,A,C,C,G` + `T` | `T` stays `T` (A and C both valid) | [1] ambiguity rule |
| M4 | No valid correction | k=3, cut-off=2; reads where erroneous position has no trusted alternative | erroneous base unchanged | [2] no correcting set |
| M5 | Count and length preserved | Any inputs | output count = input count; per-read length unchanged | INV-1/INV-2 [1][2] |
| M6 | Null reads → exception | `ErrorCorrectReads(null)` | `ArgumentNullException` | Documented failure mode |
| M7 | kmerSize < 1 → exception | `ErrorCorrectReads(reads, 0)` | `ArgumentOutOfRangeException` | Documented failure mode |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Correct reads pass through | No errors present | output equals input (upper-cased) | Trusted reads not altered |
| S2 | Case-insensitive | lowercase erroneous read | corrected and upper-cased | Matches analyzer convention |
| S3 | Determinism | run twice | identical outputs | INV-5 |

| P1 | Property (INV-2/INV-3) | fixed-seed random reads + injected substitutions | length preserved; trusted-in-input positions unchanged | INV-2, INV-3; O(n·r·k²) property test |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Read shorter than k | k=5, read `ACG` | returned unchanged (upper-cased) | No covering k-mers |
| C2 | Empty read list | `[]` | empty result | Boundary |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No prior test file for `ErrorCorrectReads`. Searched `tests/Seqeron/Seqeron.Genomics.Tests/` — only sibling `SequenceAssembler_*` tests exist (MergeContigs, Scaffold, CalculateCoverage, ComputeConsensus, QualityTrimReads). No `ErrorCorrect`/`CORRECT` test present.
- Pre-existing production code existed but was non-conforming (single-window middle-base heuristic, first-match accepted, no ambiguity rule, no validation).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| M6 | ❌ Missing | New unit |
| M7 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| S3 | ❌ Missing | New unit |
| P1 | ❌ Missing | New unit (property test, O(n·r·k²) DoD) |
| C1 | ❌ Missing | New unit |
| C2 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_ErrorCorrectReads_Tests.cs` — all cases for this unit.
- **Remove:** none (no prior tests existed).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| SequenceAssembler_ErrorCorrectReads_Tests.cs | Canonical | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented worked-example correction test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented trusted-base-unchanged test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented ambiguity test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented no-valid-correction test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented count/length preservation test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented null-reads exception test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented kmerSize<1 exception test | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented correct-reads passthrough test | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented case-insensitive test | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented determinism test | ✅ Done |
| 11 | P1 | ❌ Missing | Implemented fixed-seed property test (INV-2/INV-3) | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented read-shorter-than-k test | ✅ Done |
| 13 | C2 | ❌ Missing | Implemented empty-list test | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `ErrorCorrectReads_SingleSubstitution_CorrectsToTrustedSequence` |
| M2 | ✅ Covered | `ErrorCorrectReads_AllReadsTrusted_LeavesUnchanged` |
| M3 | ✅ Covered | `ErrorCorrectReads_AmbiguousAlternative_LeavesBaseUnchanged` |
| M4 | ✅ Covered | `ErrorCorrectReads_NoTrustedAlternative_LeavesBaseUnchanged` |
| M5 | ✅ Covered | `ErrorCorrectReads_Always_PreservesCountAndLength` |
| M6 | ✅ Covered | `ErrorCorrectReads_NullReads_Throws` |
| M7 | ✅ Covered | `ErrorCorrectReads_KmerSizeBelowOne_Throws` |
| S1 | ✅ Covered | `ErrorCorrectReads_ErrorFreeReads_PassThroughUnchanged` |
| S2 | ✅ Covered | `ErrorCorrectReads_LowercaseInput_CorrectsAndUpperCases` |
| S3 | ✅ Covered | `ErrorCorrectReads_RunTwice_IsDeterministic` |
| P1 | ✅ Covered | `ErrorCorrectReads_RandomReads_PreservesLengthAndTrustedPositions` |
| C1 | ✅ Covered | `ErrorCorrectReads_ReadShorterThanK_ReturnsUnchanged` |
| C2 | ✅ Covered | `ErrorCorrectReads_EmptyReadList_ReturnsEmpty` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default `kmerSize`/`minKmerFrequency` values (Musket/Quake auto-select from data); non-behavioral because every behavioral test passes them explicitly | Implementation defaults only |

---

## 7. Open Questions / Decisions

1. Quake/Musket additionally trim or discard reads that cannot be corrected and use quality-weighted (q-mer) counting; this unit implements only the substitution-correction core (counts and lengths preserved). Recorded in algorithm doc §5.3 "Not implemented".
