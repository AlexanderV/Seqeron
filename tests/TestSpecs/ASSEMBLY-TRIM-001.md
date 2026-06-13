# Test Specification: ASSEMBLY-TRIM-001

**Test Unit ID:** ASSEMBLY-TRIM-001
**Area:** Assembly
**Algorithm:** Quality Trimming (BWA / cutadapt running-sum)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Cutadapt — Algorithm details (quality trimming) | 3 | https://cutadapt.readthedocs.io/en/stable/algorithms.html | 2026-06-13 |
| 2 | BWA `bwa_trim_read` (bwaseqio.c) | 3 | https://github.com/lh3/bwa/blob/master/bwaseqio.c | 2026-06-13 |
| 3 | BWA `BWA_MIN_RDLEN` (bwtaln.h) | 3 | https://github.com/lh3/bwa/blob/master/bwtaln.h | 2026-06-13 |
| 4 | Cock et al. (2010), Sanger FASTQ format | 1 | https://doi.org/10.1093/nar/gkp1137 | 2026-06-13 |

### 1.2 Key Evidence Points

1. Algorithm: subtract cutoff from every quality, compute partial sums from each index to the 3' end, cut at the index of minimal partial sum — Cutadapt algorithm docs (#1).
2. Both ends trimmed "in turn" by repeating the procedure for the other end — Cutadapt algorithm docs (#1).
3. BWA decodes quality as `qual - 33` (Phred+33) and disables trimming when `trim_qual < 1` — BWA `bwa_trim_read` (#2).
4. Worked example: qualities 42,40,26,27,8,7,11,4,2,3 at threshold 10 → trimmed to first 4 bases (min partial sum -25 at index 4) — Cutadapt docs (#1).
5. Sanger FASTQ uses ASCII offset 33, encoding Phred 0–93 as ASCII 33–126 — Cock et al. 2010 (#4).

### 1.3 Documented Corner Cases

- Threshold < 1 (or = 0) disables trimming (BWA guard; subtracting 0 keeps all partial sums non-negative).
- All-high-quality read → unchanged (minimum partial sum at the last index).
- All-low-quality read → fully removed by the two end passes.
- A high-quality base inside a low-quality tail is retained only when the running sum has not already reached a new minimum (cutadapt refinement).

### 1.4 Known Failure Modes / Pitfalls

1. Using a naive per-base hard cutoff instead of the running sum — produces different boundaries than BWA/cutadapt (Cutadapt docs #1, refinement note).
2. Wrong ASCII offset (Phred+64 vs Phred+33) — mis-decodes qualities (Cock et al. 2010 #4).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `QualityTrimReads(reads, minQuality, minLength)` | SequenceAssembler | Canonical | Running-sum trim on both ends + min-length filter |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Trimmed sequence is a contiguous substring of the input read | Yes | Cutadapt/BWA cut at a single index per end (#1, #2) |
| INV-2 | Trimmed length ≤ original length | Yes | Trimming only removes bases (#1, #2) |
| INV-3 | threshold ≤ 0 ⇒ read returned unchanged (only min-length filter may drop it) | Yes | BWA `trim_qual < 1` guard (#2) |
| INV-4 | Every output read has length ≥ minLength | Yes | Min-length filter (Evidence Assumption 2) |
| INV-5 | Quality is decoded as Phred = ASCII − 33 | Yes | BWA `qual - 33` (#2); Cock et al. (#4) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Cutadapt worked example | seq len 10, quality `KI;<)(,%#$` (Phred 42,40,26,27,8,7,11,4,2,3), threshold 10, minLength 1 | trimmed to first 4 bases | Cutadapt docs (#1) |
| M2 | All-high-quality unchanged | quality all `I` (Phred 40), threshold 20, minLength 1 | sequence unchanged | partial-sum min at end (#1) |
| M3 | All-low-quality dropped | quality all `!` (Phred 0), threshold 20, minLength 1 | read removed (length 0) | running-sum both ends (#1) |
| M4 | Threshold ≤ 0 disables trim | mixed quality, threshold 0, minLength 1 | sequence unchanged | BWA `trim_qual < 1` (#2) |
| M5 | Min-length filter drops short | M1 read but minLength 5 | read removed (4 < 5) | min-length filter (Assumption 2) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | 5'-end trimming | low-quality prefix, high-quality suffix | leading low-quality bases removed | cutadapt "repeat for other end" |
| S2 | Min-length keeps survivor | trimmed length == minLength | read kept | boundary of INV-4 |
| S3 | Multiple reads | list with one trimmable, one droppable, one clean | per-read independent result | batch semantics |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty reads list | no reads | empty result | defensive contract |
| C2 | Empty sequence read | sequence "" / quality "" | dropped (length 0 < minLength≥1) | defensive contract |
| C3 | Null reads argument | null | ArgumentNullException | defensive contract |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssemblerTests.cs` — region `QualityTrimReads Tests`: `QualityTrimReads_TrimsLowQualityEnds`, `QualityTrimReads_RemovesTooShort`, `QualityTrimReads_KeepsHighQuality`. These test the prior naive per-base hard-cutoff implementation and do not encode the BWA/cutadapt running-sum.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 worked example | ❌ Missing | not tested previously |
| M2 all-high unchanged | ⚠ Weak | `KeepsHighQuality` exists but tests naive cutoff, no exact running-sum semantics |
| M3 all-low dropped | ❌ Missing | not tested |
| M4 threshold ≤ 0 | ❌ Missing | not tested |
| M5 min-length drop | ⚠ Weak | `RemovesTooShort` exists but built on naive cutoff |
| S1 5' trimming | ⚠ Weak | `TrimsLowQualityEnds` exists but naive cutoff, not running-sum |
| S2 min-length keep boundary | ❌ Missing | not tested |
| S3 multiple reads | ❌ Missing | not tested |
| C1 empty list | ❌ Missing | not tested |
| C2 empty sequence | ❌ Missing | not tested |
| C3 null arg | ❌ Missing | not tested |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_QualityTrimReads_Tests.cs` — all evidence-based cases for this unit.
- **Remove:** the three `QualityTrimReads_*` tests in `SequenceAssemblerTests.cs` (built on the superseded naive cutoff; superseded by the canonical file).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_QualityTrimReads_Tests.cs` | Canonical, this unit | 11 |
| `SequenceAssemblerTests.cs` | Other assembler methods (QualityTrimReads region removed) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented worked-example test | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote as exact unchanged-output test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented all-low drop test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented threshold≤0 test | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote as exact min-length drop test | ✅ Done |
| 6 | S1 | ⚠ Weak | Rewrote as exact 5'-trim test | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented min-length keep boundary | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented multiple-reads test | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented empty-list test | ✅ Done |
| 10 | C2 | ❌ Missing | Implemented empty-sequence test | ✅ Done |
| 11 | C3 | ❌ Missing | Implemented null-arg test | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `QualityTrimReads_CutadaptWorkedExample_TrimsToFirstFourBases` |
| M2 | ✅ Covered | `QualityTrimReads_AllHighQuality_ReturnsUnchanged` |
| M3 | ✅ Covered | `QualityTrimReads_AllLowQuality_DropsRead` |
| M4 | ✅ Covered | `QualityTrimReads_ThresholdZero_ReturnsUnchanged` |
| M5 | ✅ Covered | `QualityTrimReads_TrimmedShorterThanMinLength_DropsRead` |
| S1 | ✅ Covered | `QualityTrimReads_LowQualityPrefix_TrimsFivePrimeEnd` |
| S2 | ✅ Covered | `QualityTrimReads_TrimmedEqualsMinLength_KeepsRead` |
| S3 | ✅ Covered | `QualityTrimReads_MultipleReads_TrimsEachIndependently` |
| C1 | ✅ Covered | `QualityTrimReads_EmptyList_ReturnsEmpty` |
| C2 | ✅ Covered | `QualityTrimReads_EmptySequence_DropsRead` |
| C3 | ✅ Covered | `QualityTrimReads_NullReads_Throws` |

Total in-scope cases: 11. ✅ = 11.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Both-end pass order (3' then 5'); numerically independent on disjoint ends | S1, implementation |
| 2 | `minLength` post-trim filter drops survivors with length < minLength | M5, S2, INV-4 |

---

## 7. Open Questions / Decisions

1. The prior implementation used a naive per-base hard cutoff (non-conforming to BWA/cutadapt). Decision: replaced with the running-sum algorithm; checklist behavioural intent ("quality trimming") is satisfied by the authoritative algorithm. No open correctness questions remain.
