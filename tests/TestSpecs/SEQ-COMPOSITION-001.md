# Test Specification: SEQ-COMPOSITION-001

**Test Unit ID:** SEQ-COMPOSITION-001
**Area:** Statistics
**Algorithm:** Sequence Composition (nucleotide + amino-acid composition)
**Status:** ☑ Complete (consolidated into SEQ-STATS-001 — see §7)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Biopython `Bio.SeqUtils` (`gc_fraction`, `GC_skew`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-14 |
| 2 | Wikipedia "GC skew" (cites Lobry 1996) | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-14 |
| 3 | Lobry, J. R. (1996) *Mol Biol Evol* 13:660–665 | 1 | https://doi.org/10.1093/oxfordjournals.molbev.a025626 | 2026-06-14 |
| 4 | IUPAC nucleotide & amino-acid single-letter codes | 2 | https://www.bioinformatics.org/sms/iupac.html | 2026-06-14 |

### 1.2 Key Evidence Points

1. GC content = (G+C)/(A+T+G+C+U); float in [0,1] — Biopython `gc_fraction`.
2. GC skew = (G−C)/(G+C); AT skew = (A−T)/(A+T) — Wikipedia "GC skew" / Lobry 1996.
3. Zero denominators (no G/C, or no A/T) → skew 0; empty sequence → composition 0 — Biopython.
4. Canonical alphabet {A,C,G,T,U}; N = any base; 20 amino-acid single-letter codes — IUPAC.

### 1.3 Documented Corner Cases

- Empty sequence → all-zero composition (Biopython).
- Zero-denominator skew → 0.0 (Biopython).
- Case-insensitive counting (Biopython counts lowercase).

### 1.4 Known Failure Modes / Pitfalls

1. Degenerate IUPAC codes are counted differently by Biopython (S→GC, W→denominator) than by the repository implementation, which routes them to `CountN`/`CountOther`. Agreement is exact over {A,T,G,C,U} — Biopython `gc_fraction` `ambiguous` parameter.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateNucleotideComposition(sequence)` | SequenceStatistics | **Canonical** | Identical to the canonical method of SEQ-STATS-001. |
| `CalculateAminoAcidComposition(sequence)` | SequenceStatistics | **Canonical** | Identical to the protein-composition method of SEQ-STATS-001. |

> **Both methods are already implemented and tested under SEQ-STATS-001.** SEQ-COMPOSITION-001 is a duplicate Registry entry for the same two methods of the same class (see §7).

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ GcContent ≤ 1; 0 ≤ AtContent ≤ 1 | Yes | Biopython `gc_fraction` returns a fraction in [0,1] |
| INV-2 | −1 ≤ GcSkew ≤ 1; −1 ≤ AtSkew ≤ 1 | Yes | Wikipedia "GC skew" formula bounds |
| INV-3 | CountA+T+G+C+U+N+Other = Length (counts partition the sequence) | Yes | definition of composition |
| INV-4 | Amino-acid composition `Length` = sum of `Counts` values | Yes | IUPAC residue counting |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact nucleotide counts | `AAUUGGCC` per-base counts + length | A2 T0 U2 G2 C2 N0 Other0 Len8 | IUPAC alphabet [4] |
| M2 | GC content | `ATGC` | 0.5 | Biopython `gc_fraction` [1] |
| M3 | GC skew | `GGGC` | 0.5; `GCCC` → −0.5 | Wikipedia/Lobry [2][3] |
| M4 | AT skew | `AAAT` | 0.5 | Wikipedia [2] |
| M5 | Empty/null | `""` / `null` | all-zero composition | Biopython empty handling [1] |
| M6 | Amino-acid counts | `MKVLWA` | each residue = 1, Length 6 | IUPAC amino-acid codes [4] |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitivity | upper == lower == mixed | equal GC content/counts | Biopython lowercase counting |
| S2 | Zero-denominator skews | no G/C → GcSkew 0; no A/T → AtSkew 0 | 0.0 | Biopython zero-division |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Counts partition | INV-3 holds for arbitrary sequence | sum == Length | invariant check |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Canonical tests for both methods already exist at
  `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs`,
  delivered and committed under **SEQ-STATS-001**.
- That file covers M1 (counts), M2 (GC content), M3 (GC skew both signs), M4 (AT skew), M5 (empty + null), M6 (amino-acid residue counts), S1 (case-insensitivity), S2 (zero-denominator skews) and C1 (INV-3 partition).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 nucleotide counts | ✅ Covered | `..._RnaSequence_ReturnsExactCounts` |
| M2 GC content | ✅ Covered | `..._BalancedDna_ReturnsGcContentHalf` |
| M3 GC skew (±) | ✅ Covered | `..._GRich_...`, `..._CRich_...` |
| M4 AT skew | ✅ Covered | `..._ARich_ReturnsPositiveAtSkew` |
| M5 empty/null | ✅ Covered | `..._EmptyString_...`, `..._Null_...` |
| M6 amino-acid counts | ✅ Covered | `CalculateAminoAcidComposition_Protein_ReturnsExactResidueCounts` |
| S1 case-insensitivity | ✅ Covered | `..._MixedCase_MatchesUpperCase` |
| S2 zero-denominator skews | ✅ Covered | `..._NoGc_...`, `..._NoAt_...` |
| C1 INV-3 partition | ✅ Covered | `..._ArbitrarySequence_CountsSumToLength` |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs` — already the single canonical fixture for these two methods (owned by SEQ-STATS-001).
- **Remove:** nothing. **Do NOT create** a second test file for the same methods — that would violate the duplicate-elimination rule (one canonical file per method). SEQ-COMPOSITION-001 reuses the existing fixture.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateNucleotideComposition_Tests.cs` | Canonical fixture for both composition methods (shared with SEQ-STATS-001) | 17 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | (all M/S/C) | ✅ Covered | None — already covered by the canonical fixture; no new/duplicate tests created | ✅ Done |

**Total items:** 1
**✅ Done:** 1 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Covered by canonical fixture |
| M2 | ✅ | Covered by canonical fixture |
| M3 | ✅ | Covered by canonical fixture |
| M4 | ✅ | Covered by canonical fixture |
| M5 | ✅ | Covered by canonical fixture |
| M6 | ✅ | Covered by canonical fixture |
| S1 | ✅ | Covered by canonical fixture |
| S2 | ✅ | Covered by canonical fixture |
| C1 | ✅ | Covered by canonical fixture |

All 9 in-scope cases ✅ (9 of 9). No ❌/⚠ remain.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Degenerate IUPAC codes are not counted toward GC/AT totals (routed to N/Other). Exact agreement with Biopython over {A,T,G,C,U}. | §1.4, algorithm doc 5.3 |

---

## 7. Open Questions / Decisions

1. **DECISION — duplicate Registry entry.** The Processing Registry contains two entries for the identical pair of methods on the same class: **SEQ-STATS-001** ("Sequence Composition Statistics", ☑ Complete) and **SEQ-COMPOSITION-001** ("Sequence Composition"). Both list `SequenceStatistics.CalculateNucleotideComposition(...)` as canonical and add `CalculateAminoAcidComposition(...)`. SEQ-STATS-001 already ships the implementation, the canonical test fixture, Evidence and an algorithm doc. Per the prompt's duplicate-elimination rule ("one canonical test file per unit; no duplicate tests remain") and workflow-control rule ("note the conflict in the TestSpec and update the checklist entry"), SEQ-COMPOSITION-001 is **resolved by consolidation**: no new production code and no duplicate test file are created; this unit reuses the existing implementation and canonical fixture. Evidence was independently re-retrieved this session to confirm the behavior is correct and source-backed.
