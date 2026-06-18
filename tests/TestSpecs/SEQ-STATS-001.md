# Test Specification: SEQ-STATS-001

**Test Unit ID:** SEQ-STATS-001
**Area:** Statistics
**Algorithm:** Sequence Composition Statistics (nucleotide composition, GC content, GC/AT skew)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Lobry J.R. (1996) Asymmetric substitution patterns in the two DNA strands of bacteria. Mol Biol Evol 13(5):660–665 | 1 | https://doi.org/10.1093/oxfordjournals.molbev.a025626 | 2026-06-13 |
| 2 | Biopython `Bio.SeqUtils` (`gc_fraction`, `GC_skew`) source | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-13 |
| 3 | Wikipedia "GC skew" | 4 | https://en.wikipedia.org/wiki/GC_skew | 2026-06-13 |

### 1.2 Key Evidence Points

1. GC content = (G+C)/length, float in [0,1]; empty sequence returns 0 — Biopython `gc_fraction`.
2. GC skew = (G − C)/(G + C); 0 when G+C = 0 — Biopython `GC_skew`, Wikipedia.
3. AT skew = (A − T)/(A + T) — Wikipedia "GC skew".
4. Positive GC skew = G-rich; negative = C-rich — Wikipedia (Lobry 1996).
5. Counting is case-insensitive (lowercase included) — Biopython `GC_skew`.

### 1.3 Documented Corner Cases

- Empty sequence → composition all zero, GC content 0 (Biopython `gc_fraction`).
- No G/C → GC skew 0 (Biopython catches `ZeroDivisionError`); by symmetry no A/T → AT skew 0.
- Mixed case handled (lowercase counted).

### 1.4 Known Failure Modes / Pitfalls

1. Sign convention: Lobry's original (C−G)/(C+G) vs modern (G−C)/(G+C); the implementation uses the modern (G−C)/(G+C) — Wikipedia "GC skew".
2. Degenerate IUPAC symbols (S/W) handled differently than Biopython — out of standard-alphabet scope (see §6 assumption).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateNucleotideComposition(string)` | SequenceStatistics | Canonical | Counts, GC/AT content, GC/AT skew |
| `SummarizeNucleotideSequence(string)` | SequenceStatistics | Delegate | Aggregates composition + entropy/complexity/Tm |
| `CalculateAminoAcidComposition(string)` | SequenceStatistics | Delegate | Residue counts/ratios; MW/pI/hydro belong to SEQ-MW/PI/HYDRO |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ GcContent ≤ 1 and 0 ≤ AtContent ≤ 1 | Yes | Biopython `gc_fraction` returns [0,1] |
| INV-2 | −1 ≤ GcSkew ≤ 1 and −1 ≤ AtSkew ≤ 1 | Yes | (x−y)/(x+y) with x,y ≥ 0 ∈ [−1,1] |
| INV-3 | CountA+CountT+CountG+CountC+CountU+CountN+CountOther = Length | Yes | Counts partition the sequence |
| INV-4 | GcSkew = 0 when CountG+CountC = 0; AtSkew = 0 when CountA+CountT = 0 | Yes | Biopython zero-division handling |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact counts | `AAUUGGCC` counts | A2 T0 G2 C2 U2 N0 Other0 Length8 | Composition definition |
| M2 | GC content | `ATGC` → (G+C)/total | 0.5 | Biopython `gc_fraction` |
| M3 | GC content all-GC | `GGGC` | 1.0 | Biopython `gc_fraction` |
| M4 | AT content | `ATGC` | 0.5 | (A+T+U)/total |
| M5 | GC skew positive | `GGGC` → (3−1)/4 | 0.5 | Wikipedia / Biopython `GC_skew` |
| M6 | GC skew negative | `GCCC` → (1−3)/4 | −0.5 | Wikipedia / Biopython `GC_skew` |
| M7 | AT skew | `AAAT` → (3−1)/4 | 0.5 | Wikipedia "GC skew" |
| M8 | Empty | `""` | Length 0, all counts 0, GcContent 0 | Biopython `gc_fraction` empty |
| M9 | Null | `null` | same as empty, no throw | Existing contract (null-safe) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive | `atgc` == `ATGC` | identical composition | Biopython counts lowercase |
| S2 | GC skew zero denom | `AAAT` (no G/C) | GcSkew 0 | Zero-division handling |
| S3 | AT skew zero denom | `GGGC` (no A/T) | AtSkew 0 | Zero-division handling |
| S4 | N/Other counts | `ATGCNNXX` | N2, Other2, Length8 | Non-ACGTU partition |
| S5 | Counts sum invariant | arbitrary seq | INV-3 holds | Partition invariant |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Summary delegation | `SummarizeNucleotideSequence("ATGC")` | Length 4, GcContent 0.5 | Delegate smoke |
| C2 | AA composition | `CalculateAminoAcidComposition("MKVLWA")` | exact residue counts, Length 6 | Delegate smoke |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatisticsTests.cs` contains pre-template tests for these methods. Several use permissive assertions (GcSkew `GreaterThan(0)`/`LessThan(0)`; Summary `GreaterThan(0)`).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 exact counts | ⚠ Weak | Existing checks some counts but not full partition |
| M2/M3 GC content | ⚠ Weak | Uses `.Within(0.001)`, not 1e-10; not in canonical file |
| M4 AT content | ❌ Missing | Not asserted |
| M5/M6 GC skew | ⚠ Weak | `GreaterThan(0)`/`LessThan(0)` — no exact value |
| M7 AT skew | ❌ Missing | Not asserted |
| M8 empty | ⚠ Weak | Partial |
| M9 null | ⚠ Weak | `DoesNotThrow` only |
| S1 case | ⚠ Weak | Checks GcContent only |
| S2/S3 zero-denom skew | ❌ Missing | Not asserted |
| S4 N/Other | ⚠ Weak | Partial |
| S5 invariant | ❌ Missing | Not asserted |
| C1 summary | ⚠ Weak | `GreaterThan(0)` |
| C2 AA composition | ⚠ Weak | Not exact |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateNucleotideComposition_Tests.cs` — all SEQ-STATS-001 cases with exact evidence values.
- **Remove:** none. Legacy `SequenceStatisticsTests.cs` is retained because it also covers other units (profiles, MW, entropy); its weak SEQ-STATS-001 cases are superseded by the canonical file. No duplication of canonical-file intent is added.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_CalculateNucleotideComposition_Tests.cs` | Canonical SEQ-STATS-001 | 16 |
| `SequenceStatisticsTests.cs` | Legacy (other units; pre-existing) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | New exact-count test | ✅ Done |
| 2 | M2 | ⚠ Weak | New exact GC content | ✅ Done |
| 3 | M3 | ⚠ Weak | New all-GC test | ✅ Done |
| 4 | M4 | ❌ Missing | New AT content test | ✅ Done |
| 5 | M5 | ⚠ Weak | Exact +0.5 | ✅ Done |
| 6 | M6 | ⚠ Weak | Exact −0.5 | ✅ Done |
| 7 | M7 | ❌ Missing | Exact AT skew 0.5 | ✅ Done |
| 8 | M8 | ⚠ Weak | Exact empty zeros | ✅ Done |
| 9 | M9 | ⚠ Weak | Null = empty | ✅ Done |
| 10 | S1 | ⚠ Weak | Full composition equality | ✅ Done |
| 11 | S2 | ❌ Missing | GC skew 0 no G/C | ✅ Done |
| 12 | S3 | ❌ Missing | AT skew 0 no A/T | ✅ Done |
| 13 | S4 | ⚠ Weak | Exact N/Other | ✅ Done |
| 14 | S5 | ❌ Missing | Partition invariant | ✅ Done |
| 15 | C1 | ⚠ Weak | Summary delegation exact | ✅ Done |
| 16 | C2 | ⚠ Weak | AA composition exact | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact A/T/G/C/U/N/Other + Length |
| M2 | ✅ Covered | GC content 0.5 within 1e-10 |
| M3 | ✅ Covered | GC content 1.0 |
| M4 | ✅ Covered | AT content 0.5 |
| M5 | ✅ Covered | GC skew 0.5 |
| M6 | ✅ Covered | GC skew −0.5 |
| M7 | ✅ Covered | AT skew 0.5 |
| M8 | ✅ Covered | Empty all-zero |
| M9 | ✅ Covered | Null = empty |
| S1 | ✅ Covered | Case-insensitive full equality |
| S2 | ✅ Covered | GC skew 0 |
| S3 | ✅ Covered | AT skew 0 |
| S4 | ✅ Covered | N2/Other2 |
| S5 | ✅ Covered | Counts sum = Length |
| C1 | ✅ Covered | Summary delegation |
| C2 | ✅ Covered | AA residue counts |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Degenerate IUPAC codes (S/W/…) not counted toward composition totals (standard {A,T,G,C,U} scope). Non-correctness-affecting within scope; differs from Biopython only on degenerate symbols. | §1.4, algorithm doc §5.3 |

---

## 7. Open Questions / Decisions

1. Decision: AT content denominator includes U (A+T+U) but AT skew uses (A−T)/(A+T) — AT skew is DNA-specific per the Lobry/Wikipedia formula; resolved, no open question.
