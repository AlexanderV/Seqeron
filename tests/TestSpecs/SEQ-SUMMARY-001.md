# Test Specification: SEQ-SUMMARY-001

**Test Unit ID:** SEQ-SUMMARY-001
**Area:** Statistics
**Algorithm:** Sequence Summary (aggregation of length, GC content, Shannon entropy, linguistic complexity, melting temperature, and nucleotide composition into one record)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Biopython `Bio.SeqUtils` (`gc_fraction`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py | 2026-06-14 |
| 2 | Biopython `Bio.SeqUtils.MeltingTemp` (`Tm_Wallace`, `Tm_GC`) | 3 | https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py | 2026-06-14 |
| 3 | Wikipedia "Entropy (information theory)" (citing Shannon 1948) | 4 | https://en.wikipedia.org/wiki/Entropy_(information_theory) | 2026-06-14 |
| 4 | Wikipedia "Linguistic sequence complexity" (citing Trifonov 1990) | 4 | https://en.wikipedia.org/wiki/Linguistic_sequence_complexity | 2026-06-14 |

### 1.2 Key Evidence Points

1. GC fraction = GC count / counted-base length, case-insensitive, returns 0 for empty — Source 1.
2. Wallace Tm = 2(A+T) + 4(G+C) for short oligos; GC/Marmur-Doty Tm for longer sequences — Source 2.
3. Shannon entropy H = −Σ p·log₂ p in bits; uniform over k symbols gives H = log₂ k — Source 3.
4. Linguistic complexity = vocabulary usage (observed/possible words) combined across word sizes, 0<C<1 — Source 4.
5. Aggregation correctness: every summary field must equal its canonical per-metric method's value on the same input — aggregation-consistency requirement (Evidence §Documented Corner Cases).

### 1.3 Documented Corner Cases

- Empty sequence → GcContent 0 (Source 1), and by aggregation Length 0 / Entropy 0 / Complexity 0 / Tm 0.
- Length-dependent Tm formula selection: Wallace for length < 14, GC/Marmur-Doty otherwise (Source 2 + SEQ-TM-001 convention).

### 1.4 Known Failure Modes / Pitfalls

1. Divergence between a summary field and its canonical method (different rounding/alphabet/formula selection) — defect of the aggregation. Evidence §Aggregation-specific.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `SummarizeNucleotideSequence(string)` | SequenceStatistics | **Canonical** | Aggregates 4 metrics + composition into `SequenceSummary` |
| `CalculateNucleotideComposition`, `CalculateShannonEntropy`, `CalculateLinguisticComplexity`, `CalculateMeltingTemperature` | SequenceStatistics | **Internal** | Per-metric canonical methods; tested as the reference each summary field must equal (own conformance covered by sibling units) |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `summary.Length == sequence.Length` | Yes | Aggregation contract (record field copies `comp.Length`) |
| INV-2 | `summary.GcContent == CalculateNucleotideComposition(seq).GcContent` | Yes | Source 1; aggregation-consistency |
| INV-3 | `summary.Entropy == CalculateShannonEntropy(seq)` | Yes | Source 3; aggregation-consistency |
| INV-4 | `summary.Complexity == CalculateLinguisticComplexity(seq)` | Yes | Source 4; aggregation-consistency |
| INV-5 | `summary.MeltingTemperature == CalculateMeltingTemperature(seq, useWallaceRule: true)` (switch on A+C+G+T+U < 14, U read as T; 2026-09 B03 F12/F20) | Yes | Source 2; aggregation-consistency |
| INV-6 | Composition dict counts A,T,G,C,U,N equal `CalculateNucleotideComposition` counts | Yes | Source 1; aggregation-consistency |
| INV-7 | 0 ≤ GcContent ≤ 1 and 0 ≤ Complexity ≤ 1 over ≤ 4 distinct symbols (ACGTN → 15/14) | Yes | Source 1, Source 4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Worked example exact values | "ATGCATGC" full summary | Length 8; GcContent 0.5; Entropy 2.0; Tm 24.0; Complexity = canonical LC | Sources 1–4 derivations |
| M2 | All fields equal canonical methods | "ATGCATGCATGCATGC" each field == its per-metric method | exact equality on every field | INV-2..INV-6 |
| M3 | Tm GC/Marmur-Doty branch (len≥14) | 16-mer GC=8 Tm via GC formula | Tm = 43.375 | Source 2 + SEQ-TM-001 constants |
| M4 | Tm Wallace branch (len<14) | "ATGC" Tm = 2(A+T)+4(G+C) | Tm = 12.0; equals canonical with useWallaceRule:true | Source 2 |
| M5 | Composition dict matches counts | RNA+N input "AUGCNNA" | counts A,U,G,C,N match composition; T=0 | Source 1; INV-6 |
| M6 | Empty sequence | "" → degenerate summary | Length 0, GcContent 0, Entropy 0, Complexity 0, Tm 0 | Source 1 empty handling |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null input | `null` argument | same degenerate summary as empty (no throw) | per-metric methods guard `IsNullOrEmpty` |
| S2 | Case-insensitivity | "atgcatgc" == "ATGCATGC" summary | identical fields | per-metric methods uppercase internally |

### 4.2b Review additions (2026-09 B03)

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| R1 | `SummarizeNucleotideSequence_MatchesPythonReferences` (5) | every field = Biopython `gc_fraction` / scipy entropy / Biopython `Tm_Wallace`/`Tm_GC` / exact-fraction mean U_k | executed Python references |
| R2 | `SummarizeNucleotideSequence_RnaSpelling_EqualsDnaSpelling` (3) | RNA == DNA spelling for GC, entropy, complexity, Tm | Biopython `_check` back-transcription (F20) |
| R3 | `CalculateLinguisticComplexity_MatchesSumFormReference` (10) | Σ V_k / Σ V_max,k: ATTTGGATT m=9 → 7/8 (Rosalind LING 0.875), m=6 → 29/34; ACGTN 15/14 | Orlov & Potapov 2004 / Troyanskaya 2002; executed Python reference (B03 F21) |
| R4 | `CalculateLinguisticComplexity_EqualsCanonicalSequenceComplexity` | 1000 random inputs, m ∈ [−1,19]: identical to `SequenceComplexity.CalculateLinguisticComplexity` | delegation lock (B03 F21) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Bounds invariant | random-ish DNA fragment | 0≤GcContent≤1, 0≤Complexity<1 | INV-7 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `SummarizeNucleotideSequence` and `SequenceSummary`. No existing test exercises the summary aggregation (`SequenceStatisticsTests.cs` covers other methods). No prior `SequenceStatistics_SummarizeNucleotideSequence_Tests.cs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_SummarizeNucleotideSequence_Tests.cs` — all SEQ-SUMMARY-001 cases.
- **Remove:** none (no pre-existing summary tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceStatistics_SummarizeNucleotideSequence_Tests.cs` | Canonical unit test file | 9 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented field==canonical test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented Marmur-Doty branch test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented Wallace branch test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented composition-dict test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented empty-input test | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented null-input test | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented bounds-invariant test | ✅ Done |

**Total items:** 9
**✅ Done:** 9 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact worked-example values |
| M2 | ✅ Covered | Every field == canonical method |
| M3 | ✅ Covered | Tm 43.375 (GC branch) |
| M4 | ✅ Covered | Tm 12.0 (Wallace branch) |
| M5 | ✅ Covered | Composition dict A,U,G,C,N |
| M6 | ✅ Covered | Degenerate summary on empty |
| S1 | ✅ Covered | Null → degenerate summary |
| S2 | ✅ Covered | Lowercase identical |
| C1 | ✅ Covered | Bounds invariant |

In-scope cases: 9; ✅ = 9.

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Tm formula-selection threshold length<14 (sibling SEQ-TM-001 convention; summary tested for equality with the canonical method, so non-correctness-affecting for this unit) | INV-5, M3, M4 |

---

## 7. Open Questions / Decisions

1. ~~The `Complexity` field was the unsourced mean of per-k vocabulary usages.~~ Resolved 2026-09 (B03 F21): `CalculateLinguisticComplexity` delegates to the canonical `SequenceComplexity` sum form Σ V_k / Σ min(4^k, N−k+1) (Orlov & Potapov 2004; Troyanskaya 2002 / Rosalind LING at m ≥ N). Trifonov's product Π U_k (Source 4) is a different published measure and is not used.
