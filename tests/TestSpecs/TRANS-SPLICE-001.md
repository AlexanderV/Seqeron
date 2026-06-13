# Test Specification: TRANS-SPLICE-001

**Test Unit ID:** TRANS-SPLICE-001
**Area:** Transcriptome
**Algorithm:** Alternative Splicing — Event Classification and Percent-Spliced-In (PSI)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Wang et al. (2008), Nature 456(7221):470–476 | 1 | https://doi.org/10.1038/nature07509 | 2026-06-13 |
| 2 | BMC Bioinformatics 13(Suppl 6):S11 (PSI definition) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3330053/ | 2026-06-13 |
| 3 | Shen et al. (2014) rMATS, PNAS 111(51):E5593 | 1 / 3 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4280593/ | 2026-06-13 |
| 4 | rMATS project docs (event types) | 3 | https://rmats.sourceforge.io/ | 2026-06-13 |
| 5 | Trincado et al. (2018) SUPPA2 | 3 | https://pubmed.ncbi.nlm.nih.gov/29571299/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. PSI = inclusion / (inclusion + exclusion); μ̃ = γᵢ/(γᵢ+γₑ) — Source 2.
2. rMATS length-normalized PSI: ψ̂ = (I/lᵢ)/(I/lᵢ + S/lₛ) — Source 3.
3. Five canonical AS event classes: skipped exon (SE), intron retention (RI), alternative 5′ splice site (A5SS), alternative 3′ splice site (A3SS), mutually exclusive exons (MXE) — Sources 1, 4.
4. Inclusion reads = exon body / junctions to adjacent exons; exclusion reads = junction between the two adjacent constitutive exons skipping the exon — Sources 2, 5.
5. An AS event is defined relative to two isoforms of one gene that differ in exon structure — Source 1.

### 1.3 Documented Corner Cases

- 0/0 (no supporting reads) → PSI undefined (Source 2, pseudo-count rationale).
- S=0 → PSI=1; I=0 → PSI=0 (direct from Ψ=I/(I+S)).
- lᵢ ≠ lₛ → normalized and unnormalized PSI differ (Source 3).
- Fewer than two isoforms for a gene → no event possible (Source 1).

### 1.4 Known Failure Modes / Pitfalls

1. Dividing by zero when an event has no supporting reads — Source 2.
2. Ignoring effective-length normalization biases PSI toward the longer isoform — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculatePSI(inclusionReads, exclusionReads, inclusionEffectiveLength?, exclusionEffectiveLength?)` | TranscriptomeAnalyzer | Canonical | Ψ=I/(I+S); rMATS length-normalized when both lengths > 0 |
| `DetectAlternativeSplicing(isoforms)` | TranscriptomeAnalyzer | Canonical | classifies isoform pairs of a gene into the 5 AS classes |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ PSI ≤ 1 for non-negative I,S with I+S>0 | Yes | Source 2 (part/whole ratio) |
| INV-2 | I+S=0 ⇒ PSI is NaN (undefined) | Yes | Source 2 |
| INV-3 | S=0,I>0 ⇒ PSI=1 ; I=0,S>0 ⇒ PSI=0 | Yes | Source 2 |
| INV-4 | length-normalized PSI uses ψ̂=(I/lᵢ)/(I/lᵢ+S/lₛ) when lᵢ,lₛ>0 | Yes | Source 3 |
| INV-5 | a detected event references two isoforms of the same gene differing in structure | Yes | Source 1 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | PSI read-count ratio | I=80,S=20 | 0.80 | Source 2; Source 5 |
| M2 | PSI length-normalized | I=80,S=20,lᵢ=200,lₛ=100 | 0.6666666667 | Source 3 |
| M3 | PSI full inclusion | I=50,S=0 | 1.0 | Source 2 (INV-3) |
| M4 | PSI full exclusion | I=0,S=40 | 0.0 | Source 2 (INV-3) |
| M5 | PSI undefined | I=0,S=0 | NaN | Source 2 (INV-2) |
| M6 | Classify SkippedExon | A=(1,100),(200,300),(400,500); B=(1,100),(400,500) | one event, SkippedExon | Source 1 |
| M7 | Classify RetainedIntron | A=(1,100),(200,300); B=(1,300) | one event, RetainedIntron | Source 1 |
| M8 | Classify AlternativeThreePrimeSS | A=(1,100),(200,300); B=(1,100),(200,350) | one event, AlternativeThreePrimeSS | Source 1 |
| M9 | Classify AlternativeFivePrimeSS | A=(1,100),(200,300); B=(1,150),(200,300) | one event, AlternativeFivePrimeSS | Source 1 |
| M10 | Classify MutuallyExclusiveExons | A=(1,100),(200,300),(500,600); B=(1,100),(350,400),(500,600) | one event, MutuallyExclusiveExons | Source 1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | PSI bounds invariant | varied I,S>0 | 0 ≤ PSI ≤ 1 | INV-1 |
| S2 | fewer than two isoforms | single isoform for a gene | no event | INV-5 (Source 1) |
| S3 | identical isoforms | A==B exons | no event | no structural difference |
| S4 | null/empty isoforms input | null or empty enumerable | empty result | contract |
| S5 | negative reads | I<0 or S<0 | ArgumentOutOfRangeException | contract (counts non-negative) |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | partial length params | only one effective length > 0 | falls back to unnormalized ratio | ASSUMPTION 1 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing tests for `CalculatePSI` or `DetectAlternativeSplicing` (neither method existed before this unit). Sibling file `TranscriptomeAnalyzer_ExpressionQuantification_Tests.cs` confirms conventions. Pre-existing `FindSkippedExonEvents`/`DetectDifferentialSplicing` are out of scope for this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ✅ Covered | implemented |
| M2 | ✅ Covered | implemented |
| M3 | ✅ Covered | implemented |
| M4 | ✅ Covered | implemented |
| M5 | ✅ Covered | implemented |
| M6 | ✅ Covered | implemented |
| M7 | ✅ Covered | implemented |
| M8 | ✅ Covered | implemented |
| M9 | ✅ Covered | implemented |
| M10 | ✅ Covered | implemented |
| S1 | ✅ Covered | implemented |
| S2 | ✅ Covered | implemented |
| S3 | ✅ Covered | implemented |
| S4 | ✅ Covered | implemented |
| S5 | ✅ Covered | implemented |
| C1 | ✅ Covered | implemented |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs` — all M/S/C cases for both methods.
- **Remove:** none.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| TranscriptomeAnalyzer_AlternativeSplicing_Tests.cs | canonical TRANS-SPLICE-001 fixture | 16 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | implemented | ✅ Done |
| 2 | M2 | ❌ Missing | implemented | ✅ Done |
| 3 | M3 | ❌ Missing | implemented | ✅ Done |
| 4 | M4 | ❌ Missing | implemented | ✅ Done |
| 5 | M5 | ❌ Missing | implemented | ✅ Done |
| 6 | M6 | ❌ Missing | implemented | ✅ Done |
| 7 | M7 | ❌ Missing | implemented | ✅ Done |
| 8 | M8 | ❌ Missing | implemented | ✅ Done |
| 9 | M9 | ❌ Missing | implemented | ✅ Done |
| 10 | M10 | ❌ Missing | implemented | ✅ Done |
| 11 | S1 | ❌ Missing | implemented | ✅ Done |
| 12 | S2 | ❌ Missing | implemented | ✅ Done |
| 13 | S3 | ❌ Missing | implemented | ✅ Done |
| 14 | S4 | ❌ Missing | implemented | ✅ Done |
| 15 | S5 | ❌ Missing | implemented | ✅ Done |
| 16 | C1 | ❌ Missing | implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact value 0.80 |
| M2 | ✅ | exact value 0.6666666667 (rMATS) |
| M3 | ✅ | 1.0 |
| M4 | ✅ | 0.0 |
| M5 | ✅ | NaN |
| M6 | ✅ | SkippedExon |
| M7 | ✅ | RetainedIntron |
| M8 | ✅ | AlternativeThreePrimeSS |
| M9 | ✅ | AlternativeFivePrimeSS |
| M10 | ✅ | MutuallyExclusiveExons |
| S1 | ✅ | bounds invariant |
| S2 | ✅ | no event for single isoform |
| S3 | ✅ | no event for identical isoforms |
| S4 | ✅ | null/empty → empty |
| S5 | ✅ | ArgumentOutOfRangeException |
| C1 | ✅ | partial length → unnormalized |

**Total in-scope cases:** 16 | **✅:** 16

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Length normalization is opt-in (both effective lengths > 0 switch to rMATS form); both forms source-backed | M2, C1 |
| 2 | Exons are ordered 5′→3′ on one strand (ascending coordinates) | M6–M10 |

---

## 7. Open Questions / Decisions

1. None. Both PSI forms and the five-class taxonomy are fully source-backed; the default (unnormalized) form matches the most widely cited definition.
