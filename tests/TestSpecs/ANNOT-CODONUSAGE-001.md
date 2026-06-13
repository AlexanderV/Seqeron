# Test Specification: ANNOT-CODONUSAGE-001

**Test Unit ID:** ANNOT-CODONUSAGE-001
**Area:** Annotation
**Algorithm:** Relative Synonymous Codon Usage (RSCU)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sharp & Li (1986), Nucleic Acids Research 14(19):7737–7749 | 1 | https://doi.org/10.1093/nar/14.19.7737 | 2026-06-13 |
| 2 | PMC2528880 — synonymous codon usage of begomoviruses | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2528880/ | 2026-06-13 |
| 3 | LIRMM "RSCU RS" methods page (Rivals et al.) | 3 | https://www.lirmm.fr/~rivals/rscu/ | 2026-06-13 |
| 4 | CodonU `internal_comp.py` `rscu` (reference impl) | 3 | https://github.com/SouradiptoC/CodonU/blob/master/CodonU/analyzer/internal_comp.py | 2026-06-13 |
| 5 | NCBI Genetic Codes — Standard (transl_table=1) | 2 | https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi | 2026-06-13 |

### 1.2 Key Evidence Points

1. RSCU_{i,j} = n_i · x_{i,j} / Σ_j x_{i,j}, where n_i = number of synonymous codons for amino acid i and x_{i,j} = observed count of codon j — Source 3 (verbatim), Source 1 (origin).
2. RSCU = 1.0 means no bias; > 1 preferred, < 1 under-represented — Source 2.
3. RSCU is bounded in [0, n_i]; single-codon amino acids (Met, Trp) always have RSCU = 1.0 — Source 3, Source 5.
4. Counts are pooled across all reference (coding) sequences before computing RSCU — Source 4.
5. Only the 61 sense codons are included; stop codons are excluded (forward_table) — Source 4, Source 5.
6. Σ RSCU over a synonymous family = n_i when the family is observed (algebraic identity of the formula) — Source 3.

### 1.3 Documented Corner Cases

- Single-codon amino acids (Met/Trp): RSCU always 1.0 (Source 3/5).
- Aggregation over a list of sequences (Source 4).
- Stop codons excluded (Source 4/5).
- Whole-family zero count: base RSCU undefined (denominator 0); CAI 0.5 pseudocount is a separate CAI convention, not plain RSCU (Source 4).

### 1.4 Known Failure Modes / Pitfalls

1. Applying the CAI 0.5 pseudocount to plain RSCU would change values — only valid for CAI (Sharp & Li 1987), not RSCU — Source 4.
2. Reading frame must step by 3 from position 0; a partial trailing codon is ignored — Source 4 (frame stepping).
3. Using RNA vs DNA codon keys must be reconciled against the genetic-code table — Source 5.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetCodonUsage(IEnumerable<string> codingSequences)` | GenomeAnnotator | Canonical | RSCU over Standard code; deep evidence-based testing |
| `GetCodonUsage(IEnumerable<string> codingSequences, GeneticCode code)` | GenomeAnnotator | Delegate | Overload that lets caller pass a non-standard table; smoke only |
| `GetCodonUsage(string dnaSequence)` | GenomeAnnotator | Internal | Pre-existing raw codon-count method; unchanged, not in scope for RSCU |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | For every observed synonymous family, Σ RSCU over its codons = n_i (size of the family) | Yes | Source 3 (formula identity) |
| INV-2 | Single-codon amino acids (Met=ATG, Trp=TGG) have RSCU exactly 1.0 whenever observed | Yes | Source 3, Source 5 |
| INV-3 | RSCU value of any codon is in [0, n_i] | Yes | Source 3 |
| INV-4 | Stop codons (TAA, TAG, TGA) never appear in the output | Yes | Source 4, Source 5 |
| INV-5 | Uniform usage within a family ⇒ each member RSCU = 1.0 | Yes | Source 2, Source 3 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Leucine worked example | CDS `CTTCTTCTGTTA`; Leu family counts CTT=2,CTG=1,TTA=1, Σ=4, n_i=6 | RSCU CTT=3.0, CTG=1.5, TTA=1.5, TTG=0, CTC=0, CTA=0 | Source 3 formula + Source 5 family |
| M2 | Uniform usage = no bias | CDS `TTTTTC`; Phe TTT=1,TTC=1 | RSCU TTT=1.0, TTC=1.0 | Source 2, Source 3 |
| M3 | Single-codon amino acid | CDS `ATGATG`; Met n_i=1 | RSCU ATG=1.0 | Source 3, Source 5 |
| M4 | Pooling across sequences | Two CDS `["CTTCTT","CTGTTA"]` pooled = M1 counts | Same as M1 (CTT=3.0,CTG=1.5,TTA=1.5) | Source 4 |
| M5 | Stop codons excluded | CDS `ATGTAA` (Met + stop) | output contains ATG=1.0, contains no TAA/TAG/TGA key | Source 4, Source 5 |
| M6 | Family-sum invariant (INV-1) | M1 dataset | Σ RSCU over Leu family = 6.0 | Source 3 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Case-insensitive input | lower-case `cttcttctgtta` | identical to M1 | sibling-method convention |
| S2 | Partial trailing codon ignored | `ATGAT` (5 nt) | RSCU ATG=1.0; trailing `AT` ignored | Source 4 frame stepping |
| S3 | Non-standard table overload | `GetCodonUsage([...], GeneticCode.Standard)` equals default overload | delegation proven | Source 5 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Null input | `GetCodonUsage((IEnumerable<string>)null!)` | throws ArgumentNullException | edge completeness |
| C2 | Empty input | empty list (or all empty strings) | empty dictionary | edge completeness |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotatorTests.cs` lines 69–112 tests the pre-existing `GetCodonUsage(string)` raw-count method (not RSCU). Left untouched — it covers a different method.
- No tests exist for an RSCU `GetCodonUsage(codingSequences)` overload.

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
| S3 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |
| C2 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAnnotator_GetCodonUsage_Tests.cs` — all RSCU tests for the new overload.
- **Remove:** nothing. The legacy `GetCodonUsage(string)` raw-count tests in `GenomeAnnotatorTests.cs` remain (different method, out of scope).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomeAnnotator_GetCodonUsage_Tests.cs` | Canonical RSCU tests | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |
| 11 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact RSCU values asserted |
| M2 | ✅ Covered | 1.0 both codons |
| M3 | ✅ Covered | Met = 1.0 |
| M4 | ✅ Covered | pooling matches M1 |
| M5 | ✅ Covered | stops absent |
| M6 | ✅ Covered | family sum = 6.0 |
| S1 | ✅ Covered | case-insensitive |
| S2 | ✅ Covered | partial codon ignored |
| S3 | ✅ Covered | overload delegation |
| C1 | ✅ Covered | null throws |
| C2 | ✅ Covered | empty → empty |

**✅ count:** 11 of 11 in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Default genetic code = Standard (NCBI table 1) when none supplied (API default; overload accepts a `GeneticCode`) | default overload |
| 2 | Whole-family zero count ⇒ RSCU 0.0 for each member (avoids /0; CAI 0.5 pseudocount NOT applied — CAI-only convention) | zero-family handling |

---

## 7. Open Questions / Decisions

1. Decision: keep the legacy `GetCodonUsage(string)` raw-count method unchanged for backward compatibility (used by MCP tools); add RSCU as a new overload taking the coding-sequence collection named in the Registry.
2. Decision: do not apply the CAI 0.5 pseudocount — it belongs to CAI (Sharp & Li 1987), not to plain RSCU (Source 4 distinguishes them).
