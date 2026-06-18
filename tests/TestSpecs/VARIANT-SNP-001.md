# Test Specification: VARIANT-SNP-001

**Test Unit ID:** VARIANT-SNP-001
**Area:** Variants
**Algorithm:** SNP Detection (single-nucleotide substitution identification: `FindSnps`, `FindSnpsDirect`)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | VCFv4.3 specification (samtools/hts-specs) | 2 | https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex | 2026-06-13 |
| 2 | Wikipedia — Transversion (primary: Futuyma 2013) | 4 | https://en.wikipedia.org/wiki/Transversion | 2026-06-13 |
| 3 | Hamming Distance as a Concept in DNA Molecular Recognition (PMC5410656) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5410656/ | 2026-06-13 |
| 4 | Collins & Jukes (1994), *Genomics* 20(3):386–396 | 1 | https://doi.org/10.1006/geno.1994.1192 | 2026-06-13 |

### 1.2 Key Evidence Points

1. A SNP is a single-base substitution: single-base REF replaced by single-base ALT (VCFv4.3 §1.1 "a good simple SNP" `G`→`A`); REF == ALT is not a variant — Source 1.
2. POS is 1-based in VCF; the in-memory `Variant.Position` is 0-based (sibling contract) — Source 1.
3. REF/ALT bases are A,C,G,T,N case-insensitive, so SNP comparison and Ti/Tv classification are case-insensitive — Source 1.
4. Positional SNP detection over two equal-length sequences = enumerating Hamming mismatch positions; each mismatch is one substitution (one SNP) — Source 3.
5. Transversions are A↔C, A↔T, G↔C, G↔T (purine↔pyrimidine); transitions are A↔G, C↔T (same ring class) — Source 2.

### 1.3 Documented Corner Cases

- REF == ALT (matching column) is not a SNP and must be skipped — Source 1.
- Hamming distance is defined for equal-length strings only; unequal-length inputs to `FindSnpsDirect` are compared over the common prefix, the trailing region of the longer sequence is out of scope (indels) — Source 3 (ASSUMPTION-1).
- Lowercase bases must classify identically to uppercase — Source 1.

### 1.4 Known Failure Modes / Pitfalls

1. Treating an indel region as a SNP — a SNP is a substitution only; insertions/deletions are distinct classes (VARIANT-INDEL-001) — Source 1.
2. Comparing beyond the common length for unequal inputs — undefined for substitution semantics — Source 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindSnpsDirect(string reference, string query)` | VariantCaller | **Canonical** | Positional (Hamming-style) substitution detection, no alignment. Deep evidence-based testing. |
| `FindSnps(DnaSequence reference, DnaSequence query)` | VariantCaller | **Delegate** | Aligns then filters `CallVariants` to `VariantType.SNP`. Smoke verification of delegation + SNP-only filtering. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Identical equal-length sequences yield zero SNPs (Hamming distance 0). | Yes | Source 3 |
| INV-02 | Every emitted variant has `Type == VariantType.SNP`. | Yes | Source 1 |
| INV-03 | Every emitted SNP has `ReferenceAllele != AlternateAllele` (a substitution). | Yes | Source 1 |
| INV-04 | `FindSnpsDirect` reports a SNP at each 0-based mismatch index `i` with `Position == i`, `ReferenceAllele == ref[i]`, `AlternateAllele == query[i]`. | Yes | Sources 1, 3 |
| INV-05 | For two equal-length sequences, the SNP count from `FindSnpsDirect` equals their Hamming distance. | Yes | Source 3 |
| INV-06 | `FindSnpsDirect` compares only the common prefix `min(len(ref), len(query))` of unequal-length inputs. | Yes | Source 3 (ASSUMPTION-1) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Direct identical → no SNPs | `FindSnpsDirect("ATGC","ATGC")` | empty (Hamming distance 0) | Source 3; INV-01 |
| M2 | Direct single substitution | `FindSnpsDirect("ATGC","ATTC")` | 1 SNP: Position 2, REF "G", ALT "T", QueryPosition 2 | Sources 1, 3; INV-04 |
| M3 | Direct multiple substitutions | `FindSnpsDirect("AAAA","TGTA")` | 3 SNPs at positions {0,1,2}; ALTs {"T","G","T"}; all REF "A" | Source 3; INV-04 |
| M4 | Direct all SNP type & distinct alleles | `FindSnpsDirect("AAAA","TGTA")` | every variant Type==SNP and REF≠ALT | Source 1; INV-02, INV-03 |
| M5 | Direct unequal length → common prefix only | `FindSnpsDirect("ATGCAA","ATTC")` | only prefix indices 0–3 compared; 1 SNP at Position 2 (G→T); index 3 matches (C==C); trailing "AA" ignored | Source 3; INV-06 |
| M6 | `FindSnps` substitution-only input → SNPs only | `FindSnps(ref "ATGCATGC", query "ATGAATGC")` | exactly 1 variant, Type SNP, REF "C", ALT "A", Position 3; no insertions/deletions | Source 1; INV-02 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Direct empty input | `FindSnpsDirect("","")` | empty | documented input contract |
| S2 | Direct one empty operand | `FindSnpsDirect("ATGC","")` | empty | nothing to compare |
| S3 | `FindSnps` null reference | `FindSnps(null, query)` | throws `ArgumentNullException` | propagated from `CallVariants` validation |
| S4 | `FindSnps` null query | `FindSnps(ref, null)` | throws `ArgumentNullException` | propagated from `CallVariants` validation |
| S5 | `FindSnps` identical sequences | `FindSnps(ref, ref-equal)` | empty | delegation: no differences → no SNPs |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property: count == Hamming distance | over a deterministic equal-length pair, `FindSnpsDirect` count equals number of differing positions | counts equal; all SNP; REF≠ALT | INV-05; O(n) but a structural-invariant property test |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for SNP-specific tests of `FindSnps`/`FindSnpsDirect`.
- `VariantCaller_CallVariants_Tests.cs` (VARIANT-CALL-001) exercises `CallVariants`/`CallVariantsFromAlignment`/`ClassifyMutation`/`CalculateTiTvRatio` but **not** `FindSnps` or `FindSnpsDirect`.
- Legacy `VariantCallerTests.cs` exists; it does not provide an evidence-based, canonical test of `FindSnpsDirect` with exact positions/alleles.
- New canonical file for this unit: `VariantCaller_FindSnps_Tests.cs` (created in Phase 7).

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
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_FindSnps_Tests.cs` — all MUST/SHOULD/COULD cases for `FindSnps` and `FindSnpsDirect`.
- **Remove:** nothing. Legacy `VariantCallerTests.cs` and `VariantCaller_CallVariants_Tests.cs` do not cover these two methods, so no duplication is introduced.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `VariantCaller_FindSnps_Tests.cs` | Canonical VARIANT-SNP-001 fixture | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented `FindSnpsDirect_IdenticalSequences_ReturnsNoSnps` | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented `FindSnpsDirect_SingleSubstitution_ReturnsExactSnp` | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented `FindSnpsDirect_MultipleSubstitutions_ReturnsSnpsAtExactPositions` | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented `FindSnpsDirect_AllResults_AreSnpsWithDistinctAlleles` | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented `FindSnpsDirect_UnequalLengths_ComparesCommonPrefixOnly` | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented `FindSnps_SubstitutionOnlyInput_ReturnsSnpsOnly` | ✅ Done |
| 7 | S1 | ❌ Missing | Implemented `FindSnpsDirect_EmptyInput_ReturnsEmpty` | ✅ Done |
| 8 | S2 | ❌ Missing | Implemented `FindSnpsDirect_OneEmptyOperand_ReturnsEmpty` | ✅ Done |
| 9 | S3 | ❌ Missing | Implemented `FindSnps_NullReference_ThrowsArgumentNullException` | ✅ Done |
| 10 | S4 | ❌ Missing | Implemented `FindSnps_NullQuery_ThrowsArgumentNullException` | ✅ Done |
| 11 | S5 | ❌ Missing | Implemented `FindSnps_IdenticalSequences_ReturnsNoSnps` | ✅ Done |
| 12 | C1 | ❌ Missing | Implemented `FindSnpsDirect_AnyEqualLengthPair_CountEqualsHammingDistance` | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | exact empty assertion |
| M2 | ✅ Covered | exact position/REF/ALT/QueryPosition |
| M3 | ✅ Covered | exact positions + ALT sequence |
| M4 | ✅ Covered | type + distinct-allele assertions |
| M5 | ✅ Covered | common-prefix-only behavior asserted |
| M6 | ✅ Covered | SNP-only + exact allele/position |
| S1 | ✅ Covered | empty → empty |
| S2 | ✅ Covered | one empty operand → empty |
| S3 | ✅ Covered | null reference throws |
| S4 | ✅ Covered | null query throws |
| S5 | ✅ Covered | identical → empty (delegation) |
| C1 | ✅ Covered | count == Hamming distance property |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | `FindSnpsDirect` compares only the common prefix `min(len(ref),len(query))` of unequal-length inputs (Hamming defined for equal length only). | M5, INV-06 |
| 2 | In-memory `Variant.Position` is 0-based (VCF 1-based POS applies only to serialized `ToVcfLines`). | M2, M3, M5, M6, INV-04 |

---

## 7. Open Questions / Decisions

1. None. Both methods conform to retrieved sources (VCF SNP = single-base substitution; positional detection = Hamming mismatch enumeration); no correctness-affecting constants are involved.
