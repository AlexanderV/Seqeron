# Test Specification: VARIANT-INDEL-001

**Test Unit ID:** VARIANT-INDEL-001
**Area:** Variants
**Algorithm:** Indel Detection (`VariantCaller.FindInsertions` / `FindDeletions`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | VCFv4.3 specification (samtools/hts-specs), REF field & §1.1 examples | 2 | https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex | 2026-06-13 |
| 2 | Tan A, Abecasis GR, Kang HM (2015). Unified representation of genetic variants. *Bioinformatics* 31(13):2202–2204 | 1 | https://doi.org/10.1093/bioinformatics/btv112 | 2026-06-13 |
| 3 | ericminikel/minimal_representation `normalize.py` (implements Tan et al. 2015 Algorithm 1) | 3 | https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py | 2026-06-13 |
| 4 | PharmCAT — Variant Normalization worked examples | 3 | https://pharmcat.clinpgx.org/using/Variant-Normalization/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. A single-base insertion is encoded with ALT longer than REF: "A single base insertion of A after position 3 becomes REF=C, ALT=CA" — Source 1 (VCFv4.3 REF field).
2. A single-base deletion is encoded with REF longer than ALT: "A single base deletion of C at position 3 becomes REF=TC, ALT=T" — Source 1.
3. Insertion, deletion and SNP are distinct variant classes; an indel changes allele length relative to the reference (microsatellite example: REF `GTC` vs ALT `G` = 2-bp deletion; vs ALT `GTCT` = 1-bp insertion) — Source 1 (§1.1).
4. The same biological indel has multiple valid representations; "A VCF entry is *left aligned* if and only if its base position is smallest among all potential VCF entries having the same allele length and representing the same variant" — Source 2. Position is alignment/normalization dependent in repeats.
5. Reference trimming confirms the directional length invariant: deletion ⇒ len(REF) > len(ALT), e.g. `(7,117199646,CTT,-)` → `(7,117199644,ATCT,A)` — Source 3.

### 1.3 Documented Corner Cases

- **Pure indel allele padding (Source 1):** a pure insertion/deletion would leave REF or ALT empty in serialized VCF, so the spec requires a left anchor base. The repository's in-memory model instead uses the `"-"` gap sentinel for the absent allele (ASM-01); padded VCF is produced only by `ToVcfLines` (out of scope).
- **Non-unique indel placement in repeats (Source 2):** in a repeated region the indel column produced by the alignment is not unique; position is implementation-dependent unless normalized. Tests assert exact position only on inputs whose global alignment is provably unique (ASM-02).

### 1.4 Known Failure Modes / Pitfalls

1. Asserting an exact indel position on a repeat-region input where the alignment is ambiguous — would couple the test to a non-canonical placement (Source 2). Mitigated by using unique-alignment inputs for position assertions.
2. Treating a substitution as an indel — insertion/deletion are length-changing; a substitution keeps length (Source 1, distinct classes).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindInsertions(DnaSequence reference, DnaSequence query)` | `VariantCaller` | Canonical | Filters `CallVariants` to `VariantType.Insertion`; verifies insertion detection semantics (RefAllele = gap, AltAllele = inserted base, ALT longer than REF). |
| `FindDeletions(DnaSequence reference, DnaSequence query)` | `VariantCaller` | Canonical | Filters `CallVariants` to `VariantType.Deletion`; verifies deletion detection semantics (RefAllele = deleted base, AltAllele = gap, REF longer than ALT). |
| `FindIndels(DnaSequence reference, DnaSequence query)` | `VariantCaller` | Delegate | Union of insertions and deletions; smoke verification only. |

<!-- The underlying detection (CallVariantsFromAlignment) was validated under VARIANT-CALL-001;
     this unit tests the indel-specific semantics surfaced through the Find* filters. -->

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Identical sequences produce no insertions and no deletions. | Yes | Source 1 (an indel is a length-changing difference; none exist for identical input). |
| INV-02 | Every variant from `FindInsertions` has `Type == Insertion`; every variant from `FindDeletions` has `Type == Deletion`. | Yes | Source 1 (insertion/deletion are distinct classes); filter contract. |
| INV-03 | An Insertion column has `ReferenceAllele == "-"` and a non-gap `AlternateAllele` (ALT longer than REF). | Yes | Source 1 (insertion ⇒ ALT longer than REF); Source 3 (directional length). |
| INV-04 | A Deletion column has `AlternateAllele == "-"` and a non-gap `ReferenceAllele` (REF longer than ALT). | Yes | Source 1 (deletion ⇒ REF longer than ALT); Source 3. |
| INV-05 | For an input where a contiguous block of `k` bases is inserted (resp. deleted) within an otherwise unique alignment, the number of insertion (resp. deletion) columns equals `k`. | Yes | Source 1 (each extra/absent base is one indel event). |
| INV-06 | Every reported indel `Position` is within `[0, reference.Length]`. | Yes | Structural invariant (reference-coordinate bookkeeping). |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | `FindInsertions_SingleInsertedBase_ReturnsExactInsertion` | ref `ATGCAT` vs query `ATGTCAT`: a `T` inserted after ref index 2 (unique alignment). | 1 Insertion; `ReferenceAllele == "-"`, `AlternateAllele == "T"`, `Type == Insertion`; `Position == 3`. | Source 1 (insertion REF=C/ALT=CA, ALT longer than REF). |
| M2 | `FindDeletions_SingleDeletedBase_ReturnsExactDeletion` | ref `ATGTCAT` vs query `ATGCAT`: the `T` at ref index 3 deleted (unique alignment). | 1 Deletion; `ReferenceAllele == "T"`, `AlternateAllele == "-"`, `Type == Deletion`; `Position == 3`. | Source 1 (deletion REF=TC/ALT=T, REF longer than ALT). |
| M3 | `FindInsertions_ReturnsInsertionsOnly` | Input with one insertion and one substitution; `FindInsertions` filters. | Only Insertion(s) returned; no SNP, no Deletion. | Source 1 (distinct classes); filter contract (INV-02). |
| M4 | `FindDeletions_ReturnsDeletionsOnly` | Input with one deletion and one substitution; `FindDeletions` filters. | Only Deletion(s) returned; no SNP, no Insertion. | Source 1; INV-02. |
| M5 | `FindInsertions_MultiBaseInsertion_ReturnsConsecutiveColumns` | ref `ATGCAT` vs query `ATGTTCAT`: two `T` bases inserted after ref index 2 (unique alignment). | 2 Insertions; both `Type == Insertion`, `AlternateAllele == "T"`, `ReferenceAllele == "-"`; both at `Position == 3`. | Source 1 (microsatellite multi-base indel; INV-05). |
| M6 | `FindDeletions_MultiBaseDeletion_ReturnsConsecutiveColumns` | ref `ATGTTCAT` vs query `ATGCAT`: two `T` bases deleted at ref indices 3,4 (unique alignment). | 2 Deletions; both `Type == Deletion`, `ReferenceAllele == "T"`, `AlternateAllele == "-"`; positions `{3,4}`. | Source 1; INV-05. |
| M7 | `FindInsertions_IdenticalSequences_ReturnsEmpty` | ref `ATGCATGC` vs identical query. | No insertions. | Source 1 (no length-changing difference); INV-01. |
| M8 | `FindDeletions_IdenticalSequences_ReturnsEmpty` | ref `ATGCATGC` vs identical query. | No deletions. | Source 1; INV-01. |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | `FindInsertions_NullReference_Throws` | Null reference to `FindInsertions`. | `ArgumentNullException`. | Propagated from `CallVariants`. |
| S2 | `FindInsertions_NullQuery_Throws` | Null query to `FindInsertions`. | `ArgumentNullException`. | Propagated from `CallVariants`. |
| S3 | `FindDeletions_NullReference_Throws` | Null reference to `FindDeletions`. | `ArgumentNullException`. | Propagated from `CallVariants`. |
| S4 | `FindDeletions_NullQuery_Throws` | Null query to `FindDeletions`. | `ArgumentNullException`. | Propagated from `CallVariants`. |
| S5 | `FindInsertions_SubstitutionOnlyInput_ReturnsEmpty` | Equal-length input with one substitution, no length change. | No insertions. | Insertion distinct from substitution (Source 1). |
| S6 | `FindDeletions_SubstitutionOnlyInput_ReturnsEmpty` | Equal-length input with one substitution. | No deletions. | Deletion distinct from substitution (Source 1). |
| S7 | `FindIndels_InsertionAndDeletionInput_ReturnsBoth` | Input containing one insertion and one deletion; `FindIndels` delegate. | Both an Insertion and a Deletion returned; no SNP. | Delegate smoke test (union filter). |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | `FindInsertions_ContiguousBlock_CountEqualsBlockLength` | Property: insert a contiguous block of `k` bases in a unique alignment; count insertion columns. | Exactly `k` insertions, all `Type == Insertion`, all `Position` within reference bounds. | INV-05, INV-06; O(n×m) algorithm property test. |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/VariantCallerTests.cs` contained three pre-template indel tests: `FindInsertions_ReturnsOnlyInsertions`, `FindDeletions_ReturnsOnlyDeletions`, `FindIndels_ReturnsBothTypes`. Each used a permissive `All(v => v.Type == ...)` predicate with no count, position, allele, or message assertions — ⚠ Weak and duplicating this unit's scope.
- No `{Class}_{Method}_Tests.cs` file existed for indel detection prior to this unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 single insertion exact | ❌ Missing | No prior exact-allele/position test. |
| M2 single deletion exact | ❌ Missing | No prior exact-allele/position test. |
| M3 FindInsertions filters | ⚠ Weak | Old `FindInsertions_ReturnsOnlyInsertions` (permissive All, no message). |
| M4 FindDeletions filters | ⚠ Weak | Old `FindDeletions_ReturnsOnlyDeletions` (permissive All, no message). |
| M5 multi-base insertion | ❌ Missing | — |
| M6 multi-base deletion | ❌ Missing | — |
| M7 identical → no insertions | ❌ Missing | — |
| M8 identical → no deletions | ❌ Missing | — |
| S1–S4 null inputs | ❌ Missing | — |
| S5/S6 substitution-only → empty | ❌ Missing | — |
| S7 FindIndels both types | ⚠ Weak | Old `FindIndels_ReturnsBothTypes` (identical input, asserts nothing meaningful). |
| C1 contiguous-block property | ❌ Missing | — |

<!-- Status values: ✅ Covered, ⚠ Weak, ❌ Missing, 🔁 Duplicate -->

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_FindIndels_Tests.cs` — all VARIANT-INDEL-001 cases.
- **Remove:** the three ⚠ Weak tests from `VariantCallerTests.cs` (replaced with a pointer comment to the canonical fixture).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `VariantCaller_FindIndels_Tests.cs` | Canonical VARIANT-INDEL-001 fixture | 16 |
| `VariantCallerTests.cs` | Legacy fixture; indel tests removed (pointer comment) | 0 (indel) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented `FindInsertions_SingleInsertedBase_ReturnsExactInsertion`. | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented `FindDeletions_SingleDeletedBase_ReturnsExactDeletion`. | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote as `FindInsertions_InsertionAndSubstitutionInput_ReturnsInsertionsOnly`; deleted old test. | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewrote as `FindDeletions_DeletionAndSubstitutionInput_ReturnsDeletionsOnly`; deleted old test. | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented `FindInsertions_MultiBaseInsertion_ReturnsConsecutiveColumns`. | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented `FindDeletions_MultiBaseDeletion_ReturnsConsecutiveColumns`. | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented `FindInsertions_IdenticalSequences_ReturnsEmpty`. | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented `FindDeletions_IdenticalSequences_ReturnsEmpty`. | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented `FindInsertions_NullReference_ThrowsArgumentNullException`. | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented `FindInsertions_NullQuery_ThrowsArgumentNullException`. | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented `FindDeletions_NullReference_ThrowsArgumentNullException`. | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented `FindDeletions_NullQuery_ThrowsArgumentNullException`. | ✅ Done |
| 13 | S5 | ❌ Missing | Implemented `FindInsertions_SubstitutionOnlyInput_ReturnsEmpty`. | ✅ Done |
| 14 | S6 | ❌ Missing | Implemented `FindDeletions_SubstitutionOnlyInput_ReturnsEmpty`. | ✅ Done |
| 15 | S7 | ⚠ Weak | Rewrote as `FindIndels_InsertionAndDeletionInput_ReturnsBothTypes`; deleted old test. | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented `FindInsertions_ContiguousBlock_CountEqualsBlockLengthAndPositionsInBounds`. | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact allele/position, unique alignment. |
| M2 | ✅ Covered | Exact allele/position, unique alignment. |
| M3 | ✅ Covered | Filter verified with mixed insertion+substitution input. |
| M4 | ✅ Covered | Filter verified with mixed deletion+substitution input. |
| M5 | ✅ Covered | 2 consecutive insertion columns, exact positions. |
| M6 | ✅ Covered | 2 consecutive deletion columns, exact positions. |
| M7 | ✅ Covered | Identical → no insertions. |
| M8 | ✅ Covered | Identical → no deletions. |
| S1–S4 | ✅ Covered | Null reference/query throw on both methods. |
| S5/S6 | ✅ Covered | Substitution-only → empty for both methods. |
| S7 | ✅ Covered | FindIndels returns both types, no SNP. |
| C1 | ✅ Covered | Property: k=3 block → 3 insertion columns, positions in bounds. |

**In-scope cases:** 16 · **✅:** 16 · **❌/⚠ remaining:** 0

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| ASM-01 | In-memory indels use the `"-"` gap sentinel for the absent allele (not the VCF padded-allele form). Not correctness-affecting for detection (counts/types unchanged); representation only. | INV-03, INV-04, M1, M2, M5, M6 |
| ASM-02 | Indels are not left-aligned / parsimony-normalized; position is alignment-dependent in repeats. | M1, M2, M5, M6 (mitigated by unique-alignment inputs) |
| ASM-03 | `Variant.Position` is 0-based in the in-memory model (1-based only in `ToVcfLines`). | M1, M2, M5, M6 |

---

## 7. Open Questions / Decisions

1. **Decision:** VARIANT-INDEL-001 covers the indel-specific semantics surfaced by `FindInsertions`/`FindDeletions`; the shared detection core (`CallVariantsFromAlignment`) was already validated under VARIANT-CALL-001 and is not re-tested here for substitution behavior.
2. **Decision:** Exact position is asserted only on inputs whose Needleman–Wunsch global alignment is unique (no internal repeat permitting a left-shift), consistent with ASM-02 (Tan et al. 2015 non-uniqueness in repeats).
