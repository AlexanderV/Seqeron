# Test Specification: VARIANT-CALL-001

**Test Unit ID:** VARIANT-CALL-001
**Area:** Variants
**Algorithm:** Variant Detection (SNP / insertion / deletion calling + transition/transversion classification)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Danecek P, et al. (2011) The variant call format and VCFtools. *Bioinformatics* 27(15):2156–2158 | 1 | https://doi.org/10.1093/bioinformatics/btr330 | 2026-06-13 |
| 2 | VCFv4.3 specification, samtools/hts-specs | 2 | https://samtools.github.io/hts-specs/VCFv4.3.pdf | 2026-06-13 |
| 3 | Collins DW, Jukes TH (1994) *Genomics* 20(3):386–396 | 1 | https://doi.org/10.1006/geno.1994.1192 | 2026-06-13 |
| 4 | Tan A, Abecasis GR, Kang HM (2015) *Bioinformatics* 31(13):2202–2204 | 1 | https://doi.org/10.1093/bioinformatics/btv112 | 2026-06-13 |
| 5 | Wikipedia — Transition (genetics) (primary: Collins & Jukes 1994) | 4 | https://en.wikipedia.org/wiki/Transition_(genetics) | 2026-06-13 |
| 6 | Wikipedia — Transversion (primary: Futuyma 2013) | 4 | https://en.wikipedia.org/wiki/Transversion | 2026-06-13 |

### 1.2 Key Evidence Points

1. A variant is any difference (SNP, insertion, deletion) of a query/sample from a reference — source 1.
2. POS is 1-based; REF/ALT bases are A,C,G,T,N (case-insensitive) — source 2 (fields 2, 4, 5).
3. VCF indel padding base: pure ins/del must include the base before the event in REF and ALT — source 2 (field 4). Applies to serialized VCF (`ToVcfLines`), not the in-memory `Variant` model.
4. Transition = A↔G or C↔T (purine↔purine / pyrimidine↔pyrimidine) — sources 5, 3.
5. Transversion = purine↔pyrimidine; the four are A↔C, A↔T, G↔C, G↔T — source 6.
6. Transitional bias (rates 1.71 vs 1.22 ×10⁻⁹) — source 3; basis for Ti/Tv being a meaningful ratio.
7. Canonical (normalized) indel = left-aligned + parsimonious — source 4; bounds the position-correctness claim of the alignment-based caller.

### 1.3 Documented Corner Cases

- Empty REF/ALT is forbidden by VCF (source 2, field 4); the in-memory model uses a `"-"` gap sentinel instead.
- Indel placement is non-unique in repeated regions unless normalized (source 4).
- REF/ALT classification is case-insensitive (source 2, field 4).

### 1.4 Known Failure Modes / Pitfalls

1. Reporting an indel at a non-left-aligned position in a repeat — source 4. Mitigated in tests by asserting position only on unambiguous inputs.
2. Treating Ti/Tv with zero transversions as +∞ — undefined; the repo contract returns 0 (ASSUMPTION ASM-03).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CallVariants(DnaSequence reference, DnaSequence query)` | VariantCaller | Canonical | Aligns then delegates to `CallVariantsFromAlignment`; the unit's headline entry point. |
| `CallVariantsFromAlignment(string alignedReference, string alignedQuery)` | VariantCaller | Canonical | Column scan over a gapped alignment; SNP/insertion/deletion classification + position bookkeeping. |
| `ClassifyMutation(Variant variant)` | VariantCaller | Canonical | Transition/transversion/Other classification (checklist names this `ClassifyVariant`; reconciled below). |
| `CalculateTiTvRatio(IEnumerable<Variant>)` | VariantCaller | Delegate | Counts via `ClassifyMutation`; thin aggregation. |

> Naming reconciliation: the Processing Registry lists `ClassifyVariant(variant)`. No such method exists; the implemented and source-faithful classification method is `ClassifyMutation`, which performs the transition/transversion classification the checklist intends. The checklist is workflow-control only; the behavioral contract follows the sources. No rename performed (would break sibling tests); recorded in §7.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | Identical reference and query produce zero variants. | Yes | Source 1 (variant = difference from reference) |
| INV-02 | Every emitted SNP has `ReferenceAllele != AlternateAllele` (single base each). | Yes | Source 2 (a SNP is a single-base substitution) |
| INV-03 | Every emitted variant's 0-based `Position` lies in `[0, reference.Length]`. | Yes | Structural (position is a reference coordinate) |
| INV-04 | A column with ref-gap is an Insertion; a column with query-gap is a Deletion; a mismatched pair is a SNP. | Yes | Source 1 (the three variant classes) |
| INV-05 | `ClassifyMutation` returns Transition iff {ref,alt}⊆{A,G} or ⊆{C,T}; Transversion for any other SNP pair; Other for non-SNP. | Yes | Sources 5, 6 |
| INV-06 | `CalculateTiTvRatio` = (#transitions)/(#transversions) over SNPs, or 0 when #transversions = 0. | Yes | Definition + source 3 (ASM-03) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Identical sequences | `CallVariants("ATGCATGC","ATGCATGC")` | empty | Source 1; INV-01 |
| M2 | Single SNP via alignment | `CallVariantsFromAlignment("ATGC","ATTC")` | 1 SNP, Position=2, REF="G", ALT="T" | Source 2 (simple SNP); INV-02/04 |
| M3 | Single SNP via `CallVariants` | `CallVariants("ATGC","ATTC")` | 1 SNP, Position=2, REF="G", ALT="T" | Source 2; INV-04 |
| M4 | Insertion column | `CallVariantsFromAlignment("AT-GC","ATTGC")` | 1 Insertion, Position=2, REF="-", ALT="T", QueryPosition=2 | Source 2 (microsat ins T); INV-04 |
| M5 | Deletion column | `CallVariantsFromAlignment("ATTGC","AT-GC")` | 1 Deletion, Position=2, REF="T", ALT="-" | Source 2 (microsat del); INV-04 |
| M6 | Two-base deletion (microsatellite) | `CallVariantsFromAlignment("GTCAA","G--AA")` | 2 Deletions of T then C at consecutive ref positions 1,2 | Source 2 §1.1 (deletion of 2 bases TC); INV-04 |
| M7 | Transition A→G | `ClassifyMutation(SNP A→G)` | Transition | Source 5; INV-05 |
| M8 | Transition C→T | `ClassifyMutation(SNP C→T)` | Transition | Source 5; INV-05 |
| M9 | Transition G→A and T→C | both | Transition | Source 5; INV-05 |
| M10 | Transversion A→C / A→T / G→C / G→T | each | Transversion | Source 6; INV-05 |
| M11 | Classification case-insensitive | `ClassifyMutation(SNP a→g)` | Transition | Source 2 (case-insensitive); INV-05 |
| M12 | Non-SNP classification | `ClassifyMutation(Deletion)` | Other | INV-05 |
| M13 | Ti/Tv equal counts | {A→G, A→C} | 1.0 | definition; INV-06 |
| M14 | Ti/Tv 2 Ti 1 Tv | {A→G, C→T, A→C} | 2.0 | definition; INV-06 |
| M15 | Ti/Tv no transversions | {A→G, C→T} | 0 | ASM-03; INV-06 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Mismatched aligned lengths | `CallVariantsFromAlignment("ATGC","ATG")` | throws `ArgumentException` | documented contract |
| S2 | Empty aligned input | `CallVariantsFromAlignment("","")` | empty | contract |
| S3 | Null reference | `CallVariants(null, q)` | throws `ArgumentNullException` | input validation |
| S4 | Null query | `CallVariants(r, null)` | throws `ArgumentNullException` | input validation |
| S5 | Multiple SNPs | `CallVariantsFromAlignment("AAAA","TGTA")` | 3 SNPs at positions 0,1,2 | INV-04 |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Property: positions in bounds & SNP REF≠ALT | random-free constructed alignment | every variant Position ∈ [0,len]; every SNP REF≠ALT | INV-02/03 (O(n×m) invariant) |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/VariantCallerTests.cs` — pre-template fixture covering the whole `VariantCaller` class (SNP/indel/classification/stats/effect/VCF). Many assertions are permissive (`Is.GreaterThan(0)`, `Any(...)`, no exact positions/alleles, no assertion messages). It is the legacy broad fixture, not the canonical per-unit file for VARIANT-CALL-001.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 identical → empty | ⚠ Weak | exists but no message; acceptable shape, rewrite for exactness/message |
| M2 single SNP exact | ⚠ Weak | `CallVariants_SingleSnp_DetectsIt` checks REF/ALT but not Position; no message |
| M3 single SNP via CallVariants | ⚠ Weak | same as M2 |
| M4 insertion exact | ⚠ Weak | `...WithInsertion_DetectsIt` uses `Any(...)`, no position/allele |
| M5 deletion exact | ⚠ Weak | `...WithDeletion_DetectsIt` uses `Any(...)` |
| M6 two-base deletion | ❌ Missing | not covered |
| M7–M12 classification | ⚠ Weak | classifications exist but no assertion messages; case-insensitive (M11) missing |
| M11 case-insensitive classify | ❌ Missing | not covered |
| M13–M15 Ti/Tv | ⚠ Weak | exist but M14 (2:1) missing; no messages |
| M14 Ti/Tv 2:1 | ❌ Missing | not covered |
| S1 length mismatch throws | ✅ Covered | exists; will re-encode in canonical file |
| S2 empty → empty | ✅ Covered | exists |
| S3 null reference | ✅ Covered | exists |
| S4 null query | ✅ Covered | exists |
| S5 multiple SNPs exact | ⚠ Weak | `MultipleSnps` counts only, no positions |
| C1 property | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_CallVariants_Tests.cs` — all VARIANT-CALL-001 cases (M1–M15, S1–S5, C1) with exact evidence-based values and messages, `#region` per method.
- **Remove:** nothing from `VariantCallerTests.cs` — it also exercises out-of-scope methods (`PredictEffect`, `ToVcfLines`, `CalculateStatistics`) owned by other (future) units; leaving it avoids dropping their coverage. The canonical file is the authoritative source for this unit's cases; no duplicate of an in-scope canonical assertion is added to the legacy file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `VariantCaller_CallVariants_Tests.cs` | Canonical VARIANT-CALL-001 fixture | 21 |
| `VariantCallerTests.cs` | Legacy broad fixture (out-of-scope methods retained) | unchanged |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | Rewrote exact+message in canonical file | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote with exact Position/REF/ALT | ✅ Done |
| 3 | M3 | ⚠ Weak | Rewrote with exact values | ✅ Done |
| 4 | M4 | ⚠ Weak | Rewrote exact insertion alleles+position | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote exact deletion alleles+position | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented two-base deletion | ✅ Done |
| 7 | M7–M10 | ⚠ Weak | Rewrote all transition/transversion with messages | ✅ Done |
| 8 | M11 | ❌ Missing | Implemented case-insensitive classify | ✅ Done |
| 9 | M12 | ⚠ Weak | Rewrote non-SNP→Other | ✅ Done |
| 10 | M13 | ⚠ Weak | Rewrote exact 1.0 | ✅ Done |
| 11 | M14 | ❌ Missing | Implemented 2:1 → 2.0 | ✅ Done |
| 12 | M15 | ⚠ Weak | Rewrote exact 0 | ✅ Done |
| 13 | S1–S4 | ✅ Covered | Re-encoded in canonical file | ✅ Done |
| 14 | S5 | ⚠ Weak | Rewrote with exact positions | ✅ Done |
| 15 | C1 | ❌ Missing | Implemented property test | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact, message |
| M2 | ✅ | exact Position/REF/ALT |
| M3 | ✅ | exact via CallVariants |
| M4 | ✅ | exact insertion |
| M5 | ✅ | exact deletion |
| M6 | ✅ | two-base deletion |
| M7 | ✅ | A→G transition |
| M8 | ✅ | C→T transition |
| M9 | ✅ | G→A, T→C transition |
| M10 | ✅ | four transversions |
| M11 | ✅ | case-insensitive |
| M12 | ✅ | non-SNP → Other |
| M13 | ✅ | Ti/Tv 1.0 |
| M14 | ✅ | Ti/Tv 2.0 |
| M15 | ✅ | Ti/Tv 0 |
| S1 | ✅ | ArgumentException |
| S2 | ✅ | empty |
| S3 | ✅ | null ref |
| S4 | ✅ | null query |
| S5 | ✅ | multiple SNPs exact positions |
| C1 | ✅ | property test |

In-scope cases: 21. ✅ count: 21.

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| ASM-01 | In-memory `Variant` uses `"-"` gap sentinel + 0-based Position (VCF padding/1-based applies only to serialized output, out of scope) | M4, M5, M6 |
| ASM-02 | Indels not left-aligned/parsimony-normalized; position asserted only on unambiguous alignments | M4, M5, M6, C1 |
| ASM-03 | Ti/Tv with zero transversions returns 0 (undefined case) | M15 |

---

## 7. Open Questions / Decisions

1. Registry method name `ClassifyVariant` does not exist; behavioral contract maps to `ClassifyMutation`. Decision: keep `ClassifyMutation` (source-faithful; renaming breaks sibling tests), document the mapping. The Registry is workflow-control only.
2. Left-alignment/parsimony normalization (Tan 2015) is out of scope for this detection unit; a future normalization unit could add it. Recorded as ASM-02.
