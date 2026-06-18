# Indel Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Variants |
| Test Unit ID | VARIANT-INDEL-001 |
| Related Projects | Seqeron.Genomics.Annotation, Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Indel detection identifies insertions and deletions between a reference DNA sequence and a query sequence. The reference and query are aligned with Needleman–Wunsch global alignment, then each alignment column where exactly one side carries a gap is reported as a length-changing event: a reference-side gap is an insertion in the query, a query-side gap is a deletion in the query [1]. The repository exposes `VariantCaller.FindInsertions` and `VariantCaller.FindDeletions`, which filter the column-based variant caller to one indel class. Detection is specification-driven (VCF variant classes [1]) over a deterministic optimal alignment; it reports per-base indel columns rather than a single normalized multi-base allele.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The Variant Call Format (VCF) distinguishes three simple variant classes: substitutions (SNPs), insertions, and deletions [1]. An insertion adds bases present in the query but absent from the reference; a deletion removes bases present in the reference but absent from the query. In VCF these are length-changing events, in contrast to a SNP which preserves length [1].

### 2.2 Core Model

A pairwise global alignment places the reference and query into equal-length gapped strings. For each column `i`:

- reference gap (`-`) opposite a query base ⇒ **insertion** of that base in the query;
- query gap (`-`) opposite a reference base ⇒ **deletion** of that reference base in the query;
- two non-gap bases that differ ⇒ SNP (out of scope for this unit);
- equal bases ⇒ match [1].

In VCF the directional length convention is: an insertion has ALT longer than REF ("A single base insertion of A after position 3 becomes REF=C, ALT=CA") and a deletion has REF longer than ALT ("A single base deletion of C at position 3 becomes REF=TC, ALT=T") [1]. For a *serialized* VCF record a pure indel must carry a left anchor (padding) base "which must be reflected in the POS field" [1]; the in-memory model in this repository instead marks the absent allele with the `-` gap sentinel (see ASM-01).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | In-memory indels use the `-` gap sentinel for the absent allele rather than the VCF left-anchored padded allele. | The in-memory `Variant` is not a literal VCF record; padded, 1-based VCF is produced only by `ToVcfLines`. Detection (counts/types) is unaffected. |
| ASM-02 | The reported indel position is the alignment column produced by global alignment; it is not left-aligned / parsimony-normalized per Tan et al. (2015) [2]. | In a tandem-repeat / low-complexity region the same indel can be placed at several positions; the reported position may not be the canonical leftmost one [2]. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Identical sequences produce no insertions and no deletions. | An indel is a length-changing difference; identical input has none [1]. |
| INV-02 | Every variant from `FindInsertions` has `Type == Insertion`; every variant from `FindDeletions` has `Type == Deletion`. | The methods filter `CallVariants` by `VariantType` (insertion/deletion are distinct classes [1]). |
| INV-03 | An insertion column has `ReferenceAllele == "-"` and a one-base non-gap `AlternateAllele` (ALT longer than REF). | Insertion ⇒ ALT longer than REF [1]; minimal-representation directional length [3]. |
| INV-04 | A deletion column has `AlternateAllele == "-"` and a one-base non-gap `ReferenceAllele` (REF longer than ALT). | Deletion ⇒ REF longer than ALT [1], [3]. |
| INV-05 | A contiguous block of `k` inserted (resp. deleted) bases in an otherwise unique alignment yields exactly `k` insertion (resp. deletion) columns. | Each extra/absent base is one indel event [1]. |
| INV-06 | Every reported indel `Position` lies in `[0, reference.Length]`. | Reference-coordinate bookkeeping advances `refPos` only on reference-consuming columns. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `reference` | `DnaSequence` | required | Reference DNA sequence. | Non-null; A/C/G/T (case-normalized to uppercase by the aligner). |
| `query` | `DnaSequence` | required | Query DNA sequence to compare. | Non-null. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IEnumerable<Variant>` | Insertions only (`FindInsertions`) or deletions only (`FindDeletions`). |
| `Variant.Position` | `int` | 0-based reference position of the event (ASM-02; 0-based per the in-memory model). |
| `Variant.ReferenceAllele` | `string` | Deleted base, or `"-"` for an insertion (ASM-01). |
| `Variant.AlternateAllele` | `string` | Inserted base, or `"-"` for a deletion (ASM-01). |
| `Variant.Type` | `VariantType` | `Insertion` or `Deletion`. |
| `Variant.QueryPosition` | `int` | 0-based query position of the event. |

### 3.3 Preconditions and Validation

A null `reference` or `query` throws `ArgumentNullException` (propagated from `CallVariants`). Empty sequences yield no variants. Bases are case-normalized to uppercase by `SequenceAligner.GlobalAlign`. Positions are 0-based (ASM-02); 1-based VCF POS is produced only by `ToVcfLines`.

## 4. Algorithm

### 4.1 High-Level Steps

1. Globally align `reference` and `query` (Needleman–Wunsch) to obtain two equal-length gapped strings.
2. Walk the aligned columns, maintaining `refPos` (reference coordinate) and `queryPos` (query coordinate).
3. For a reference-gap column emit an `Insertion` (`ReferenceAllele = "-"`, `AlternateAllele = query base`) at `refPos`; advance `queryPos`.
4. For a query-gap column emit a `Deletion` (`ReferenceAllele = ref base`, `AlternateAllele = "-"`) at `refPos`; advance `refPos`.
5. `FindInsertions` keeps `Type == Insertion`; `FindDeletions` keeps `Type == Deletion` [1].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Global alignment uses the `SimpleDna` scoring matrix (match +1, mismatch −1, linear gap −1). On inputs without internal repeats this yields a unique optimal alignment, so indel positions are deterministic; in repeat regions the placement is not canonical (ASM-02, [2]).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindInsertions` / `FindDeletions` | O(n × m) | O(n × m) | Dominated by Needleman–Wunsch DP over reference length `n` and query length `m`; the column walk and filter are O(n + m). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [VariantCaller.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs)

- `VariantCaller.FindInsertions(reference, query)`: `CallVariants` filtered to `VariantType.Insertion`.
- `VariantCaller.FindDeletions(reference, query)`: `CallVariants` filtered to `VariantType.Deletion`.
- `VariantCaller.FindIndels(reference, query)`: union of insertions and deletions (delegate).
- `VariantCaller.CallVariantsFromAlignment(...)`: the shared column walk that classifies each gap column (validated under VARIANT-CALL-001).

### 5.2 Current Behavior

Detection is alignment-based: `CallVariants` calls `SequenceAligner.GlobalAlign` and then `CallVariantsFromAlignment`. A multi-base indel is reported as several consecutive single-base indel columns, not one merged allele. The absent allele uses the `-` sentinel (ASM-01). Indels are not normalized (ASM-02).

**Search-reuse decision (suffix tree):** not applicable. Indel detection is scoring-based optimal alignment (dynamic programming with gaps), not exact substring matching; the repository suffix tree (exact-occurrence enumeration) does not fit this problem, so it is not used. The shared `SequenceAligner` is reused instead.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Insertion vs deletion classification by gap side, with the VCF directional length convention (insertion ⇒ ALT longer than REF; deletion ⇒ REF longer than ALT) [1].
- Insertion/deletion are distinct variant classes, separated from SNPs [1].

**Intentionally simplified:**

- In-memory allele representation uses the `-` gap sentinel instead of the VCF left-anchored padding base; **consequence:** an in-memory indel `Variant` is not a literal VCF record (the padded, 1-based form is produced by `ToVcfLines`) (ASM-01).

**Not implemented:**

- Left-alignment / parsimony normalization of indels (Tan et al. 2015 Algorithm 1) [2]; **users should rely on:** an external normalizer (e.g. `vt normalize`, `bcftools norm`) when canonical placement in repeat regions is required (ASM-02).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Gap-sentinel in-memory alleles | Assumption | In-memory `Variant` ≠ literal VCF record | accepted | ASM-01; VCF padding produced by `ToVcfLines` |
| 2 | Indels not normalized | Assumption | Non-canonical position in repeats | accepted | ASM-02; cite Tan et al. 2015 [2] |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical sequences | No insertions, no deletions | No length-changing difference [1] (INV-01). |
| Null reference / null query | `ArgumentNullException` | Input validation propagated from `CallVariants`. |
| Substitution-only input (equal length) | No insertions / deletions | A substitution is length-preserving, not an indel [1]. |
| Multi-base indel | One indel column per base, consecutive | Each extra/absent base is one event [1] (INV-05). |

### 6.2 Limitations

Indel positions are not normalized (left-aligned / parsimonious); in tandem-repeat or low-complexity regions the reported position is one of several valid placements [2]. Detection depends on the global-alignment placement, which is deterministic but unique only when no internal repeat permits an equal-score left shift.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reference = new DnaSequence("ATGCAT");
var query = new DnaSequence("ATGTCAT"); // a T inserted after index 2

var insertions = VariantCaller.FindInsertions(reference, query).ToList();
// insertions[0]: Type = Insertion, ReferenceAllele = "-", AlternateAllele = "T", Position = 3
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [VariantCaller_FindIndels_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/VariantCaller_FindIndels_Tests.cs) — covers INV-01..INV-06
- Evidence: [VARIANT-INDEL-001-Evidence.md](../../../docs/Evidence/VARIANT-INDEL-001-Evidence.md)
- Related algorithms: [Variant_Detection](../Variants/Variant_Detection.md), [SNP_Detection](../Variants/SNP_Detection.md)

### 7.4 Change History

| Date | Version | Changes |
|------|---------|---------|
| 2026-06-13 | 1.0 | Initial documentation (VARIANT-INDEL-001). |

## 8. References

1. The SAM/BCF/VCF working group. 2024. The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex
2. Tan A, Abecasis GR, Kang HM. 2015. Unified representation of genetic variants. *Bioinformatics* 31(13):2202–2204. https://doi.org/10.1093/bioinformatics/btv112
3. Minikel E. minimal_representation: find the minimal representation of a variant (implements Tan et al. 2015 normalization). https://raw.githubusercontent.com/ericminikel/minimal_representation/master/normalize.py
