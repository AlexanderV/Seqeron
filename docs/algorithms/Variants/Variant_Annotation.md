# Variant Annotation (Functional Impact / Sequence Ontology Consequences)

| Field | Value |
|-------|-------|
| Algorithm Group | Variants |
| Test Unit ID | VARIANT-ANNOT-001 |
| Related Projects | Seqeron.Genomics.Annotation, Seqeron.Genomics.Core |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Variant annotation assigns to each genetic variant the functional consequence it has on overlapping transcripts, expressed as a Sequence Ontology (SO) consequence term (e.g. `missense_variant`, `stop_gained`, `frameshift_variant`) and a coarse IMPACT rating (HIGH / MODERATE / LOW / MODIFIER). This unit implements the consequence-determination core of the Ensembl Variant Effect Predictor (VEP) [1]: it locates the variant relative to a transcript's exons, introns, UTRs and coding sequence, and — for coding substitutions — translates the affected reference and alternate codons with the NCBI Standard genetic code [4] to distinguish synonymous, missense, stop-gained, stop-lost and start-lost changes. It is specification-driven (the SO terms, IMPACT ratings and severity ranks are taken verbatim from the VEP reference implementation [2]) and deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A variant overlaps zero or more transcripts. VEP describes each variant–transcript effect with a standardized SO consequence term and reports, per variant, the most severe consequence across all overlaps [1]. Each consequence term carries a fixed IMPACT class and an integer severity rank (1 = most severe) defined in the VEP `Constants.pm` `OverlapConsequence` table [2].

### 2.2 Core Model

For a coding single-/multi-nucleotide substitution, classification follows the VEP `VariationEffect.pm` predicates over the translated reference peptide `ref_pep` and alternate peptide `alt_pep` [3]:

- `stop_gained` ⇔ `alt_pep` contains `*` and `ref_pep` does not [3].
- `stop_lost` ⇔ `ref_pep` is `*` and `alt_pep` is not [3].
- `start_lost` ⇔ the substitution disrupts the canonical start codon (precedence over missense) [3].
- `synonymous_variant` ⇔ `alt_pep == ref_pep`, neither peptide unknown (`X`) [3].
- `missense_variant` ⇔ `ref_pep ≠ alt_pep`, equal length, and not start/stop_lost/stop_gained [3].

For a coding indel, classification is length-based [3]:

- `frameshift_variant` ⇔ `|len(alt) − len(ref)| mod 3 ≠ 0` [3].
- otherwise `inframe_insertion` (alt longer) or `inframe_deletion` (alt shorter) [3].

Peptides are obtained by translating codons with the NCBI Standard code (`transl_table 1`): AAs `FFLLSSSS…GGGG`, stops TAA/TAG/TGA, start ATG [4].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Codon translation uses the NCBI Standard code (table 1) | Organisms with alternate codes (e.g. vertebrate mitochondria) misclassify stop/sense codons |
| ASM-02 | The codon affected by a substitution is contiguous on the forward strand reference window | A codon split across an intron boundary yields an incorrect codon and consequence |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | IMPACT(term) equals the `impact` field of that term in VEP `Constants.pm` | Mapping copied verbatim from [2] |
| INV-02 | A coding SNV is `synonymous_variant` ⇔ `ref_pep == alt_pep` | VEP predicate [3] over Standard code [4] |
| INV-03 | A coding indel is `frameshift_variant` ⇔ length difference not divisible by 3 | VEP frameshift predicate [3] |
| INV-04 | `Annotate` returns the lowest-rank (most severe) consequence per variant | VEP most-severe reporting [1] + ranks [2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `variant` | `Variant` | required | Variant; 1-based genomic `Position`, forward-strand `Reference`/`Alternate` | non-empty alleles |
| `transcript` | `Transcript` | required | Transcript model (exons, coding exons, CDS bounds) | CDS bounds for coding calls |
| `referenceSequence` | `string` | required | Forward-strand reference window covering the variant codon context | non-null/non-empty |
| `sequenceStart` | `int` | `1` | 1-based genomic coordinate of `referenceSequence[0]` | ≥ 1 |
| `variants` | `IEnumerable<Variant>` | required (`Annotate`) | Variants to annotate | non-null |
| `annotations` | `IEnumerable<Transcript>` | required (`Annotate`) | Transcript models | non-null |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `FunctionalImpact.Consequence` | `ConsequenceType` | SO consequence term |
| `FunctionalImpact.Impact` | `ImpactLevel` | HIGH/MODERATE/LOW/MODIFIER per [2] |
| `FunctionalImpact.CodonChange` | `string?` | HGVS-style `c.` change for coding substitutions |
| `FunctionalImpact.AminoAcidChange` | `string?` | HGVS-style `p.` change for coding substitutions |
| `Annotate → VariantAnnotation` | record | One most-severe annotation per variant |

### 3.3 Preconditions and Validation

Positions are 1-based; coordinates inclusive. Alleles are upper-cased; T/U handled by the translator. `PredictFunctionalImpact` throws `ArgumentException` on null/empty `referenceSequence`. `Annotate` throws `ArgumentNullException` on null `variants`/`annotations`. Codons containing IUPAC ambiguity translate to `X` and are excluded from `synonymous_variant` (reported as `coding_sequence_variant`).

## 4. Algorithm

### 4.1 High-Level Steps

1. Locate the variant relative to the transcript (upstream/downstream/intron/UTR/splice/coding) via `DetermineConsequence`.
2. If coding and a substitution: compute CDS position, extract the affected reference codon, build the alternate codon, translate both with the Standard code.
3. Apply the VEP peptide predicates to assign `stop_gained`/`stop_lost`/`start_lost`/`synonymous`/`missense` [3].
4. Map the consequence term to its IMPACT via `GetImpactLevel` [2].
5. For `Annotate`, repeat over all overlapping transcripts and return the lowest-rank consequence [1][2].

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Consequence severity ranks** (`ConsequenceRank`): copied verbatim from VEP `Constants.pm` `OverlapConsequence` `rank` fields (release/110) [2].
- **NCBI Standard genetic code (table 1)**: codon→amino-acid translation via `GeneticCode.Standard` [4].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `PredictFunctionalImpact` | O(E) | O(1) | E = exons (CDS position scan); codon translation is O(1) |
| `Annotate` | O(v × g) | O(g) | v variants, g transcripts |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [VariantAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantAnnotator.cs)

- `VariantAnnotator.PredictFunctionalImpact(Variant, Transcript, string, int)`: codon-translation consequence + IMPACT + HGVS-style change.
- `VariantAnnotator.Annotate(IEnumerable<Variant>, IEnumerable<Transcript>, string?, int)`: most-severe annotation per variant.
- `VariantAnnotator.GetImpactLevel(ConsequenceType)`: term → IMPACT mapping.
- `VariantAnnotator.GetConsequenceRank(ConsequenceType)`: term → severity rank.

### 5.2 Current Behavior

Coding-substitution refinement (synonymous/missense/stop/start) requires a reference sequence window; without one, `AnnotateVariant`'s coordinate-only routing returns the coarse coding term. Codon extraction assumes the codon is contiguous on the forward strand (ASM-02), satisfied by single-exon coding tests and typical exonic SNVs. No suffix-tree reuse: consequence determination is coordinate/codon arithmetic, not substring search, so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- SO consequence terms, IMPACT classes and severity ranks copied from VEP `Constants.pm` [2].
- Coding-substitution predicates (synonymous/missense/stop_gained/stop_lost/start_lost) per VEP `VariationEffect.pm` [3].
- Frameshift vs inframe insertion/deletion length rule [3].
- Codon translation via NCBI Standard code (table 1) [4].
- Most-severe-consequence selection per variant [1].

**Intentionally simplified:**

- Forward-strand codon extraction only; **consequence:** minus-strand coding SNV consequences are not refined by translation in this unit (ASM-02).
- Single genetic code (table 1); **consequence:** non-standard organism codes are not modelled (ASM-01).

**Not implemented:**

- NMD, downstream re-initiation, splice_donor_5th_base and polypyrimidine sub-terms; **users should rely on:** Ensembl VEP for full transcript-level effects.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Pre-existing `PredictPathogenicity` / SIFT / PolyPhen / conservation methods carry invented constants | Deviation | Those outputs are not source-traceable | open | Out of scope for VARIANT-ANNOT-001; flagged for a future Test Unit |
| 2 | Forward-strand codon extraction | Assumption (ASM-02) | Minus-strand coding SNV refinement skipped | accepted | Documented limitation |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Premature stop substitution | `StopGained`, HIGH | `alt_pep` contains `*` not in `ref_pep` [3] |
| Stop→sense substitution | `StopLost`, HIGH | `ref_pep` is `*`, `alt_pep` is not [3] |
| Same-amino-acid substitution | `SynonymousVariant`, LOW | `ref_pep == alt_pep` [3] |
| Ambiguous codon (contains N) | not `SynonymousVariant` (→ `CodingSequenceVariant`) | `X` peptide excluded from synonymous [3] |
| Coding indel Δ not multiple of 3 | `FrameshiftVariant`, HIGH | length rule [3] |
| Empty `referenceSequence` | `ArgumentException` | contract |
| Null `variants`/`annotations` | `ArgumentNullException` | contract |

### 6.2 Limitations

Minus-strand coding SNVs are not codon-refined (ASM-02); only the Standard genetic code is supported (ASM-01); transcript-level effects such as NMD and alternative re-initiation are not modelled. Codons split across intron boundaries are not handled.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var variant = new VariantAnnotator.Variant("chr1", 101, "A", "T", VariantAnnotator.VariantType.SNV);
var fi = VariantAnnotator.PredictFunctionalImpact(variant, codingTranscript, referenceWindow, sequenceStart: 100);
// codon GAA (Glu) -> GTA (Val): fi.Consequence == MissenseVariant, fi.Impact == Moderate
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [VariantAnnotator_FunctionalImpact_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/VariantAnnotator_FunctionalImpact_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [VARIANT-ANNOT-001-Evidence.md](../../../docs/Evidence/VARIANT-ANNOT-001-Evidence.md)

## 8. References

1. McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, Flicek P, Cunningham F. 2016. The Ensembl Variant Effect Predictor. Genome Biology 17:122. https://doi.org/10.1186/s13059-016-0974-4
2. Ensembl. ensembl-variation (release/110). Bio::EnsEMBL::Variation::Utils::Constants — OverlapConsequence impact/rank table. https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/Constants.pm
3. Ensembl. ensembl-variation (release/110). Bio::EnsEMBL::Variation::Utils::VariationEffect — consequence predicates. https://raw.githubusercontent.com/Ensembl/ensembl-variation/release/110/modules/Bio/EnsEMBL/Variation/Utils/VariationEffect.pm
4. NCBI. The Genetic Codes (gc.prt), Standard code (transl_table 1). https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
5. Eilbeck K, Lewis SE, Mungall CJ, et al. 2005. The Sequence Ontology: a tool for the unification of genome annotations. Genome Biology 6:R44. https://doi.org/10.1186/gb-2005-6-5-r44
