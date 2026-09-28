# Variant Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Variants |
| Test Unit ID | VARIANT-CALL-001 |
| Related Projects | Seqeron.Genomics.Annotation, Seqeron.Genomics.Alignment |
| Implementation Status | Simplified (linear-gap default scoring; see §5.3) |
| Last Reviewed | 2026-09-28 (B21 review) |

## 1. Overview

Variant detection identifies the differences (SNPs, insertions, deletions) of a query DNA sequence relative to a reference. The library performs a pairwise global alignment of reference and query, then scans the resulting gapped alignment column by column: a substituted column is a single-nucleotide polymorphism (SNP), a reference-side gap is an insertion in the query, and a query-side gap is a deletion. SNPs are further classified as transitions or transversions. The procedure is deterministic and specification-driven (the variant classes follow the Variant Call Format [1]); Indels found by `CallVariants` are left-aligned in the sense of Tan et al. [4] (verified against `bcftools norm`, see §5.3); multi-base indels are reported one event per gap column.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A *variant* is a difference between a sample (query) sequence and a reference. The Variant Call Format (VCF) stores "DNA polymorphism data such as SNPs, insertions, deletions and structural variants" [1]. The three small-variant classes handled here are: substitution (SNP), insertion, and deletion [1][2].

### 2.2 Core Model

Given a global alignment of reference `R` and query `Q` (both strings over `{A,C,G,T,…}` plus the gap symbol `-`), scan positions `i = 0 … L-1` of the aligned strings:

- `R[i] = -`, `Q[i] ≠ -` → **insertion** of `Q[i]` in the query at reference coordinate `refPos`.
- `R[i] ≠ -`, `Q[i] = -` → **deletion** of `R[i]` from the query at reference coordinate `refPos`.
- `R[i] ≠ -`, `Q[i] ≠ -`, `R[i] ≠ Q[i]` → **SNP** `R[i]→Q[i]` at `refPos`.
- otherwise → match (no variant).

Reference and query coordinates `refPos`, `queryPos` advance by one for each non-gap base consumed on their respective side.

**SNP sub-classification** [5][6]:
- **Transition** — purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T) [5].
- **Transversion** — purine↔pyrimidine: A↔C, A↔T, G↔C, G↔T [6].
- The **transition/transversion (Ti/Tv) ratio** = (#transitions)/(#transversions). Transitions exceed transversions in real genomes (transitional bias) [3], so a meaningful Ti/Tv is typically > 0.5.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Identical reference and query yield zero variants. | No alignment column differs [1]. |
| INV-02 | Every emitted SNP has distinct single-base REF and ALT. | A SNP column is, by definition, a mismatch of two bases [2]. |
| INV-03 | Every emitted variant's 0-based `Position` lies in `[0, reference.Length]`. | `Position` is the reference coordinate advanced only by consumed reference bases. |
| INV-04 | Ref-gap → Insertion; query-gap → Deletion; mismatch → SNP. | Column-classification rule of §2.2 [1]. |
| INV-05 | `ClassifyMutation` = Transition iff {ref,alt}⊆{A,G} or ⊆{C,T}; Transversion for purine↔pyrimidine over A/C/G/T; Other for non-SNPs and for SNPs involving a non-ACGT base such as N (bcftools stats counts A>N as neither ts nor tv). | Definitions of transition/transversion [5][6]. |
| INV-06 | `CalculateTiTvRatio` = #Ti / #Tv over SNPs, or 0 when #Tv = 0. | Ratio definition; undefined denominator mapped to 0 (see 5.4). |
| INV-07 | `CallVariants` reports every indel at its leftmost equivalent column (left-aligned [4]). | NW traceback prefers the diagonal move, pushing gap runs left; cross-checked vs `bcftools norm` (§5.3). |
| INV-08 | Base comparison in the column scan is case-insensitive (a case-only difference is not a variant). | VCF REF/ALT bases are case insensitive [2]; `bcftools norm` rejects REF=G/ALT=g as duplicate alleles. |

### 2.5 Comparison with Related Methods

| Aspect | This (alignment column scan) | Direct positional scan (`FindSnpsDirect`) |
|--------|------------------------------|--------------------------------------------|
| Indels detected | Yes (via gaps) | No (length-equal comparison only) |
| Cost | O(n×m) alignment + O(L) scan | O(min(n,m)) |
| Indel normalization | Left-aligned via aligner traceback (`CallVariants`); literal columns (`CallVariantsFromAlignment`) [4] | N/A |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reference | DnaSequence | required | Reference sequence | non-null |
| query | DnaSequence | required | Query sequence to compare | non-null |
| alignedReference | string | required | Gapped reference (`-` = gap) | equal length to alignedQuery |
| alignedQuery | string | required | Gapped query (`-` = gap) | equal length to alignedReference |
| variant | Variant | required | Variant to classify | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Variant.Position | int | 0-based reference coordinate of the variant |
| Variant.ReferenceAllele | string | Reference base, or `"-"` for an insertion |
| Variant.AlternateAllele | string | Query base, or `"-"` for a deletion |
| Variant.Type | VariantType | SNP / Insertion / Deletion |
| Variant.QueryPosition | int | 0-based query coordinate |
| ClassifyMutation → MutationType | enum | Transition / Transversion / Other |
| CalculateTiTvRatio → double | double | #Ti / #Tv, or 0 when #Tv = 0 |

### 3.3 Preconditions and Validation

- `CallVariants` / `CalculateStatistics`: null reference or query → `ArgumentNullException`.
- `CallVariantsFromAlignment`: null/empty input → empty result; unequal aligned lengths → `ArgumentException`.
- Indexing is 0-based for the in-memory `Variant`. Serialized VCF (`ToVcfLines`, out of scope for this unit) shifts to 1-based and applies VCF padding rules [2].
- Classification is case-insensitive (REF/ALT bases are case-insensitive in VCF [2]).

## 4. Algorithm

### 4.1 High-Level Steps

1. `CallVariants`: globally align reference and query (`SequenceAligner.GlobalAlign`).
2. Scan the aligned columns; emit SNP / insertion / deletion per §2.2, advancing reference and query coordinates.
3. (Classification) `ClassifyMutation`: for a SNP, compare purine/pyrimidine membership of REF vs ALT.
4. (Aggregation) `CalculateTiTvRatio`: count transitions and transversions over the SNPs.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Transition/transversion lookup [5][6]:

| REF→ALT examples | Class |
|------------------|-------|
| A↔G, C↔T | Transition |
| A↔C, A↔T, G↔C, G↔T | Transversion |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CallVariants` | O(n×m) | O(n×m) | Needleman-Wunsch global alignment dominates |
| `CallVariantsFromAlignment` | O(L) | O(1) extra | L = aligned length; single pass |
| `ClassifyMutation` | O(1) | O(1) | constant lookup |
| `CalculateTiTvRatio` | O(k) | O(1) | k = number of variants |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [VariantCaller.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs)

- `VariantCaller.CallVariants(reference, query)`: aligns then delegates to the column scan.
- `VariantCaller.CallVariantsFromAlignment(alignedReference, alignedQuery)`: column scan producing `Variant` records.
- `VariantCaller.ClassifyMutation(variant)`: transition/transversion/Other classification.
- `VariantCaller.CalculateTiTvRatio(variants)`: Ti/Tv aggregation.

### 5.2 Current Behavior

The in-memory `Variant` uses a `"-"` gap sentinel for the absent allele of an indel and a 0-based `Position` — this is the internal model, distinct from serialized VCF (which requires a padding base and 1-based POS [2]). `CallVariants` indel positions are the columns produced by `SequenceAligner.GlobalAlign`; its traceback prefers the diagonal move, so every gap run sits at its leftmost equivalent column, i.e. the indel is left-aligned [4]. `CallVariantsFromAlignment` reports a caller-supplied alignment literally (no re-positioning). Column comparison is case-insensitive [2].

**Search reuse (suffix tree):** Not applicable. Variant detection is a scoring-based pairwise *alignment* (edit-distance / Needleman-Wunsch), not exact substring search; the repository suffix tree (`Contains`/`FindAllOccurrences`) addresses exact-match occurrence enumeration and does not fit edit-distance variant calling. No suffix-tree usage.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- The three small-variant classes SNP / insertion / deletion [1][2].
- Transition = A↔G / C↔T; transversion = purine↔pyrimidine [5][6].
- Ti/Tv ratio = #transitions / #transversions [3].
- Case-insensitive base handling for classification and for the column scan [2] (B21 fix, 2026-09: `CallVariantsFromAlignment` previously reported `a`/`A` columns as SNPs).
- Left-alignment [4] of `CallVariants` indels. Reference cross-check: `bcftools norm -f` (htslib/bcftools via pysam 0.24.1) left 1000/1000 randomly generated single-indel records (repeat-rich AC/ACGT sequences, indel length 1–3) unchanged; worked cases `ACGTTTTACG→ACGTTTACG` (norm POS 3 GT>G → Position 3), `CAGAGAGT→CAGAGT` (POS 1 CAG>C → Positions 1,2), `GCACAT→GCACACAT` (POS 1 G>GCA → Position 1 ×2), `ATTG→ATTTG` (POS 1 A>AT → Position 1).

**Intentionally simplified:**

- Indel representation: a multi-base indel is reported as one event per gap column rather than one VCF record; `CallVariantsFromAlignment` does not re-position gaps of a caller-supplied alignment.
- Default scoring: `SequenceAligner.SimpleDna` (linear gap −1 = mismatch −1); **consequence:** an adjacent two-base swap in equal-length sequences (`GGACGG→GGCAGG`) is reported as insertion + deletion instead of two substitutions. Read aligners penalise gaps above mismatches with affine costs (BWA-MEM defaults: mismatch 4, gap open 6, extend 1 — `bwa.1`); the library has no public affine (Gotoh) pairwise aligner yet (cross-batch request, B21 2026-09).
- In-memory allele model uses a `"-"` sentinel rather than the VCF padded-allele representation; **consequence:** users serializing to VCF must add the padding base (handled by `ToVcfLines`, a separate unit).

**Not implemented:**

- MNP / complex-substitution detection; **users should rely on:** no current alternative in this class (out of scope for VARIANT-CALL-001).
- Normalization of caller-supplied alignments (`CallVariantsFromAlignment`) and merging of per-column events into single VCF records; **users should rely on:** `CallVariants` or external tools (e.g. `bcftools norm`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Ti/Tv with #Tv = 0 returns 0 | Assumption | undefined ratio mapped to a sentinel rather than +∞/throw | accepted | INV-06; no source mandates a sentinel |
| 2 | `CallVariantsFromAlignment` reports gaps literally | Design | a non-left-aligned input alignment yields non-left-aligned positions | accepted | `CallVariants` output is left-aligned (INV-07) |
| 3 | Linear gap cost = mismatch cost (SimpleDna) | Limitation | adjacent swaps become ins+del | declared | needs affine global aligner (cross-batch) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical sequences | empty result | INV-01 [1] |
| Empty aligned input | empty result | contract |
| Unequal aligned lengths | `ArgumentException` | columns must align |
| Null reference/query | `ArgumentNullException` | input validation |
| Non-SNP passed to `ClassifyMutation` | `Other` | classification defined only for SNPs |
| Lowercase SNP bases | classified identically | case-insensitive [2] |
| Case-only column difference (`a`/`A`) | no variant | case-insensitive [2] |

### 6.2 Limitations

Detection quality is bounded by the upstream alignment and its linear scoring (gap = mismatch cost, §5.3). Indel placement from `CallVariants` is left-aligned [4]; caller-supplied alignments are not re-positioned. Multi-nucleotide and complex variants are reported as adjacent SNPs/indels rather than a single combined event. The aligner is O(n×m) in time and space, limiting practical sequence length.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reference = new DnaSequence("ATGC");
var query = new DnaSequence("ATTC");
var variants = VariantCaller.CallVariants(reference, query).ToList();
// 1 variant: { Position = 2, ReferenceAllele = "G", AlternateAllele = "T", Type = SNP }
var cls = VariantCaller.ClassifyMutation(variants[0]); // Transversion (G→T)
```

**Numerical / biological walk-through:**

VCFv4.3 §1.1 microsatellite `GTC → G` is "a deletion of 2 bases (TC)" [2]. Aligning `GTCAA` (ref) to `G--AA` (query) yields two deletion columns at reference positions 1 (`T`) and 2 (`C`) — matching the spec's TC deletion (here split per base in the in-memory model).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [VariantCaller_CallVariants_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/VariantCaller_CallVariants_Tests.cs) — covers INV-01…INV-06.
- Evidence: [VARIANT-CALL-001-Evidence.md](../../../docs/Evidence/VARIANT-CALL-001-Evidence.md)

## 8. References

1. Danecek P, Auton A, Abecasis G, et al. 2011. The variant call format and VCFtools. *Bioinformatics* 27(15):2156–2158. https://doi.org/10.1093/bioinformatics/btr330
2. The SAM/BCF/VCF working group. 2024. The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://samtools.github.io/hts-specs/VCFv4.3.pdf
3. Collins DW, Jukes TH. 1994. Rates of transition and transversion in coding sequences since the human-rodent divergence. *Genomics* 20(3):386–396. https://doi.org/10.1006/geno.1994.1192
4. Tan A, Abecasis GR, Kang HM. 2015. Unified representation of genetic variants. *Bioinformatics* 31(13):2202–2204. https://doi.org/10.1093/bioinformatics/btv112
5. Wikipedia. Transition (genetics). https://en.wikipedia.org/wiki/Transition_(genetics) (accessed 2026-06-13; primary: Collins & Jukes 1994)
6. Wikipedia. Transversion. https://en.wikipedia.org/wiki/Transversion (accessed 2026-06-13; primary: Futuyma DJ, *Evolution*, 3rd ed., 2013)
