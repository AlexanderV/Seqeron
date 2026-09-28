# SNP Detection

| Field | Value |
|-------|-------|
| Algorithm Group | Variants |
| Test Unit ID | VARIANT-SNP-001 |
| Related Projects | Seqeron.Genomics.Annotation, Seqeron.Genomics.Core, Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

SNP detection identifies single-nucleotide substitutions between a reference and a query DNA sequence. A SNP (single-nucleotide polymorphism) is a position where one base in the reference is replaced by a different single base in the query [1]. Two entry points are provided: `FindSnpsDirect`, an exact positional comparison (the enumeration of Hamming mismatches between two aligned/equal-length strings) [3]; and `FindSnps`, which first globally aligns the inputs and then filters the resulting variants to substitutions only, so that indels are excluded [1]. The detection is exact and deterministic, not heuristic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A variant is a difference between a query sequence and a reference. The Variant Call Format (VCF) distinguishes SNPs, insertions, deletions and structural variants by the properties of the REF and ALT alleles [1]. A SNP is the substitution case: a single-base REF allele and a single-base ALT allele that differ at one position. When the two sequences are already aligned and of equal length, the set of substitutions is exactly the set of positions at which the sequences differ — the Hamming-distance mismatch set [3].

### 2.2 Core Model

For two strings `r` (reference) and `q` (query):

- **Positional model (`FindSnpsDirect`):** the inputs are pre-aligned and of equal length `n`; for each column `i` in `[0, n)` where neither side is a gap and `upper(r[i]) != upper(q[i])` (VCF REF/ALT bases are case-insensitive [1]), a SNP is emitted with reference allele `r[i]` and alternate allele `q[i]` at the ungapped reference/query coordinates. For gap-free inputs the positions are the column indices and their number is the Hamming distance, defined as "the number of positions that two codewords of the same length differ" [3]. Unequal lengths are rejected (`ArgumentException`), because the Hamming mismatch set is undefined for them [3]. Implemented by delegating to `CallVariantsFromAlignment` and keeping SNP columns.
- **Alignment model (`FindSnps`):** the inputs are globally aligned; in each alignment column where both sides hold a non-gap base and the bases differ, a SNP is emitted; insertion and deletion columns are excluded by filtering on `VariantType.SNP` [1].

A position where `r[i] == q[i]` is a match, not a variant: a SNP is by definition a substitution [1].

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | `FindSnpsDirect` treats the inputs as positionally corresponding (already aligned, no indels). | If the inputs are mutually shifted by an indel, downstream mismatches are reported as SNPs rather than as the indel that caused them — use `FindSnps` (which aligns first) for un-aligned inputs. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Identical (case-insensitively) equal-length sequences yield zero SNPs. | Hamming distance of equal strings is 0 [3]. |
| INV-02 | Every emitted variant has `Type == VariantType.SNP`. | `FindSnpsDirect` only constructs `SNP` variants; `FindSnps` filters to `SNP` [1]. |
| INV-03 | Every emitted SNP has `ReferenceAllele != AlternateAllele`. | A SNP is a substitution; equal bases are skipped [1]. |
| INV-04 | For gap-free inputs `FindSnpsDirect` reports a SNP at each 0-based case-insensitive mismatch index `i` with `Position == QueryPosition == i`, `ReferenceAllele == r[i]`, `AlternateAllele == q[i]`; gap columns are never SNPs. | Direct positional comparison [1][3]. |
| INV-05 | For two equal-length sequences, the SNP count from `FindSnpsDirect` equals their Hamming distance. | The mismatch set is the Hamming mismatch set [3]. |
| INV-06 | `FindSnpsDirect` throws `ArgumentException` when both inputs are non-empty and their lengths differ. | Hamming distance is defined for equal-length strings only [3]; SciPy `hamming` and scikit-bio `mismatches` also reject unequal lengths. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| reference (`FindSnps`) | `DnaSequence` | required | Reference sequence. | Non-null. |
| query (`FindSnps`) | `DnaSequence` | required | Query sequence to compare. | Non-null. |
| reference (`FindSnpsDirect`) | `string` | required | Aligned reference bases (may contain `-`). | A,C,G,T,N case-insensitive [1]; same length as `query`. |
| query (`FindSnpsDirect`) | `string` | required | Aligned query bases (may contain `-`). | Same length as `reference`. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IEnumerable<Variant>` | One `Variant` per detected substitution. |
| `Variant.Position` | `int` | 0-based reference index of the substituted base. (VCF POS is 1-based [1]; the 1-based form is produced only by `ToVcfLines`.) |
| `Variant.ReferenceAllele` | `string` | Single reference base at the SNP position. |
| `Variant.AlternateAllele` | `string` | Single query base at the SNP position. |
| `Variant.Type` | `VariantType` | Always `SNP` for these methods. |
| `Variant.QueryPosition` | `int` | 0-based query index of the substituted base. |

### 3.3 Preconditions and Validation

- `FindSnps`: throws `ArgumentNullException` for a null reference or null query (validated in `CallVariants`).
- `FindSnpsDirect`: returns an empty sequence when either input is null or empty; throws `ArgumentException` when both are non-empty and their lengths differ (INV-06).
- Indexing is 0-based for the in-memory `Variant`. Inputs are treated as DNA over A,C,G,T(,N); base comparison is case-insensitive per VCF [1] and alleles are reported as given.

## 4. Algorithm

### 4.1 High-Level Steps

`FindSnpsDirect`:
1. If either input is null/empty → return empty; if lengths differ → `ArgumentException`.
2. Scan columns via `CallVariantsFromAlignment` (tracks ungapped reference/query coordinates).
3. Keep columns with two non-gap bases that differ case-insensitively → `Variant(Position=refPos, REF=reference[i], ALT=query[i], Type=SNP, QueryPosition=queryPos)`.

`FindSnps`:
1. Globally align reference and query (`CallVariants` → `SequenceAligner.GlobalAlign`).
2. Emit one variant per differing alignment column.
3. Filter the result to `VariantType.SNP` (drop insertions/deletions).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

No scoring tables or thresholds. The only decision rule is the per-column case-insensitive equality test `upper(reference[i]) != upper(query[i])` on non-gap columns; there are no tunable numeric constants in either method.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `FindSnpsDirect` | O(n) | O(1) lazy / O(k) materialized | `n` = aligned length; `k` = number of SNPs. |
| `FindSnps` | O(n·m) | O(n·m) | Dominated by the global alignment (`CallVariants`); the SNP filter is O(n). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [VariantCaller.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/VariantCaller.cs)

- `VariantCaller.FindSnpsDirect(string reference, string query)`: positional Hamming-mismatch SNP enumeration over equal-length aligned inputs; delegates to `CallVariantsFromAlignment` (single implementation of the column scan) and keeps SNP columns. Canonical entry point for this unit.
- `VariantCaller.FindSnps(DnaSequence reference, DnaSequence query)`: aligns then filters `CallVariants` to `VariantType.SNP`. Delegate.

### 5.2 Current Behavior

`FindSnpsDirect` is a single forward scan (the shared `CallVariantsFromAlignment` column scan) and does not use the repository suffix tree (see §7 below): SNP detection is a positional equality test between two corresponding strings, not an occurrence/substring search, so a suffix tree is not applicable. `FindSnps` delegates to `CallVariants`, which performs a Needleman–Wunsch-style global alignment and then a LINQ `Where` filter; because alignment of an arbitrary pair is not unique in repeated regions, the *positions* of SNPs adjacent to indels follow the chosen alignment (substitution-only inputs are unaffected).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- SNP = single-base substitution (REF≠ALT at one position); equal columns are not variants [1].
- Positional substitution set = Hamming mismatch set over equal-length inputs [3]; unequal lengths rejected [3].
- REF/ALT bases compared case-insensitively; a gap is not a base, so gap columns are not SNPs [1].

**Not implemented:**

- Quality/genotype/confidence scoring and VCF normalization (left-alignment, parsimony); **users should rely on:** `ToVcfLines` for serialization and a dedicated normalizer for canonical VCF positions — no in-repo alternative for normalization.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | 0-based in-memory `Variant.Position` vs 1-based VCF POS | Assumption | Position interpretation for downstream consumers | accepted | ASM consistent with sibling VARIANT-CALL-001; VCF 1-based POS emitted only by `ToVcfLines`. |
| 2 | Unequal-length `FindSnpsDirect` input | Contract | `ArgumentException` | fixed 2026-09-28 | Previously silently compared the common prefix (hid the tail); Hamming defined for equal length only [3]. |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical equal-length inputs | empty | Hamming distance 0 [3] (INV-01). |
| Gap column in `FindSnpsDirect` input | not a SNP; positions ungapped | a gap is not a base [1]. |
| `FindSnpsDirect` empty / null input | empty | nothing to compare (documented contract). |
| `FindSnpsDirect` unequal lengths | `ArgumentException` | Hamming defined for equal length only [3] (INV-06). |
| `FindSnps` null reference/query | `ArgumentNullException` | input validation in `CallVariants`. |
| Lowercase bases | compared and classified case-insensitively (`acgt` vs `ACGT` → empty) | VCF REF/ALT are case-insensitive [1]. |

### 6.2 Limitations

`FindSnpsDirect` assumes its inputs are positionally corresponding (no internal indels); a frameshift between the two strings will surface as a run of spurious SNPs — use `FindSnps` (which aligns first) for un-aligned inputs. SNP *positions* near indels reported by `FindSnps` depend on the alignment chosen and are not VCF-normalized.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
// Direct positional SNP detection (inputs already aligned / equal length).
var snps = VariantCaller.FindSnpsDirect("ATGC", "ATTC").ToList();
// snps[0] => Position 2, ReferenceAllele "G", AlternateAllele "T", Type SNP, QueryPosition 2.

// Alignment-based SNP-only detection.
var aligned = VariantCaller.FindSnps(new DnaSequence("ATGCATGC"), new DnaSequence("ATGAATGC"));
// one SNP at Position 3, REF "C", ALT "A"; no indels.
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [VariantCaller_FindSnps_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/VariantCaller_FindSnps_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [VARIANT-SNP-001-Evidence.md](../../../docs/Evidence/VARIANT-SNP-001-Evidence.md)
- Related algorithms: [Variant_Detection](./Variant_Detection.md)

## 8. References

1. The SAM/BCF/VCF working group. 2024. The Variant Call Format Specification VCFv4.3. samtools/hts-specs. https://samtools.github.io/hts-specs/VCFv4.3.pdf (source: https://raw.githubusercontent.com/samtools/hts-specs/master/VCFv4.3.tex)
2. Wikipedia. Transversion. https://en.wikipedia.org/wiki/Transversion (primary: Futuyma DJ, *Evolution*, 3rd ed., 2013, ISBN 978-1605351155)
3. Acharya T, et al. 2017. Hamming Distance as a Concept in DNA Molecular Recognition. *ACS Omega*. PMC5410656. https://pmc.ncbi.nlm.nih.gov/articles/PMC5410656/
4. Collins DW, Jukes TH. 1994. Rates of transition and transversion in coding sequences since the human-rodent divergence. *Genomics* 20(3):386–396. https://doi.org/10.1006/geno.1994.1192
