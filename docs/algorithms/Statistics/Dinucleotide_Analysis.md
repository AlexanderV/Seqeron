# Dinucleotide Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-DINUC-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Dinucleotide Analysis computes three exact, specification-driven compositional statistics over a nucleotide sequence: normalized dinucleotide frequencies, the dinucleotide relative abundance (observed/expected "odds ratio", the basis of Karlin's genomic signature), and codon usage frequencies for a chosen reading frame. The dinucleotide odds ratio quantifies whether neighbouring bases co-occur more or less than expected under positional independence; a value of 1 indicates no bias [1]. All three are exact O(n) counting computations rather than heuristics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Neighbouring bases in genomes are not independent: dinucleotides such as TpA are pervasively under-represented and others over-represented, and the resulting set of 16 relative-abundance values forms a stable "genomic signature" that discriminates organisms [1]. Codon usage describes the non-uniform use of synonymous codons in protein-coding sequences and is tabulated per organism by counting non-overlapping triplets across all CDS [4].

### 2.2 Core Model

**Dinucleotide frequency.** For a sequence of length N, f_XY = count(XY) / (N − 1), the normalized frequency over the N−1 adjacent dinucleotide positions [1].

**Dinucleotide relative abundance (odds ratio).** ρ_XY = f_XY / (f_X · f_Y), where f_X is the normalized single-base frequency [1]. ρ = 1 means the dinucleotide occurs exactly as expected under independence (no bias); ρ > 1 over-representation, ρ < 1 under-representation [1]. The widely used Karlin & Burge (1995) interpretive criterion classifies a dinucleotide as under-represented when ρ ≤ 0.78 and over-represented when ρ ≥ 1.23 [2]; this library returns the raw ρ and leaves classification to the caller.

The same odds-ratio shape underlies the CpG O/E ratio of Gardiner-Garden & Frommer (1987), O/E = (#CpG/N) / ((#C/N)·(#G/N)), which differs only by normalizing the dinucleotide count by N rather than N−1 [3]; this library follows the Karlin N−1 convention.

**Codon frequency.** Reading consecutive non-overlapping triplets from a frame offset, frequency(codon) = count(codon) / (total counted codons) [4].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Dinucleotide frequencies sum to 1.0 (≥1 valid dinucleotide) | each valid position contributes 1/total [1] |
| INV-02 | Codon frequencies sum to 1.0 (≥1 valid codon) | count/total over all counted codons [4] |
| INV-03 | ρ_XY = 1.0 ⇔ f_XY = f_X·f_Y (no bias); ρ ≥ 0 | definition of relative abundance [1] |
| INV-04 | All returned frequencies and ratios are finite and ≥ 0 | non-negative ratios of non-negative counts |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | string | required | Nucleotide sequence | DNA/RNA; case-insensitive; alphabet {A,T,G,C,U} for dinucleotides, {A,T,G,C} for codons |
| readingFrame | int | 0 | Codon frame offset | 0, 1, or 2 (only `CalculateCodonFrequencies`) |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IReadOnlyDictionary<string,double>` | dinucleotide/codon → frequency, or dinucleotide → ρ ratio |

### 3.3 Preconditions and Validation

Input is normalized to upper case (`ToUpperInvariant`). `CalculateDinucleotideFrequencies`/`CalculateDinucleotideRatios` return an empty dictionary for null/empty input or length < 2; `CalculateCodonFrequencies` returns empty for null/empty or length < 3. Dinucleotides outside {A,T,G,C,U} and triplets outside {A,T,G,C} are excluded from counts. When a constituent base of a dinucleotide is absent (expected = 0), its ρ is reported as 0 to avoid division by zero.

## 4. Algorithm

### 4.1 High-Level Steps

1. Upper-case the sequence.
2. Scan adjacent pairs (or frame-stepped triplets), counting only alphabet-valid units; track the total.
3. Divide each count by the total to obtain frequencies; for ρ, divide the dinucleotide frequency by the product of its two base frequencies.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Dinucleotide alphabet filter: `{A,T,G,C,U}`; codon alphabet filter: `{A,T,G,C}`.
- Karlin (N−1) normalization for dinucleotide frequency [1]; interpretive thresholds 0.78 / 1.23 are documentation-only and not applied in code [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| any of the three methods | O(n) | O(k) | k = number of distinct units observed (≤ 25 dinucleotides, ≤ 64 codons) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateDinucleotideFrequencies(string)`: normalized dinucleotide frequencies, count/(N−1).
- `SequenceStatistics.CalculateDinucleotideRatios(string)`: ρ_XY = f_XY/(f_X·f_Y).
- `SequenceStatistics.CalculateCodonFrequencies(string,int)`: non-overlapping triplet frequencies for a frame.

### 5.2 Current Behavior

`CalculateDinucleotideRatios` derives single-base frequencies from `CalculateNucleotideComposition` over {A,T,G,C,U} and dinucleotide frequencies from `CalculateDinucleotideFrequencies`, so RNA (U) is supported. No suffix-tree reuse: these are single linear counting passes (no substring search / occurrence enumeration), so the repository suffix tree does not apply.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- ρ_XY = f_XY/(f_X·f_Y) with ρ=1 as the no-bias baseline [1].
- Normalized dinucleotide frequency over N−1 positions [1].
- Codon frequency as count/total over non-overlapping triplets per frame; non-ACGT triplets excluded; trailing bases ignored [4].

**Intentionally simplified:**

- Dinucleotide odds ratio uses single-strand frequencies, not the strand-symmetrized ρ* (concatenation with the reverse complement) used by Karlin for genomic signatures [1]; **consequence:** values are single-strand relative abundances, appropriate for per-sequence O/E (e.g. CpG) but not strand-symmetrized signatures.

**Not implemented:**

- Over-/under-representation classification (0.78 / 1.23 thresholds) [2]; **users should rely on:** comparing the returned ρ against those thresholds themselves.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | (N−1) vs N dinucleotide normalization | Assumption | numeric ratio differs by N/(N−1) from the Gardiner-Garden CpG form | accepted | Karlin convention [1]; both authoritative |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty / length < 2 (dinuc) | empty dictionary | input guard |
| length < 3 (codon) | empty dictionary | input guard |
| constituent base absent (expected = 0) | ρ = 0 for that dinucleotide | division-by-zero guard [1] |
| non-ACGT triplet | excluded from codon counts | [4] |
| trailing 1–2 bases | ignored | non-overlapping triplets [4] |

### 6.2 Limitations

Single-strand only (no ρ* symmetrization); no classification output; codon analysis assumes the caller supplies the correct frame and treats the input as a single contiguous coding region.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** For `ATGCGCGT` (A=1,T=2,G=3,C=2; N=8; dinucleotide positions = 7): f_GC = 2/7, f_G = 3/8, f_C = 2/8 ⇒ ρ_GC = (2/7)/((3/8)(2/8)) = 64/21 ≈ 3.0476. ρ_AT = (1/7)/((1/8)(2/8)) = 32/7 ≈ 4.5714. For codons of `ATGATGAAA` in frame 0: ATG, ATG, AAA ⇒ ATG = 2/3, AAA = 1/3.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateDinucleotide_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateDinucleotide_Tests.cs) — covers INV-01..INV-04
- Evidence: [SEQ-DINUC-001-Evidence.md](../../../docs/Evidence/SEQ-DINUC-001-Evidence.md)

## 8. References

1. Karlin S. 1998. Pervasive properties of the genomic signature. PMC. https://pmc.ncbi.nlm.nih.gov/articles/PMC126251/
2. Karlin S., Burge C. 1995. Dinucleotide relative abundance extremes: a genomic signature. Trends in Genetics 11(7):283-290. https://doi.org/10.1016/S0168-9525(00)89076-9 (criterion ρ≤0.78 / ρ≥1.23 retrieved from https://academic.oup.com/mbe/article/19/6/964/1095097)
3. Gardiner-Garden M., Frommer M. 1987. CpG islands in vertebrate genomes. J Mol Biol 196(2):261-282. https://doi.org/10.1016/0022-2836(87)90689-9
4. Nakamura Y., Gojobori T., Ikemura T. Codon Usage Database (CUTG), Kazusa. https://www.kazusa.or.jp/codon/readme_codon.html
