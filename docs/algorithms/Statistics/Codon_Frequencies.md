# Codon Frequencies

| Field | Value |
|-------|-------|
| Algorithm Group | Statistics |
| Test Unit ID | SEQ-CODON-FREQ-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Computes the usage frequency of each codon (nucleotide triplet) in a DNA coding sequence by reading consecutive, non-overlapping in-frame triplets and dividing each codon's count by the total number of counted codons. This is the count/total convention of the Kazusa Codon Usage Database (CUTG) [1][2], the standard way to summarise codon usage of a gene or genome. The result is exact (not heuristic): it is a deterministic frequency table over the observed codons.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code reads coding DNA in non-overlapping triplets (codons) within a reading frame. Because the code is degenerate (synonymous codons encode the same amino acid), organisms differ in how often each codon is used — codon usage bias [4]. Quantifying usage starts from raw codon frequencies: the fraction of all codons in the input that are a given triplet [2].

### 2.2 Core Model

For a sequence read from frame offset `f`, let the non-overlapping triplets be `c_1, c_2, …` where `c_k` covers bases `[f + 3(k-1), f + 3k)`. Let `V` be the set of triplets composed only of A, C, G, T (ambiguous codons are excluded [2]). Then for codon `x`:

- count(x) = number of `c_k ∈ V` equal to `x`
- total = |{ c_k : c_k ∈ V }|
- frequency(x) = count(x) / total  [2]

Kazusa CUTG reports this scaled per thousand: "the frequency (per thousand) of codon use in each organism was calculated by summing up the numbers of codons used" [2]; hence frequency(x) = (CUTG per-thousand value) / 1000. This count/total fraction is distinct from the per-amino-acid "Fraction" column of EMBOSS cusp, which divides by the synonymous-codon group total rather than by all codons [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Each frequency ∈ (0, 1]; only observed codons are keys | count(x) ≥ 1 for every key and count(x) ≤ total [2] |
| INV-02 | Σ frequency(x) = 1 over all keys (when total ≥ 1) | Σ count(x) = total by definition [2] |
| INV-03 | Codons with any non-ACGT base never appear and never change total | ambiguous codons excluded from count [2] |
| INV-04 | Output is independent of input letter case | input is upper-cased before counting |
| INV-05 | frequency(x) = CUTG per-thousand value ÷ 1000 | cusp dataset: 22/386×1000 = 56.995 [3] |

### 2.5 Comparison with Related Methods

| Aspect | This method (count / total) | EMBOSS cusp "Fraction" | `CodonOptimizer.CalculateCodonUsage` |
|--------|-----------------------------|------------------------|--------------------------------------|
| Denominator | all counted codons | synonymous-codon group of the amino acid [3] | n/a (returns raw counts) |
| Output | frequencies summing to 1 | per-amino-acid proportions | integer counts |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dnaSequence | string | required | DNA coding sequence | case-insensitive; non-ACGT bases excluded from counting |
| readingFrame | int | 0 | 0-based offset of the first codon | typically 0, 1, or 2 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | IReadOnlyDictionary&lt;string, double&gt; | codon → frequency (count / total counted codons); empty if no valid codon |

### 3.3 Preconditions and Validation

`null`, empty, or length &lt; 3 returns an empty dictionary. Input is upper-cased (T↔U is not performed; RNA `U` is treated as a non-ACGT base and excluded). Counting is 0-based starting at `readingFrame`; only complete non-overlapping triplets are read, so trailing 1–2 bases are ignored. If no triplet is composed solely of ACGT (`total = 0`), the result is empty — there is no division by zero.

## 4. Algorithm

### 4.1 High-Level Steps

1. Guard: null / empty / length &lt; 3 → empty table.
2. Upper-case the sequence.
3. Step `i` from `readingFrame` to `length − 3` in increments of 3; take the triplet at `i`.
4. If the triplet is all ACGT, increment its count and the running total.
5. Divide each codon count by the total to produce frequencies.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Codon length = 3 (fixed by the genetic code; non-overlapping reading) [2].
- Valid-base alphabet = {A, C, G, T}; any other character makes the whole triplet ineligible [2].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateCodonFrequencies | O(n) | O(k) | n = sequence length; k = number of distinct observed codons (≤ 64) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceStatistics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs)

- `SequenceStatistics.CalculateCodonFrequencies(string dnaSequence, int readingFrame = 0)`: returns the codon → frequency map per the count/total definition.

### 5.2 Current Behavior

Single linear pass over the sequence. The method is a frequency tabulation, not a search/matching operation, so the repository suffix tree is not applicable (no occurrence enumeration; one O(n) scan with a small hash map). Behavior matches the contract: case-insensitive, non-overlapping in-frame triplets, ambiguous codons and trailing bases ignored, empty table when no valid codon exists.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- count/total codon frequency over non-overlapping in-frame triplets [2].
- exclusion of ambiguous (non-ACGT) codons from the count [2].
- frequencies sum to 1 (INV-02), equal to CUTG per-thousand ÷ 1000 (INV-05) [2][3].

**Intentionally simplified:**

- Per-thousand scaling and the per-amino-acid Fraction column reported by Kazusa/cusp are not produced here; **consequence:** callers needing per-thousand multiply by 1000, and per-amino-acid proportions are out of scope for this method.

**Not implemented:**

- RNA `U`-aware counting; **users should rely on:** converting U→T before calling, since CUTG tabulates DNA CDS [2].

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null / empty / length &lt; 3 | empty dictionary | no full codon to count |
| trailing 1–2 bases | ignored | non-overlapping triplets only [2] |
| triplet with non-ACGT base | excluded from count and total | ambiguous codons excluded [2] |
| all triplets ambiguous (total = 0) | empty dictionary | only well-defined count/total result; no division by zero |
| lowercase input | same as uppercase | input upper-cased (INV-04) |

### 6.2 Limitations

Computes raw codon usage only; does not derive codon-usage indices (CAI, Fop, Nc) [4], does not interpret reading frames biologically (no ORF detection), and treats RNA `U` as ambiguous unless converted to `T` first.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var freq = SequenceStatistics.CalculateCodonFrequencies("ATGATGAAA", readingFrame: 0);
// freq["ATG"] == 2.0/3.0, freq["AAA"] == 1.0/3.0
```

**Numerical walk-through:** `ATGATGAAA` from frame 0 yields triplets ATG, ATG, AAA. total = 3, count(ATG) = 2, count(AAA) = 1, so frequency(ATG) = 2/3, frequency(AAA) = 1/3, summing to 1 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceStatistics_CalculateCodonFrequencies_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/SequenceStatistics_CalculateCodonFrequencies_Tests.cs) — covers INV-01..INV-05
- Evidence: [SEQ-CODON-FREQ-001-Evidence.md](../../../docs/Evidence/SEQ-CODON-FREQ-001-Evidence.md)
- Related algorithms: [Dinucleotide_Analysis](Dinucleotide_Analysis.md)

## 8. References

1. Nakamura Y, Gojobori T, Ikemura T. 2000. Codon usage tabulated from international DNA sequence databases: status for the year 2000. Nucleic Acids Research 28(1):292. https://doi.org/10.1093/nar/28.1.292
2. Kazusa DNA Research Institute. Codon Usage Database (CUTG) — README. https://www.kazusa.or.jp/codon/readme_codon.html
3. Rice P, Longden I, Bleasby A. 2000. EMBOSS — `cusp` application documentation. https://emboss.sourceforge.net/apps/cvs/emboss/apps/cusp.html
4. Wikipedia. Codon usage bias. https://en.wikipedia.org/wiki/Codon_usage_bias
