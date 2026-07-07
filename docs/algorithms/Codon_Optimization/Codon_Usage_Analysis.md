# Codon Usage Analysis

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Optimization |
| Test Unit ID | CODON-USAGE-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Codon usage analysis measures how often codons appear in a coding sequence and compares codon distributions between sequences.[1][2] In this repository, `CodonOptimizer.CalculateCodonUsage` returns raw codon counts, while `CodonOptimizer.CompareCodonUsage` compares normalized codon-frequency distributions using a total-variation-distance similarity score. The implementation is case-insensitive, normalizes DNA to RNA notation, and ignores incomplete trailing bases. These methods are intended for direct sequence-level codon-profile analysis rather than organism-wide bias modeling by themselves.[2][4]

## 2. Scientific / Formal Basis

### 2.1 Domain Context

The genetic code contains 64 codons encoding 20 amino acids plus stop signals, and most amino acids have multiple synonymous codons. Codon usage bias reflects biological factors such as tRNA abundance, genome GC content, translational selection, and mutational bias.[2][3]

### 2.2 Core Model

Given a coding sequence, codon usage is computed by splitting the sequence into non-overlapping triplets and counting occurrences of each codon:

$$
\mathrm{Count}(c) = |\{ i : \text{seq}[3i:3i+3] = c \}|
$$

Normalized codon frequencies are then:

$$
f(c) = \frac{\text{Count}(c)}{\sum_{c'} \text{Count}(c')}
$$

The repository compares two sequences using the normalized absolute-difference metric documented in the original file and the test specification:

$$
\mathrm{Similarity} = 1 - \frac{\sum_c |f_1(c) - f_2(c)|}{2}
$$

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The input is interpreted in-frame as coding triplets. | Trailing bases are ignored and counts may not represent the intended codon stream. |
| ASM-02 | Comparing codon-frequency distributions is a meaningful proxy for codon-usage similarity. | Similarity can be mathematically correct yet biologically uninformative for the use case. |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `sum(counts.Values) == floor(sequence.Length / 3)` after normalization to complete codons. | The counter increments exactly once per extracted codon. |
| INV-02 | `0 <= similarity <= 1`. | The method uses a normalized total-variation-distance formula. |
| INV-03 | `CompareCodonUsage(a, b) == CompareCodonUsage(b, a)`. | Absolute differences are symmetric. |
| INV-04 | `CompareCodonUsage(s, s) == 1` for non-empty `s`. | The two normalized distributions are identical. |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `[CalculateCodonUsage] codingSequence` | `string` | required | DNA or RNA sequence to count by codon. | `null` or empty input returns an empty dictionary. |
| `[CompareCodonUsage] sequence1` | `string` | required | First sequence in the comparison. | Empty input contributes no codons. |
| `[CompareCodonUsage] sequence2` | `string` | required | Second sequence in the comparison. | Empty input contributes no codons. |

### 3.2 Output / Return Value

| Name | Type | Description |
|------|------|-------------|
| `CalculateCodonUsage` result | `Dictionary<string, int>` | Raw counts for the codons observed in the normalized input. |
| `CompareCodonUsage` result | `double` | Similarity in `[0, 1]` based on normalized codon-frequency differences. |

### 3.3 Preconditions and Validation

Both methods uppercase the input and convert `T` to `U`. Codons are extracted only from complete triplets, so trailing one or two bases are ignored. `CompareCodonUsage` returns `0` when either sequence yields zero codons after preprocessing.

## 4. Algorithm

### 4.1 High-Level Steps

1. Normalize the input to uppercase RNA notation.
2. Split the sequence into complete codons.
3. For `CalculateCodonUsage`, increment a count for each observed codon.
4. For `CompareCodonUsage`, compute counts for both sequences.
5. Form the union of codons observed in either sequence.
6. Convert counts to frequencies using the codon totals for each sequence.
7. Sum the absolute frequency differences and return `1 - sum / 2`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The comparison metric is the same total-variation-distance similarity documented in the spec.[5]

| Situation | Expected Similarity |
|-----------|---------------------|
| Identical non-empty codon distributions | `1.0` |
| Completely disjoint codon distributions | `0.0` |
| Partially overlapping distributions | Strictly between `0` and `1` |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CalculateCodonUsage` | `O(n)` | `O(64)` | Raw codon counting over the sequence. |
| `CompareCodonUsage` | `O(n + m)` | `O(64)` | Counts both sequences and compares the resulting distributions. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonOptimizer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs)

- `CodonOptimizer.CalculateCodonUsage(string)`
- `CodonOptimizer.CompareCodonUsage(string, string)`

### 5.2 Current Behavior

`CalculateCodonUsage` returns counts only for codons observed in the normalized input; it does not pre-populate all 64 codons. `CompareCodonUsage` calls `CalculateCodonUsage` for both sequences, unions the observed codon keys, and computes `1 - (sum(abs(freq1 - freq2)) / 2)`. If both sequences are empty, or if either sequence has zero complete codons, the method returns `0`.[5]

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Codon usage is computed by counting non-overlapping triplets.[2]
- Comparison uses the normalized absolute-difference similarity `1 - Σ|f₁-f₂|/2` documented in the repository test specification.[5]

**Intentionally simplified:**

- The count output includes only codons present in the sequence; **consequence:** absent codons are represented implicitly rather than as explicit zero entries.
- The comparison metric is a single scalar similarity rather than a richer codon-bias profile; **consequence:** different distribution shapes can collapse to the same summary score.

**Not implemented:**

- Relative synonymous codon usage (RSCU), codon-pair bias, or position-specific codon statistics; **users should rely on:** no current alternative.
- Organism-table-aware normalization inside these methods; **users should rely on:** [CAI_Calculation.md](CAI_Calculation.md) or caller-side interpretation.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty sequence | `CalculateCodonUsage` returns `{}`. | Explicit empty-input handling. |
| One empty sequence in a comparison | Similarity is `0`. | One distribution has zero total codons. |
| Both sequences empty | Similarity is `0`. | There is no data to compare. |
| Incomplete trailing bases | Ignored. | Codon splitting requires three bases. |
| DNA input | Converted to RNA notation before counting. | Internal normalization uses `T -> U`. |

### 6.2 Limitations

These methods operate only on direct codon counts and normalized frequency differences. They do not attach biological weighting, do not validate whether a sequence is a true CDS, and do not infer organism-specific codon bias without external context.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonOptimizer_CodonUsage_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonOptimizer_CodonUsage_Tests.cs) — covers `INV-01`, `INV-02`, `INV-03`, `INV-04`
- Test specification: [CODON-USAGE-001.md](../../../tests/TestSpecs/CODON-USAGE-001.md)
- Related algorithms: [CAI_Calculation.md](CAI_Calculation.md), [Rare_Codon_Detection.md](Rare_Codon_Detection.md), [Sequence_Optimization.md](Sequence_Optimization.md)

## 8. References

1. Sharp PM, Li WH. 1987. The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research. N/A
2. Plotkin JB, Kudla G. 2011. Synonymous but not the same: the causes and consequences of codon bias. Nature Reviews Genetics. N/A
3. Wikipedia contributors. 2026. Codon usage bias. Wikipedia. https://en.wikipedia.org/wiki/Codon_usage_bias
4. Kazusa Codon Usage Database. 2026. https://www.kazusa.or.jp/codon/
5. Test specification: [CODON-USAGE-001.md](../../../tests/TestSpecs/CODON-USAGE-001.md)
