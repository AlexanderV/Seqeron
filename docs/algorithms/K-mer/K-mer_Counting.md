# K-mer Counting

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer Analysis |
| Test Unit ID | KMER-COUNT-001 |
| Related Projects | N/A |
| Implementation Status | N/A |
| Last Reviewed | 2026-04-30 |

## 1. Overview

K-mer counting extracts and counts all overlapping substrings of length `k` from a biological sequence. It is a foundational operation for genome assembly, metagenomics binning, sequence comparison, and repeat analysis. In this repository, the implementation provides string, `DnaSequence`, span-based, cancellation-aware, async, and both-strand counting surfaces. The counting logic is DNA-oriented in its framing, but the raw-string and span-based surfaces treat any uppercase symbols as literal k-mer characters rather than filtering to `A/C/G/T`.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A k-mer is a substring of fixed length `k` contained within a sequence. For a sequence of length `L`, there are `L - k + 1` overlapping k-mers when `k <= L`. For a DNA alphabet of size 4, the number of possible distinct k-mers is `4^k`. Sources: Wikipedia (K-mer), Rosalind (K-mer Composition), Compeau et al. (2011), Marçais & Kingsford (2011).

### 2.2 Core Model

The repository uses the standard sliding-window count model:

$$
Count(kmer) = \sum_{i=0}^{L-k} \mathbf{1}(sequence[i..i+k-1] = kmer)
$$

and stores counts in a dictionary keyed by the uppercased k-mer string.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The sum of all counts equals `L - k + 1` whenever `k <= L` | The sliding window emits exactly one k-mer per valid start position |
| INV-02 | For DNA-alphabet inputs, the number of unique k-mers is at most `min(4^k, L - k + 1)` | There are at most `4^k` DNA words and at most one emitted k-mer per window |
| INV-03 | The implementation is case-insensitive | Input is normalized to uppercase before counting |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string`, `DnaSequence`, or `ReadOnlySpan<char>` | required | Sequence to analyze | Null or empty string returns an empty dictionary |
| `k` | `int` | required | K-mer length | For non-empty string input and span-based counting, `k <= 0` throws `ArgumentOutOfRangeException`; `k > sequence.Length` returns an empty dictionary |
| `cancellationToken` | `CancellationToken` | optional | Cancellation support for long-running counting | Used in cancellation-aware overloads |
| `progress` | `IProgress<double>?` | optional | Progress reporter for long-running counting | Reports values between `0.0` and `1.0` |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `counts` | `Dictionary<string, int>` | Maps each observed k-mer to its count |

### 3.3 Preconditions and Validation

String-based counting returns an empty dictionary for null or empty input before validating `k`. The `DnaSequence` overloads delegate to that string-based path, so an empty `DnaSequence.Sequence` also returns an empty dictionary before `k` is checked. For non-empty string input, the cancellation-aware and synchronous string overloads throw `ArgumentOutOfRangeException` when `k <= 0`, and the span-based path validates `k` before checking sequence length. When `k` is greater than the sequence length, counting returns an empty dictionary. The string and cancellation-aware overloads uppercase the input before emitting dictionary keys.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input sequence and `k`.
2. Normalize the sequence to uppercase when using string-based overloads.
3. Slide a window of length `k` across every valid start position.
4. Increment the count for the observed k-mer in the output dictionary.
5. Return the completed count map.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Entry points documented in the original file and confirmed in source:

| Method | Class | Description |
|--------|-------|-------------|
| `CountKmers(string, k)` | `KmerAnalyzer` | Canonical string-based counting |
| `CountKmers(DnaSequence, k)` | `KmerAnalyzer` | Wrapper for `DnaSequence` input |
| `CountKmersSpan(ReadOnlySpan<char>, k)` | `KmerAnalyzer` / `SequenceExtensions` | Span-based counting variant |
| `CountKmersBothStrands(DnaSequence, k)` | `KmerAnalyzer` | Forward plus reverse-complement counting |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `CountKmers` | `O(n × k)` effective | `O(u)` | `n` is sequence length, `k` is k-mer length, and each window creates and hashes a length-`k` string key |
| `CountKmersBothStrands` | `O(n × k)` effective | `O(u)` | Counts forward and reverse-complement sequences separately, with the same per-window string allocation and hashing cost |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs), [SequenceExtensions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs)

- `KmerAnalyzer.CountKmers(...)`: Canonical counting overloads for strings and `DnaSequence` values.
- `KmerAnalyzer.CountKmersAsync(...)`: Async wrapper over cancellation-aware counting.
- `KmerAnalyzer.CountKmersSpan(ReadOnlySpan<char>, int)`: Span-based counting entry point.
- `KmerAnalyzer.CountKmersBothStrands(DnaSequence, int)`: Counts forward and reverse-complement k-mers and combines the totals.

### 5.2 Current Behavior

All string-based entry points uppercase the input before counting. The synchronous `CountKmers(...)` methods are stateless, and the cancellation-aware overload checks cancellation periodically while optionally reporting progress. `CountKmersSpan(...)` delegates to the span-based helper in `SequenceExtensions`, and `CountKmersBothStrands(...)` counts forward and reverse-complement sequences independently before summing counts. The implementation materializes string keys for observed windows, so runtime includes per-window length-`k` string allocation and hashing work in addition to the sliding-window scan. The raw-string and span-based paths do not restrict the alphabet, so ambiguous or non-ACGT symbols are preserved as literal k-mer keys.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Overlapping sliding-window k-mer extraction.
- Exact count accumulation in a dictionary keyed by the observed k-mer.
- Case-insensitive normalization for string-based input.

**Intentionally simplified:**

- Both-strand counting sums forward and reverse-complement counts rather than canonicalizing each k-mer to a single representative key; **consequence:** forward and reverse-complement words remain separate dictionary entries unless they are identical strings.
- Raw-string and span-based counting do not enforce a DNA alphabet; **consequence:** the usual `4^k` combinatorial bound applies only when the input is restricted to DNA symbols.

**Not implemented:**

- Canonical reverse-complement collapsing as the default count representation; **users should rely on:** downstream post-processing if canonical k-mer keys are required.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty or null string input | Returns an empty dictionary | Explicit source guard |
| `k <= 0` with non-empty string input or span input | Throws `ArgumentOutOfRangeException` | Invalid k-mer length after the relevant input guard |
| Empty or null string input with `k <= 0` | Returns an empty dictionary | The string-based overload short-circuits before validating `k` |
| `k > sequence.Length` | Returns an empty dictionary | No valid windows exist |
| Homopolymer such as `AAAA`, `k = 2` | Returns one key with count `L - k + 1` | Every window is identical |
| String input with ambiguous or non-ACGT symbols | Counts those symbols literally after uppercasing | The string/span counting logic does not filter the alphabet |

### 6.2 Limitations

The current implementation uses string keys for observed k-mers and does not canonicalize reverse complements automatically. It is exact and general-purpose, but large genomes or large `k` values can still create substantial runtime and dictionary pressure because each window materializes and hashes a length-`k` key. For raw-string and span-based inputs, ambiguous symbols are retained rather than normalized to a DNA alphabet.

## 7. Examples and Related Material

### 7.2 Applications and Use Cases (Optional)

- Genome assembly via de Bruijn graph construction.
- Metagenomics binning based on k-mer signatures.
- Alignment-free sequence comparison.
- Error detection through k-mer spectra.
- Repeat analysis in repetitive regions.

## 8. References

1. Wikipedia. "K-mer." https://en.wikipedia.org/wiki/K-mer
2. Rosalind. "K-mer Composition." https://rosalind.info/problems/kmer/
3. Compeau, P.E.C., Pevzner, P.A., Tesler, G. (2011). "How to apply de Bruijn graphs to genome assembly." Nature Biotechnology, 29(11), 987–991.
4. Marçais, G., Kingsford, C. (2011). "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers." Bioinformatics, 27(6), 764–770.
