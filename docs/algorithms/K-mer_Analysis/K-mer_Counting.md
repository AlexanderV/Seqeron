# K-mer Counting

## Overview

K-mer counting is a fundamental operation in computational genomics that extracts and counts all substrings of length *k* from a biological sequence. This technique is widely used for sequence assembly, genome comparison, metagenomics binning, and alignment-free sequence analysis.

## Definition

A **k-mer** is a substring of length *k* contained within a biological sequence. For a sequence of length *L*, there are exactly **L − k + 1** k-mers (overlapping substrings).

**Key Properties:**
- Total possible k-mers for DNA alphabet (ACGT): 4^k
- Sum of all k-mer counts = L − k + 1 (the **invariant**)
- K-mers are typically case-insensitive (converted to uppercase)

## Algorithm

### Sliding Window Approach

```
procedure CountKmers(sequence, k):
    L ← length(sequence)
    counts ← empty dictionary

    for i ← 0 to L − k:
        kmer ← substring(sequence, i, k)
        counts[kmer] ← counts[kmer] + 1

    return counts
```

**Complexity:** O(n) where n = sequence length

### K-mer Spectrum

The **k-mer spectrum** shows the frequency distribution: how many k-mers appear 1 time, 2 times, etc. This is computed by inverting the k-mer count dictionary.

## Edge Cases

| Scenario | Expected Result |
|----------|-----------------|
| Empty sequence | Empty dictionary |
| k ≤ 0 | ArgumentOutOfRangeException |
| k > sequence length | Empty dictionary |
| null sequence | Empty dictionary |
| Homopolymer (e.g., "AAAA", k=2) | Single k-mer with count L-k+1 |

## Invariants

1. **Total count invariant:** Sum of all k-mer counts = L − k + 1
2. **Unique k-mers bound:** Number of unique k-mers ≤ min(4^k, L − k + 1)
3. **Symmetry (DNA):** For canonical k-mers, forward and reverse complement can be combined

## Applications

- **Genome assembly:** De Bruijn graph construction
- **Metagenomics:** Species binning based on k-mer signatures
- **Sequence comparison:** Alignment-free distance metrics
- **Error detection:** Identifying sequencing errors via k-mer spectra
- **Repeat analysis:** Detecting repetitive regions

## Implementation Notes

### Current Implementation

The library provides multiple entry points:

| Method | Class | Description |
|--------|-------|-------------|
| `CountKmers(string, k)` | KmerAnalyzer | Canonical string-based counting |
| `CountKmers(DnaSequence, k)` | KmerAnalyzer | Wrapper for DnaSequence type |
| `CountKmersSpan(ReadOnlySpan<char>, k)` | SequenceExtensions | Memory-efficient span-based variant |
| `CountKmersBothStrands(DnaSequence, k)` | KmerAnalyzer | Forward + reverse complement |

### Case Handling

All implementations normalize to uppercase before counting, ensuring case-insensitive matching.

### Thread Safety

The synchronous `CountKmers` methods are stateless and thread-safe. The async variants support `CancellationToken` for cancellation.

## References

1. Wikipedia. "K-mer." https://en.wikipedia.org/wiki/K-mer
2. Rosalind. "K-mer Composition." https://rosalind.info/problems/kmer/
3. Compeau, P.E.C., Pevzner, P.A., Tesler, G. (2011). "How to apply de Bruijn graphs to genome assembly." Nature Biotechnology, 29(11), 987–991.
4. Marçais, G., Kingsford, C. (2011). "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers." Bioinformatics, 27(6), 764–770. (Jellyfish)

## Test Unit

- **ID:** KMER-COUNT-001
- **Methods:** CountKmers, CountKmersSpan, CountKmersBothStrands
