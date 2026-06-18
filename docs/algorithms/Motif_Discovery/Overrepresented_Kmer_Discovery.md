# Overrepresented k-mer Motif Discovery

| Field | Value |
|-------|-------|
| Algorithm Group | Matching / Motif Discovery |
| Test Unit ID | MOTIF-DISCOVER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Discovers candidate motifs in a single DNA sequence by enumerating every length-`k` substring (k-mer), counting how often each occurs, and ranking them by how much their observed count exceeds the count expected by chance. Overrepresentation is the observed/expected (O/E) ratio under a zero-order i.i.d. uniform background where each nucleotide is equally likely [1]. The method is deterministic and exact (no sampling): it returns, for each k-mer meeting a minimum-count cutoff, its count, its 0-based occurrence positions, and its enrichment.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Regulatory motifs (transcription-factor binding sites, etc.) tend to recur within a sequence more often than random words of the same length. A simple, well-defined way to surface candidates is to compare each k-mer's observed frequency against its expectation under a null model of random DNA [1].

### 2.2 Core Model

For a sequence of length `N`, there are `N − k + 1` length-`k` windows. Under the zero-order background model in which each of the four nucleotides is drawn independently with probability `1/4`, the expected number of occurrences of any specific k-mer is

```
E = (N − k + 1) / 4^k
```

where `4^k` is the number of distinct DNA k-mers [1]. The overrepresentation (enrichment) of a k-mer with observed count `c` is the observed/expected ratio

```
enrichment = c / E
```

A value > 1 indicates the k-mer occurs more often than chance predicts [1][2].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Zero-order i.i.d. uniform background (each base p = 1/4) | If real base composition is skewed (e.g. GC-rich), expected counts are biased and the O/E ratio over/under-states enrichment [1] |
| ASM-02 | Occurrences are counted with overlap allowed at every window | The published probability statistic warns its approximation ignores self-overlap [1]; the deterministic count used here is exact and does count overlaps |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Count equals the number of occurrences; Positions are the 0-based window starts of those occurrences | Direct window enumeration |
| INV-02 | enrichment = Count / ((N − k + 1) / 4^k) | Definition in 2.2 [1] |
| INV-03 | Every returned motif has Count ≥ minCount | Filter applied before yielding |
| INV-04 | enrichment > 0 for every returned motif | `E > 0` since `N − k + 1 ≥ 1` whenever a k-mer was counted [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | DnaSequence | required | DNA sequence to analyse | non-null |
| k | int | 6 | k-mer length | k ≥ 1 |
| minCount | int | 2 | minimum occurrence count for a k-mer to be returned | ≥ 1 in practice |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Sequence | string | the k-mer |
| Count | int | observed number of occurrences |
| Positions | IReadOnlyList&lt;int&gt; | 0-based start positions of every occurrence |
| Enrichment | double | observed/expected ratio `Count / ((N − k + 1) / 4^k)` |

### 3.3 Preconditions and Validation

Null `sequence` raises `ArgumentNullException`; `k < 1` raises `ArgumentOutOfRangeException`. Positions are 0-based window starts. The k-mer text is taken verbatim from the sequence (no normalization beyond what the `DnaSequence` already holds). When `k > N` there are no windows and the result is empty.

## 4. Algorithm

### 4.1 High-Level Steps

1. Slide a length-`k` window across the sequence; record each k-mer and the positions where it starts.
2. Compute the expected count `E = (N − k + 1) / 4^k`.
3. For each k-mer with `Count ≥ minCount`, emit a record with its count, positions, and `enrichment = Count / E`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Alphabet size `4` (A, C, G, T) is the base of `4^k`; defined as the named constant `DnaAlphabetSize`.
- A hash map from k-mer string to its list of start positions accumulates counts in one pass.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DiscoverMotifs | O(N · k) | O(N · k) | N − k + 1 substrings of length k hashed; positions stored per distinct k-mer |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.DiscoverMotifs(DnaSequence, int, int)`: enumerates k-mers, computes O/E enrichment, returns overrepresented motifs.

### 5.2 Current Behavior

The expected count is computed once per call from the closed-form `(N − k + 1) / 4^k`; there is no clamp/floor on the denominator (a previous `max(E, 0.1)` floor was an untraceable value and was removed — `E` is always strictly positive when any k-mer exists, INV-04). Counting allows overlapping occurrences. Results are yielded in hash-map enumeration order (the caller sorts if a specific order is needed).

**Search reuse:** The suffix tree was evaluated. Motif *discovery* here requires counting *all distinct k-mers and their positions in one pass* — not searching for a known pattern — so a single linear scan with a hash map is the appropriate structure; the suffix tree (best for many queries of known patterns against one text) is not used for this method.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Expected count `E = (N − k + 1) / 4^k` under the i.i.d. uniform background [1].
- Overrepresentation as the observed/expected ratio [1][2].
- Overlapping window enumeration of k-mers [1].

**Intentionally simplified:**

- Background model: zero-order uniform only; higher-order Markov backgrounds (as in monaLisa [2]) are not modelled; **consequence:** on compositionally biased sequences the enrichment may over/under-state true overrepresentation.

**Not implemented:**

- The closed-form probability/E-value `Pr(N,4,k,t)` for ≥ t occurrences [1]; **users should rely on:** the deterministic Count and O/E Enrichment fields for ranking, or an external statistical-significance tool for p-values.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null sequence | ArgumentNullException | Validation contract |
| k < 1 | ArgumentOutOfRangeException | Validation contract |
| k > N | empty result | No length-k windows exist [1] |
| Homopolymer "AAAA…" | the single k-mer dominates with high enrichment | All windows are identical |

### 6.2 Limitations

Zero-order background only (ASM-01); no statistical p-value/E-value; single-sequence (cross-sequence shared motifs are a separate unit, `FindSharedMotifs`); DNA alphabet only.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** Sequence `ATGCATGCATGC` (N=12), k=4. Windows = 12 − 4 + 1 = 9. Expected count `E = 9 / 4^4 = 9/256 = 0.03515625`. The k-mer `ATGC` occurs at positions 0, 4, 8 (Count = 3). Enrichment = `3 / (9/256) = 768/9 ≈ 85.333` [1].

**API usage example:**

```csharp
var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), k: 4, minCount: 2);
var atgc = motifs.First(m => m.Sequence == "ATGC");
// atgc.Count == 3, atgc.Positions == [0,4,8], atgc.Enrichment == 768.0/9
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_DiscoverMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/MotifFinder_DiscoverMotifs_Tests.cs) — covers `INV-01`..`INV-04`
- Evidence: [MOTIF-DISCOVER-001-Evidence.md](../../../docs/Evidence/MOTIF-DISCOVER-001-Evidence.md)

## 8. References

1. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach*, 2nd ed., Ch. 2 (Finding Regulatory Motifs). Active Learning Publishers. Formula and worked example reproduced at https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string
2. fmicompbio. monaLisa `getKmerFreq` — observed vs expected k-mer frequencies and log2 enrichment. https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html
