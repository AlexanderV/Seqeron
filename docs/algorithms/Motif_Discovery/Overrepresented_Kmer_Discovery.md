# Overrepresented k-mer Motif Discovery

| Field | Value |
|-------|-------|
| Algorithm Group | Matching / Motif Discovery |
| Test Unit ID | MOTIF-DISCOVER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-09-29 |

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

**Bernoulli (independent, non-uniform) background** — RSAT `oligo-analysis` [3][4]: with residue
probabilities `q = (q_A, q_C, q_G, q_T)` (normalised to sum 1), the expected frequency of a word
`w = w_1…w_k` is `p(w) = ∏ q[w_i]` and its expected number of occurrences is
`E(w) = p(w) · (N − k + 1)` (RSAT: `exp_occ = exp_freq × sum_occurrences`, overlapping windows);
`enrichment = c / E(w)` (RSAT "ratio"). With `q = (¼,¼,¼,¼)` this is exactly the uniform model above
(RSAT "equiprobable": `exp_freq = 1/4^k`).

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Zero-order background: i.i.d. uniform (each base p = 1/4) by default, or a caller-supplied Bernoulli composition | With the uniform default on skewed composition (e.g. GC-rich) expected counts are biased [1]; pass the composition via the background overload [3] |
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
| background (overload) | IReadOnlyList&lt;double&gt; | — | Bernoulli residue probabilities A, C, G, T | exactly 4 finite, strictly positive values; normalised to sum 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| Sequence | string | the k-mer |
| Count | int | observed number of occurrences |
| Positions | IReadOnlyList&lt;int&gt; | 0-based start positions of every occurrence |
| Enrichment | double | observed/expected ratio `Count / ((N − k + 1) / 4^k)` |

### 3.3 Preconditions and Validation

Null `sequence` (or `background`) raises `ArgumentNullException`; `k < 1` or a non-finite/non-positive background value raises `ArgumentOutOfRangeException`; a background without exactly 4 values raises `ArgumentException`. Validation is eager; enumeration is lazy. Positions are 0-based window starts. The k-mer text is taken verbatim from the sequence (no normalization beyond what the `DnaSequence` already holds). When `k > N` there are no windows and the result is empty.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count every overlapping k-mer with the canonical `SequenceExtensions.CountKmersSpan` (the counter `KmerAnalyzer.CountKmersSpan` delegates to).
2. For the k-mers with `Count ≥ minCount` only, collect their 0-based start positions in a second window pass (span-keyed dictionary lookup, no per-window string allocation).
3. Emit, in order of first occurrence, a record with count, ascending positions and `enrichment = Count / (W · ∏ q[w_i])`, `W = N − k + 1`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- The uniform default is the background `(¼,¼,¼,¼)`, so both entry points share one code path.
- The product `∏ q[w_i]` is carried as mantissa · 2^exponent with exact power-of-two rescaling (`Math.ScaleB`), so the ratio stays finite whenever its true value is representable (4^k alone overflows a double for k ≥ 512); for `q = ¼` every step is exact and the result is the correctly rounded `Count · 4^k / W`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| DiscoverMotifs | O(N · k) | O(N · k) | two passes of N − k + 1 hashed windows; positions stored only for k-mers with Count ≥ minCount |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [MotifFinder.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.cs)

- `MotifFinder.DiscoverMotifs(DnaSequence, int, int)`: uniform background (Compeau & Pevzner / RSAT equiprobable).
- `MotifFinder.DiscoverMotifs(DnaSequence, int, int, IReadOnlyList<double>)`: Bernoulli background (RSAT oligo-analysis).
- MCP `discover_motifs` (Seqeron.Mcp.Analysis `AnalysisTools.DiscoverMotifs`) delegates to the uniform overload.

### 5.2 Current Behavior

The expected count is computed once per call from the closed-form `(N − k + 1) / 4^k`; there is no clamp/floor on the denominator (a previous `max(E, 0.1)` floor was an untraceable value and was removed — `E` is always strictly positive when any k-mer exists, INV-04). Counting allows overlapping occurrences. Results are yielded in order of each k-mer's first occurrence; positions are ascending.

**Search reuse:** The suffix tree was evaluated. Motif *discovery* here requires counting *all distinct k-mers and their positions in one pass* — not searching for a known pattern — so a single linear scan with a hash map is the appropriate structure; the suffix tree (best for many queries of known patterns against one text) is not used for this method.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Expected count `E = (N − k + 1) / 4^k` under the i.i.d. uniform background [1].
- Overrepresentation as the observed/expected ratio [1][2].
- Overlapping window enumeration of k-mers [1].

- Bernoulli (independent, non-uniform) background, `E = (N − k + 1) · ∏ q[w_i]` [3][4] (overload; uniform default unchanged).

**Not implemented (out of the O/E contract; declared):**

- Markov-chain backgrounds of order ≥ 1 (RSAT `-bg` Markov models, monaLisa [2]).
- Significance statistics — RSAT binomial `occ_P = P(X ≥ occ)`, `occ_E = occ_P × number of tested words`, `occ_sig = −log10 occ_E` [3][4] — and the textbook `Pr(N,4,k,t)` [1]; **users should rely on:** the deterministic Count and O/E Enrichment for ranking, or an external significance tool.
- Reverse-complement grouping / both-strand counting (RSAT `-2str`).

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

With the Bernoulli background `q = (0.3, 0.2, 0.2, 0.3)`: `p(ATGC) = 0.3·0.3·0.2·0.2 = 0.0036`, `E = 9 · 0.0036 = 0.0324`, enrichment = `3 / 0.0324 = 92.5926` (TGCA/GCAT/CATG: `2 / 0.0324 = 61.7284`) [3][4].

**API usage example:**

```csharp
var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), k: 4, minCount: 2);
var atgc = motifs.First(m => m.Sequence == "ATGC");
// atgc.Count == 3, atgc.Positions == [0,4,8], atgc.Enrichment == 768.0/9
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_DiscoverMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_DiscoverMotifs_Tests.cs) — covers `INV-01`..`INV-04`
- Evidence: [MOTIF-DISCOVER-001-Evidence.md](../../../docs/Evidence/MOTIF-DISCOVER-001-Evidence.md)

## 8. References

1. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach*, 2nd ed., Ch. 2 (Finding Regulatory Motifs). Active Learning Publishers. Formula and worked example reproduced at https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string
2. fmicompbio. monaLisa `getKmerFreq` — observed vs expected k-mer frequencies and log2 enrichment. https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html
3. van Helden J, André B, Collado-Vides J. 1998. Extracting regulatory sites from the upstream region of yeast genes by computational analysis of oligonucleotide frequencies. *J Mol Biol* 281:827–842.
4. RSAT `oligo-analysis` source — https://raw.githubusercontent.com/rsa-tools/rsat-code/master/perl-scripts/oligo-analysis (Bernoulli `exp_freq *= residue_proba`, `exp_occ = exp_freq * sum_occurrences`, `sum_of_binomials` occ_P) and `perl-scripts/lib/RSA.disco.lib` `MultiTestCorrections` (occ_E, occ_sig); opened 2026-09-29.
