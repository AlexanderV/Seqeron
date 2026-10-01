# K-mer Euclidean Distance

| Field | Value |
|-------|-------|
| Algorithm Group | K-mer / Alignment-free sequence comparison |
| Test Unit ID | KMER-DIST-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 |

## 1. Overview

K-mer Euclidean distance is an alignment-free dissimilarity measure between two biological
sequences. Each sequence is summarized by the frequencies of its length-*k* substrings
(k-mers); the distance is the Euclidean (L2) distance between the two frequency vectors over
the union of observed k-mers [1][2]. It is exact (deterministic), requires no alignment, and
is used for fast whole-genome phylogeny, clustering, and database screening where alignment
is impractical [1][3]. Identical sequences have distance 0 and more dissimilar sequences have
larger values [1].

Since audit round 1 (WP2) the same unit also provides the other word-vector metrics of the
alignment-free literature through `KmerDistance(seq1, seq2, k, KmerDistanceMetric)` (count-based
squared Euclidean d_E [3][5], Manhattan, Chebyshev, Canberra, cosine, D2 [6]), the exact k-mer
**Jaccard** index and **Mash distance** [7][8] (§2.7), and **spaced-word** counts [4] (§2.8).

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Alignment-free methods map each sequence into a fixed-length vector of word (k-mer)
statistics and compare those vectors, sidestepping the cost and the recombination/shuffling
assumptions of alignment [3]. Word-composition methods map a sequence into a 4^k-dimensional
vector of k-word frequencies (for a 4-letter nucleotide alphabet), and a dissimilarity score
is obtained by a vector metric such as the Euclidean distance, Pearson correlation,
Kullback–Leibler discrepancy, or cosine distance [3].

### 2.2 Core Model

For a sequence *s* of length *L*, the number of overlapping k-mer windows is *L − k + 1*. The
k-mer **count** vector counts how many times each word *w* occurs; the k-mer **frequency**
vector normalizes each count by the total number of k-mers in the sequence (i.e. the sequence
length minus the k-mer length) [2]:

```
f_s(w) = count_s(w) / (L_s − k + 1)
```

Given the union *W* of k-mers occurring in either sequence *x* or *y*, the Euclidean distance is

```
d(x, y) = sqrt( Σ_{w ∈ W} ( f_x(w) − f_y(w) )² )
```

where a word absent from a sequence contributes a 0 component [1]. The difference between the
two word vectors "is very commonly computed by the Euclidean distance" [1], applied to the
relative-frequency vectors [4].

**Variant note (2026-09 review).** The classical word-count distance of Blaisdell (1986), as
reviewed by Vinga & Almeida (2003) [3], is the *squared* Euclidean distance on raw counts,
d_E = Σ (c_x(w) − c_y(w))²; for the Fig. 1 example this is 3 (√3 ≈ 1.7320508 unsquared). This
method implements the **non-squared Euclidean on relative frequencies** [2][4], giving √0.11.
The Lau et al. [2] parenthetical "(i.e. the sequence length minus the k-mer length)" is loose
wording: the total number of overlapping k-mers is L − k + 1 (scikit-bio
`Sequence.kmer_frequencies(relative=True, overlap=True)` uses `len(self) - k + 1`), which is
what the implementation uses. MUSCLE's k-mer distance (Edgar 2004) is a different measure
(fractional common k-mer count) and is not what this method computes.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | d(x, x) = 0 | Equal frequency vectors give a zero sum of squares; "identical sequences yield a distance of 0" [1] |
| INV-02 | d(x, y) = d(y, x) | (f_x − f_y)² = (f_y − f_x)²; Euclidean distance is a metric [3] |
| INV-03 | d(x, y) ≥ 0 | Square root of a sum of squares is non-negative [1][3] |
| INV-05 | `Euclidean` metric overload = `KmerDistance(seq1, seq2, k)` bit for bit; every dissimilarity member is 0 for identical inputs, symmetric | same loop, same union order |
| INV-06 | Jaccard ∈ [0, 1], symmetric, 1 for identical non-empty sets; Mash D ∈ [0, 1] | [7][8]; Mash boundary rules |
| INV-04 | Two sequences each consisting of a single distinct k-mer (frequency 1) with disjoint word sets give d = √2 | Frequency vectors are (1,0) and (0,1); √((1−0)²+(0−1)²)=√2 [2] |

### 2.5 Comparison with Related Methods

| Aspect | K-mer Euclidean (frequency) | D2 / cosine |
|--------|-----------------------------|-------------|
| Operates on | normalized frequency vectors [2][4] | raw counts (D2 is an uncentered correlation of counts) [3] |
| Range | ≥ 0, length-normalized | unbounded (D2) / [−1,1]-derived (cosine) |
| Metric | yes (Euclidean) [3] | D2 is a similarity statistic, not a metric [3] |

### 2.6 Metric variants (audit round 1, WP2)

`KmerDistance(seq1, seq2, k, KmerDistanceMetric)` and `KmerDistance(counts1, counts2, KmerDistanceMetric)`
(any two count tables: literal, Jellyfish `-C` canonical, or spaced-word) evaluate over the union of words
(absent = 0; words absent from both contribute 0 to every metric, so this equals the full |Σ|^k vector):

| Member | Definition | Vector | Source / reference implementation |
|--------|------------|--------|-----------------------------------|
| `Euclidean` (default; = `KmerDistance(seq1, seq2, k)` bit for bit) | √Σ(f_x − f_y)² | frequencies f = c/(L−k+1) | [1][2][3]; alfpy `euclid_norm` on `Freqs`; scipy `euclidean` |
| `SquaredEuclideanCounts` | Σ(c_x − c_y)² | raw counts | Blaisdell 1986 d_E [5], reviewed in [3]; alfpy `euclid_squared` on `Counts`; scipy `sqeuclidean` |
| `Manhattan` | Σ\|f_x − f_y\| | frequencies | alfpy `manhattan` on `Freqs`; scipy `cityblock` |
| `Chebyshev` | max\|f_x − f_y\| | frequencies | alfpy `chebyshev`; scipy `chebyshev` |
| `Canberra` | Σ\|f_x − f_y\|/(f_x + f_y), 0/0 terms omitted | frequencies | alfpy `canberra`; scipy `canberra` |
| `Cosine` | 1 − c_x·c_y/(‖c_x‖‖c_y‖), clipped to [0, 2] | counts (scale-invariant) | [3]; scipy `cosine` (clip as scipy); zero vector → similarity 0 → distance 1 |
| `D2` | Σ c_x(w)·c_y(w) | raw counts | Torney et al. 1990; Lippert et al. 2005; Reinert et al. 2009 [6]; [3] — a **similarity**, not a metric |

Counts vs frequencies follow the sources: d_E and D2 are defined on counts [3][5][6]; the L1/L∞/Canberra
members use the same relative-frequency vectors as the default (alfpy `Freqs`, total = L − k + 1), so they are
length-normalized like it. alfpy's own `word_d2` module is a squared Euclidean summed over k, *not* Torney's D2.

### 2.7 Exact k-mer Jaccard index and Mash distance

`JaccardSimilarity(a, b, k[, KmerCountingOptions])` = |K(a) ∩ K(b)| / |K(a) ∪ K(b)| over the **distinct** k-mer
sets (Jaccard 1901/1912 [7]; the quantity Mash estimates by MinHash, Ondov et al. 2016 eq. 1 [8]); fraction in
[0, 1], case-insensitive, both sets empty → 0 (the convention of `GenomicAnalyzer.CalculateSimilarity` and
`ComparativeGenomics`; DUP_MAP §16 canonical). With `Canonical = true` the k-mers are those of Mash and sourmash
(upper-cased, windows with a non-ACGT base skipped, key min(w, RC(w)) — Mash `Sketch.cpp` `addMinHashes`), so
the value is the exact Jaccard that `mash dist` / sourmash estimate; `AcgtOnly = true` matches `mash sketch -n`.

`MashDistance(a, b, k, options)` = −(1/k)·ln(2J/(1+J)) (Ondov 2016 eq. 4 [8]) with the boundary rules of Mash
2.x `CommandDistance.cpp`: shared = union (incl. both empty) → 0; no shared k-mer → 1; values above 1 capped
at 1. `MashDistanceFromJaccard(J, k)` is the bare conversion (J = 1 → 0, J = 0 → 1, cap 1). The method uses the
exact sets, so it equals `mash dist -s s` whenever s ≥ |K(a) ∪ K(b)|.

### 2.8 Spaced words

`CountSpacedWords(sequence, pattern)` (Leimeister et al. 2014 [4]): a pattern P ∈ {0,1}^ℓ with P[1] = P[ℓ] = 1
('1' = match, '0' = don't care; weight k = number of '1's). The spaced word at window i is the string of
sequence[i + j] over the match positions j; counts sum to L − ℓ + 1; the all-'1' pattern equals `CountKmers`.
Spaced-word frequency vectors are compared with `KmerDistance(counts1, counts2, metric)` (Leimeister et al. use
the Euclidean distance of the relative frequencies). Multi-pattern averaging (the `spaced` program) is a
caller-side loop over patterns.

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| seq1 | string | required | First sequence | nucleotide/protein text; null/empty allowed (empty vector) |
| seq2 | string | required | Second sequence | same as seq1 |
| k | int | required | K-mer length | must be > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | Non-negative Euclidean distance between the two frequency vectors; 0 for identical or both-empty inputs |

### 3.3 Preconditions and Validation

- `k ≤ 0` throws `ArgumentOutOfRangeException`.
- Inputs are upper-cased before counting, so comparison is case-insensitive (ASM-01).
- A null/empty sequence, or one shorter than *k*, produces an empty frequency vector, treated
  as the zero vector (ASM-02); the distance then equals the L2 norm of the other sequence's
  frequency vector, and 0 when both are empty.
- No alphabet restriction is enforced; any character may form a k-mer.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate `k > 0`.
2. Compute the normalized k-mer frequency vector of each sequence (count ÷ total windows) [2].
3. Form the union of k-mers present in either vector.
4. Sum the squared per-word frequency differences over the union (absent words = 0).
5. Return the square root of that sum [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| KmerDistance | O(n + m) | O(u) | n, m = sequence lengths; u = number of distinct k-mers in the union; one linear pass per sequence to count, one pass over the union |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.KmerDistance(string seq1, string seq2, int k)`: returns the Euclidean distance between the two normalized k-mer frequency vectors.
- `KmerAnalyzer.GetKmerFrequencies(string, int)`: builds the per-sequence frequency vector (count ÷ sum of counts).
- `KmerAnalyzer.CountKmers(string, int)`: underlying k-mer counter (upper-cases input; throws for k ≤ 0).
- `KmerAnalyzer.KmerDistance(string, string, int, KmerDistanceMetric)` / `KmerDistance(IReadOnlyDictionary<string,int>, IReadOnlyDictionary<string,int>, KmerDistanceMetric)`: the metric variants of §2.6 (single word-vector loop; the legacy method delegates to it).
- `KmerAnalyzer.JaccardSimilarity(string, string, int[, KmerCountingOptions])`, `MashDistance(string, string, int, KmerCountingOptions)`, `MashDistanceFromJaccard(double, int)`: §2.7 (over `DistinctKmers`).
- `KmerAnalyzer.CountSpacedWords(string, string)`: §2.8.
- MCP: `kmer_distance` (Analysis) optional `metric`; `kmer_jaccard` (Analysis) with optional `canonical` / `acgtOnly`.

### 5.2 Current Behavior

- Uses the **frequency** variant: `GetKmerFrequencies` divides each count by the sum of counts,
  which equals the number of k-mer windows (L − k + 1) for inputs whose every position forms a
  k-mer, matching the normalization in [2].
- Distance is taken over the union of observed k-mers only (sparse representation), which is
  equivalent to the full 4^k-dimensional vector since unobserved words have frequency 0 on
  both sides and contribute nothing.
- **Suffix tree not used:** this is not a substring-search / occurrence-enumeration task; it is
  a single linear k-mer counting pass per sequence followed by a vector difference. A suffix
  tree would add construction overhead without changing the O(n + m) cost, so the dictionary
  counter is the correct structure.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Word-vector model and union over both sequences with 0 components for absent words [1].
- Frequency normalization: count ÷ (L − k + 1) [2].
- Euclidean distance √(Σ (f_x − f_y)²) over the frequency vectors [1][4].
- Identity property d(x, x) = 0 [1].

**Intentionally simplified:**

- (none)

**Implemented in audit round 1 (WP2):**

- Count-based squared d_E (Blaisdell 1986 [5]), Manhattan, Chebyshev, Canberra, cosine and D2 [6] — `KmerDistanceMetric` (§2.6).
- Exact k-mer Jaccard and Mash distance [7][8] (§2.7); spaced-word counts [4] (§2.8).

**Not implemented:**

- Background-corrected D2* / D2S (Reinert et al. 2009 [6]; Wan et al. 2010): they need a Markov background model
  of word probabilities; **users should rely on:** `D2` (raw) and the distances above.
- MinHash *sketching* (Mash/sourmash estimate J from a bottom-s sketch); this unit computes the exact sets, which
  is the quantity the sketch estimates.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Case-insensitive (upper-casing) | Assumption | Mixed-case inputs treated as identical k-mers | accepted | ASM-01; benign for canonical upper-case input |
| 2 | L < k ⇒ empty (zero) vector | Assumption | Defines distance for too-short inputs not covered by sources | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical sequences | 0 | INV-01 [1] |
| k ≤ 0 | ArgumentOutOfRangeException | Input validation |
| One sequence shorter than k | distance = L2 norm of the other's frequency vector | ASM-02 (empty = zero vector) |
| Both sequences empty | 0 | ASM-02 |
| Disjoint single-k-mer sequences | √2 | INV-04 |
| Lower-case vs upper-case | same as upper-case result | ASM-01 |

### 6.2 Limitations

- Frequency normalization removes overall length information, so two sequences with the same
  relative composition but different lengths can have distance 0.
- The default uses contiguous literal k-mers; spaced words (§2.8) and canonical (strand-collapsed)
  count tables (`KmerDistance(counts1, counts2, metric)`) are available. No statistical background
  correction (D2*, D2S). The raw Euclidean value is a
  dissimilarity, not a calibrated phylogenetic distance [2].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
double d = KmerAnalyzer.KmerDistance("ATGTGTG", "CATGTG", k: 3);
// d ≈ 0.33166247903553997
```

**Numerical walk-through (Zielezinski et al. 2017, Fig. 1):**

For x = "ATGTGTG" and y = "CATGTG", k = 3, the union of words is {ATG, CAT, GTG, TGT} with
counts c_x = (1, 0, 2, 2) and c_y = (1, 1, 1, 1) [1]. Normalizing by the window counts (5 for
x, 4 for y) gives frequencies f_x = (0.2, 0, 0.4, 0.4) and f_y = (0.25, 0.25, 0.25, 0.25).
The squared differences are (0.0025, 0.0625, 0.0225, 0.0225), summing to 0.11, so
d = √0.11 ≈ 0.33166247903553997.

### 7.2 Reference cross-check (audit round 1, WP2)

Reference implementations executed 2026-10-01: scipy 1.17.1 `scipy.spatial.distance` on scikit-bio 0.7.4
`Sequence.kmer_frequencies(overlap=True, relative=False/True)` over the union; alfpy 1.0.6 (`word_distance`,
source tarball from PyPI); Python set Jaccard; sourmash 4.9.4 `MinHash(n=0, ksize=k, scaled=1)` + `add_sequence(force=True)`;
the **Mash 2.3 binary** (`mash dist -k k -s 100000 [-n]`). C# equals every value (tolerance 1e-12; Mash output is printed to 6 digits).

| Pair (k) | Euclid (f) | sqEuclid (c) | Manhattan | Chebyshev | Canberra | Cosine | D2 |
|----------|-----------|--------------|-----------|-----------|----------|--------|----|
| ATGTGTG / CATGTG (3) | 0.33166247903554 | 3 | 0.6 | 0.25 | 1.5726495726495728 | 0.16666666666666663 | 5 |
| ACGTTGCAACGGT / ACGTAGCATCGGTA (2) | 0.26600633232367216 | 11 | 0.7692307692307692 | 0.15384615384615385 | 7.568421052631578 | 0.29704050843336227 | 13 |
| GATTACAGATTACA / GATTACCGATTTCA (3) | 0.31180478223116176 | 14 | 1.0 | 0.16666666666666666 | 9.666666666666666 | 0.3603978509331687 | 12 |
| acgtNNacgtacgRtTTGCAnA / ACGTACGTTTGCAAA (3) | 0.2179788817739193 | 13 | 0.9 | 0.07692307692307693 | 13.49750671269659 | 0.26664120237743094 | 16 |
| R1 / R2 120 nt (5) | 0.09829098492233948 | 130 | 1.0517241379310347 | 0.017241379310344827 | 114.66666666666666 | 0.5199846392626807 | 60 |

| Pair (k) | J literal | J ACGT-only | J canonical (sourmash) | Mash D (`mash dist`) | Mash D `-n` |
|----------|-----------|-------------|------------------------|----------------------|-------------|
| ATGTGTG / CATGTG (3) | 0.75 | 0.75 | 1.0 | 0 (3/3) | 0.05138355994241945 (0.0513836, 3/4) |
| ACGTTGCAACGGT / ACGTAGCATCGGTA (2) | 0.46153846153846156 | 0.46153846153846156 | 0.5 | 0.20273255405408222 (0.202733, 5/10) | 0.2297661646892201 (6/13) |
| GATTACAGATTACA / GATTACCGATTTCA (3) | 0.3076923076923077 | = | = | 0.2512572674587934 (4/13) | = |
| acgtNNacgtacgRtTTGCAnA / ACGTACGTTTGCAAA (3) | 0.4 | 0.7272727272727273 | 0.8333333333333334 | 0.03177005993477496 (5/6) | 0.05728341897555305 (8/11) |
| R1 / R2 (5) | 0.3273809523809524 | = | 0.3670886075949367 | 0.12433764331556005 (0.124338, 58/158) | 0.1413382811335405 (55/168) |
| R1 / R2 (11) | 0 | 0 | 0 | 1 (0/220) | 1 |

R1/R2 are the 120-nt pair in the test file (Python `random.seed(2026)`, every 9th base substituted). Spaced words
(Python replica of the [4] definition): ATGTGTG/`101` → AG:1 GG:2 TT:2; ACGTTGCAACGGT/`1101` → 9 words, CGT:2;
GATTACAGATTACA/`11011` → ATAC:2 GATA:2 TTCA:2 ACGA/AGTT/CAAT/TAAG:1; R1/`1100111` → 108 distinct, 114 total;
`1101` frequency-Euclidean GATTACAGATTACA vs GATTACCGATTTCA 0.3149183286488868 (sqEuclid counts 12).
The `spaced` reference program (spaced.gobics.de) was not reachable (host not allow-listed; no GitHub mirror found).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_KmerDistance_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_KmerDistance_Tests.cs) — covers `INV-01`–`INV-04`
- Tests: [KmerAnalyzer_DistanceMetrics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_DistanceMetrics_Tests.cs) — §7.2 values, `INV-05`/`INV-06`, conventions, spaced words
- Evidence: [KMER-DIST-001-Evidence.md](../../../docs/Evidence/KMER-DIST-001-Evidence.md)

## 8. References

1. Zielezinski A, Vinga S, Almeida J, Karlowski WM. 2017. Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology 18:186. https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ (DOI: 10.1186/s13059-017-1319-7)
2. Lau AK, et al. 2022. Interpreting alignment-free sequence comparison: what makes a score a good score? NAR Genomics and Bioinformatics. https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
3. Vinga S, Almeida J. 2003. Alignment-free sequence comparison—a review. Bioinformatics 19(4):513–523. https://academic.oup.com/bioinformatics/article/19/4/513/218529 (DOI: 10.1093/bioinformatics/btg005)
4. Leimeister C-A, Boden M, Horwege S, Lindner S, Morgenstern B. 2014. Fast alignment-free sequence comparison using spaced-word frequencies. Bioinformatics 30(14):1991–1999. https://academic.oup.com/bioinformatics/article/30/14/1991/2391234 (earlier revisions of this doc listed the second author, Boden, as first author)
5. Blaisdell BE. 1986. A measure of the similarity of sets of sequences not requiring sequence alignment. PNAS 83:5155–5159. doi:10.1073/pnas.83.14.5155 (cited by alfpy `pwdist_euclid_squared`)
6. Torney DC, Burks C, Davison D, Sirotkin KM. 1990. Computation of d2: a measure of sequence dissimilarity. In: Computers and DNA, SFI Studies XII:109–125; Lippert RA, Huang H, Waterman MS. 2005. PNAS 102:13980–13989; Reinert G, Chew D, Sun F, Waterman MS. 2009. Alignment-free sequence comparison (I): statistics and power. J Comput Biol 16:1615–1634.
7. Jaccard P. 1901. Bull Soc Vaudoise Sci Nat 37:547–579; Jaccard P. 1912. The distribution of the flora in the alpine zone. New Phytol 11:37–50.
8. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. 2016. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology 17:132. Source: github.com/marbl/Mash `src/mash/CommandDistance.cpp`, `src/mash/Sketch.cpp` (raw.githubusercontent.com, opened 2026-10-01).
9. Zielezinski A, Girgis HZ, Bernard G, et al. 2019. Benchmarking of alignment-free sequence comparison methods. Genome Biology 20:144 (alfpy reference implementation; PyPI alfpy 1.0.6).
