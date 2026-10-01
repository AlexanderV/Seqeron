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
| `JensenShannon` (audit round 2) | ½Σ f_x log₂(f_x/m) + ½Σ f_y log₂(f_y/m), m = ½(f_x + f_y) | frequencies | Lin 1991 [12]; the JS measure of [4] and `spaced -d JS`; scipy `jensenshannon(p, q, base=2)²`; zero vector vs non-empty → ½ |
| `EuclideanCounts` (audit round 2) | √Σ(c_x − c_y)² | raw counts | the per-pattern value of `spaced -d EU` (§7.5) |

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

`ContainmentIndex(a, b, k[, options])` (audit round 2, WP7) = |K(a) ∩ K(b)| / |K(a)| — the containment index of
Koslicki & Zabeti 2019 [14] and of sourmash `compare --containment` (exact for `scaled=1` sketches); asymmetric,
C(a, b)·|K(a)| = C(b, a)·|K(b)| = |K(a) ∩ K(b)|; K(a) empty → 0. It shares the set-intersection helper with
`JaccardSimilarity`/`MashDistance`.

### 2.8 Spaced words

`CountSpacedWords(sequence, pattern)` (Leimeister et al. 2014 [4]): a pattern P ∈ {0,1}^ℓ with P[1] = P[ℓ] = 1
('1' = match, '0' = don't care; weight k = number of '1's). The spaced word at window i is the string of
sequence[i + j] over the match positions j; counts sum to L − ℓ + 1; the all-'1' pattern equals `CountKmers`.
Spaced-word frequency vectors are compared with `KmerDistance(counts1, counts2, metric)` (Leimeister et al. use
the Euclidean distance of the relative frequencies).

`SpacedWordDistance(seq1, seq2, patterns, metric = Euclidean)` (audit round 2) is the multiple-pattern distance of
[4]: "the distance between two sequences is the average of the distances based on the individual patterns",
d_P = (1/m) Σᵢ d(N_{Pᵢ}(S₁), N_{Pᵢ}(S₂)), with all patterns of the same weight (as the pattern sets of [4] and of
the `spaced` program). [4] applies the Euclidean distance and the Jensen–Shannon distance to relative spaced-word
frequencies (`Euclidean`, `JensenShannon`); every count-table metric of §2.6 except D2*/D2S is accepted. Pattern
sets are passed by the caller: the `spaced` program generates them by randomised optimisation (`variance::Improve`,
rasbhari-style), which is not deterministic and is not reproduced.

`CountSpacedWords(sequence, pattern, KmerCountingOptions)` and
`SpacedWordDistance(seq1, seq2, patterns, metric, KmerCountingOptions, bothStrands = false)` (audit round 2, WP8) add
the two remaining `spaced` 1.2.0 conventions as options (§7.7): `AcgtOnly` drops a window whose symbol at a **match**
position is not A/C/G/T (don't-care positions are ignored), with frequencies count ÷ W, W = L − ℓ + 1 windows (the
dropped windows stay in W, as in `spaced`); `bothStrands` is `spaced`'s default reverse-complement mode, in which
`seq1` (the first input record) is counted on both strands (forward + reverse-complement strand, total 2·W₁) and
`seq2` on its forward strand only. That mode is asymmetric in the argument order and is a convention of the tool, not
of [4]. `Canonical` is rejected for spaced words: the reverse strand reads a window with the mirrored pattern, so
min(word, RC(word)) is not strand-independent for an asymmetric pattern, and no source defines it.

**Spaced-faithful reader (audit round 3, WP9).** `spaced`'s FASTA reader (`spacedDNA`) deletes every character for
which `isalpha` is false (gap `-`, `*`, digits, blanks) before windowing; only letters are stored (non-ACGT letters as
N). `SpacedWordDistance(…, AcgtOnly, …)` — the spaced-faithful mode, with or without `bothStrands` — now reads both
sequences that way, so L and W = L − ℓ + 1 are those of the letters only: `ACG-TACGT` is windowed as `ACGTACGT`
(`spaced -r -d JS -f {1011, 1101}` against `ACGTTACGA` prints 0.57013316426; keeping the gap as an N window, as WP8
did, gave 0.25). The literal default and `CountSpacedWords` window the string as given.

**Evolutionary distance `spaced -d EV` (audit round 3, WP9).** `KmerDistanceMetric.SpacedEvolutionary`
(Morgenstern, Zhu, Horwege & Leimeister 2015 [16]) in `SpacedWordDistance`, implemented exactly as `spaced` 1.2.0
`sort.h` (EV branch) [13]: N = Σ_P Σ_w min(c₂(w), c₁(w)) is the number of spaced-word matches summed over all
patterns (c₂ = forward counts of the second record, c₁ = the first record's forward, or forward + reverse-strand,
counts; words reading N excluded); with read lengths L₁, L₂ (N letters included), ℓ the common pattern length and w the
weight, m = min(L₁, L₂) − ℓ + 1, M = max(L₁, L₂) − ℓ + 1, f(a) = count(a)/L, q = Σ_a f₁(a)f₂(a) (strand-averaged
f(a) = ½(f(a) + f(ā)) in reverse-complement mode), V = N/(|P|·m) − s·M·q^w with s = 2 in reverse-complement mode and
1 with `-r`. V ≥ 0: p = V^(1/w) (estimated match probability per site) and d = −¾ ln(4p/3 − 1/3) (Jukes–Cantor);
V < 0: `spaced` prints 1.2 (`SpacedEvolutionarySaturationDistance`). As in the program, p < ¼ gives NaN (`-nan`)
and p = ¼ gives +∞. EV is always read in the spaced-faithful mode; all patterns must have the same length (spaced keeps
one weight and one don't-care count, `ell = dontCare + weight − 1`) and both read sequences at least ℓ letters
(shorter input makes `spaced` index outside the record). `KmerDistance(seq1, seq2, k, SpacedEvolutionary, 0,
bothStrands)` is the same with the contiguous pattern 1^k. The count-table overload rejects it (it needs lengths and
base composition).

### 2.9 Background-adjusted D2* and D2S (audit round 1, WP4)

`BackgroundAdjustedD2(seq1, seq2, k, markovOrder = 0)` and the metrics `KmerDistanceMetric.D2Star` /
`D2Shepherd` (Reinert, Chew, Sun & Waterman 2009 [6]; Wan et al. 2010 [10]; dissimilarities as Song et al. 2014
[10]; CAFE `D2star` / `D2shepp` [11]). X_w, Y_w are the single-strand counts over the ACGT windows (Jellyfish
convention, as CAFE), n̄ = Σ X_w, m̄ = Σ Y_w. Each sequence gets its own order-r Markov background, fitted by
maximum likelihood on its ACGT r-mer and (r+1)-mer counts:
p̂(w) = N(w₁..w_r)/Σ N(r-mers) · Π_{i>r} N(w_{i−r}..w_i)/Σ_a N(w_{i−r}..w_{i−1}a) (r = 0: the product of letter
frequencies). Expected counts E_X = n̄·p̂_X(w), E_Y = m̄·p̂_Y(w); centred counts X̃ = X − E_X, Ỹ = Y − E_Y. Sums run
over all 4^k words (absent words contribute through −E). `markovOrder = −1` selects each sequence's order by
BIC(r) = −2 ln L̂_r + 3·4^r·ln N_r (N_r = number of ACGT (r+1)-mers; Schwarz 1978, Katz 1981; the criterion of
CAFE `-M -1`, orders 0..min(k − 1, 10)):

- D2* = Σ X̃Ỹ/√(E_X E_Y) (E_X·E_Y = 0 omitted); d2* = ½(1 − D2*/√(Σ X̃²/E_X · Σ Ỹ²/E_Y)).
- D2S = Σ X̃Ỹ/√(X̃² + Ỹ²) (X̃ = Ỹ = 0 omitted); d2S = ½(1 − D2S/√(Σ X̃²/√(X̃²+Ỹ²) · Σ Ỹ²/√(X̃²+Ỹ²))).

With a common background p (Reinert 2009's setting) √(E_X E_Y) = √(n̄ m̄)·p_w, the published D2* denominator.
d2*, d2S ∈ [0, 1]; identical sequences → 0. k ≤ 12 (`MaxBackgroundAdjustedK`: 4^12 words are enumerated), 0 ≤ r < k.
NaN when a normaliser is 0 (e.g. a homopolymer at order 0, whose counts equal their expectation). Audit round 3 (WP9):
d2* and d2S are clamped to [0, 1] — Cauchy–Schwarz bounds them there, and rounding gave −1.1102230246251565e-16 for
identical sequences (S1/S1, k = 2, r = 0); NaN passes through, the raw D2*/D2S are not clamped. The counts
overload `KmerDistance(counts1, counts2, metric)` rejects D2*/D2S (a background needs the sequence).

**Both strands (audit round 2, CAFE `-R`).** `BackgroundAdjustedD2(seq1, seq2, k, markovOrder, bothStrands)` and
`KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands)`. CAFE's semantics, read from the source
(`kmer.cpp`): `KmerModel::load` adds every word's count to its reverse complement's entry as well, so
X^R(w) = X(w) + X(RC(w)) (palindromes count twice; Σ_w X^R = 2n̄), and `KmerProbEnsembDelegate::getKmerlogProb`
returns log(½(p(w) + p(RC(w)))), where p(RC(w)) is the probability of the reverse-complement word under the chain
fitted to the given strand (the reverse-complement `KmerProbDelegate` walks w left to right with complemented,
reversed contexts). With total 2n̄ the expected count is E^R(w) = n̄·(p̂(w) + p̂(RC(w))) = E[X(w) + X(RC(w))]. The sums
still run over all 4^k words (a non-palindromic pair contributes twice). The Markov order chosen by BIC uses the
given strand (CAFE `getEstMarkovOrder` builds single-strand models). As in single-strand mode, p̂ is the maximum-likelihood
chain on the sequence, not CAFE's prefix-marginal estimator (§7.4). The incremental odometer carries a second prefix
product for p̂(RC(w)): at depth d ≥ r the window w_{d−r..d} reverse-complemented gives the transition
P(c(w_{d−r}) | c(w_d)…c(w_{d−r+1})), and at d = k − 1 the initial r-mer c(w_{k−1})…c(w_{k−r}) is multiplied in.

**Memory (audit round 2).** The Markov tables are dense arrays only for r + 1 ≤ 8 (4^8 doubles); above that
they are dictionaries holding the r-mers / (r+1)-mers present, so k = 12, r = 11 no longer allocates a 4^12-double
(134 MB) array per sequence (locked by an allocation test).

### 2.10 MinHash sketches (Mash) (audit round 2, WP7)

`CreateMinHashSketch(sequence | records, k, sketchSize = 1000, canonical = true, seed = 42)` builds the bottom-s
sketch of `mash sketch` 2.3 [8]: the k-mer set is `DistinctKmers` with `Canonical` (or `AcgtOnly` for `-n`) —
upper-cased, non-ACGT windows skipped, min(w, RC(w)) by `memcmp` — so no second canonicalisation exists; each k-mer's
k ASCII bytes are hashed with MurmurHash3_x64_128 (`MurmurHash3X64_128`, Appleby; Mash's bundled `MurmurHash3.cpp`,
same as `mmh3.hash64`), keeping h1 when 4^k > 2^32 (k ≥ 17, Mash `use64`) and its low 32 bits otherwise
(`hash.cpp`); the sketch holds the s smallest **distinct** hash values ascending (`MinHashHeap`). With several
records (contigs of one FASTA), no k-mer spans two records, and records shorter than k are skipped; the sketch length
is the summed length of the remaining records (all symbols, N included; Mash `sketchFile`). k ∈ 1..32 (`mash -k`).

`CompareMinHashSketches(ref, query)` is `mash dist`'s `compareSketches`: merge the sorted sketches until
s = min(s_ref, s_query) union hashes have been visited (if one sketch runs out first, the rest of the other is
added and the total capped at s); x = hashes in both; J = x / denominator (Mash's "x/s" column); distance by
`MashDistanceFromJaccard` (x = denominator → 0, x = 0 → 1). The p-value (`MashPValue`, Mash `pValue`) is
P(X ≥ x) for X ~ Binomial(denominator, r), r = p₁p₂/(p₁ + p₂ − p₁p₂), pᵢ = 1/(1 + 4^k/lengthᵢ) — Mash's
`gsl_cdf_binomial_Q(x − 1, r, denominator)`, here `StatisticsHelper.BinomialUpperTail` (log-space, Loader 2000);
x = 0 → 1. Both sketches empty: denominator 0, J = 0 (exact-Jaccard convention; Mash prints nan), distance 0,
p-value 1. Sketches with different k, seed or canonical mode are rejected (Mash skips them). Audit round 3 (WP9):
`CompareMinHashSketches` also rejects malformed sketches (hashes not strictly ascending, more hashes than
`SketchSize`, `SketchSize` < 1, negative length) with `ArgumentException`; `MinHashSketch.FromHashes` builds a
well-formed sketch from arbitrary hash values. `MashPValue` rejects inputs no comparison can produce: x > s (Mash's
`gsl_cdf_binomial_Q(x − 1, r, s)` silently returns 0 for x − 1 ≥ s, GSL `cdf/binomial.c`) and x ≥ 1 with a length
of 0 (Mash would compute r = 0/0 = NaN; an empty set has no hashes); earlier, the first case returned 0 and the second
threw `ArgumentOutOfRangeException` from inside the binomial tail. Audit round 4 (WP10): k outside Mash's 1..32
(`Command.cpp`: `Option(Option::Integer, "k", …, "21", 1, 32)`) is rejected with `ArgumentOutOfRangeException("k")`;
before, k = 600 gave 4^600 = ∞, r = NaN and an undocumented `ArgumentOutOfRangeException("p")` from the binomial tail.

### 2.11 FracMinHash (sourmash `scaled`) sketches (audit round 3, WP9)

`CreateFracMinHashSketch(sequence, k, scaled, canonical = true, seed = 42)` is sourmash 4.9.4
`MinHash(n=0, ksize=k, scaled=S).add_sequence(seq, force=True)` [17][18][19]: canonical k-mers (Rust `SeqToHashes`:
`std::cmp::min(kmer, krc)` on the upper-cased bytes, the same set as `DistinctKmers` with `Canonical`; windows with a
non-ACGT base are skipped, which is `force=True` — without it sourmash raises), hashed with `MurmurHash3X64_128`
(h1, seed 42; `_hash_murmur`), and **every** distinct hash h ≤ max_hash kept (Rust `add_hash`). max_hash =
`FracMinHashMaxHash(S)` = Rust `max_hash_for_scaled`: S = 1 → 2^64 − 1, else `(u64::MAX as f64 / S as f64) as u64`
(double division by 2^64, truncated). sourmash's Python helper `_get_max_hash_for_scaled` rounds instead; the two
agree while 2^64/S ≥ 2^53 and differ by one above (S = 7919: Rust 2329428472497733 = `MinHash(0, 21,
scaled=7919)._max_hash`, Python helper …734); the Rust value is the one applied. `canonical = false` hashes forward
k-mers (Mash `-n` style; sourmash has no such DNA mode).

`CompareFracMinHashSketches(a, b)` → `FracMinHashComparison`: x = |A ∩ B|, u = |A ∪ B|, Jaccard = x / max(1, u)
(Rust `KmerMinHash::jaccard`); `contained_by` (A in B) = x / (|A|·b), b = 1 − (1 − 1/S)^(|A|·S) (the bias factor of
Hera et al. 2023 [18], sourmash `minhash.py`), clamped to [0, 1], 0 for an empty A; `max_containment` uses
min(|A|, |B|). S = 1 gives b = 1, i.e. the exact `JaccardSimilarity`/`ContainmentIndex` (canonical). Sketches must
share k, seed, canonical mode and scaled unless `downsample` is set (§2.12; sourmash refuses otherwise unless asked to downsample; the downsampled
sketch at S′ ≥ S is the subset h ≤ max_hash(S′), i.e. the sketch built at S′). Malformed sketches (hashes not strictly
ascending, a hash above `MaxHash`, `MaxHash` ≠ max_hash(`Scaled`)) are rejected.

### 2.12 FracMinHash downsampling, abundance tracking, u32 `scaled` (audit round 4, WP10)

`scaled` is a `long` restricted to 1..4294967295 (sourmash `ScaledType` = u32; `MaxSourmashScaled`): S = 4294967295
gives max_hash 4294967297 = sourmash `MinHash(0, 21, scaled=4294967295)._max_hash`.

`DownsampleFracMinHash(sketch, S′)` = sourmash `MinHash.downsample(scaled=S′)` / Rust `downsample_scaled`: the hashes
h ≤ max_hash(S′) (same truncating `FracMinHashMaxHash`; Python's `downsample` converts its rounded max_hash back to S′
and the Rust constructor recomputes the truncated value — `downsample(scaled=7919)._max_hash` = 2329428472497733),
abundances kept; S′ < S is rejected (sourmash "new scaled … is lower than current sample scaled", Rust
`CannotUpsampleScaled`). It equals sketching the sequence at S′ (checked on sequence C: 1011 hashes at S = 10 → 15,
identical to the S = 1000 sketch).

`CompareFracMinHashSketches(a, b, downsample: true)` = sourmash `downsample=True`: sketches of different scaled are both
downsampled to max(S_a, S_b) and compared there. Jaccard (and angular similarity) equal `a.jaccard(b, downsample=True)`
(Rust `similarity` downsamples the smaller-scaled sketch). The containments equal `sourmash compare --containment` /
`--max-containment`, which downsamples all signatures to the common maximum scaled first (`commands.py`), i.e.
`a.downsample(scaled=S).contained_by(b.downsample(scaled=S))`. The Python methods `contained_by(…, downsample=True)` /
`max_containment(…, downsample=True)` only downsample inside `count_common` and keep `len(self)` and `self.scaled` of
the undownsampled sketch in the denominator, so when `self` has the smaller scaled they differ (sourmash 4.9.4, A at
S = 10 in B at S = 100, k = 21: method 0.029064039408866996, CLI / downsampled 0.3155080213903743); that mixed-scale
artefact is not reproduced.

`CreateFracMinHashSketch(…, trackAbundance: true)` = sourmash `track_abundance=True`: every k-mer occurrence adds 1 to
its hash (Rust `add_hash_with_abundance`), so `FracMinHashSketch.Abundances` holds the canonical counts of
`CountKmers(…, Canonical)` (summed if two k-mers collide on one hash). With both sketches tracking abundance,
`FracMinHashComparison.AngularSimilarity` = sourmash `angular_similarity` = `similarity(ignore_abundance=False)` (Rust
`KmerMinHash::angular_similarity`): cos = min(1, Σ_{h∈A∩B} a_h·b_h / (‖a‖·‖b‖)) with ‖a‖² = Σ over all of A's
abundances (u64), 0 when a norm is 0, similarity = 1 − 2·acos(cos)/π. `WeightedContainmentAInB` = sourmash
`contained_by_weighted` = Σ_{h∈A∩B} a_h / Σ_{h∈A} a_h (not bias-corrected; needs only A's abundances). With one flat
sketch, `AngularSimilarity` is null (sourmash raises `TypeError`) and `similarity` falls back to Jaccard.

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
| BackgroundAdjustedD2 | O(n + m + 4^k·k) | O(u + distinct (r+1)-mers) | counts + Markov fit per sequence, then one odometer pass over all 4^k words (prefix probabilities updated incrementally; O(k) per leaf for the code, O(r) per factor) |
| SpacedWordDistance | O(m·(n + n')·k) | O(distinct spaced words) | m patterns of weight k |
| CreateMinHashSketch | O(n·k + d log d) | O(d) | d distinct (canonical) k-mers, all hashed and sorted; the sketch keeps s |
| CompareMinHashSketches | O(s) | O(1) | one merge of two sorted sketches |
| ContainmentIndex | O(n + m) | O(u) | same set pass as `JaccardSimilarity` |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [KmerAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs)

- `KmerAnalyzer.KmerDistance(string seq1, string seq2, int k)`: returns the Euclidean distance between the two normalized k-mer frequency vectors.
- `KmerAnalyzer.GetKmerFrequencies(string, int)`: builds the per-sequence frequency vector (count ÷ sum of counts).
- `KmerAnalyzer.CountKmers(string, int)`: underlying k-mer counter (upper-cases input; throws for k ≤ 0).
- `KmerAnalyzer.KmerDistance(string, string, int, KmerDistanceMetric)` / `KmerDistance(IReadOnlyDictionary<string,int>, IReadOnlyDictionary<string,int>, KmerDistanceMetric)`: the metric variants of §2.6 (single word-vector loop; the legacy method delegates to it).
- `KmerAnalyzer.JaccardSimilarity(string, string, int[, KmerCountingOptions])`, `MashDistance(string, string, int, KmerCountingOptions)`, `MashDistanceFromJaccard(double, int)`: §2.7 (over `DistinctKmers`).
- `KmerAnalyzer.CountSpacedWords(string, string[, KmerCountingOptions])`: §2.8.
- `KmerAnalyzer.SpacedWordDistance(string, string, IReadOnlyList<string>, KmerDistanceMetric, KmerCountingOptions, bool bothStrands = false)`: `spaced` N-word rule and reverse-complement mode (§2.8, §7.7; audit round 2, WP8).
- `KmerAnalyzer.BackgroundAdjustedD2(string, string, int, int)` → `D2StarStatistics(D2Star, D2Shepherd, D2StarDistance, D2ShepherdDistance)`; `KmerDistance(seq1, seq2, k, metric, markovOrder)`; metrics `D2Star` / `D2Shepherd`; `ParseDistanceMetric(string)` (MCP metric names): §2.9.
- `KmerAnalyzer.BackgroundAdjustedD2(string, string, int, int, bool bothStrands)`, `KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands)` (CAFE `-R`), `SpacedWordDistance(string, string, IReadOnlyList<string>, KmerDistanceMetric)`, metrics `JensenShannon` / `EuclideanCounts` (audit round 2, WP6).
- `KmerAnalyzer.ContainmentIndex(string, string, int[, KmerCountingOptions])` (§2.7), `CreateMinHashSketch` → `MinHashSketch`, `CompareMinHashSketches` → `MashComparison(SharedHashes, Denominator, Jaccard, Distance, PValue)`, `MashPValue`, `MurmurHash3X64_128` (§2.10; audit round 2, WP7).
- Audit round 3 (WP9): metric `KmerDistanceMetric.SpacedEvolutionary` (`spaced -d EV`, §2.8) in `SpacedWordDistance` and `KmerDistance(seq1, seq2, k, metric, markovOrder, bothStrands)`, `SpacedEvolutionarySaturationDistance`; spaced reader in the `AcgtOnly` spaced path (§2.8); `FracMinHashMaxHash`, `CreateFracMinHashSketch` → `FracMinHashSketch`, `CompareFracMinHashSketches` → `FracMinHashComparison` (§2.11); `MinHashSketch.FromHashes`.
- Audit round 4 (WP10): `DownsampleFracMinHash`, `CompareFracMinHashSketches(a, b, downsample)`, `CreateFracMinHashSketch(…, trackAbundance)` with `FracMinHashSketch.Abundances`, `FracMinHashComparison.AngularSimilarity` / `WeightedContainmentAInB` / `WeightedContainmentBInA`, `long` scaled in 1..4294967295 (`MaxSourmashScaled`) (§2.12); `MashPValue` k limited to Mash's 1..32.
- MCP: `kmer_distance` (Analysis and Sequence servers) optional `metric` (incl. `d2star`, `d2shepherd`, `jensen_shannon`, `euclidean_counts`), `markovOrder` and `bothStrands`; `kmer_jaccard` (Analysis) with optional `canonical` / `acgtOnly` and `sketchSize` (0 = exact; s > 0 = Mash sketch estimate + `sharedHashes`/`sketchDenominator`/`pValue`), always returning both exact containment indices; `kmer_d2_statistics` (raw D2*/D2S, d2*/d2S, orders, BIC) and `spaced_word_distance` (Analysis; optional `acgtOnly` and `bothStrands`, WP8). WP9: metric `ev` on `spaced_word_distance` and both `kmer_distance` tools; `kmer_jaccard` optional `scaled` (sourmash FracMinHash; exclusive with `sketchSize`) and output `maxContainment`. WP10: `kmer_jaccard` `scaled` accepts 1..4294967295 and requires `canonical = true` (sourmash DNA hashing is canonical; `canonical = false` was silently ignored), optional `trackAbundance` → output `angularSimilarity`.

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

**Implemented in audit round 1 (WP4):**

- Background-adjusted D2* / D2S and d2* / d2S [6][10][11] with an order-r Markov background per sequence (§2.9). `markovOrder = -1` (`AutoMarkovOrder`) chooses each sequence's order in [0, min(k − 1, 10)] by BIC (`MarkovOrderBic`, `SelectMarkovOrder`; Schwarz 1978, Katz 1981 — CAFE `-M -1`). CAFE's both-strand mode (`-R`) was added in audit round 2 (§2.9).

**Implemented in audit round 2 (WP6):**

- CAFE `-R` both-strand D2*/D2S (§2.9, §7.5); sparse Markov tables for high orders; null = empty for every `KmerDistance` metric.
- Multiple-pattern spaced-word distance [4] (`SpacedWordDistance`, §2.8) with the Jensen–Shannon (`JensenShannon`) and count-Euclidean (`EuclideanCounts`) metrics; cross-checked against the `spaced` 1.2.0 program (§7.5).

**Implemented in audit round 2 (WP7):**

- Mash bottom-s MinHash sketches, `mash dist` comparison (x/s, Jaccard estimate, distance, binomial p-value) and
  MurmurHash3_x64_128 (§2.10), = the Mash 2.3 binary on 63 runs (§7.6); exact containment index [14] (§2.7),
  = sourmash 4.9.4 `scaled=1` `contained_by`.

**Implemented in audit round 2 (WP8):**

- `spaced` 1.2.0's N-word rule (`AcgtOnly`) and default reverse-complement mode (`bothStrands`) for
  `SpacedWordDistance`; = the `spaced` binary with and without `-r` on 36 runs (§7.7).

**Implemented in audit round 3 (WP9):**

- `spaced -d EV` evolutionary distance [16] (`SpacedEvolutionary`, §2.8) and `spaced`'s reader (non-letters deleted)
  in the `AcgtOnly` spaced path; = the `spaced` binary on 248 runs with and without `-r` (§7.8).
- sourmash FracMinHash (`scaled`) sketches with Jaccard, bias-corrected `contained_by` and `max_containment`
  [17][18] (§2.11); = sourmash 4.9.4 on 40 comparisons (k = 21/31 × scaled 1/10/100/1000 × 5 pairs, §7.8).
- Validation: `MashPValue` (x > s, x ≥ 1 with length 0), malformed sketches in `CompareMinHashSketches`; d2*/d2S
  clamped to [0, 1] (§2.9, §2.10).

**Not implemented:**

- Mash `screen` (containment score with multiplicities, a different tool) and Mash's 32-bit `ARCH_32` hash variant
  (two MurmurHash3_x86_32 calls; only in 32-bit builds); sourmash abundance tracking (angular similarity) and
  protein/Dayhoff/HP hash functions (not k-mer DNA set measures).
- `spaced`'s randomised pattern-set generation (`variance::Improve`; not deterministic, so not reproducible).

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
  count tables (`KmerDistance(counts1, counts2, metric)`) are available. Background-corrected d2* / d2S are
  available (§2.9) for k ≤ 12. The raw Euclidean value is a
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
The `spaced` reference program (spaced.gobics.de) was not reachable in round 1; round 2 obtained it from the
Ubuntu archive (§7.5).

### 7.4 D2* / D2S reference cross-check (audit round 1, WP4)

Sources opened: CAFE source (github.com/younglululu/CAFE, cloned 2026-10-01): `dist_model.cpp` (D2starStrategy /
D2sheppStrategy `dealWithQuad`, `getDist` = 0.5·(1 − num/(√Σ·√Σ))), `kmer.cpp` (`getCntExpDist`: expected count =
exp(log total + log p); with lower count 0 the traverse iterator visits all 4^k words; `getMarkovModel`,
`saveFromLargerK`), `seq_model.cpp` (Markov model, log-probabilities); Song et al. 2014 definition of d2S (search
snippet: Ñ_w = N_w − E N_w estimated from each sequence, sum over w ∈ A^k). The Reinert 2009 / Song 2014 full texts
(PMC, arXiv) are blocked by the proxy.

Reference programs executed:

1. **CAFE binary** (built with `g++ -O2`, run as `cafe -M r -K k -D D2star,D2shepp -J jellyfish` with Jellyfish 2.3.1,
   single strand). CAFE derives the r- and (r+1)-mer counts by marginalising the k-mer table on the prefix, and it
   stores log-probabilities and tests `== 0` for "missing", so a transition probability of exactly 1 (log 0) prunes the
   word (E = 0).
2. An independent **Python replica** of the formulas (`itertools.product` over ACGT^k). With CAFE's estimator
   (prefix marginal + the log-zero pruning) it reproduces every CAFE output to the 6 printed digits: random 300/250-nt
   pair (k, r) = (3,0) 0.414504/0.420488, (3,1) 0.463955/0.417225, (4,2) 0.456430/0.498431, (5,0) 0.471348/0.419045,
   (5,1) 0.480534/0.428895 (d2*/d2S); fixtures S1/S2 k=3 r=0 0.446849/0.510267, k=4 r=2 0.336470/0.488786 (needs the
   log-zero emulation), k=5 r=1 0.531238/0.451976; S1/S3 0.565131/0.525654, 0.488625/0.488019, 0.470944/0.422149;
   S2/S3 0.568492/0.528325, 0.691869/0.530888, 0.521633/0.445835.
3. The same replica with the **maximum-likelihood estimator on the sequence** (the definition implemented here) gives
   the values locked in the tests; C# equals all 18 rows (3 pairs × (k, r) ∈ {(2,0), (3,0), (3,1), (4,0), (4,2), (5,1)})
   to 1e-12, e.g. S1/S2 k=3 r=0: D2* 7.136981012975184, D2S −0.6884013403313263, d2* 0.44457941706964565,
   d2S 0.5084031847096179; S1/S3 k=5 r=1: 27.25084993741092, 15.703707350850248, 0.47723568878974626, 0.41877711459909545.

BIC (Python replica `bic_seq` = C# to 1e-9; S1 231.88984329117494 / 255.74056400499072 / 370.90681674887685 /
914.1030179340269 for r = 0..3; a periodic 100-nt sequence picks r = 1). CAFE's printed BIC uses its prefix-marginal
tables and ln(n − r + 1), so the values differ slightly (300-nt sequence A: CAFE 848.929 / 900.92 / 1156.14 /
1816.02 / 4655.87 vs 847.339 / 885.089 / 1047.615 / 1692.768 / 4633.977 here), but the selected order agrees
(0 for both test sequences).

The estimator difference is deliberate: the published statistics estimate the background from the sequence;
CAFE's prefix marginal drops the last k − 1 letters (e.g. T 76 vs 78 in the 300-nt sequence), and its log-zero
sentinel is an implementation artefact (a probability of 1 is valid).

### 7.5 Both-strand D2* / D2S and multiple-pattern spaced words (audit round 2, WP6)

**CAFE `-R`.** Source read: `kmer.cpp` `KmerModel::load` (`(*kmerCntUnorderMap)[index2revCompleIdx(idx, k)] += cnt`
when not single-strand), `KmerProbDelegate::push` (reverse-complement branch), `KmerProbEnsembDelegate::getKmerlogProb`
(log_sum(log p, log p_RC) − log 2), `main.cpp` (`-R` sets `singleStrain = false`; Jellyfish is called without `-C`).
The WP4 CAFE binary (Jellyfish 2.3.1) was run with `-R` on the WP4 pairs; the Python replica extended with
X^R = X(w) + X(RC w), p^R = ½(p + p_RC), total 2n̄ and CAFE's estimator (prefix marginals + log-zero pruning in both
delegates) reproduces all **20 runs to the 6 printed digits** (d2* / d2S): A/B (k, r) = (3,0) 0.352575/0.401636,
(3,1) 0.416406/0.435342, (4,2) 0.500353/0.530720, (5,0) 0.453000/0.454436, (5,1) 0.465944/0.460541; S1/S2
0.381753/0.432310, 0.582753/0.531700, 0.475315/0.495437, 0.513293/0.399786, 0.535129/0.454515; S1/S3
0.555598/0.509182, 0.342658/0.394448, 0.355213/0.452993, 0.487593/0.359372, 0.454912/0.404730; S2/S3
0.660850/0.627895, 0.644076/0.570558, 0.523529/0.500030, 0.524609/0.373953, 0.513226/0.430231 (single-strand A/B
k=3 r=0 = 0.414504/0.420488, the WP4 value). With the sequence maximum-likelihood estimator the replica gives the
locked values; C# equals all 18 rows to 1e-12 (S1/S2 k=3 r=0: D2* 15.34415418761105, D2S 8.071125717307769,
d2* 0.3838876158581438, d2S 0.4345248995381073), and two order-8 rows (k = 10, sparse tables) to 1e-9.

**Spaced words.** Sources: [4] (search snippets: "the distance between two sequences is the average of the distances
based on the individual patterns"; Euclidean and Jensen–Shannon distances of relative spaced-word frequencies);
the `spaced` 1.2.0 source (Ubuntu archive `spaced_1.2.0-201605+dfsg.orig.tar.xz`, `src/sort.h` `spacedDNA`) and
binary (`apt-get install spaced`; earlier tries: spaced.gobics.de not allow-listed, no GitHub repository found).
The source shows: per pattern, `-d EU` accumulates (c₁ − c₂)² on **raw counts** and takes the square root; `-d JS`
divides counts by the number of word positions and sums ½ f log₂(f/m) on both sides; the matrix is divided by the
number of patterns; `-r` disables the reverse complement. Runs (`spaced -r -t 1 -f patterns -d JS|EU`, 12 digits)
= Python replica = scipy (`jensenshannon(base=2)**2`, `euclidean` on frequencies) = C# (1e-12):

| Pair | Patterns | JS (`spaced -d JS`) | EU (`spaced -d EU`) | Euclidean on frequencies |
|------|----------|---------------------|---------------------|--------------------------|
| S1/S2 | 11011, 10111, 11101 | 0.816322854161 (0.8163228541607376) | 12.4089758166 | 0.17567404832368613 |
| S1/S2 | 1101011, 1011101, 1110011 | 0.946781915619 | 11.6604179806 | 0.16946832495600062 |
| S1/S2 | 1111 | 0.800398953366 | 12.4899959968 | 0.17553282636413806 (= `KmerDistance(k=4)`) |
| A/B | 11011, 10111, 11101 | 0.414517255412 | 25.7044109754 | 0.09756721735575496 |
| A/B | 1101011, 1011101, 1110011 | 0.757979428497 | 23.8455999144 | 0.0898354125393724 |
| A/B | 1111 | 0.418942116525 | 25.8069758011 | 0.09773005238795092 |

### 7.6 MinHash sketches and containment (audit round 2, WP7)

Sources read: Mash v2.3 `src/mash/Sketch.cpp` (`addMinHashes`, `sketchFile`, `use64 = pow(alphabetSize, k) > 2^32`),
`hash.cpp`, `MinHashHeap.cpp`, `CommandDistance.cpp` (`compareSketches`, `pValue`), `MurmurHash3.cpp`, `Command.cpp`
(`-k` 1..32, `-s` 1000) — raw.githubusercontent.com/marbl/Mash/v2.3. Executed: the Mash 2.3 binary (`/usr/bin/mash`,
`mash dist -k K -s S [-n]`, `mash info -d`), sourmash 4.9.4 (`MinHash(n=…)`, `MinHash(scaled=1).contained_by` / CLI
`compare --containment`), mmh3 (`hash64(key, 42, signed=False)`). Sequences: 64-bit LCG (`gen.py`, reproduced in the
tests): A random 10 kb; B = A with ~1 % substitutions; C ~5 %; D unrelated random 10 kb; E = A[2000, 7000) lower-case
with "NNNNN" inserted (5005 nt); AA2 = A as two records (6000 + 4000).

- Hashes: `mash info -d` of `mash sketch -k 21 -s 5 A.fa` = 6460448764372083, 7468640311819670, 7856285671510867,
  9659045342622411, 10133432196212684 (= sourmash `MinHash(n=5, ksize=21)`); `-k 16 -s 5` (32-bit) = 357156, 675973,
  681892, 948453, 1205556. Eight `mmh3.hash64` vectors (lengths 0–32) equal `MurmurHash3X64_128`.
- `mash dist` output (distance, p-value, x/s), all 63 runs reproduced **exactly at the printed 6 digits**, e.g.:

| Pair | k | s | Mash output | Note |
|------|---|---|-------------|------|
| A/B | 21 | 1000 | 0.0101459  0  678/1000 | sourmash `num=1000` Jaccard 0.678 (same hashes) |
| A/E, E/A | 21 | 1000 | 0.0196919  0  494/1000 | symmetric; E has N + lower case |
| AA2/A | 21 | 1000 | 2.38274e-05  0  999/1000 | the junction k-mers are missing |
| A/B | 21 | 50 | 0.0118017  1.09468e-273  32/50 | |
| A/B | 21 | 100000 | 0.0101216  0  8069/11891 | = exact Jaccard (sourmash `scaled=1` 0.678580438987470) |
| A/B | 16 | 1000 | 0.00949197  0  753/1000 | 32-bit hashes (sourmash 64-bit `num=1000`: 0.738) |
| A/B | 16 | 50 | 0.0111051  2.23086e-202  36/50 | |
| A/B -n | 21 | 1000 | 0.00997914  0  682/1000 | non-canonical |
| A/D | 8 | 1000 | 0.163715  3.51356e-20  156/1000 | unrelated, small k |
| C/D | 9 | 200 | 0.335603  0.320676  5/200 | |
| A/D | 11 | 1000 | 0.502133  0.33404  2/1000 | |
| A/D | 21 | 1000 | 1  1  0/1000 | nothing shared |

- Containment (canonical, exact) = sourmash `scaled=1` `contained_by`: k=21 A in B 0.808517034068136, A in C
  0.347795591182365, A in E 0.496993987975952, E in A 1; k=16 A in B 0.850575863795693, A in C 0.446169253880821,
  A in E 0.497746619929895. C(A,B)·|K(A)| = 8069 = the `mash dist -s 100000` shared count.

### 7.7 `spaced` N-word rule and reverse-complement mode (audit round 2, WP8)

Source: `spaced` 1.2.0 `src/sort.h` `spacedDNA` (Ubuntu archive `spaced_1.2.0-201605+dfsg.orig.tar.xz`). While
reading, every letter other than A/C/G/T (after `toupper`) is stored as `N`; non-letters are skipped. The reverse
complement of each record is stored after the forward data (N ↦ N). A word is kept only while `correctWord` holds,
i.e. no match position reads `N`. Word positions per record are `seqWordEnd − seqStart` = L − ℓ + 1, including the
dropped ones. With `revComp` (no `-r`), the entry d[i][j] (i > j, i.e. j is the earlier record) uses
`row[i]` = forward counts of record i and `row[j] + row[j + seqNum]` = forward + reverse-strand counts of record j:
EU sums |row_i − (row_j + row_j′)|², JS uses row_i ÷ W_i against (row_j + row_j′) ÷ (2·W_j).

Runs: the `spaced` binary (`spaced [-r] -t 1 -f patterns -d JS|EU`, 12 printed digits) on S1/S2, A/B and N1/N2 =
S1/S2 with N, R, Y and lower case inserted, × the three pattern sets of §7.5, with and without `-r`. All 36 values
= the Python replica = C# (`AcgtOnly`, `bothStrands`; 1e-12, EU 1e-10). Without `-r`:

| Pair (first, second) | Patterns | JS | EU |
|---|---|---|---|
| S1, S2 | 11011, 10111, 11101 | 0.741748311452 | 16.165781203 |
| S1, S2 | 1101011, 1011101, 1110011 | 0.899930296951 | 14.5571570182 |
| S1, S2 | 1111 | 0.705385434836 | 15.7797338381 |
| A, B | 11011, 10111, 11101 | 0.347532568204 | 35.6349133067 |
| A, B | 1101011, 1011101, 1110011 | 0.672311954762 | 29.6025878705 |
| A, B | 1111 | 0.355631754711 | 36.6196668472 |
| N1, N2 | 11011, 10111, 11101 | 0.644779801994 | 14.4441651361 |
| N1, N2 | 1101011, 1011101, 1110011 | 0.724905982831 | 12.8296572574 |
| N1, N2 | 1111 | 0.618160592821 | 14.2478068488 |
| N2, N1 (records swapped) | 11011, 10111, 11101 | 0.632325490347 | 14.604611969 |

With `-r`, N1/N2 (N words dropped): JS 0.706137367852 / 0.767473374082 / 0.693975547992, EU 11.2827316031 /
10.59658835 / 11.3578166916; S1/S2 and A/B as in §7.5. The swapped row shows the order dependence of the mode.

Correction (audit round 3, WP9): the sentence "non-letters are skipped" above describes `spaced`, but until WP9 the C#
`AcgtOnly` path kept non-letters as window positions (none of the 36 inputs contained one). Now the spaced path
deletes them too (§2.8, §7.8).

### 7.8 `spaced -d EV`, the spaced reader on gapped input, and sourmash FracMinHash (audit round 3, WP9)

**spaced.** The binary (`spaced [-r] -t 1 -f patterns -d EV|JS|EU`, 12 printed digits) on 11 pairs: m10/m25/dirty/ab
(WP9 scratch), N1/N2, a gapped pair G1 = `ACG-TACGTACGGT-ACCA*TTG` / G2 = `ACGTTACGTAC1GGTAACCATT`, three LCG pairs
(A = 1500-nt LCG sequence with seed 11; B, C = A with 10 % / 25 % substitutions; D = A with 12 % substitutions, then
gaps `-`, `*`, digits, N, R and lower case inserted, against E = A[200:1400]), `X = ACG-TACGT` / `Y = ACGTTACGA` and
S1/S2; pattern sets {11011, 10111, 11101}, {1101011, 1011101, 1110011}, {1111}, {110101100111, 101110011011,
111001010111}, {11011011, 10111011, 11101101}, {1011, 1101}, {11111111}. **All 248 runs** equal the Python replica
(`ref.py`: reader + N-word rule + the `sort.h` EV formula) to the printed digits, and C# equals all of them to one unit
of the 12th digit (one JS value differs in the 12th digit by summation order). Small weights saturate EV (V < 0 →
1.2). Selected EV values:

| Pair (first, second) | Patterns | `-r` | default (both strands) |
|---|---|---|---|
| A, B (10 %) | weight 8, ℓ 12 | 0.117255304438 | 0.119414426485 |
| A, B (10 %) | 11111111 (`KmerDistance`, k = 8) | 0.119340289783 | 0.122709335653 |
| A, B (10 %) | weight 6, ℓ 8 | 0.27326159765 | 1.2 |
| A, C (25 %) | weight 8, ℓ 12 | 0.308204981926 | 0.311479608973 |
| A, C (25 %) | weight 6, ℓ 8 | −nan (p < ¼) | 1.2 |
| D (dirty), E | weight 8, ℓ 12 | 0.154793408617 | 0.156900799006 |
| N1, N2 | 1101011, 1011101, 1110011 | −nan | 0.767483449137 |
| G1, G2 | 11011, 10111, 11101 | 0.159471857578 | 0.194616371339 |
| X, Y | 1011, 1101 | 0.427666596474 | 0.578620038162 |

Gapped input (JS / EU, `-r`): G1/G2 {11011, 10111, 11101} JS 0.330976739232, EU 3.77806939218 (WP8's C# gave
0.305757130732 / 3.40370085031 with the gaps kept as N windows); X/Y JS 0.57013316426 (was 0.25).

**sourmash 4.9.4** (`pip`; `MinHash(n=0, ksize=k, scaled=S)`, `add_sequence(seq, force=True)`): A = 20 000-nt LCG
sequence (seed 1), B = A with 5 % substitutions, C = A[5000:15000], D = an unrelated 8000-nt sequence, E = B with
every 997th base N and the first 100 lower case; pairs A/B, A/C, A/D, A/E, C/B × k ∈ {21, 31} × S ∈ {1, 10, 100,
1000}. All 40 rows (sketch sizes, `count_common`, `jaccard`, both `contained_by`, `max_containment`) equal C#
exactly (max |Δ| = 0). Examples (k = 21): S = 1 A/B J 0.21773579155873837, C(A,B) 0.3576076076076076; S = 1000 A/B
len 23 / 14, common 5, J 0.15625, C(A,B) 0.21739130436987927, C(B,A) = max 0.35714315204470304; A/C S = 10 C(C,A)
clamped to 1. The five smallest hashes of C at k = 21, S = 1000: 4371404454895749, 5919670475642754,
6049817491049527, 9659045342622411, 10788914448772484 (15 hashes). `_max_hash` for S = 1, 3, 10, 100, 1000, 7919:
18446744073709551615, 6148914691236516864, 1844674407370955264, 184467440737095520, 18446744073709552,
2329428472497733.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [KmerAnalyzer_KmerDistance_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_KmerDistance_Tests.cs) — covers `INV-01`–`INV-04`
- Tests: [KmerAnalyzer_DistanceMetrics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_DistanceMetrics_Tests.cs) — §7.2 values, `INV-05`/`INV-06`, conventions, spaced words
- Tests: [KmerAnalyzer_ParallelAndBackgroundD2_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_ParallelAndBackgroundD2_Tests.cs) — §7.4 values, D2*/D2S conventions and validation
- Tests: [KmerAnalyzer_BothStrandD2AndSpacedWords_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_BothStrandD2AndSpacedWords_Tests.cs) — §7.5 values, null-as-empty, sparse tables, JS metric
- Tests: [KmerAnalyzer_MinHashContainment_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_MinHashContainment_Tests.cs) — §7.6 values (63 `mash dist` rows, hashes, containment)
- Tests: [KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests.cs) — §7.7 values (36 `spaced` runs), N-word rule, order dependence
- Tests: [KmerAnalyzer_SpacedEvAndFracMinHash_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/KmerAnalyzer_SpacedEvAndFracMinHash_Tests.cs) — §7.8 values (120 `spaced` rows, 40 sourmash rows), EV contract, spaced reader, p-value / sketch validation, d2* clamp
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
10. Wan L, Reinert G, Sun F, Waterman MS. 2010. Alignment-free sequence comparison (II): theoretical power of comparison statistics. J Comput Biol 17(11):1467–1490; Song K, Ren J, Reinert G, Deng M, Waterman MS, Sun F. 2014. New developments of alignment-free sequence comparison: measures, statistics and next-generation sequencing. Brief Bioinform 15(3):343–353 (d2* / d2S dissimilarities; full text blocked, definition from search snippets).
11. Lu YY, Tang K, Ren J, Fuhrman JA, Waterman MS, Sun F. 2017. CAFE: aCcelerated Alignment-FrEe sequence analysis. Nucleic Acids Res 45(W1):W554–W559. Source: github.com/younglululu/CAFE `code/dist_model.cpp`, `code/kmer.cpp`, `code/seq_model.cpp` (cloned and built 2026-10-01).
12. Lin J. 1991. Divergence measures based on the Shannon entropy. IEEE Trans Inf Theory 37(1):145–151.
13. `spaced` 1.2.0 (Leimeister, Hahn, Morgenstern), Debian Med package source `spaced_1.2.0-201605+dfsg` (archive.ubuntu.com, `src/sort.h`, `src/spaced.cc`), opened and run 2026-10-01.
14. Koslicki D, Zabeti H. 2019. Improving MinHash via the containment index with applications to metagenomic analysis. Applied Mathematics and Computation 354:206–215 (containment index C(A, B) = |A ∩ B| / |A|).
15. Appleby A. MurmurHash3 (public domain), `MurmurHash3_x64_128`, as bundled in Mash `src/mash/MurmurHash3.cpp`.
16. Morgenstern B, Zhu B, Horwege S, Leimeister CA. 2015. Estimating evolutionary distances between genomic sequences from spaced-word matches. Algorithms for Molecular Biology 10:5. doi:10.1186/s13015-015-0032-x (implemented as `spaced` 1.2.0 `sort.h` EV branch [13]).
17. Irber L, Brooks PT, Reiter T, Pierce-Ward NT, Hera MR, Koslicki D, Brown CT. 2022. Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers. bioRxiv 2022.01.11.475838.
18. Hera MR, Pierce-Ward NT, Koslicki D. 2023. Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash. Genome Research 33:1061–1068 (containment bias factor 1 − (1 − 1/S)^(|A|·S)).
19. sourmash 4.9.4 (PyPI; `sourmash/minhash.py` installed) and its Rust core at tag v4.9.4 (raw.githubusercontent.com `src/core/src/sketch/minhash.rs` `max_hash_for_scaled`, `jaccard`, `add_hash`; `src/core/src/signature.rs` `SeqToHashes`; `src/core/src/lib.rs` `_hash_murmur`), opened 2026-10-01.
