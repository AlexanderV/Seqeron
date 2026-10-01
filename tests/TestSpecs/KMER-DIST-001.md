# Test Specification: KMER-DIST-001

**Test Unit ID:** KMER-DIST-001
**Area:** K-mer
**Algorithm:** K-mer Euclidean Distance (alignment-free word-frequency distance)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-10-01 (audit round 1, WP2)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Zielezinski, Vinga, Almeida & Karlowski (2017). Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology 18:186. | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ (10.1186/s13059-017-1319-7) | 2026-06-13 |
| 2 | Lau AK et al. (2022). Interpreting alignment-free sequence comparison: what makes a score a good score? NAR Genom Bioinform. | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/ | 2026-06-13 |
| 3 | Vinga S, Almeida J (2003). Alignment-free sequence comparison—a review. Bioinformatics 19(4):513–523. | 1 | https://academic.oup.com/bioinformatics/article/19/4/513/218529 (10.1093/bioinformatics/btg005) | 2026-06-13 |
| 4 | Leimeister C-A, Boden M, Horwege S, Lindner S, Morgenstern B (2014). Fast alignment-free sequence comparison using spaced-word frequencies. Bioinformatics 30(14):1991. | 1 | https://academic.oup.com/bioinformatics/article/30/14/1991/2391234 | 2026-10-01 (definition via search snippets; PDF 403) |
| 5 | Blaisdell BE (1986). PNAS 83:5155 — squared count Euclidean d_E. | 1 | 10.1073/pnas.83.14.5155 | 2026-10-01 (via alfpy citation + Vinga & Almeida 2003) |
| 6 | Torney et al. 1990; Lippert, Huang & Waterman 2005; Reinert et al. 2009 — D2 = Σ X_w Y_w. | 1 | J Comput Biol 16:1615 | 2026-10-01 |
| 7 | Ondov BD et al. (2016). Mash. Genome Biol 17:132 (eq. 1 Jaccard, eq. 4 Mash distance); Mash source `CommandDistance.cpp`, `Sketch.cpp`. | 1 | raw.githubusercontent.com/marbl/Mash/master/src/mash/ | 2026-10-01 |
| 8 | Reference implementations: scipy 1.17.1 `spatial.distance`; scikit-bio 0.7.4; alfpy 1.0.6 (PyPI sdist); sourmash 4.9.4; Mash 2.3 binary. | 2 | PyPI / apt | 2026-10-01 |

### 1.2 Key Evidence Points

1. Each sequence is mapped to a vector of word (k-mer) counts over the union of words; identical sequences give distance 0 — Source 1 (Figure 1).
2. Worked example: x="ATGTGTG", y="CATGTG", k=3 ⇒ count vectors c_X=(1,0,2,2), c_Y=(1,1,1,1) over (ATG,CAT,GTG,TGT) — Source 1 (Figure 1).
3. "This difference is very commonly computed by the Euclidean distance" — Source 1 (Figure 1).
4. k-mer frequency = count ÷ (sequence length − k) [= number of k-mer windows L − k + 1]; Euclidean is applied to the frequency vectors — Source 2.
5. The standard Euclidean alignment-free distance is taken over relative-frequency (normalized) vectors — Source 4 ("Euclidean distance to the relative-frequency vectors").
6. Word-composition methods map a sequence into a 4^k-dimensional k-word frequency vector; Euclidean distance is one of the standard dissimilarity measures — Source 3.

### 1.3 Documented Corner Cases

- Identical sequences ⇒ distance 0 (Source 1).
- A word absent from a sequence contributes a 0 vector component (Source 1: c_X has 0 for CAT).
- Frequencies require L ≥ k (Source 2).

### 1.4 Known Failure Modes / Pitfalls

1. Mixing count-based and frequency-based variants yields different numbers; the implementation under test uses the **frequency** variant (Source 2/4). — Sources 2, 4.
2. No authoritative source defines the distance when a sequence has zero k-mer windows (L < k); this is an ASSUMPTION (see §6).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `KmerDistance(string seq1, string seq2, int k)` | KmerAnalyzer | Canonical | Euclidean distance over normalized k-mer frequency vectors (delegates to the metric overload) |
| `KmerDistance(string, string, int, KmerDistanceMetric)` / `KmerDistance(IReadOnlyDictionary<string,int>, IReadOnlyDictionary<string,int>, KmerDistanceMetric)` | KmerAnalyzer | Canonical (WP2) | Euclidean, SquaredEuclideanCounts, Manhattan, Chebyshev, Canberra, Cosine, D2 |
| `JaccardSimilarity(string, string, int[, KmerCountingOptions])` | KmerAnalyzer | Canonical (WP2, DUP_MAP §16) | Exact k-mer set Jaccard, fraction, both-empty → 0 |
| `MashDistance(string, string, int, KmerCountingOptions)`, `MashDistanceFromJaccard(double, int)` | KmerAnalyzer | Canonical (WP2) | Ondov 2016 eq. 4 + Mash boundary rules |
| `CountSpacedWords(string, string)` | KmerAnalyzer | Canonical (WP2) | Leimeister 2014 spaced words |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | d(x, x) = 0 for any x with at least one k-mer | Yes | Source 1 (Fig. 1) |
| INV-2 | d(x, y) = d(y, x) (symmetry) | Yes | Source 3 (Euclidean is a metric) |
| INV-3 | d(x, y) ≥ 0 | Yes | Source 1/3 (Euclidean norm is non-negative) |
| INV-4 | For two sequences each consisting of a single distinct k-mer (frequency 1) with disjoint word sets, d = √2 | Yes | Derived from Source 2 frequency definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Zielezinski Fig.1 example | x="ATGTGTG", y="CATGTG", k=3 | 0.33166247903553997 (√0.11) | Source 1 (Fig.1 vectors) + Source 2 (freq def) |
| M2 | Identical sequences | x="ATGTGTG", y="ATGTGTG", k=3 | 0.0 exactly | Source 1 (Fig.1: identical ⇒ 0) |
| M3 | Single substitution, k=1 | "AAAA" vs "AAAT", k=1 | 0.3535533905932738 (√0.125) | Source 2 (freq def, derivation) |
| M4 | Symmetry | d("ATGTGTG","CATGTG",3) == d("CATGTG","ATGTGTG",3) | equal (= √0.11) | Source 3 (metric symmetry) / INV-2 |
| M5 | Non-negativity | d("AAAA","TTTT",2) ≥ 0 and equals √2 | 1.4142135623730951 | INV-3, INV-4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Disjoint single-kmer sequences | "AAAA" vs "TTTT", k=2 | √2 = 1.4142135623730951 | INV-4; both vectors are single 1.0 components |
| S2 | One sequence too short for k | "ACGT" (k=5 ⇒ empty) vs "AAAAAA" (k=5 ⇒ "AAAAA"=1.0) | 1.0 | ASSUMPTION A2 (empty vector = zero vector) |
| S3 | Case-insensitivity | "atgtgtg" vs "CATGTG", k=3 | 0.33166247903553997 | ASSUMPTION A1 (upper-casing); equals M1 |
| S4 | skbio/scipy reference cross-check | ACGTTGCAACGGT/ACGTAGCATCGGTA k=2; GATTACAGATTACA/GATTACCGATTTCA k=3; Fig.1 | 0.26600633232367216; 0.31180478223116176; 0.33166247903554 | scikit-bio kmer_frequencies(relative) + scipy euclidean (2026-09 review) |
| S5 | Frequency, not count, variant | Fig.1 inputs | ≠ √3 and ≠ 3 (count Euclidean / Blaisdell squared d_E) | Source 1 count vectors; Vinga & Almeida 2003 d_E |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Invalid k | k = 0 | throws ArgumentOutOfRangeException | Validation inherited from CountKmers |
| C2 | Both sequences empty | "" vs "", k=3 | 0.0 | ASSUMPTION A2 (both empty ⇒ 0) |

### 4.4 Audit round 1 (WP2) — metrics, Jaccard/Mash, spaced words (`KmerAnalyzer_DistanceMetrics_Tests.cs`)

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| W1 | 7 metrics × 6 pairs (Fig.1, ACGTTGCAACGGT k=2, GATTACA k=3, lower-case/IUPAC k=3, R1/R2 k=5, k=11), both orders | values of algorithm doc §7.2 (tol 1e-12) | scipy + alfpy (Sources 5, 6, 8) |
| W2 | Fig.1 SquaredEuclideanCounts | exactly 3 | Source 5 via Vinga & Almeida 2003 |
| W3 | Euclidean metric = legacy KmerDistance | bit-identical (4 inputs) | refactor invariant |
| W4 | identical inputs → 0 for every dissimilarity; D2(x,x) = Σc² (ATGTGTG k=3 → 9) | exact | definitions |
| W5 | zero-vector conventions (additive metrics 0; cosine 1; one empty: Manhattan 1, Canberra = #words, d_E = Σc²) | exact | scipy clip / similarity-0 convention (A3) |
| W6 | contracts: k ≤ 0, undefined metric, null table, negative count | throw | — |
| W7 | metrics on canonical (-C) count tables | d_E 3, D2 6 | count-table overload |
| W8 | Jaccard literal / ACGT-only / canonical and Mash D (canonical, -n) × 6 pairs | doc §7.2 table | Python sets, sourmash scaled=1, Mash 2.3 binary |
| W9 | Mash binary shared-hashes = exact counts (58/158, 55/168) and printed distances 0.124338 / 0.141338 | exact / 5e-7 | Mash 2.3 |
| W10 | Jaccard conventions (both empty 0, null, case-insensitive, Mash both-empty 0, nothing shared 1, k ≤ 0 throws) | exact | Mash `CommandDistance.cpp`; GenomicAnalyzer convention (A4) |
| W11 | Jaccard = GenomicAnalyzer.CalculateSimilarity/100 (k = 2, 3, 5, 8) | equal | DUP_MAP §16 |
| W12 | MashDistanceFromJaccard (1→0, 0→1, 0.5/k2, 0.75/k3, cap) and out-of-range throws | exact | Ondov eq. 4 + Mash cap |
| W13 | spaced-word tables (101, 1101, 11011, lower-case), R1/1100111 108 distinct/114 total, all-ones = CountKmers, spaced distance 0.3149183286488868 / 12 | exact | Source 4 definition, Python replica |
| W14 | invalid patterns ("", 0110, 1100, 1021), null pattern, short/empty sequence | throw / empty | Source 4 (P[1] = P[ℓ] = 1) |

---

### 4.5 Audit round 1 (WP4) — background-adjusted D2* / D2S (`KmerAnalyzer_ParallelAndBackgroundD2_Tests.cs`)

| ID | Test | Evidence |
|----|------|----------|
| D1 | 18 rows (3 pairs × (k, r) ∈ {(2,0),(3,0),(3,1),(4,0),(4,2),(5,1)}): D2*, D2S, d2*, d2S equal the Python replica to 1e-12 | Reinert et al. 2009; Song et al. 2014; CAFE formulas (replica engine = CAFE binary to 6 digits) |
| D2 | `KmerDistance(.., D2Star/D2Shepherd)` = order-0 statistics; symmetric | definition |
| D3 | Identical sequences → 0; values in [0, 1] | definition |
| D4 | Case-insensitive; non-ACGT windows skipped | Jellyfish convention (CAFE counts with Jellyfish) |
| D5 | Homopolymer at order 0 → NaN distance (0/0), D2* = 0 | definition (all centred counts 0) |
| D6 | Validation: null, k ∉ [1, 12], r ∉ [−1, k), no ACGT window, counts overload rejects D2*/D2S | contract |
| D7 | `MarkovOrderBic` = replica (S1, S2, S3, periodic; r = 0..3, 1e-9); `markovOrder = −1` picks (1, 0) for periodic/S1 and gives the replica values | Schwarz 1978; Katz 1981; CAFE `-M -1` |
| D8 | `ParseDistanceMetric` names incl. `d2star`, `d2shepherd`, `d2s`; unknown → `ArgumentException` | MCP contract |

### 4.6 Audit round 2 (WP6) — both strands, sparse tables, null-as-empty, JS, spaced-word distance (`KmerAnalyzer_BothStrandD2AndSpacedWords_Tests.cs`)

| ID | Test | Evidence |
|----|------|----------|
| R1 | 18 rows (S1/S2, S1/S3, S2/S3 × (k, r) ∈ {(2,0),(3,0),(3,1),(4,2),(5,1)} + A/B × {(3,0),(4,2),(5,1)}): `bothStrands` D2*, D2S, d2*, d2S = replica (1e-12) | CAFE `-R` (`KmerModel::load`, `getKmerlogProb`); replica engine = CAFE `-R` binary on 20 runs (6 digits) |
| R2 | `bothStrands = false` = single-strand statistic (WP4 value; 4-argument overload equal) | backward compatibility |
| R3 | `bothStrands` invariant under reverse-complementing a sequence (r = 0) | X(w)+X(RC w) and p(w)+p(RC w) are strand-symmetric |
| R4 | `KmerDistance(.., bothStrands)` = statistics; auto order (−1) reports the BIC minimisers; `bothStrands` with a plain metric → `ArgumentException` (`bothStrands`) | contract |
| R5 | Order 8, k = 10 (sparse Markov tables), single and both strands = replica (1e-9) | ML Markov chain |
| R6 | k = 12, r = 11 allocates < 16 MB (no 4^12-double table) | audit item (134 MB per sequence before) |
| R7 | `KmerDistance(null, ..)` = `KmerDistance("", ..)` for every metric (D2*/D2S: same `ArgumentException`, not `ArgumentNullException`) | consistency |
| R8 | `JensenShannon` = scipy `jensenshannon(base=2)²` (0.1522040934665307), symmetric, identity 0, disjoint 1, zero vector ½; `EuclideanCounts` = √3 | Lin 1991; scipy 1.17.1 |
| R9 | `SpacedWordDistance` (6 rows: S1/S2, A/B × {w4, w5, contiguous 1111}): JS = `spaced -r -d JS`, `EuclideanCounts` = `spaced -r -d EU`, `Euclidean` = scipy frequency Euclidean (1e-12 / 1e-10) | Leimeister et al. 2014; spaced 1.2.0 binary + source |
| R10 | Contiguous pattern = `KmerDistance(k)`; mean of per-pattern values; identity 0; validation (null/empty set, null pattern, unequal weight, malformed pattern, D2*) | definition / contract |

### 4.7 Audit round 2 (WP7) — Mash MinHash sketches, containment index (`KmerAnalyzer_MinHashContainment_Tests.cs`)

| ID | Test | Evidence |
|----|------|----------|
| M1 | `MurmurHash3X64_128(key, 42)` = `mmh3.hash64(key, 42, signed=False)` on 8 keys of length 0–32 (both words) | Appleby MurmurHash3; Mash `MurmurHash3.cpp`; mmh3 |
| M2 | `CreateMinHashSketch(A, 21, 5)` hashes = `mash info -d` (= sourmash `MinHash(n=5, ksize=21)`); `(A, 16, 5)` = the 32-bit `mash info -d` hashes, `Use64` false | Mash `Sketch.cpp` (`use64 = 4^k > 2^32`), `hash.cpp` |
| M3 | Sketch length counts N (E: 5005 = `mash info -d`); records shorter than k are skipped (length and k-mers) | Mash `sketchFile` (`l < kmerSize`) |
| M4 | 63 `mash dist` rows (A–E, AA2; k ∈ {8, 9, 10, 11, 12, 16, 21}; s ∈ {10, 50, 200, 1000, 100000}; canonical and `-n`): x, denominator exact; J = x/denominator; distance and p-value formatted `G6` = Mash's printed text | Mash 2.3 binary; `CommandDistance.cpp` `compareSketches`/`pValue` |
| M5 | s ≥ union → sketch J and distance = exact `JaccardSimilarity`/`MashDistance` (8069/11891 = sourmash `scaled=1` 0.678580438987470) | Mash; sourmash 4.9.4 |
| M6 | Different sketch sizes → merge capped at min(s) (32/50); mismatched k/seed/canonical → `ArgumentException`; k ∉ 1..32, s < 1, null → exceptions; both sketches empty → (0, 0, 0, 0, 1) | Mash (sketch-size reduction, parameter checks, `-k` 1..32); convention |
| M7 | `MashPValue`: x = 0 → 1; A/D k=9 s=200 x=9 → 0.0137623 (`mash dist`); negative inputs / k ≤ 0 throw | Mash `pValue` |
| M8 | `ContainmentIndex` canonical (9 rows, k = 21 and 16) = sourmash `MinHash(scaled=1).contained_by` (1e-14) | Koslicki & Zabeti 2019; sourmash `compare --containment` |
| M9 | C(A,B)·|K(A)| = C(B,A)·|K(B)| = 8069 (Mash shared count); J ≤ min(C(A,B), C(B,A)); literal GATTACA/GATTACCA k=3 → 0.8; empty/null → 0; k ≤ 0 throws | definition / contract |

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzerTests.cs` — contains three permissive KmerDistance tests (`KmerDistance_IdenticalSequences_ReturnsZero`, `KmerDistance_DifferentSequences_ReturnsPositive`, `KmerDistance_SimilarSequences_SmallDistance`) in a "K-mer Distance" region.
- No dedicated `KmerAnalyzer_KmerDistance_Tests.cs` existed before this unit.
- `tests/Seqeron/Seqeron.Mcp.Sequence.Tests/KmerDistanceTests.cs` — MCP-layer tests, out of scope for this algorithm unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 (Fig.1 example) | ❌ Missing | No exact-value test for the canonical worked example |
| M2 (identical ⇒ 0) | ⚠ Weak | Existing `KmerDistance_IdenticalSequences_ReturnsZero` uses tolerance 0.0001, not exact; in auxiliary file — rewrite |
| M3 (single sub, k=1) | ❌ Missing | — |
| M4 (symmetry) | ❌ Missing | — |
| M5 (non-negativity √2) | ⚠ Weak | Existing `KmerDistance_DifferentSequences_ReturnsPositive` uses `GreaterThan(0)` only — rewrite to exact |
| S1 (disjoint √2) | ⚠ Weak | Existing `KmerDistance_SimilarSequences_SmallDistance` uses `LessThan` ordering only — rewrite |
| S2 (short input) | ❌ Missing | — |
| S3 (case-insensitive) | ❌ Missing | — |
| C1 (invalid k) | ❌ Missing | — |
| C2 (both empty) | ❌ Missing | — |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/KmerAnalyzer_KmerDistance_Tests.cs` — all KMER-DIST-001 evidence-based tests (M/S/C).
- **Remove:** the three permissive KmerDistance tests and the "K-mer Distance" region from `KmerAnalyzerTests.cs`; update that file's header comment to point KmerDistance to the new dedicated file.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `KmerAnalyzer_KmerDistance_Tests.cs` | Canonical KMER-DIST-001 fixture | 10 |
| `KmerAnalyzerTests.cs` | Auxiliary (KmerDistance region removed) | unchanged minus 3 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact √0.11 test | ✅ Done |
| 2 | M2 | ⚠ Weak | Rewrote as exact 0.0 (Within 1e-10) in canonical file; removed weak original | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented exact √0.125 test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented symmetry test | ✅ Done |
| 5 | M5 | ⚠ Weak | Rewrote as exact √2; removed weak original | ✅ Done |
| 6 | S1 | ⚠ Weak | Rewrote as exact √2; removed weak original | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented short-input test (=1.0) | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 9 | C1 | ❌ Missing | Implemented invalid-k throw test | ✅ Done |
| 10 | C2 | ❌ Missing | Implemented both-empty test | ✅ Done |

**Total items:** 10
**✅ Done:** 10 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact √0.11 |
| M2 | ✅ Covered | Exact 0.0 |
| M3 | ✅ Covered | Exact √0.125 |
| M4 | ✅ Covered | Symmetry |
| M5 | ✅ Covered | Exact √2 |
| S1 | ✅ Covered | Exact √2 |
| S2 | ✅ Covered | Exact 1.0 |
| S3 | ✅ Covered | Equals M1 |
| C1 | ✅ Covered | Throws |
| C2 | ✅ Covered | 0.0 |

**Total in-scope cases:** 10 — **✅ Covered:** 10

---

### 5.7 Audit round 2, WP8 (B06) — `spaced` N-word rule and reverse-complement mode

`CountSpacedWords(sequence, pattern, options)`, `SpacedWordDistance(seq1, seq2, patterns, metric, options, bothStrands)`;
tests in `KmerAnalyzer_StrandOptionsAndSpacedConventions_Tests.cs`; reference values K-mer_Euclidean_Distance.md §7.7.

| ID | Case | Expected (source) |
|----|------|-------------------|
| W1 | 36 rows: S1/S2, A/B, N1/N2 × 3 pattern sets × JS/EU × with/without `-r`, `AcgtOnly` | `spaced` 1.2.0 binary (12 digits) = replica (1e-12; EU 1e-10) |
| W2 | `bothStrands`, records swapped (N2, N1) | `spaced` 0.632325490347 / 14.604611969 (order dependence) |
| W3 | `CountSpacedWords("ANGTCN", "1011", AcgtOnly)` | {AGT: 1} (N at the don't-care position kept, N at a match position dropped; `sort.h` `correctWord`) |
| W4 | all-'1' pattern + `AcgtOnly` | = `CountKmers(…, AcgtOnly)` |
| W5 | Default options, single strand | = 4-argument `SpacedWordDistance` (bit-identical); ACGT input + `AcgtOnly` = literal |
| W6 | Literal words on N input | ≠ `spaced` value (documents why `AcgtOnly` is needed) |
| W7 | `Canonical` | `ArgumentException` (CountSpacedWords and SpacedWordDistance) |
| W8 | D2*, unequal weights, malformed pattern, null sequences with `bothStrands` | `ArgumentException` ×3; null = empty → 0 |
| W9 | MCP `spaced_word_distance(acgtOnly, bothStrands)` | N1/N2 JS 0.706137367852 / 0.644779801994, EU both 14.4441651361; S1/S2 JS both 0.741748311452 |

### 5.8 Audit round 3, WP9 (B06) — `spaced -d EV`, spaced reader, FracMinHash, validation

Tests in `KmerAnalyzer_SpacedEvAndFracMinHash_Tests.cs`; reference values K-mer_Euclidean_Distance.md §7.8.

| ID | Case | Expected (source) |
|----|------|-------------------|
| V1 | 120 rows: LCG pairs A/B, A/C, D(dirty)/E, G1/G2 (gaps), N1/N2 × pattern sets × EV/JS/EU × with/without `-r`, `AcgtOnly` | `spaced` 1.2.0 binary, 12 printed digits (one unit of the 12th digit; `-nan` = NaN) |
| V2 | `ACG-TACGT` / `ACGTTACGA`, {1011, 1101} | `spaced`: JS 0.57013316426 (both modes), EU 2.64575131106 / 4.69041575982, EV 0.427666596474 / 0.578620038162; non-letters deleted = stripped input |
| V3 | `KmerDistance(…, SpacedEvolutionary, 0, bothStrands)` | = `SpacedWordDistance` with 1^k; A/B k = 8: `spaced -f 11111111` 0.119340289783 / 0.122709335653; 1111 → 1.2 |
| V4 | EV with literal options; count-table overload; unequal pattern lengths; too-short read; `Canonical`; `markovOrder` ≠ 0 | literal = `AcgtOnly`; `ArgumentException` (param metric / patterns / seq1 / seq2) |
| V5 | `ParseDistanceMetric("ev" / " EV " / "evolutionary")` | `SpacedEvolutionary` |
| V6 | `FracMinHashMaxHash(1, 3, 10, 100, 1000, 7919)` | sourmash `_max_hash` (Rust truncation) |
| V7 | 40 rows: 5 pairs × k 21/31 × scaled 1/10/100/1000 | sourmash 4.9.4 len, `count_common`, `jaccard`, `contained_by` ×2, `max_containment` (1e-15) |
| V8 | First 5 hashes of C, k = 21, scaled 1000 | sourmash hashes (15 total) |
| V9 | scaled 1 | = exact canonical `JaccardSimilarity` / `ContainmentIndex` |
| V10 | Empty sketches; k/scaled < 1; mismatched scaled/k/canonical; unsorted/duplicate hashes; wrong MaxHash; hash > MaxHash | 0s; `ArgumentOutOfRangeException` / `ArgumentException` |
| V11 | `MashPValue(1, 0, 0, 21, 10)`, `(1, 10, 0, …)`, `(5, 1000, 1000, 21, 3)` | `ArgumentException` (length1 / length2 / sharedHashes); x = 0 → 1 |
| V12 | `CompareMinHashSketches` with reversed / duplicate hashes, Count > SketchSize, null hashes, Length < 0 | `ArgumentException`; `MinHashSketch.FromHashes` sorts/dedups/keeps s |
| V13 | S1/S1 k = 2 r = 0 (raw d2* = d2S = −1.11e-16) | 0 (clamped to [0, 1]) |
| V14 | MCP `kmer_jaccard(scaled: 3)`, `spaced_word_distance(metric: "ev")`, `kmer_distance(metric: "ev")` (both servers) | sourmash S1/S2 k = 4 S = 3 values; `spaced` 0.767483449137, NaN, 0.427666596474, 1.2 |
| V15 | `FracMinHashMaxHash(2049, 123456789, 4294967295)`; S = 4294967296 / 0; sketch at S = 4294967295 (WP10) | 9002803354665472, 149418628356, 4294967297 (sourmash `_max_hash`); `ArgumentOutOfRangeException` |
| V16 | `DownsampleFracMinHash` C k = 21: 10 → 1000; 1 → 7919; S′ < S | 1011 → 15 hashes = the S = 1000 sketch, `_max_hash` 18446744073709552 / 2329428472497733 (sourmash `downsample`); `ArgumentException` (newScaled) |
| V17 | `CompareFracMinHashSketches(…, downsample: true)`, 6 pairs of mixed scaled | sourmash `jaccard(downsample=True)` and `sourmash compare --containment/--max-containment` (= downsample both, then `contained_by`) (1e-15); without `downsample` → `ArgumentException` |
| V18 | `trackAbundance` sketch of P (k = 21, S = 1) and its downsample to 10 | 2600 hashes, abundances {1, 2, 3}, Σ 3380; 230 hashes, Σ 280 (sourmash `track_abundance=True`) |
| V19 | 9 abundance comparisons (P/Q/W tandem-repeat inputs, A/B/C; k 7–31, S 1–100, mixed S) | sourmash `angular_similarity` (= `similarity`), Jaccard, `contained_by_weighted` both ways (1e-15) |
| V20 | One flat sketch; identical sketches; empty sketch; bad `Abundances` | angular null / Jaccard 0.63748031496063; 1; 0; `ArgumentException` |
| V21 | `MashPValue(1, 1000, 1000, 600 / 33 / 0, 10)`; k = 32 | `ArgumentOutOfRangeException("k")` (Mash `-k` 1..32); 2.7105054312137606e-16 (scipy `binom.sf`) |
| V22 | MCP `kmer_jaccard`: `scaled` without `canonical`; `trackAbundance` without `scaled`; S = 4294967295 / 4294967296; N1/N2 with N/IUPAC; S1/S2 `trackAbundance` S = 1/3 | `ArgumentException` (canonical / trackAbundance / scaled); sourmash `force=True` values; `angularSimilarity` 0.2363801370444173 / 0.3123095603640216 |
| V23 | MCP `kmer_jaccard` S1/S2 k = 4 `trackAbundance` S = 1 / 3; without `trackAbundance` (WP11) | `weightedContainmentSeq1InSeq2` / `Seq2InSeq1` = sourmash `contained_by_weighted` 0.37662337662337664 / 0.417910447761194 (S = 1), 0.4642857142857143 / 0.5238095238095238 (S = 3); null |
| V24 | `CompareMinHashSketches` with `MinHashSketch(40, 10, true, 42, true, 1000, {1,2,3})`, K = 0 / 33, `Use64` ≠ K > 16, 32-bit sketch with a hash > 2^32 − 1; `FromHashes(16, …, [2^32])` (WP11) | `ArgumentException` naming the sketch (Mash `-k` 1..32, `use64 = 4^k > 2^32`); K = 32 / K = 16 boundaries valid; `FromHashes` `ArgumentException("hashes")` |

## 6. Assumption Register

**Total assumptions:** 4

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Inputs are upper-cased before counting (case-insensitive) | S3 |
| A2 | A sequence with no k-mer windows (L < k, or empty) is treated as the zero frequency vector | S2, C2 |
| A3 | Cosine with a zero vector: similarity 0 → distance 1 (scipy returns NaN there) | W5 |
| A4 | Jaccard of two empty sets is 0 (undefined 0/0; GenomicAnalyzer / ComparativeGenomics convention) | W10 |

---

## 7. Open Questions / Decisions

1. **Count vs frequency variant.** Source 1's Figure 1 shows counts; Sources 2 and 4 define the distance over normalized frequencies. The implementation uses the **frequency** variant (counts ÷ total windows), which is explicitly endorsed by Sources 2 and 4. Decision: keep the frequency variant; the count-based numbers from Source 1 are recomputed as frequency-based expected values for M1. No conflict — both are documented; the unit tests the implemented (frequency) behavior with source-derived numbers.
