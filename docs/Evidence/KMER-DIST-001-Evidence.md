# Evidence Artifact: KMER-DIST-001

**Test Unit ID:** KMER-DIST-001
**Algorithm:** K-mer Euclidean Distance (alignment-free word-frequency distance)
**Date Collected:** 2026-06-13

---

## Online Sources

### Zielezinski, Vinga, Almeida & Karlowski (2017) — "Alignment-free sequence comparison: benefits, applications, and tools" (Genome Biology)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; located the worked example in Figure 1)
**Authority rank:** 1 (peer-reviewed review, Genome Biology 18:186, co-authored by Vinga & Almeida who authored the foundational 2003 review)

**Key Extracted Points:**

1. **Word-vector representation:** Each sequence is transformed into a vector by "counting the number of times each particular word (from W₃) appears within the sequences." The example uses sequences x = "ATGTGTG" and y = "CATGTG" with word size k = 3.
2. **Word set (union over both sequences):** W₃ = {ATG, CAT, GTG, TGT}; per-sequence unique words W_X_3 = {ATG, TGT, GTG}, W_Y_3 = {CAT, ATG, TGT, GTG}.
3. **Count vectors (verbatim from Figure 1):** c_X_3 = (1, 0, 2, 2) and c_Y_3 = (1, 1, 1, 1), in the column order (ATG, CAT, GTG, TGT).
4. **Distance rule:** "This difference is very commonly computed by the Euclidean distance, although any metric can be applied."
5. **Zero/identity property:** "identical sequences yield a distance of 0"; higher value ⇒ more distant sequences.

### Lau, Kläne, Leimeister, Morgenstern et al. — "Interpreting alignment-free sequence comparison: what makes a score a good score?" (NAR Genomics and Bioinformatics, 2022)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; section "Overview of alignment-free methods based on k-mer frequencies" and Table 1)
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Frequency definition (verbatim):** "The k-mer frequencies are derived from the counts by dividing each k-mer count by the total number of k-mers in the sequence (i.e. the sequence length minus the k-mer length)." For a sequence of length L over a complete alphabet the number of k-mer windows is L − k + 1.
2. **Distance choice (verbatim):** "Once the frequency vectors have been calculated, various distance measures can be used to estimate similarity, which include the Euclidian, Manhattan, Canberra or Chebyshev distances."
3. **Variable length suitability:** Table 1 lists "Euclidian distance (euclid)" as a metric suitable for sequences of variable length.

### Vinga & Almeida (2003) — "Alignment-free sequence comparison—a review" (Bioinformatics)

**URL:** https://academic.oup.com/bioinformatics/article/19/4/513/218529
**Accessed:** 2026-06-13 (abstract/metadata retrieved via WebFetch; full text paywalled — used only for the high-level statement below, not for the numeric formula)
**Authority rank:** 1 (foundational peer-reviewed review)

**Key Extracted Points:**

1. **Word-composition mapping:** word-composition methods map each sequence into a 4^k-dimensional vector according to k-word frequency; the similarity score is obtained by measures including the Euclidean distance, Pearson correlation, Kullback–Leibler discrepancy and cosine distance.

### Boden, Schöneich, Horwege, Lindner, Leimeister & Morgenstern (2014) — "Fast alignment-free sequence comparison using spaced-word frequencies" (Bioinformatics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4080745/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; "Benchmark Set-Up" section)
**Authority rank:** 1 (peer-reviewed; reference method)

**Key Extracted Points:**

1. **Euclidean on relative frequencies (verbatim):** "we applied the Euclidean distance to the relative-frequency vectors obtained with our multiple-pattern approach." Confirms that the standard alignment-free Euclidean distance is taken over relative (normalized) word-frequency vectors, not raw counts.

---

## Documented Corner Cases and Failure Modes

### From Zielezinski et al. (2017)

1. **Identical sequences:** distance = 0 (stated explicitly in Figure 1 description).
2. **Union of words:** the comparison vector spans the union of words occurring in either sequence; words absent from a sequence contribute a 0 component (c_X_3 has a 0 for CAT, which occurs only in y).

### From Lau et al. (2022)

1. **Normalization domain:** frequencies are counts divided by (L − k + 1); this requires at least one k-mer window (L ≥ k), otherwise the frequency vector is empty.

---

## Test Datasets

### Dataset: Zielezinski 2017 Figure 1 worked example

**Source:** Zielezinski et al. (2017), Genome Biology 18:186, Figure 1.

| Parameter | Value |
|-----------|-------|
| Sequence x | ATGTGTG |
| Sequence y | CATGTG |
| k | 3 |
| Word order | ATG, CAT, GTG, TGT |
| c_X (counts) | (1, 0, 2, 2) |
| c_Y (counts) | (1, 1, 1, 1) |
| Total windows x | 5 (= 7 − 3 + 1) |
| Total windows y | 4 (= 6 − 3 + 1) |

**Derived count-based Euclidean distance** (verbatim vectors from source): differences (0, −1, 1, 1) ⇒ √(0² + 1² + 1² + 1²) = √3 ≈ 1.7320508075688772.

**Derived frequency-based Euclidean distance** (counts divided by total windows, per Lau et al. 2022 definition; this is what `KmerAnalyzer.KmerDistance` computes):
f_X = (1/5, 0, 2/5, 2/5) = (0.2, 0.0, 0.4, 0.4); f_Y = (1/4, 1/4, 1/4, 1/4) = (0.25, 0.25, 0.25, 0.25).
differences (−0.05, −0.25, 0.15, 0.15) ⇒ squares (0.0025, 0.0625, 0.0225, 0.0225), sum = 0.11 ⇒ √0.11 = 0.33166247903553997.

### Dataset: Single-substitution small case (derivation, k = 1)

**Source:** Direct derivation from Lau et al. (2022) frequency definition; trivially verifiable.

| Parameter | Value |
|-----------|-------|
| Sequence 1 | AAAA |
| Sequence 2 | AAAT |
| k | 1 |
| Words | A, T |
| f_1 | A=4/4=1.0, T=0 |
| f_2 | A=3/4=0.75, T=1/4=0.25 |

Euclidean distance = √((1.0−0.75)² + (0−0.25)²) = √(0.0625 + 0.0625) = √0.125 = 0.3535533905932738.

---

## Assumptions

1. **ASSUMPTION: Case folding / alphabet handling.** Authoritative sources work over a fixed nucleotide alphabet and do not specify case sensitivity. The implementation upper-cases input before counting (via `CountKmers`), so mixed-case inputs are treated as the same k-mer. This is a benign normalization that does not change the source-defined output for canonical (already upper-case) inputs.
2. **ASSUMPTION: Empty / too-short input.** Sources state frequencies require L ≥ k but do not define the distance when a sequence has no k-mer windows. The implementation returns an empty frequency vector for such inputs, so the distance equals the Euclidean norm of the other sequence's frequency vector, and 0 when both are empty. This is the natural extension (a sequence with no words is the zero vector); marked as an assumption because no source defines it explicitly.

---

## Recommendations for Test Coverage

1. **MUST Test:** Zielezinski 2017 Figure 1 example (x="ATGTGTG", y="CATGTG", k=3) ⇒ frequency-based distance √0.11 ≈ 0.3316624790. — Evidence: Zielezinski et al. (2017) Fig. 1 + Lau et al. (2022) frequency definition.
2. **MUST Test:** Identical sequences ⇒ distance exactly 0. — Evidence: Zielezinski et al. (2017) Fig. 1.
3. **MUST Test:** Single-substitution k=1 derivation (AAAA vs AAAT) ⇒ √0.125 ≈ 0.3535533906. — Evidence: Lau et al. (2022) frequency definition.
4. **MUST Test:** Symmetry d(x,y) = d(y,x). — Rationale: Euclidean distance is a metric (symmetric).
5. **SHOULD Test:** Non-overlapping word sets (e.g. all-A vs all-T, k≥2) ⇒ √(1²+1²)=√2 (each sequence is a single distinct k-mer with frequency 1). — Rationale: maximal-disjoint case has a closed form.
6. **SHOULD Test:** k > min(length) for one sequence ⇒ that sequence has an empty vector; distance = norm of the other's frequency vector. — Rationale: documented short-input behavior (ASSUMPTION).
7. **COULD Test:** k ≤ 0 throws ArgumentOutOfRangeException. — Rationale: input validation contract inherited from CountKmers.
8. **COULD Test:** Case-insensitivity (lower vs upper) yields identical distance. — Rationale: documents the normalization assumption.

---

## References

1. Zielezinski A, Vinga S, Almeida J, Karlowski WM. (2017). Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology 18:186. https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ (DOI: 10.1186/s13059-017-1319-7)
2. Lau AK, et al. (2022). Interpreting alignment-free sequence comparison: what makes a score a good score? NAR Genomics and Bioinformatics. https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
3. Vinga S, Almeida J. (2003). Alignment-free sequence comparison—a review. Bioinformatics 19(4):513–523. https://academic.oup.com/bioinformatics/article/19/4/513/218529 (DOI: 10.1093/bioinformatics/btg005)
4. Leimeister C-A, Boden M, Horwege S, Lindner S, Morgenstern B. (2014). Fast alignment-free sequence comparison using spaced-word frequencies. Bioinformatics 30(14):1991–1999. https://academic.oup.com/bioinformatics/article/30/14/1991/2391234 (first author corrected 2026-10-01)
5. Blaisdell BE. (1986). PNAS 83:5155–5159 (squared count Euclidean d_E).
6. Torney et al. (1990); Lippert, Huang & Waterman (2005) PNAS 102:13980; Reinert, Chew, Sun & Waterman (2009) J Comput Biol 16:1615 (D2).
7. Jaccard P. (1901, 1912); Ondov BD et al. (2016) Mash, Genome Biology 17:132.

---

## Independent Cross-check (2026-09 review)

Reference: scikit-bio `Sequence.kmer_frequencies(k, overlap=True, relative=True)` (source
inspected: denominator `len(self) - k + 1`) over the union of k-mers, then
`scipy.spatial.distance.euclidean`.

| seq1 | seq2 | k | reference distance |
|------|------|---|--------------------|
| ATGTGTG | CATGTG | 3 | 0.33166247903554 (√0.11) |
| AAAA | AAAT | 1 | 0.3535533905932738 |
| AAAA | TTTT | 2 | 1.4142135623730951 |
| ACGT | AAAAAA | 5 | 1.0 |
| ACGTTGCAACGGT | ACGTAGCATCGGTA | 2 | 0.26600633232367216 |
| GATTACAGATTACA | GATTACCGATTTCA | 3 | 0.31180478223116176 |

Raw-count Counter vectors for Fig. 1 reproduce the source exactly: x = {ATG:1, GTG:2, TGT:2},
y = {ATG:1, CAT:1, GTG:1, TGT:1}; count Euclidean = √3 = 1.7320508075688772. Blaisdell (1986) /
Vinga & Almeida (2003) d_E is the squared count form (= 3) — a different variant from the
implemented frequency form (confirmed via search snippets of Höhl, Rigoutsos & Ragan 2006, which
quote d_E = Σ(c_i^X − c_i^Y)² and attribute it to Blaisdell 1986).

Note: Lau et al. (2022) write "(i.e. the sequence length minus the k-mer length)"; the number of
overlapping k-mers is L − k + 1, which is what both scikit-bio and the implementation use.

## Audit round 1 (WP2, 2026-10-01) — metrics, exact Jaccard / Mash distance, spaced words

Sources opened: Mash master `src/mash/CommandDistance.cpp` (`jaccard = common/denom`; `common == denom` → 0;
`common == 0` → 1; `distance = -log(2*jaccard/(1.+jaccard))/kmerSize`, capped at 1) and `src/mash/Sketch.cpp`
`addMinHashes` (upper-case unless `-Z`; k-mers containing a character outside the alphabet skipped;
`memcmp(kmer_fwd, kmer_rev) <= 0 ? fwd : rev` canonical unless `-n`) — via raw.githubusercontent.com. alfpy 1.0.6
sdist (PyPI) `word_distance.py` (`euclid_squared` cites Blaisdell 1986; `euclid_norm` cites Vinga & Almeida 2003;
`manhattan`, `chebyshev`, `canberra` with 0/0 = 0), `word_vector.Freqs` (total = len − k + 1), `word_sets_distance`
(set Jaccard). Leimeister et al. 2014 definition of spaced words (pattern over {0,1}, P[1] = P[ℓ] = 1, '1' = match,
'0' = don't care) from search snippets of the OUP/ResearchGate pages; the PDF (OUP, Semantic Scholar) returned 403
and spaced.gobics.de is not allow-listed, so the spaced-word reference is an independent Python replica of the
definition.

Executed references (all agree with C#; full table in docs/algorithms/K-mer/K-mer_Euclidean_Distance.md §7.2):
scipy 1.17.1 (`euclidean`, `sqeuclidean`, `cityblock`, `chebyshev`, `canberra`, `cosine`, `numpy.dot` for D2) on
scikit-bio 0.7.4 count/frequency vectors; alfpy (identical to scipy up to the last ulp); Python set Jaccard;
sourmash 4.9.4 `MinHash(n=0, ksize=k, scaled=1)` (canonical Jaccard); the Mash 2.3 binary
(`mash dist -k k -s 100000 [-n]`). Examples: Fig. 1 d_E = 3, Manhattan 0.6, Chebyshev 0.25, Canberra
1.5726495726495728, cosine 0.16666666666666663, D2 5; Jaccard literal 0.75, canonical 1.0 (Mash 0, 3/3), `-n` Mash
0.0513836 (3/4); R1/R2 k=5 canonical J = 58/158 = 0.3670886075949367, Mash 0.124338.

## Audit round 1 (WP4, 2026-10-01) — background-adjusted D2* / D2S

- Sources opened: CAFE (Lu et al. 2017, NAR 45:W554) source, cloned from github.com/younglululu/CAFE: `code/dist_model.cpp`
  (`D2starStrategy::dealWithQuad`: X̃ = X − E_X, numerator X̃Ỹ/√(E_X E_Y), normalisers X̃²/E_X; `D2sheppStrategy`:
  X̃Ỹ/√(X̃²+Ỹ²), normalisers X̃²/√(X̃²+Ỹ²); `getDist` = 0.5·(1 − num/(√Σ_X·√Σ_Y))), `code/kmer.cpp` (`getCntExpDist`:
  E = total·p, all 4^k words when lower count is 0; `getMarkovModel`; `saveFromLargerK` prefix marginal; `getEstMarkovOrder`
  BIC), `code/seq_model.cpp` (log-probability storage). Song et al. 2014 (Brief Bioinform 15:343) d2S definition via search
  snippet (centred counts estimated from each sequence, sum over A^k). Reinert et al. 2009 / Wan et al. 2010 full texts
  blocked (PMC, arXiv); the D2* form √(n̄ m̄)·p_w is the common-background special case of CAFE's √(E_X E_Y).
- alfpy 1.0.6 has no D2*/D2S (its `word_d2` is a squared Euclidean, WP2 note), so CAFE is the reference implementation.
- Executed: the CAFE binary (g++ build, Jellyfish 2.3.1) and an independent Python replica. Replica + CAFE's estimator
  (prefix-marginal counts, log-zero pruning) = CAFE output on 14 (pair, k, r) runs to 6 digits; replica + sequence
  maximum-likelihood estimator = C# on 18 rows to 1e-12 (tables: K-mer_Euclidean_Distance.md §7.4).
- Assumptions: ASM-D1 background estimated per sequence by maximum likelihood on the sequence (not CAFE's prefix marginal);
  ASM-D2 a transition probability of 1 is a valid probability (CAFE treats log 0 as missing); ASM-D3 k ≤ 12 because the
  sums run over all 4^k words; ASM-D4 NaN for a zero normaliser (0/0, as CAFE).

## Audit round 2 (WP6, 2026-10-01) — both-strand D2* / D2S, spaced-word distance

- CAFE `-R` read in `code/kmer.cpp`: `KmerModel::load` adds each count to the reverse complement's entry
  (X^R(w) = X(w) + X(RC w)); `KmerProbEnsembDelegate::getKmerlogProb` = log(½(p(w) + p(RC w))) with p(RC w) from the
  reverse-complement `KmerProbDelegate` on the same single-strand chain; `getEstMarkovOrder` is single-strand;
  `main.cpp` runs Jellyfish without `-C`. The CAFE binary with `-R` = replica (CAFE estimator) on 20 runs to 6 digits;
  replica (sequence MLE) = C# on 18 + 2 rows (K-mer_Euclidean_Distance.md §7.5).
- Spaced words: `spaced` 1.2.0 obtained from the Ubuntu archive (source tarball + binary; spaced.gobics.de not
  allow-listed, no GitHub repository found). `src/sort.h` `spacedDNA`: per-pattern Euclidean on raw counts (`-d EU`),
  JS divergence base 2 on frequencies (`-d JS`), average over patterns. Binary runs with `-r -f` = replica = scipy = C#
  on 6 (pair, pattern set) rows × 2 measures. Leimeister et al. 2014 definition via search snippets (average of the
  per-pattern distances; Euclidean / JS of relative spaced-word frequencies).
- Assumptions: ASM-D5 `-R` background = ½(p̂(w) + p̂(RC w)) with the sequence MLE chain (CAFE semantics, published
  estimator); ASM-D6 null sequence = empty for every metric; ASM-D7 spaced words are literal (no N filtering), so
  `spaced` equality is claimed on ACGT sequences only.

## Audit round 2 (WP7, 2026-10-01) — Mash MinHash sketches, containment index

- Mash v2.3 source read (raw.githubusercontent.com/marbl/Mash/v2.3): `Sketch.cpp` `addMinHashes` (upper-case,
  skip windows with a base outside the ACGT alphabet, canonical = `memcmp(fwd, rev) <= 0 ? fwd : rev`), `sketchFile`
  (records with `l < kmerSize` skipped; `reference.length += l`), `setAlphabetFromString` (`use64 = pow(alphabetSize,
  k) > 2^32`); `hash.cpp` (`MurmurHash3_x64_128`; `hash64` = first 64-bit word, `hash32` = its first 4 bytes =
  low 32 bits on little-endian); `MinHashHeap::tryInsert` (non-redundant bottom-s); `CommandDistance.cpp`
  `compareSketches` (merge until `denom = sketchSize`, then complete the union capped at s; `common == denom` → 0,
  `common == 0` → 1) and `pValue` (`gsl_cdf_binomial_Q(x − 1, r, sketchSize)`, r = pX·pY/(pX + pY − pX·pY),
  p = 1/(1 + kmerSpace/length)); `Command.cpp` (`-k` 1..32, default 21; `-s` default 1000); `MurmurHash3.cpp`.
- No MurmurHash3 existed in the repository (grep); `System.IO.Hashing` has no Murmur (XxHash/CRC only), so
  `KmerAnalyzer.MurmurHash3X64_128` was written from `MurmurHash3.cpp` and checked against mmh3 on 8 vectors.
- Reference runs: Mash 2.3 binary (`/usr/bin/mash`), 63 `mash dist` rows + 3 `mash info -d` dumps; sourmash 4.9.4
  (`MinHash(n=1000)` gives Mash's k = 21 hashes and Jaccard 0.678 = 678/1000; its 64-bit hashes differ from Mash's
  32-bit ones at k = 16: 0.738 vs 753/1000; `MinHash(scaled=1)` Jaccard 0.678580438987470 = `mash dist -s 100000`
  8069/11891; `contained_by` and CLI `compare --containment` give the containment values). All reproduced: Mash rows
  exactly at the printed 6 digits (x/s exact), containment to 1e-14 (K-mer_Euclidean_Distance.md §7.6).
- Koslicki & Zabeti 2019 (Appl Math Comput 354:206): containment index C(A, B) = |A ∩ B| / |A| (definition; the
  same quantity as sourmash `contained_by` on `scaled=1` sketches).
- Assumptions: ASM-D8 a single input string is one FASTA record (sketch length = its length incl. N); ASM-D9 both
  sketches empty → J = 0, distance 0, p-value 1 (Mash prints nan for J); ASM-D10 empty K(A) → containment 0.

## Audit round 2 (WP8, 2026-10-01) — `spaced` N-word rule and reverse-complement mode

- `spaced` 1.2.0 `src/sort.h` `spacedDNA` re-read (WP6's Ubuntu archive tarball): non-ACGT letters are stored as `N`
  (N complements to N); a word is dropped when a match position reads `N` (`correctWord`); word positions per record
  = L − ℓ + 1 incl. dropped words (JS denominator); with `revComp` (no `-r`) d[i][j], i > j, compares the forward
  counts of record i with forward + reverse-strand counts of the earlier record j (EU: |row_i − (row_j + row_j′)|;
  JS: row_i / W_i vs (row_j + row_j′) / (2 W_j)).
- Binary runs `spaced [-r] -t 1 -f patterns -d JS|EU` on S1/S2, A/B, N1/N2 (S1/S2 with N, R, Y, lower case) × 3
  pattern sets × 2 measures × 2 modes = 36 values = the Python replica (`rep.py`) = C# (`AcgtOnly`, `bothStrands`)
  to the 12 printed digits; records swapped: JS 0.632325490347 vs 0.644779801994 (order dependence).
- Supersedes ASM-D7: `spaced` equality now holds on any input with `AcgtOnly` (literal words remain the default).
  ASM-D11: `bothStrands` maps `seq1` to the first FASTA record (both strands) and `seq2` to the second (forward), the
  `spaced` matrix entry d[second][first].

## Audit round 3 (WP9, 2026-10-01) — `spaced -d EV`, spaced reader, sourmash FracMinHash, validation

- Sources opened: `spaced` 1.2.0 `src/sort.h` (Ubuntu archive `spaced_1.2.0-201605+dfsg`, Debian patch 0002 moves
  `ell = dontCare + weight − 1` out of `#ifdef _OPENMP`) — EV branch: Σ min(row_i, row_j [+ row_j′]) over words and
  patterns, `min/max = min/max(length) − ell`, q = Σ f_i f_j (strand-averaged with `revComp`), V = M/(P·min) − (2·)max·q^w,
  p = V^(1/w), d = −0.75·ln(4p/3 − 1/3), else 1.2; reader: `isalpha` letters only, `toupper`, non-ACGT → N
  (Morgenstern, Zhu, Horwege & Leimeister 2015, Algorithms Mol Biol 10:5). sourmash 4.9.4 Python `minhash.py`
  (`contained_by`, `max_containment`, `_get_max_hash_for_scaled`) and Rust core v4.9.4 (`max_hash_for_scaled`,
  `jaccard`, `add_hash`, `SeqToHashes` `std::cmp::min(kmer, krc)`, `_hash_murmur`). GSL `cdf/binomial.c`
  (`gsl_cdf_binomial_Q`: k ≥ n → 0) for Mash's p-value edge cases.
- `spaced` binary: 248 runs (EV/JS/EU × `-r`/default × 11 pairs incl. gaps, `*`, digits, N, IUPAC, lower case × up to
  7 pattern sets) = Python replica `ref.py` at the printed digits = C# within one unit of the 12th digit. E.g. A/B
  (10 %) weight-8 EV 0.117255304438 / 0.119414426485; N1/N2 P7 EV −nan / 0.767483449137; X = `ACG-TACGT`, Y =
  `ACGTTACGA` {1011, 1101} JS 0.57013316426 (pre-WP9 C# 0.25), EV 0.427666596474.
- sourmash: 40 comparisons (5 pairs × k 21/31 × scaled 1/10/100/1000; `add_sequence(force=True)`) — sketch sizes,
  `count_common`, `jaccard`, `contained_by` both ways, `max_containment` all equal C# exactly; `_max_hash` S = 7919 =
  2329428472497733 (Rust truncation; the Python helper's rounding gives …734).
- ASM-D12: EV requires equal pattern lengths and read lengths ≥ ℓ (`spaced` assumes one ℓ and reads outside the
  record otherwise). ASM-D13: `MashPValue` rejects x > s and x ≥ 1 with a zero length (unreachable from a comparison;
  Mash would return 0 / NaN). ASM-D14: d2*/d2S clamped to [0, 1] (Cauchy–Schwarz; removes −1.1e-16 rounding).

## Audit round 4 (WP10, 2026-10-01) — FracMinHash downsampling, abundance tracking, u32 scaled, Mash k range

- Sources opened: sourmash 4.9.4 Python `minhash.py` (`downsample`, `jaccard`, `similarity`, `angular_similarity`,
  `contained_by`, `max_containment`, `contained_by_weighted`, `_get_max_hash_for_scaled` / `_get_scaled_for_max_hash`),
  `commands.py` (`compare` downsamples every signature to the maximum scaled), `compare.py`; Rust core v4.9.4
  `sketch/minhash.rs` (`max_hash_for_scaled`, `scaled_for_max_hash`, `add_hash_with_abundance`, `count_common`,
  `similarity`, `angular_similarity`, `downsample_scaled`) and `signature.rs` (`SeqToHashes`: always
  `min(kmer, krc)` for DNA; `force` skips, otherwise `InvalidDNA`). Mash `src/mash/Command.cpp` (`-k`: Integer, 1..32)
  and `sketchParameterSetup.cpp`.
- sourmash (executed): `downsample(scaled=S′)._max_hash` = Rust truncation (2049 → 9002803354665472, 7919 →
  2329428472497733, 123456789 → 149418628356, 4294967295 → 4294967297); C k = 21 S 10 → 1000: 1011 → 15 hashes = the
  S = 1000 sketch. Downsample comparisons (x at S_a, y at S_b): A/B k 21 10/100 jaccard 0.18322981366459629,
  `compare --containment` 0.3155080213903743 / 0.30412371134020616 (method `contained_by(downsample=True)` 0.029064039408866996
  — denominator not downsampled), max 0.3155080213903743; C/B k 31 100/10 0.06985294117647059, 0.19791666666666666 /
  0.09743589743589744; B/A k 31 1000/10 0.06060606060606061, 0.1250000139547377 / 0.10526315847892492.
- Abundance (`track_abundance=True`): P/Q k 21 S 1 angular 0.4701473676328215 (Jaccard 0.63748031496063, weighted
  0.8236686390532545 / 0.8590686274509803); P/W k 15 0.7252694187711534; W/Q k 21 S 10/100 (downsample) 0.6670115706749962;
  A/B k 21 S 10 0.2325202182944529; S1/S2 k 4 S 1 / 3: 0.2363801370444173 / 0.3123095603640216. All equal C# within 1e-15.
- ASM-D15: with `downsample`, containments follow `sourmash compare` (downsample both); the Python methods' mixed-scale
  denominator is not reproduced. ASM-D16: MCP `kmer_jaccard` with `scaled` requires `canonical = true` (no non-canonical
  sourmash DNA mode); non-ACGT k-mers are skipped (`force=True`). ASM-D17: `MashPValue` k ∈ 1..32 (Mash).

## Audit round 5 (WP11, 2026-10-01) — MCP weighted containments, MinHash sketch K / Use64 validation

- sourmash 4.9.4 (pip, executed): `MinHash(0, 4, scaled=S, track_abundance=True).add_sequence(force=True)` on the MCP
  S1/S2 inputs: `S1.contained_by_weighted(S2)` / `S2.contained_by_weighted(S1)` = 0.37662337662337664 / 0.417910447761194
  (S = 1), 0.4642857142857143 / 0.5238095238095238 (S = 3); `angular_similarity` 0.2363801370444173 / 0.3123095603640216
  (unchanged). MCP `kmer_jaccard` returns them as `weightedContainmentSeq1InSeq2` / `weightedContainmentSeq2InSeq1`.
- Mash (raw.githubusercontent.com marbl/Mash master): `src/mash/Sketch.cpp` line 1136 `parameters.use64 =
  pow(parameters.alphabetSize, parameters.kmerSize) > pow(2, 32)` (DNA: k > 16; 32-bit sketches store uint32 hashes);
  `Command.cpp` `-k` range 1..32. ASM-D18: `CompareMinHashSketches` rejects sketches with K ∉ 1..32, `Use64` ≠ K > 16
  or a hash > 2^32 − 1 in a 32-bit sketch (`ArgumentException`); `FromHashes` rejects a 64-bit value for k ≤ 16.

## Change History

- **2026-10-01**: Audit round 5 WP11 — MCP `kmer_jaccard` weighted containments; `CompareMinHashSketches` / `FromHashes` K, `Use64` and 32-bit hash validation.

- **2026-10-01**: Audit round 4 WP10 — `DownsampleFracMinHash` / `downsample`, `trackAbundance` + angular similarity / weighted containment, u32 `scaled`, `MashPValue` k range, MCP `kmer_jaccard` canonical requirement + `trackAbundance`.

- **2026-10-01**: Audit round 3 WP9 — `spaced -d EV` (`SpacedEvolutionary`), spaced reader in the `AcgtOnly` path, sourmash FracMinHash sketches, `MashPValue` / sketch validation, d2*/d2S clamp.

- **2026-10-01**: Audit round 2 WP8 — `spaced` N-word rule (`AcgtOnly`) and reverse-complement mode (`bothStrands`) for `SpacedWordDistance`; 36 binary runs.

- **2026-10-01**: Audit round 2 WP7 — Mash MinHash sketches (MurmurHash3_x64_128, x/s, distance, p-value), exact containment index.

- **2026-10-01**: Audit round 2 WP6 — CAFE `-R` both-strand D2*/D2S, sparse Markov tables, JS / count-Euclidean metrics, multiple-pattern spaced-word distance.

- **2026-10-01**: Audit round 1 WP2 — metric variants, exact Jaccard / Mash distance, spaced words; reference [4] first author corrected.

- **2026-09-28**: 2026-09 review — added scikit-bio/scipy cross-check table and variant note.

- **2026-06-13**: Initial documentation.
