# Overrepresented k-mer Motif Discovery

| Field | Value |
|-------|-------|
| Algorithm Group | Matching / Motif Discovery |
| Test Unit ID | MOTIF-DISCOVER-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Production |
| Last Reviewed | 2026-10-01 |

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
- `MotifFinder.DiscoverMotifs(DnaSequence, int k, int minCount, OligoBackgroundModel, OligoStrandMode = Single, bool countOverlapping = true)` → `OligoAnalysisResult` ([MotifFinder.OligoAnalysis.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysis.cs)): full RSAT `oligo-analysis` occurrence statistics — exp_freq, exp_occ, ratio, binomial `occ_P`/`occ_E`/`occ_sig` — with backgrounds `OligoBackgroundModel.Equiprobable` (`-bg equi`), `.BernoulliFromInput` (default `-bg input`), `.Bernoulli(acgt)`, `.MarkovFromInput(m)` (`-markov m`), `.MarkovFromOligoFrequencies(table, ψ = 0.01, strandInsensitive)` (`-bgfile`, `RSAT::MarkovModel`), strands `-1str`/`-2str`, and `-ovlp`/`-noov`. The binomial tail is the shared `StatisticsHelper.LogBinomialUpperTail` (Loader 2000 saddle-point terms, log space) [5][6].
- MCP `discover_motifs` (Seqeron.Mcp.Analysis `AnalysisTools.DiscoverMotifs`) delegates to the uniform overload.
- `MotifFinder.AnalyzeOligos(IEnumerable<DnaSequence>, int k, OligoAnalysisOptions?)` → `OligoAnalysisReport` ([MotifFinder.OligoAnalysisOptions.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.OligoAnalysisOptions.cs)): the RSAT command-line options on a sequence set — `MinCount` (`-lth occ`), `Background` (incl. `OligoBackgroundModel.Lexicon`, `-lexicon`), `Strands`, `CountOverlapping`, `PseudoFrequency` (`-pseudo`), `Degeneracy` (`-oneN` / `-onedeg`), `Calibration` (`OligoCalibration.Parse` / `FromEntries`, `-calibN` / `-calib1`); per pattern occ, overlaps (`ovl_occ`), obs_freq, exp_freq, exp_occ, exp_var, overlap coefficient, z-score, ratio, occ_P / occ_E / occ_sig, fitted distribution and lexicon segmentation. With default options it returns exactly `DiscoverMotifs` (tested bit for bit). Calibration P-values use the shared `StatisticsHelper.LogPoissonRangeProbability` / `LogNegativeBinomialRangeProbability`.
- `MotifFinder.FindSharedMotifs(…, OligoStrandMode, double pseudoFrequency)`: matching-sequence statistics with `-pseudo`.
- `MotifFinder.AnalyzeDyads(IEnumerable<DnaSequence>, DyadAnalysisOptions?)` → `DyadAnalysisReport` ([MotifFinder.DyadAnalysis.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/MotifFinder.DyadAnalysis.cs)): RSAT `dyad-analysis` (spaced dyads, §5.4).
- MCP `oligo_analysis` (RSAT overload; with `extraSequences` / `zscore` / `expectedFrequencyPseudo` / `degenerate` / `calibrationTable` / background `lexicon` → `AnalyzeOligos`) and `dyad_analysis` (→ `AnalyzeDyads`).

### 5.2 Current Behavior

The expected count is computed once per call from the closed-form `(N − k + 1) / 4^k`; there is no clamp/floor on the denominator (a previous `max(E, 0.1)` floor was an untraceable value and was removed — `E` is always strictly positive when any k-mer exists, INV-04). Counting allows overlapping occurrences. Results are yielded in order of each k-mer's first occurrence; positions are ascending.

**Search reuse:** The suffix tree was evaluated. Motif *discovery* here requires counting *all distinct k-mers and their positions in one pass* — not searching for a known pattern — so a single linear scan with a hash map is the appropriate structure; the suffix tree (best for many queries of known patterns against one text) is not used for this method.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Expected count `E = (N − k + 1) / 4^k` under the i.i.d. uniform background [1].
- Overrepresentation as the observed/expected ratio [1][2].
- Overlapping window enumeration of k-mers [1].

- Bernoulli (independent, non-uniform) background, `E = (N − k + 1) · ∏ q[w_i]` [3][4] (overload; uniform default unchanged).
- RSAT significance (`OligoBackgroundModel` overload) [3][4][5]: with n = `sum_occurrences` = N − k + 1 windows (also with `-noov`) and p = exp_freq, `occ_P = P(X ≥ occ)`, X ~ Binomial(n, p) (`RSAT::stats::sum_of_binomials`, right tail); `occ_E = occ_P × T` with T = number of patterns tested after the `-lth occ` threshold (`MultiTestCorrections($nb_tested_patterns)` — the value RSAT prints; the manual's NPO, `nb_possible_oligos`, is exposed as `OligoAnalysisResult.PossibleOligos` and equals T with `-zeroocc`); `occ_sig = −log10 occ_E`.
- Markov backgrounds of order m (`-markov m`): exp_freq(w) = f(w[0..m]) · ∏_{o≥1} f(w[o..o+m]) / f(w[o..o+m−1]) with f the relative frequencies of all overlapping (m+1)- and m-mers of the input, single strand (`CalcSubWordFrequencies`, `CalcExpected`); order constraint m ≤ k − 2 for m > 0. From a table (`-bgfile`): `RSAT::MarkovModel` — P(b|x) = (1 − ψ)·f(xb)/S(x) + ψ/4, P(x) = (1 − ψ)·S(x)/F + ψ/4^m, ψ = 0.01; strand-insensitive tables are halved for non-palindromes (`load_from_file_oligos`).
- Both strands (`-2str`, `-grouprc`): occ(W|W') = occ(W) + occ(W'), palindromes once; exp_freq(W|W') = exp_freq(W) + exp_freq(W') (manual, "Specific treatment for double strand counts"); input Bernoulli pools complementary residues; NPO = (4^k + 4^{k/2})/2 (even k) or 4^k/2; pattern reported under the lexicographically smaller member.
- `-noov`: windows read from the 3′ end; an occurrence closer than k to the last counted occurrence of the same word (or, with `-2str`, of its reverse complement) is discarded.

**Deliberate deviations from the RSAT code (documented formula / exact arithmetic kept):**

- `-2str` with a strand-asymmetric background: current RSAT code copies the kept member's expected frequency to its partner (2·exp_freq(min(W, W'))) instead of the documented sum; on t2.fa `-markov 1` it prints `aaat|attt` 0.0166535 = 2 × 0.0083268, and with `-markov 2` exp_freq 0 (probability "NA") for observed pairs. Both coincide for strand-symmetric backgrounds (equiprobable, input Bernoulli), where the RSAT run is reproduced exactly.
- `segment_proba` returns the probability formatted with `%5g` (6 significant digits); the unrounded value is used.
- occ_P is evaluated in log space, so occ_sig stays finite and exact where RSAT's occ_P underflows (RSAT prints `$MAX_SIG` = 350).
- For k = 1 with `-markov 0`, RSAT has no 1-mer sub-word table (exp_freq undefined); the model is the residue composition here.

**RSAT options (`AnalyzeOligos`, rsat-code 10043f2, oligo-analysis 1.169) [5]:**

- `-return zscore` (`CalcZscore`): z = (occ − exp_occ)/√exp_var with exp_var = notforb · p · (2·ovlp − 1 − (2k + 1)·p) (overlapping counts; notforb = nb_pos), exp_var = exp_occ under `-noov` (Poisson-like, no self-overlap), or the calibration variance; NaN ("NA") when exp_var ≤ 0. ovlp (`OverlapCoeff`, Pevzner, Borodovsky & Mironov 1989) = 1 + Σ over the periods i of w of ∏_{j<i} q(w_j), with q the Bernoulli residue probabilities (input composition for the default, lexicon and calibration runs; ¼ for Markov models, where RSAT's `%residue_proba` is empty); reverse-complement overlaps are not added (RSAT never sets `$sum_strands`). RSAT also emits `exp_var` and `ovlp` columns.
- `-pseudo ψ`: applied in `CalcExpected` after any background (Bernoulli, Markov, lexicon, calibration) and before the reverse-complement sum: exp_freq ← (1 − ψ)·exp_freq + ψ/NPO, NPO the strand- and degeneracy-aware `NbPossibleOligos`; with a calibration only the displayed exp_freq changes (P-values use exp_occ / exp_var). Also on `FindSharedMotifs`.
- `-oneN` / `-onedeg` (`Degenerate`): occ(D) = Σ occ(w) over the words matching D at its single IUPAC position (N, or R Y W S M K H B V D N); NPO = k·|codes|·4^(k−1) (RSAT, then the `-2str` reduction with RSAT's palindrome count).
- `-lexicon`: sub-word frequencies from the (occ + overlaps) prefixes of the counted k-mers (no trailing sub-words), maximal segmentation frequencies M per length, exp_freq(w) = max_s M(prefix)·M(suffix) (XML doc of `OligoBackgroundModel.Lexicon`); the best split is reported (RSAT `segments`).
- `-calibN` / `-calib1` (`ReadCalibration`, `CalibrateSetFromSingleSequence`, `CalcProba`): the table (RSAT `fit-distribution` output: pattern [id] avg std var …) gives exp_occ and exp_var (× number of sequences for `-calib1`; a word whose reverse complement has no positive mean passes its values to it); exp_freq = exp_occ / (nb_pos − forbidden positions); occ_P = P(occ ≤ X ≤ n) under a negative binomial with p = var/mean − 1, size = mean/p (`sum_of_negbin2`) when mean < var, else Poisson(mean) (`sum_of_poisson`).

**Deliberate deviations for the options (RSAT code defects, each reproduced by a run of the unmodified script):**

- `-oneN` / `-onedeg` produce no output in RSAT 1.169: `Degenerate` reads the global `%IUPAC`, which is never filled (`&IUPAC()` returns a lexical hash), and assigns scalars to `%patterns`; `oligo-analysis -i t2.fa -l 3 -onedeg -return occ,proba` prints "oligomers tested for significance 0" and no row. The documented behaviour (manual: "successively insert one ambiguous nucleotide code at each position of each pattern") is implemented: exp_freq(D) = Σ exp_freq(w) over the matching words (RSAT's own aggregation for frequency files; identical to RSAT's IUPAC residue products for Bernoulli models), n = the number of windows (each window matches D at most once; RSAT's `CalcOccSum` would sum the degenerate counts, 7k × the windows for `-onedeg`), and the overlap coefficient uses the IUPAC residue sums. Cross-checked against RSAT with only those lines fixed (60 random runs, ≤ 1e-13).
- `sum_of_poisson` keeps `$prev_value` as a package global, so a Poisson tail stops early when the previous call ended above the current partial sum (t2.fa, `-l 3 -calibN`: tgc Poisson(1.5), occ 6 → RSAT 0.003529988861723204, exact 0.0044559807752478 = scipy / mpmath). The exact sum is used.
- Negative-binomial terms are rounded to 5 significant digits by `LogToEng` (atg: RSAT 0.0067141894 vs exact 0.0067142352); the exact sum is used.
- With `-pseudo`, an observed word absent from a calibration table gets exp_freq > 0 but exp_occ = 0, and RSAT dies (`sum_of_poisson` → undefined `&Warning`); here such a pattern is not tested.
- With `-2str` and the input Bernoulli model, a residue absent from the input gets no probability in RSAT (`CalcAlphabet` loops over the observed residues only) although its pooled value (occ(b) + occ(b̄))/2N is positive (e.g. no C but some G); the pooled value is used for exp_freq and the overlap coefficient (the F21 behaviour).
- RSAT also lists the reverse complements of observed words as empty "NA" rows in calibration runs (autovivified by `ReadCalibration`) and drops occurrence-0 calibration rows whose rounded occ_P exceeds 1 (`-uth occ_P 1`); neither affects any statistic or the tested count, and they are not reported here.

### 5.4 RSAT dyad-analysis (`AnalyzeDyads`)

Port of `perl-scripts/dyad-analysis` 1.78 (van Helden, Rios & Collado-Vides 2000) [5][7]. A dyad M₁ n{s} M₂ joins two monads of length m with s unspecified residues (s ∈ [min, max]); types any / dr (M₂ = M₁) / ir (M₂ = rc M₁) / rep (dr or ir). Per spacing, T_s = Σ max(0, L − 2m − s + 1) positions (the binomial trials), occ_sum_s / ovl_sum_s the counted / overlap-discarded single-strand occurrences; `-noov` (RSAT default) discards an occurrence closer than 2m + s to the last counted occurrence of the same dyad (or, with both strands, of its reverse complement). Monad frequencies f(M) = occ(M)/Σ over all monad windows (single strand, overlapping); exp_freq(D) = f(M₁)f(M₂) (+ f(rc M₁)f(rc M₂) with both strands unless D is a reverse palindrome), or from a dyad / monad frequency table (`-expfreq`, `-mncf`; a missing dyad falls back to monads when max spacing > 20 or occ = 0, otherwise it is not tested); exp_occ = exp_freq · occ_sum_s; obs_freq = occ / (occ_sum_s + ovl_sum_s); ovlp = `OverlapCoeff` of M₁M₂ with equiprobable residues; exp_var = T_s · p · (2·ovlp − 1 − (4m + 1)·p); occ_P = P(Bin(T_s, exp_freq) ≥ occ); occ_E = occ_P × tested; `MinCount ≤ 0` tests every (observed monad, spacing, allowed second monad) combination as RSAT does without `-lth occ`. RSAT defaults are the option defaults (m = 3, spacing 0–20, `-2str`, `-noov`).

Cross-check: RSAT dyad-analysis itself, single strand (all types, `-ovlp`/`-noov`, thresholds 0–3, 1–3 sequences, m 1–3, spacings 0–7): 200 runs, 14,746 dyads — every count, per-spacing sum, exp_freq, exp_occ, variance, ov_coef, ratio, occ_P, occ_E, occ_sig ≤ 1.3e-14 (z-score to RSAT's 2 printed decimals) except the `rep` / z-score-request artefacts below; `-expfreq` tables 40/40 (≤ 1.6e-14); both strands vs an independent Python port of the pair rules 150/150 (≤ 5.2e-15) and vs RSAT on every pair RSAT reports (100 runs). Deliberate deviations (RSAT code defects):

- With `-2str`, `SumRCDyads` sums a pair only when the lexicographically smaller member is a key, so pairs seen only as their larger member are lost (t2.fa `-l 3 -sp 0-2 -lth occ 2`: aagn{0}cat|atgn{0}ctt and aagn{1}atg|catn{1}ctt, 2 occurrences each; 9 tested instead of 11; 19 of 139 pairs at `-ovlp -lth occ 1`). Every observed pair is summed here.
- With `-type rep`, `SecondElement(M)` returns (M, rc M), i.e. M twice for a reverse-palindromic monad, so `CalcOccSum` and the tested count include those dyads twice (40/40 brute-force runs). They are counted once here.
- When z-scores are requested (`-return zscore`), RSAT deletes a dyad whose variance estimate is ≤ 0 (z = "NA", `CheckPatternThresholds`) before its P-value, so its tested set depends on the output columns (e.g. A n{3} A in an A-rich sequence: 15 tested with zscore, 16 without); such a dyad is tested here (z = NaN), as RSAT does without `-return zscore`.
- Without an occurrence threshold RSAT leaves exp_occ empty for unseen combinations (its z-score and ratio for them are 0); exp_occ = exp_freq · occ_sum_s is reported. z-scores are not rounded (RSAT `%7.2f`). Background-table monad frequencies are used for every monad (RSAT assigns them only to monads seen in the input).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null sequence | ArgumentNullException | Validation contract |
| k < 1 | ArgumentOutOfRangeException | Validation contract |
| k > N | empty result | No length-k windows exist [1] |
| Homopolymer "AAAA…" | the single k-mer dominates with high enrichment | All windows are identical |

### 6.2 Limitations

The legacy overloads report the O/E ratio only (uniform or Bernoulli background); significance, Markov / lexicon backgrounds, both-strand counting and the RSAT options (`-zscore`, `-pseudo`, `-oneN`/`-onedeg`, `-calibN`/`-calib1`, multi-sequence input) are in the `OligoBackgroundModel` overload and `AnalyzeOligos` (§5.3), spaced dyads in `AnalyzeDyads` (§5.4). DNA alphabet only. RSAT's organism-specific frequency and calibration files (`-bg upstream …`) belong to an RSAT installation's genome data; pass such a table to `MarkovFromOligoFrequencies`, `OligoCalibration` or `DyadBackgroundModel`. The RSAT code defects listed in §5.3 / §5.4 are not reproduced.

## 7. Examples and Related Material

### 7.1 Worked Example

**Numerical walk-through:** Sequence `ATGCATGCATGC` (N=12), k=4. Windows = 12 − 4 + 1 = 9. Expected count `E = 9 / 4^4 = 9/256 = 0.03515625`. The k-mer `ATGC` occurs at positions 0, 4, 8 (Count = 3). Enrichment = `3 / (9/256) = 768/9 ≈ 85.333` [1].

With the Bernoulli background `q = (0.3, 0.2, 0.2, 0.3)`: `p(ATGC) = 0.3·0.3·0.2·0.2 = 0.0036`, `E = 9 · 0.0036 = 0.0324`, enrichment = `3 / 0.0324 = 92.5926` (TGCA/GCAT/CATG: `2 / 0.0324 = 61.7284`) [3][4].

**RSAT significance (run of RSAT oligo-analysis 1.169):** `ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC`, k = 4, `-1str -bg equi -lth occ 2`: n = 60, 10 tested; `ATGC` occ 6, exp_freq 1/256, `occ_P = P(Bin(60, 1/256) ≥ 6) = 1.4844942103072834e-07`, `occ_E = 1.48449e-06`, `occ_sig = 5.8284214918614703`. `-2str` (input Bernoulli): `atgc|gcat` occ 9, exp_freq 0.0077771093320170084, occ_P 1.0760647330191343e-09, occ_sig 7.7920703428970466 (15 tested, NPO 136).

**API usage example:**

```csharp
var motifs = MotifFinder.DiscoverMotifs(new DnaSequence("ATGCATGCATGC"), k: 4, minCount: 2);
var atgc = motifs.First(m => m.Sequence == "ATGC");
// atgc.Count == 3, atgc.Positions == [0,4,8], atgc.Enrichment == 768.0/9
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [MotifFinder_DiscoverMotifs_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_DiscoverMotifs_Tests.cs) — covers `INV-01`..`INV-04`; [MotifFinder_OligoAnalysis_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/MotifFinder_OligoAnalysis_Tests.cs) — RSAT-run occ_P/occ_E/occ_sig, Markov, `-2str`, `-noov`; property/metamorphic/fuzz: `MotifOligoAnalysisProperties`, `MotifOligoAnalysisMetamorphicTests`, `MotifOligoAnalysisFuzzTests`; binomial tail: `StatisticsHelper_BinomialUpperTail_Tests`
- Evidence: [MOTIF-DISCOVER-001-Evidence.md](../../../docs/Evidence/MOTIF-DISCOVER-001-Evidence.md)

## 8. References

1. Compeau P, Pevzner P. 2015. *Bioinformatics Algorithms: An Active Learning Approach*, 2nd ed., Ch. 2 (Finding Regulatory Motifs). Active Learning Publishers. Formula and worked example reproduced at https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string
2. fmicompbio. monaLisa `getKmerFreq` — observed vs expected k-mer frequencies and log2 enrichment. https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html
3. van Helden J, André B, Collado-Vides J. 1998. Extracting regulatory sites from the upstream region of yeast genes by computational analysis of oligonucleotide frequencies. *J Mol Biol* 281:827–842.
4. RSAT `oligo-analysis` source — https://raw.githubusercontent.com/rsa-tools/rsat-code/master/perl-scripts/oligo-analysis (Bernoulli `exp_freq *= residue_proba`, `exp_occ = exp_freq * sum_occurrences`, `sum_of_binomials` occ_P) and `perl-scripts/lib/RSA.disco.lib` `MultiTestCorrections` (occ_E, occ_sig); opened 2026-09-29.
5. RSAT source, rsa-tools/rsat-code master 10043f2 (2026-09-23), git clone run with perl 5.38: `perl-scripts/oligo-analysis` (v1.169: `CountOligos`, `SumReverseComplements`, `CalcSubWordFrequencies`, `CalcExpected`, `CalcProba`), `perl-scripts/lib/RSA.disco.lib` (`NbPossibleOligos`, `MultiTestCorrections`, `GroupRC`), `perl-scripts/lib/RSA.lib` (`SumExpectedFrequencies`, `ReadPatternFrequencies`), `perl-scripts/lib/RSAT/stats.pm` (`sum_of_binomials`, `binomial_boe`), `perl-scripts/lib/RSAT/MarkovModel.pm` (`load_from_file_oligos`, `add_pseudo_freq`, `normalize_transition_frequencies`, `segment_proba`); opened 2026-09-30.
6. Loader C. 2000. Fast and accurate computation of binomial probabilities (R nmath `dbinom_raw`, `stirlerr`, `bd0`).
7. van Helden J, Rios AF, Collado-Vides J. 2000. Discovering regulatory elements in non-coding sequences by analysis of spaced dyads. *Nucleic Acids Res* 28:1808–1818. RSAT source rsa-tools/rsat-code master 10043f2: `perl-scripts/dyad-analysis` v1.78, `perl-scripts/oligo-analysis` v1.169 (`Degenerate`, `CalcSubWordFrequencies`, `CalcExpected` lexicon and `-pseudo` blocks, `CalcOverlapCoefficient`, `CalcZscore`, `CalibrateSetFromSingleSequence`, `CalcProba`), `perl-scripts/calibrate-oligos`, `perl-scripts/fit-distribution`, `perl-scripts/lib/RSA.lib` (`ReadCalibration`, `ReadExpectedFrequencies`), `perl-scripts/lib/RSA.seq.lib` (`OverlapCoeff`, `IUPAC`, `ReverseComplement`), `perl-scripts/lib/RSA.disco.lib` (`NbPossibleOligos`, `GroupRC`), `perl-scripts/lib/RSAT/stats.pm` (`sum_of_poisson`, `sum_of_negbin2`, `LogToEng`); opened and run (perl 5.38) 2026-10-01.
8. Pevzner PA, Borodovsky MY, Mironov AA. 1989. Linguistics of nucleotide sequences I: the significance of deviations from mean statistical characteristics and prediction of the frequencies of occurrence of words. *J Biomol Struct Dyn* 6:1013–1026 (overlap coefficient).
