# Evidence Artifact: MOTIF-DISCOVER-001

**Test Unit ID:** MOTIF-DISCOVER-001
**Algorithm:** Motif Discovery via Overrepresented k-mers (observed/expected enrichment)
**Date Collected:** 2026-06-14

---

## Online Sources

### Compeau & Pevzner, *Bioinformatics Algorithms* — "Expected number of occurrences of a k-mer in a DNA string"

**URL:** https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string
**Accessed:** 2026-06-14 (retrieved with WebFetch; this wiki page reproduces the formula and worked example from Compeau & Pevzner, *Bioinformatics Algorithms: An Active Learning Approach*, "Finding Regulatory Motifs in DNA Sequences")
**Authority rank:** 1 (peer-reviewed textbook content, reproduced verbatim)

**Key Extracted Points:**

1. **Background model:** Retrieved text states "the sequences are formed by selecting each nucleotide (A, C, G, T) with the same probability (0.25)." This is the zero-order i.i.d. uniform background.
2. **Probability formula:** Verbatim from the page: `Pr(N,A,Pattern,t) ≈ ( N−t·(k−1) | t )/A^(t·k)`, where N = sequence length, A = alphabet size (=4 for DNA), k = pattern length, t = number of occurrences, and `( m | n )` is a binomial coefficient.
3. **Expected count (t = 1):** For a single expected occurrence the binomial `( N−(k−1) | 1 ) = N−k+1`, so the expected number of occurrences of a specific k-mer is `(N − k + 1) / 4^k`: the number of length-k windows `N − k + 1` divided by the number `4^k` of distinct DNA k-mers.
4. **Worked example:** Retrieved text: `Pr(1000, 4, 9, 1)·500 ≈ 1.9` — a random 9-mer is expected ≈ 2 times across 500 sequences of length 1000.

### RSAT oligo-analysis source (review 2026-09-29)

**URL:** https://raw.githubusercontent.com/rsa-tools/rsat-code/master/perl-scripts/oligo-analysis, `perl-scripts/lib/RSA.disco.lib`, `perl-scripts/lib/RSAT/stats.pm` (opened with curl, 2026-09-29)
**Authority rank:** 1 (original code of van Helden, André & Collado-Vides 1998, J Mol Biol 281:827)

**Key Extracted Points:**

1. Equiprobable model: `$common_exp_freq = 1/($alphabet_size**$oligo_length)`.
2. Bernoulli model: `exp_freq = 1; foreach residue: exp_freq *= $residue_proba{$nt}`.
3. `exp_occ = exp_freq * $sum_occurrences`, where `sum_occurrences = sum_overlaps + sum_noov` (all overlapping windows = N − k + 1 for one strand).
4. Significance (implemented 2026-09-30, B05 F21 — `DiscoverMotifs(seq, k, minCount, OligoBackgroundModel, OligoStrandMode, countOverlapping)` → `OligoAnalysisResult` occ_P/occ_E/occ_sig; RSAT options `-zscore`, `-pseudo`, `-oneN`/`-onedeg`, `-lexicon`, `-calibN`/`-calib1` added by F31 via `AnalyzeOligos(sequences, k, OligoAnalysisOptions)`; originally declared not implemented): `occ_P = sum_of_binomials(exp_freq, sum_occurrences, occ, sum_occurrences)`; `occ_E = occ_P × nb_tested_patterns`; `occ_sig = −log10(occ_E)`.

Reference values (exact Python `Fraction` re-computation of the RSAT formulas): uniform "ATGCATGCATGC", k=4 → ATGC 85.33333333333333; Bernoulli (0.3,0.2,0.2,0.3) → ATGC 92.5925925925926, TGCA/GCAT/CATG 61.72839506172839; "A"×10, k=3 → 37.03703703703704; LCG 512-mer X, X+X, k=512 → 7.008550233381348e+305 (uniform), 3.2383751592973165e+306 (bg 26/24/24/26).

### O/E ratio confirmation (search corroboration)

**URL:** WebSearch query "expected number of occurrences k-mer in random sequence overrepresented motif observed expected ratio formula" (2026-06-14); top results include monaLisa `getKmerFreq` and the PeerJ supplemental "Ratio of Observed and Expected k-mer counts".
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation / methods documentation)

**Key Extracted Points:**

1. **O/E ratio:** The retrieved search summary states the standard overrepresentation statistic is the "ratio of observed to expected frequency of a k-mer (O/E ratio), where O is the observed frequency." Enrichment > 1 ⇒ the k-mer is overrepresented.
2. **Expected = length/4^k:** The retrieved summary states "The expected number of occurrences of a k-mer is N/4^k" for the zero-order uniform model (consistent with `(N−k+1)/4^k` once the exact window count is used).

---

## Documented Corner Cases and Failure Modes

### From Compeau & Pevzner

1. **Self-overlap approximation:** The probability formula "is only an approximation because it assumes that the pattern cannot overlap with itself, which is not always the case." This affects the *probability* statistic, not the deterministic observed count or the expected-count denominator used for the O/E ratio.
2. **k > N (no windows):** When `k > N` there are zero length-k windows (`N − k + 1 ≤ 0`); no k-mer can be counted, so the discovery returns nothing.

---

## Test Datasets

### Dataset: Tandem-repeat string "ATGCATGCATGC"

**Source:** Derived directly from the Compeau & Pevzner expected-count formula `(N−k+1)/4^k`.

| Parameter | Value |
|-----------|-------|
| Sequence | ATGCATGCATGC |
| N (length) | 12 |
| k | 4 |
| Window count N−k+1 | 9 |
| Expected count E = 9/4^4 | 9/256 = 0.03515625 |
| Observed count of "ATGC" | 3 (positions 0, 4, 8) |
| Enrichment = 3 / (9/256) | 768/9 = 85.3333… |

### Dataset: Homopolymer "AAAAAAAAAA"

**Source:** Derived from `(N−k+1)/4^k`.

| Parameter | Value |
|-----------|-------|
| Sequence | AAAAAAAAAA |
| N (length) | 10 |
| k | 3 |
| Window count N−k+1 | 8 |
| Expected count E = 8/4^3 | 8/64 = 0.125 |
| Observed count of "AAA" | 8 (positions 0..7) |
| Enrichment = 8 / 0.125 | 64.0 |

---

## Assumptions

1. **ASSUMPTION: minCount filter** — The `minCount` parameter (default 2) that filters which k-mers are returned is a presentation threshold, not part of the published statistic; the formal O/E value is defined for every k-mer regardless of the cutoff. Changing it only includes/excludes rows, never alters a returned k-mer's Count, Positions, or Enrichment, so it is not correctness-affecting for any returned record.

---

## Recommendations for Test Coverage

1. **MUST Test:** Exact O/E enrichment for "ATGC" in "ATGCATGCATGC" (k=4) equals 768/9 — Evidence: Compeau & Pevzner expected-count formula.
2. **MUST Test:** Exact O/E enrichment for "AAA" in "AAAAAAAAAA" (k=3) equals 64.0 — Evidence: same formula.
3. **MUST Test:** Observed count and positions of a repeated k-mer (e.g. "ATGC" at 0/4/8) — Evidence: deterministic window enumeration.
4. **SHOULD Test:** minCount filter excludes k-mers below the threshold — Rationale: documented filter semantics.
5. **SHOULD Test:** null sequence → ArgumentNullException; k < 1 → ArgumentOutOfRangeException — Rationale: documented validation.
6. **COULD Test:** k > N returns no motifs — Rationale: no windows exist (corner case).

---

## RSAT oligo-analysis run (review 2026-09 follow-up, 2026-09-30)

**Source opened:** git clone of rsa-tools/rsat-code master 10043f2 (2026-09-23): `perl-scripts/oligo-analysis` (v1.169),
`perl-scripts/lib/RSA.disco.lib` (`NbPossibleOligos`, `MultiTestCorrections`, `GroupRC`), `perl-scripts/lib/RSA.lib`
(`SumExpectedFrequencies`, `ReadPatternFrequencies`), `perl-scripts/lib/RSA.seq.lib` (`SmartRC`),
`perl-scripts/lib/RSAT/stats.pm` (`sum_of_binomials`, `binomial_boe`), `perl-scripts/lib/RSAT/MarkovModel.pm`
(`load_from_file_oligos`, `add_pseudo_freq`, `normalize_transition_frequencies`, `segment_proba`).

**How it was run:** `RSAT=<clone> perl -I perl-scripts/lib perl-scripts/oligo-analysis -i x.fa -l k ...` (an empty
`RSAT_config.props`; no other installation needed), plus a copy with one debug line after `MultiTestCorrections`
printing every pattern's exp_freq, occ_P, occ_E, occ_sig, exp_ms, ms_P, ms_E, ms_sig with `%.17g`.

**Source facts established (code, confirmed by the runs):**
- Trials n = `sum_occurrences` = all N − k + 1 windows (also with `-noov`, where overlaps are added back); `-2str` pairs
  are grouped before summing, so n is not doubled.
- occ_P = `sum_of_binomials(exp_freq, sum_occurrences, occ, sum_occurrences)` (right tail, exact summation).
- occ_E = occ_P × `$nb_tested_patterns` — the number of patterns reaching `CalcProba`, i.e. after `-lth occ`
  (the manual's "NPO" equals it only with `-zeroocc`); occ_sig = −log10 occ_E (350 when occ_E = 0).
- ms_P = `binomial_boe(exp_ms/S, S, mseq)`, exp_ms = S·(1 − (1 − exp_freq)^(nb_pos/S)), S = all sequences,
  nb_pos over sequences with L ≥ k; ms_E = ms_P × `nb_possible_oligos` (4^k; with `-2str` (4^k + 4^{k/2})/2 for even k, 4^k/2 for odd k).
- `-markov m` (slow mode): sub-word frequencies of all overlapping (m+1)- and m-mers, single strand;
  exp_freq = f(w[0..m]) ∏ f(w[o..o+m]) / f(w[o..o+m−1]); m ≤ k − 2 for m > 0.
- `-2str`: occ summed with the reverse complement (palindromes once); input Bernoulli pools complementary residues.
  The code copies the kept member's exp_freq to its partner (2·p(min)), contradicting the manual's
  exp_freq(W) + exp_freq(W') for asymmetric backgrounds (e.g. `-markov 2` gives exp_freq 0 for observed pairs);
  the documented sum is implemented.
- `-bgfile` probabilities are returned as `%5g` strings by `segment_proba` (6 significant digits).

**Independent cross-checks:** Python port of the above (scipy.stats.binom.sf) = RSAT to ≤ 2e-12 relative on 68 runs
(t2/long sequences, k = 3/4/6, equi/input/`-markov 0..4`/`-noov`, 1str/2str, mseq on 4 sequences) and to %5g on 18
`-bgfile` runs, tables of order 0–2 (strand-sensitive and -insensitive); C# = port to ≤ 3e-13 on 400 random configurations.
Binomial tail vs mpmath exact sums ≤ 1e-13 (scipy ≤ 5e-14).

**Locked values (see TestSpec §9):** t2 = `ATGCATGCATGCAAATTTGGGCCCATGCTTAGCGGATCCATGCATGCTTTAAACGTACGTAGC`, k = 4,
`-1str -bg equi -lth occ 2`: ATGC occ 6, occ_P 1.4844942103072834e-07, occ_E 1.4844942103072835e-06 (10 tested),
occ_sig 5.8284214918614703; `-2str` input: ATGC|GCAT occ 9, occ_P 1.0760647330191343e-09, occ_sig 7.7920703428970466;
`-markov 1`: ATGC exp_freq 0.032915191520534237, occ_P 0.013954892377198637; `-noov` t3 k = 2: AA occ 5 at {1,3,5,27,29}.


**Audit group E (2026-10-01) — RSAT options and dyad-analysis.** Sources opened and run: rsa-tools/rsat-code 10043f2
`perl-scripts/oligo-analysis` 1.169 (`CountOligos`, `Degenerate`, `CalcSubWordFrequencies`, `CalcExpected` lexicon / `-pseudo`
blocks, `CalcOverlapCoefficient`, `CalcZscore`, `CalibrateSetFromSingleSequence`, `CalcProba`), `perl-scripts/dyad-analysis` 1.78,
`perl-scripts/calibrate-oligos`, `perl-scripts/fit-distribution` (calibration file layout), `lib/RSA.lib` (`ReadCalibration`,
`ReadExpectedFrequencies`), `lib/RSA.seq.lib` (`OverlapCoeff`, `IUPAC`), `lib/RSA.disco.lib` (`NbPossibleOligos`, `GroupRC`,
`CheckPatternThresholds`), `lib/RSAT/stats.pm` (`sum_of_poisson`, `sum_of_negbin2`, `LogToEng`). RSAT code defects found by
running the scripts (each bypassed with a documented rule, see the algorithm doc §5.3/§5.4): `-oneN`/`-onedeg` return nothing
(global `%IUPAC` never filled, scalars stored in `%patterns`); `sum_of_poisson` keeps `$prev_value` global between calls
(Poisson tails truncated: 0.003529988861723204 vs exact 0.0044559807752478468); negative-binomial terms rounded to 5 digits;
dyad `-2str` loses pairs seen only as their larger member (9 vs 11 tested on t2); `-type rep` double-counts palindromic
monads; with `-return zscore` dyads with variance ≤ 0 are dropped; with `-pseudo` an uncalibrated observed word crashes
RSAT; with `-2str` a residue absent from the input gets no Bernoulli probability (`CalcAlphabet` loops over observed
residues only), although its pooled value is positive. Cross-check counts in TestSpec §10.

**Audit round 2, groups G2 + G4 (2026-10-01) — `-seqtype dna|prot|other` and degenerate matching sequences.** Sources
opened and run: rsat-code 10043f2 `perl-scripts/oligo-analysis` 1.169 (`ReadArguments` `-seqtype` — prot/other set
`$sum_rc = 0`; the alphabet / `%accepted_residue` / `%accepted_oligo` block; `CountOligos` residue counts, pattern / residue
deletion, `$nb_possible_pos -= $discarded_occurrences`, the mseq block; `sub alphabet`; `CalcAlphabet`; `Degenerate`
incl. `$deg_mseq{$deg} += …{mseq}`; `CalcExpected` equiprobable / Markov / `-pseudo` / exp_ms; `CalcProba` ms_P),
`lib/RSA.disco.lib` `NbPossibleOligos` (`lc($seq_type) eq "dna"` guards the `-2str` / degenerate corrections),
`lib/RSA.seq.lib` `FoldSequence`, `OverlapCoeff`, `ReadNextSequence`, `perl-scripts/calibrate-oligos` (no sequence type),
`lib/RSA.lib` `ReadCalibration` (unconditional reverse-complement inference), `lib/RSAT/stats.pm` `binomial_boe`. RSAT code
defects found by running the unmodified script (oracle copies `oligo-analysis-g2` / `-g4` fix only these lines, env-guarded):
Markov `OverlapCoeff` fallback a,c,g,t = ¼ / other residues 0 for protein and text (LLL ovlp 1, AAA 1.3125; fixed: 1/|A|);
`OverlapCoeff` regex built from the unquoted word (`.` `?` … act as operators for `-seqtype other`; fixed: `\Q…\E`); `-markov`
dies with "Illegal division by zero" when a sub-word before a discarded residue was never counted (fixed: exp_freq stays 0);
`-onedeg` / `-oneN` with `-return mseq,proba` print no row ("oligomers tested for significance 0"), and the summed
`$deg_mseq` exceeds the number of sequences (AYG: 4 > 3 → "Successes (4) cannot be higher than trials (3)"; fixed: the union
of the words' matching-sequence sets); k = 1 `-2str` degenerate pairs with exp_freq > 1 die in `sum_of_binomials` (not
compared). Comparisons: `-seqtype` 450 runs / 15,066 patterns, every column ≤ 1.7e-13; degenerate mseq 300 runs / 68,997
patterns, mseq exact, ms columns ≤ 1.3e-13. Locked values in TestSpec §11.

---

## References

1. Compeau P, Pevzner P (2015). *Bioinformatics Algorithms: An Active Learning Approach*, 2nd ed., Chapter 2 ("Which DNA Patterns Play the Role of Molecular Clocks?"). Active Learning Publishers. Formula reproduced at https://github.com/wikiselev/bioinformatics-algorithms/wiki/Kmer-expected-number-of-occurrences-in-a-DNA-string (accessed 2026-06-14).
2. monaLisa (fmicompbio), `getKmerFreq` — observed vs expected k-mer frequencies and log2 enrichment. https://fmicompbio.github.io/monaLisa/reference/getKmerFreq.html (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-29**: RSAT oligo-analysis source + Bernoulli-background and long-k reference values (review 2026-09, B05).
- **2026-09-30**: RSAT oligo-analysis significance / Markov / `-2str` implemented; RSAT run + Python port + mpmath references (review 2026-09 follow-up, B05).
- **2026-10-01**: RSAT options (`-zscore`, `-pseudo`, `-oneN`/`-onedeg`, `-lexicon`, `-calibN`/`-calib1`) and `dyad-analysis` (review 2026-09 audit group E, B05 F31).
- **2026-10-01**: RSAT `-seqtype dna|prot|other` and degenerate matching-sequence statistics (review 2026-09 audit round 2, B05 F33).
- **2026-10-01**: Doc sync — key point 4 (significance) marked implemented (B05 F21, F31 options).
