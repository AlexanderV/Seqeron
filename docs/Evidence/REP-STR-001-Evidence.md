# Evidence Artifact: REP-STR-001

> **2026-09-30:** the approximate-detector (TRF) parts of this artifact are superseded by `docs/Evidence/REP-APPROX-001-Evidence.md` (TRF 4.10.0 compiled and used as the oracle; % matches / % indels are between adjacent copies, not vs the consensus).

**Test Unit ID:** REP-STR-001
**Algorithm:** Microsatellite / Short Tandem Repeat (STR) detection — perfect (default) and approximate / imperfect / interrupted (opt-in, Tandem Repeats Finder model)
**Date Collected:** 2026-06-24

---

## Online Sources

### Wikipedia: Microsatellite

**URL:** https://en.wikipedia.org/wiki/Microsatellite
**Accessed:** 2026-06-24 (prior validation session)
**Authority rank:** 4 (encyclopedia citing primaries)

**Key Extracted Points:**

1. **Motif length & copies:** repeat unit "typically ten nucleotides or less" (library uses 1–6 bp), "repeated 5–50 times".
2. **Worked examples (verbatim):** "TATATATATA is a dinucleotide microsatellite" and "GTCGTCGTCGTCGTC is a trinucleotide microsatellite".
3. **History:** first microsatellite was "a polymorphic GGAT repeat in the human myoglobin gene" (Weller, Jeffreys et al. 1984).

These cover the **perfect** STR detector (unchanged default) and are recorded in full in the per-unit
validation report `docs/Validation/reports/REP-STR-001.md`.

### Tandem Repeats Finder — Benson (1999), Nucleic Acids Research 27(2):573–580

**URL:** https://academic.oup.com/nar/article/27/2/573/1061099 (DOI https://doi.org/10.1093/nar/27.2.573)
**Accessed:** 2026-06-24 (fetched via WebFetch in this session)
**Authority rank:** 1 (peer-reviewed paper — the canonical TRF model)

**Key Extracted Points (verbatim where quoted):**

1. **Approximate-repeat definition:** "A tandem repeat in DNA is two or more contiguous, *approximate* copies of a pattern of nucleotides." The model treats these as the "alignment of two tandem copies of a pattern of length *n* by a sequence of *n*-independent Bernoulli trials."
2. **Reported statistics (output, verbatim):** "indices of the repeat in the sequence; (ii) period size; (iii) number of copies aligned with the consensus pattern; (iv) size of the consensus pattern (may differ from the period size); (v) percent of matches between adjacent copies overall; (vi) percent of indels between adjacent copies overall; (vii) alignment score".
3. **Alignment scoring (verbatim):** "We used one of two sets of alignment parameters (match, mismatch, gap), either (+2,−7,−7) or (+2,−5,−7). Only those repeats scoring at least 50 with these parameters are reported." The alignment is "Smith-Waterman style local alignment using wraparound dynamic programming."
4. **How a repeat is scored (verbatim):** "a candidate pattern … is selected from the nucleotide sequence and aligned with the surrounding sequence using wraparound dynamic programming (WDP)."
5. **Consensus pattern (verbatim):** "An initial candidate pattern P is drawn from the sequence, but this is usually not the best pattern to align with the tandem repeat. To improve the alignment, we determine a consensus pattern by majority rule from the alignment."

### Tandem Repeats Finder — definitions page

**URL:** https://tandem.bu.edu/trf/trf.definitions.html
**Accessed:** 2026-06-24 (fetched via WebFetch in this session)
**Authority rank:** 2 (official tool documentation)

**Key Extracted Points (verbatim):**

1. **Statistic definitions:** "Period Size: Period size of the repeat"; "Copy Number: Number of copies aligned with the consensus pattern"; "Consensus Size: Size of consensus pattern (may differ slightly from the period size)"; "Percent Matches: Percent of matches between adjacent copies overall"; "Percent Indels: Percent of indels between adjacent copies overall"; "Alignment score".
2. **Scoring (verbatim):** match weight is "+2 in all options here. Mismatch and indel weights (interpreted as negative numbers) are either 3, 5, or 7."

### Tandem Repeats Finder — detailed description / Bernoulli model (trf.desc.html, trf.definitions.html)

**URL:** https://tandem.bu.edu/trf/trf.desc.html and https://tandem.bu.edu/trf/trf.definitions.html
**Accessed:** 2026-06-24 (fetched via WebFetch in this session)
**Authority rank:** 2 (official tool documentation accompanying the Benson 1999 paper)

**Key Extracted Points (verbatim):**

1. **Bernoulli model (verbatim):** "We model alignment of two tandem copies of a pattern of length n by a sequence of n independent Bernoulli trials (coin-tosses)."
2. **Matching probability PM (verbatim):** "The probability of success, P(Heads), which we also call PM or matching probability, represents the average percent identity between the copies."
3. **Indel probability PI (verbatim):** "A second probability, PI or indel probability, specifies the average percentage of insertions and deletions between the copies."
4. **Statistics are between ADJACENT copies (verbatim):** the statistics refer to "the matches, mismatches and indels overall between adjacent copies in the sequence, not between the sequence and the consensus pattern."
5. **Sum-of-heads random variable (verbatim):** "Let the random variable R(d,k,pM) = the total number of heads in head runs of length k or longer in an iid Bernoulli sequence of length d with success probability pM." "The distribution of R(d,k,pM) is well approximated by the normal distribution and its exact mean and variance can be calculated in constant time." The criterion uses "the largest number, x, such that 95% of the time R(d,k,pM) ≥ x".
6. **Random walk variable (verbatim):** "Let the random variable W(d,pI) = the maximum displacement from the origin of a one dimensional random walk with expected number of steps equal to pI·d."
7. **Default probabilities (verbatim):** "PM = .80 and PI = .10 by default"; "The best performance can be achieved with values of PM=80 and PI=10."

> Note on reproducibility: points 1–4 and 7 are the *reported probabilistic measures* and are faithfully
> reproducible (PM/PI estimated between adjacent copies). Points 5–6 (R(d,k,pM) and W(d,pI)) are the
> *k-tuple SEEDING* machinery; the actual 95% percentile cut-offs are obtained by TRF from
> non-redistributable simulation tables and the closed-form mean/variance of R(d,k,pM) is not stated
> verbatim on these pages — so the percentile-based seeding is the honest, non-reproducible residual.

### Tandem Repeats Finder — reference implementation (Benson-Genomics-Lab/TRF)

**URL:** https://github.com/Benson-Genomics-Lab/TRF
**Accessed:** 2026-06-24 (fetched via WebFetch in this session)
**Authority rank:** 3 (reference implementation / official usage)

**Key Extracted Points (verbatim):**

1. **Usage:** "Please use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]".
2. **Recommended parameters (verbatim):** "The recomended values for Match Mismatch and Delta are 2, 7, and 7 respectively." Example command "trf yoursequence.txt 2 7 7 80 10 50 500".
3. **Minscore meaning (verbatim):** "if we set the matching weight to 2 and the minimun score to 50, assuming perfect alignment, we will need to align at least 25 characters to meet the minimum score (for example 5 copies with a period of size 5)."

---

### MISA — Thiel et al. (2003) Theor Appl Genet 106:411-422, `misa.pl` v1.0 source

**URL:** https://raw.githubusercontent.com/cfljam/SSR_marker_design/master/misa.pl (mirror of the IPK script; the IPK host is blocked)
**Accessed:** 2026-09-29 (downloaded and read in this session)

1. Per motif size: `my $search = "(([acgt]{$motiflen})\\2{$minreps,})"; while ( $seq =~ /$search/ig )` — leftmost, greedy, non-overlapping matches; only `a/c/g/t` motifs.
2. `#reject false type motifs [e.g. (TT)6 or (ACAC)5]` — non-primitive motifs are dropped.
3. `$repeats{$nr} = length($ssr) / $motiflen` over the matched complete copies.

### pytrf 1.5.0 (Krait engine) — Du et al. (2018) Bioinformatics 34(4):681-683

**URL:** PyPI `pytrf==1.5.0` sdist, `src/str.c` (installed and read in this session)

1. Run found by `while ((i < b) && (self->seq[i] == self->seq[i+j])) ++i; rl = i + j - cs;` — the period-j run from the seed `cs`.
2. `ssr->repeat = rl/j; ssr->length = ssr->repeat * j;` — complete copies only; `next_start = end`.
3. `if (self->seq[i] == 78) continue;` — `N` is skipped. No primitivity check (smaller sizes are tried first).

### Kolpakov R, Kucherov G (1999) "Finding maximal repetitions in a word in linear time", FOCS

Maximal repetition (run): a periodic factor `S[a..e)` with minimal period p, exponent ≥ 2, not extendable left or right with the same period. (Definition as used here; paper not opened.)

### Numerical cross-check (review 2026-09)

Harness: C# `RepeatFinder.FindMicrosatellites(string, p, p, k)` vs (a) Python brute-force maximal-repetition reference, (b) MISA regex rule, (c) pytrf `STRFinder` with only size p enabled. 3,132 cases (11 crafted × p=1..6 × k∈{2,3} + 3,000 random, seed 20260929, alphabets incl. `N`, planted repeats).

| Reference | Agreement | All disagreements explained by |
|-----------|-----------|--------------------------------|
| Brute-force maximal runs | 3,132 / 3,132 | — |
| MISA rule | 3,099 / 3,132 | 17 same-period runs overlapping by < p (MISA restarts after the first run's copies); 16 non-primitive region consumed by the regex first |
| pytrf single size | 2,752 / 3,132 | 17 same overlap case; 363 non-primitive/`N` units (pytrf has no primitivity check) |

Before the fix the code additionally reported every rotation of a run that reached past the run-start's complete copies (e.g. `ATATATA` → `AT×3@0` **and** `TA×3@1`; `AAACACACACACAAA` → `AC×5@2` **and** `CA×5@3`) — one locus reported up to p times, disagreeing with all three references — and reported runs of `N` as mononucleotide microsatellites.

### MISA per-unit-size thresholds, compound SSRs and repeat-type classes (review 2026-09-30, audit WP3)

**Source:** the same `misa.pl` v1.0 (read in full and executed with `perl` 5; an instrumented copy that additionally prints
the per-sequence SSR list in `@order` and the rejected non-primitive matches was used only as a harness).

1. `misa.ini` header example: `definition(unit_size,min_repeats): 1-10 2-6 3-5 4-5 5-5 6-5`, `interruptions(max_difference_for_2_SSRs): 100`.
2. Compound loop: `@order = sort { $start{$a} <=> $start{$b} } keys %start`; `$space = $amb + 1`; two SSRs join when
   `$start{next} - $end{prev} <= $space` (1-based inclusive coordinates ⇒ ≤ `amb` bases in between); `< 1` ⇒ overlap,
   type `c*`, notation `($motif)$repeats*`; otherwise `$interssr = lc substr($seq, $end, $start − $end − 1)` and type `c`;
   the inner `while` compares with `$end{$order[$i]}` of the previous SSR and sets `$end = $end{$order[$i+1]}`.
3. `.statistics` "Frequency of classified repeat types (considering sequence complementary)": for each motif, the
   smallest rotation of the motif and of its reverse complement (`tr/ACGT/TGCA/`, `reverse`), joined `A/B` with the
   smaller first; counts summed over the group.
4. Krait `src/motif.py` (`StandardMotif`, raw.githubusercontent lmdu/krait): `similar_motif` = rotations,
   `reverse_complete_motif`, `complete_motif`, `reverse_motif`; `motif_sorted` by `motif_to_number` (A=1, T=2, C=3, G=4);
   levels 0–4; `src/widgets.py` default `ssr/level` = 3. `_motifs` is a class attribute, so its cache is shared by all
   levels until `setLevel` — the reference runs used a fresh cache per call.

**Numerical cross-check (C# harness vs `perl misa.pl`, PERL_HASH_SEED=0):** 6 048 sequences (6 × 8 crafted + 6 000 random
SSR-rich, with interruptions of 0–130 bp, N runs, lowercase) in six configurations — `1-10 2-6 3-5 4-5 5-5 6-5` with
interruptions 100 and 0; `1-5 2-3 3-3 4-3 5-3 6-3` / 20; `1-3 2-2 3-2 4-2 5-2 6-2` / 5; `2-4 3-3 5-2` / 50;
defaults + `7-4 8-3 10-3` / 100 — 72 974 misa.pl SSRs, 12 330 compounds.

| Check | Result |
|-------|--------|
| `FindMicrosatellites(seq, map)` vs brute-force maximal primitive runs with per-size thresholds | 0 mismatching sequences |
| `AssembleCompoundMicrosatellites` fed misa.pl's SSR list in misa.pl order vs misa.pl `.misa` (type, notation, size, start, end) | 0 mismatching sequences |
| End to end vs misa.pl, sequences with identical SSR lists | `.misa` rows identical in all of them |
| SSR-list differences (971 sequences) | 1 220 SSRs: same-size run overlapping the previous match by < p (misa.pl truncates); 1 365 SSRs: primitive run inside a rejected non-primitive match (misa.pl consumed it); 0 unexplained |
| Same SSR set, equal-start SSRs ordered differently (53 sequences, configurations with minimum copies 2–3) | Perl hash order in misa.pl; stable input order here |
| `GetCanonicalMotifFrequencies` on misa.pl's SSR list vs misa.pl classified table | identical in all 6 configurations |
| `GetCanonicalMotifClass` vs misa.pl `.statistics` (one run per motif, motif × 12) | 5 356 / 5 356 primitive motifs of 1–6 bp |
| `GetStandardMotif(m, level)` vs Krait `StandardMotif(level).standard(m)` | 5 460 motifs × 5 levels, 0 mismatches |

**Progress contract (`IProgress<double>` on the cancellable overloads):** values `(k·n + i)/(K·n)` for the k-th of K
unit lengths every 1000 visited run starts — non-decreasing, in [0, 1) — then exactly 1.0; the token is checked at
the same points and before the final report. (The 2026-06 TestSpec entry "progress reporting not implemented" was wrong.)

### misa.pl-parity scan (review 2026-10-01, audit WP8)

**Source:** `misa.pl` v1.0 (raw GitHub `cfljam/SSR_marker_design`, re-downloaded; perl 5.38.2). Scan loop (lines
101–125): `for` unit sizes in `sort { $a <=> $b } keys %typrep`; `$search = "(([acgt]{$motiflen})\\2{$minreps,})"`;
`while ($seq =~ /$search/ig)`; redundancy test `([ACGT]{$j})\\1{($motiflen/$j-1)}` for `$j = $motiflen-1 … 1` —
for j ∤ p the count is fractional, perl warns "Unescaped left brace in regex is passed through" and the brace is a
literal, so only divisors can match; `next if $redundant` after `pos()` advanced; `$end = pos($seq)`,
`$start = $end - length($ssr) + 1`. Order (line 130): `sort { $start{$a} <=> $start{$b} } keys %start`.

**Hash order (task check):** misa.pl run with `PERL_HASH_SEED` = 1, 2, 3 on the 6 000-sequence set: outputs differ in
the settings `1-3 2-2 3-2 4-2 5-2 6-2`/10 (933 `.misa` rows between seeds 1 and 2), `1-5 2-3 3-3 4-3 5-3 6-3`/0 (115),
`1-2 3-2 5-2`/5 (1 667); identical for `1-10 2-6 3-5 4-5 5-5 6-5`/100, `1-12 2-4 3-4 4-3 5-3 6-3`/50, `2-3 4-2 6-2`/20
(no start ties). Same seed twice: identical; two unseeded runs: differ. Example `TTTGTTTGTTTGTTTGTTTGTTT`
(`1-3 2-2 …`/10): seed 1 `c* (TTTG)5(T)3*g(T)3g(T)3g(T)3g(T)3g(T)3 23 1 23`, seed 2
`c* (T)3(TTTG)5*(T)3*g(T)3g(T)3g(T)3g(T)3 23 1 23`. Ties need two primitive runs of periods p ≠ q starting at one
position; with run lengths ≥ p + q − gcd(p, q) they would share period gcd (Fine–Wilf), which the MISA default
thresholds always guarantee — so default-setting output never depends on the hash order.

**Cross-check (`scratchpad/wp8/misa`, C# harness `xcm`):** 6 000 SSR-rich sequences (20–600 bp; runs of 1–6-bp units
incl. non-primitive ones, partial copies, AC/AT/A-rich spacers, N runs, 20 % mixed case) × 6 settings above; misa.pl
instrumented to dump its SSR list, and a copy with the tie broken by SSR number (`|| $a <=> $b`).

| Setting | misa.pl SSRs | `.misa` rows | `MisaRegex` vs misa.pl (tie by SSR nr): SSR lists / rows | vs stock misa.pl seed 1 / seed 2 (sequences) | all explained by start ties with equal SSR set | `MaximalRuns` SSR lists differing |
|---|---|---|---|---|---|---|
| 1-10 2-6 3-5 4-5 5-5 6-5 / 100 | 37 299 | 5 944 | 0 / 0 | 0 / 0 | — | 534 |
| 1-3 2-2 3-2 4-2 5-2 6-2 / 10 | 174 866 | 8 772 | 0 / 0 | 908 / 900 | yes | 3 048 |
| 1-5 2-3 3-3 4-3 5-3 6-3 / 0 | 58 909 | 41 103 | 0 / 0 | 14 / 10 | yes | 1 037 |
| 1-12 2-4 3-4 4-3 5-3 6-3 / 50 | 41 252 | 6 965 | 0 / 0 | 0 / 0 | — | 868 |
| 2-3 4-2 6-2 / 20 | 25 885 | 14 938 | 0 / 0 | 0 / 0 | — | 1 038 |
| 1-2 3-2 5-2 / 5 | 324 737 | 21 578 | 0 / 0 | 1 501 / 1 556 | yes | 1 425 |

## Documented Corner Cases and Failure Modes

### From Benson (1999)

1. **Imperfect copies are intrinsic:** copies are *approximate* — substitutions and indels within the repeat are expected, not errors; the percent-matches / percent-indels statistics quantify them.
2. **Consensus ≠ period:** the consensus pattern size "may differ from the period size" (indels can shift the consensus length). In the implemented subset the consensus is exactly `period` bases by construction (majority rule over period-aligned columns), so `ConsensusSize == Period`.
3. **Minimum two copies:** a tandem repeat requires "two or more contiguous … copies".

### From the perfect detector (prior validation)

1. **Interrupted repeat fragmentation:** the perfect `FindMicrosatellites` detector breaks an interrupted tract into separate short perfect runs (or misses copies below `minRepeats`). This is the precise gap the approximate detector closes.

---

## Test Datasets

### Dataset: Perfect dinucleotide (control)

**Source:** derived; Benson (1999) statistics definitions.

| Parameter | Value |
|-----------|-------|
| Sequence | `CACACACACA` (10 bp) |
| Period | 2, consensus `CA`, ref `CA`×5 |
| Alignment | 10 match columns, 0 mismatch, 0 indel |
| Alignment score | 10 × 2 = **20** |
| Percent matches | 10/10 = **100 %** |
| Percent indels | 0/10 = **0 %** |
| Copy number | 10/2 = **5** |

### Dataset: Interrupted trinucleotide (one substitution)

**Source:** derived from a perfect `(CAG)×6` tract with one mid-tract substitution.

| Parameter | Value |
|-----------|-------|
| Sequence | `CAGCAGCAGTAGCAGCAG` (18 bp) — copy 4 `C`→`T` at index 9 |
| Period | 3, consensus `CAG`, ref `CAG`×6 (18) |
| Alignment | 17 match, 1 mismatch, 0 indel (18 columns) |
| Alignment score | 17 × 2 + 1 × (−7) = **27** |
| Percent matches | 17/18 = **94.4̄ %** |
| Percent indels | 0/18 = **0 %** |
| Copy number | 18/3 = **6** |
| Perfect detector | reports only `CAG`×3 at pos 0 (breaks at the `T`) |

### Dataset: Interrupted dinucleotide (one substitution)

**Source:** derived from `(CA)×6` with one mid-tract substitution.

| Parameter | Value |
|-----------|-------|
| Sequence | `CACACATACACA` (12 bp) — index 6 is `T` |
| Period | 2, consensus `CA`, ref `CA`×6 (12) |
| Alignment | 11 match, 1 mismatch, 0 indel (12 columns) |
| Alignment score | 11 × 2 + 1 × (−7) = **15** |
| Percent matches | 11/12 = **91.6̄ %** |
| Percent indels | 0/12 = **0 %** |
| Copy number | 12/2 = **6** |
| Perfect detector | reports only `CA`×3 at pos 0 |

### Dataset: Single-base deletion (indel)

**Source:** derived from a perfect `(CAG)×10` tract (30 bp) with one base deleted at index 15.

| Parameter | Value |
|-----------|-------|
| Sequence | `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (29 bp) |
| Period | 3, consensus `CAG`, ref `CAG`×10 (30) |
| Alignment | 29 match, 0 mismatch, 1 indel (30 columns) |
| Alignment score | 29 × 2 + 1 × (−7) = **51** |
| Percent matches | 29/30 = **96.6̄ %** |
| Percent indels | 1/30 = **3.3̄ %** |
| Copy number | 29/3 = **9.6̄** |
| Minscore | 51 ≥ 50 → reported at the default `minScore = 50` |
| Perfect detector | fragments into `CAG`×5 + several 4-copy frame rotations |

### Dataset: Bernoulli PM/PI estimates between adjacent copies (TRF statistical model)

**Source:** Benson (1999) Bernoulli model — PM = average percent identity between adjacent copies,
PI = average percentage of indels between adjacent copies; statistics "between adjacent copies … not
between the sequence and the consensus pattern". Each Bernoulli trial = one alignment column between two
adjacent copies (heads = match). Expected matches = Bernoulli mean E[heads] = PM·d over d trials.

| Sequence | Period | Adjacent copy pairs (cols) | Matches / Mismatches / Indels | PM | PI | E[matches] = PM·d |
|----------|--------|----------------------------|-------------------------------|----|----|-------------------|
| `CACACACACA` (perfect CA×5) | 2 | 4 pairs (8 cols) | 8 / 0 / 0 | 8/8 = **1.00** | **0** | 1.00·8 = **8** |
| `CAGCAGCAGTAGCAGCAG` (CAG×6, copy 4 = TAG) | 3 | 5 pairs (15 cols) | 13 / 2 / 0 | 13/15 ≈ **0.8667** | **0** | (13/15)·15 = **13** |
| `CACACATACACA` (CA×6, index 6 = T) | 2 | 5 pairs (10 cols) | 8 / 2 / 0 | 8/10 = **0.80** | **0** | 0.80·10 = **8** |
| `ACACTGTG` (two unrelated period-4 copies) | 4 | 1 pair (4 cols) | 0 / 4 / 0 | 0/4 = **0.00** | **0** | **0** |

The single-substitution tracts show the adjacent-copy PM (`13/15`, `0.80`) is distinct from the
consensus-based percent matches (`17/18`, `11/12`) computed by `FindApproximateTandemRepeats`, which is
exactly the distinction Benson draws ("between adjacent copies … not between the sequence and the
consensus pattern"). A tract is "significant" under the model when its estimated PM is at least Benson's
default PM = 0.80; the `CACACATACACA` tract sits exactly on that threshold (PM = 0.80, inclusive).

---

## Assumptions

1. **ASSUMPTION: Deterministic exhaustive period scan in place of TRF probabilistic k-tuple seeding.**
   The implemented subset enumerates every (start, period) window and scores it by alignment, instead of
   reproducing TRF's probabilistic k-tuple distance-list seeding and the sum-of-Bernoulli statistical
   significance test. This is the honest residual recorded in `LIMITATIONS.md` and §5.3 of the algorithm
   doc. It does not affect the *statistics* of a reported repeat (those follow Benson 1999 exactly); it
   only changes *which candidate windows are examined* (the subset examines all of them deterministically).
2. **ASSUMPTION: percent-matches / percent-indels denominator = total alignment columns.** Benson (1999)
   defines these as "percent of matches/indels between adjacent copies overall" without a verbatim formula;
   the implementation expresses each as a fraction of the total alignment-column count, which is internally
   consistent and reproduces the worked numbers above. (Match weight, mismatch/indel penalties, Minscore,
   and the statistic *names* are all source-verbatim; only the percentage denominator convention is fixed
   here.)
3. **ASSUMPTION: adjacent-copy segmentation for the Bernoulli PM/PI estimate.** `ComputeBernoulliStatistics`
   estimates PM/PI "between adjacent copies" (Benson 1999, verbatim) by segmenting the tract into
   period-length copies and aligning each adjacent pair with the recommended TRF scoring. For pure-
   substitution and perfect tracts this segmentation is unambiguous and the PM/PI values are exact (see
   the dataset above). For tracts containing indels the per-pair alignment frame is alignment-dependent,
   so only the *qualitative* Bernoulli property (PI > 0; PM + mismatch-fraction + PI = 1) is asserted for
   those cases, not a fragile exact PI. The Bernoulli mean E[heads] = PM·d is the only moment stated
   exactly by the source ("average percent identity"); the closed-form mean/variance of the sum-of-heads
   R(d,k,pM) and its 95% percentile are NOT reproduced (non-redistributable simulation tables — the
   genome-scale seeding residual).

---

## Recommendations for Test Coverage

1. **MUST Test:** perfect dinucleotide → approximate detector reports 100 % matches, 0 % indels, period 2, copy number 5 — Evidence: Benson (1999) statistics; perfect-alignment control.
2. **MUST Test:** interrupted (one substitution) tract → approximate detector reports it as ONE repeat with the exact percent-matches, while the perfect detector fragments it — Evidence: Benson (1999) approximate definition.
3. **MUST Test:** single-base-deletion tract → reported with the exact percent-indels and copy number — Evidence: Benson (1999) percent-indels statistic.
4. **MUST Test:** Minscore gate — a tract below the threshold is not reported; default 50 — Evidence: Benson (1999) "Only those repeats scoring at least 50 … are reported".
5. **SHOULD Test:** scoring constants are the recommended TRF set (+2, −7, −7) — Rationale: traceability of the score to the source.
6. **COULD Test:** determinism / null / empty / invalid-parameter guards — Rationale: API contract.
7. **MUST Test (Bernoulli):** perfect tract → PM = 1.0, PI = 0, E[matches] = PM·d — Evidence: Benson (1999) PM = average percent identity.
8. **MUST Test (Bernoulli):** one-substitution tract → exact adjacent-copy PM (13/15; 8/10), distinct from the consensus percent matches — Evidence: "between adjacent copies … not … the consensus pattern".
9. **MUST Test (Bernoulli):** indel tract → PI > 0 and match + mismatch + indel fractions partition the trials — Evidence: Benson (1999) PI = average percentage of indels.
10. **SHOULD Test (Bernoulli):** PM threshold (default 0.80, custom) and exposed defaults PM = .80 / PI = .10 — Evidence: TRF desc "PM = .80 and PI = .10 by default".

---

## References

1. Benson G (1999). Tandem repeats finder: a program to analyze DNA sequences. Nucleic Acids Research 27(2):573–580. https://doi.org/10.1093/nar/27.2.573
2. Tandem Repeats Finder — definitions. https://tandem.bu.edu/trf/trf.definitions.html (accessed 2026-06-24)
3. Tandem Repeats Finder — reference implementation. https://github.com/Benson-Genomics-Lab/TRF (accessed 2026-06-24)
4. Wikipedia: Microsatellite. https://en.wikipedia.org/wiki/Microsatellite (accessed 2026-06-24)

---

## Change History

- **2026-06-24**: Initial documentation — Benson (1999) TRF approximate-repeat model added to support the opt-in `FindApproximateTandemRepeats` detector (REP-STR-001 limitation fix). Perfect-repeat detector evidence (Wikipedia / MISA) carried from the prior validation.
- **2026-06-24**: Added the TRF Bernoulli statistical-significance model (Benson 1999) — verbatim PM/PI/Bernoulli-trial definitions from the TRF desc/definitions pages, the adjacent-copy PM/PI dataset, and the supporting assumption — for the new opt-in `ComputeBernoulliStatistics`. The R(d,k,pM)/W(d,pI) k-tuple seeding remains the documented genome-scale-performance residual.
- **2026-09-30** (review 2026-09, B04 audit WP3): MISA per-unit-size thresholds, compound SSRs (types c / c*), MISA
  repeat-type classes and Krait standard motifs — sources, misa.pl / Krait cross-checks; progress-reporting contract.
- **2026-10-01** (B04 audit WP8): misa.pl-parity scan (`MicrosatelliteScanMode.MisaRegex`) — scan-loop source, Perl hash-order
  check (PERL_HASH_SEED), 6 000 × 6 cross-check (0 differences with deterministic tie order).
