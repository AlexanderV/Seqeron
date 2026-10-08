# Evidence Artifact: REP-APPROX-001

**Test Unit ID:** REP-APPROX-001
**Algorithm:** Approximate (imperfect / interrupted) tandem-repeat detection with the Tandem Repeats Finder model (`RepeatFinder.FindApproximateTandemRepeats`) and the TRF Bernoulli statistics (`RepeatFinder.ComputeBernoulliStatistics`)
**Date Collected:** 2026-09-30 (campaign 2026-09, batch B04); revised 2026-10-01 (WP6, WP7)

---

## Sources opened

### Benson G (1999) Nucleic Acids Res 27(2):573–580 — via the TRF README

The paper PDF (OUP / PMC) is blocked by the proxy; its content is reproduced in the TRF 4.10.0 README sections
"How does Tandem Repeats Finder work?" (Probabilistic Model, Detection Component, Sum of Heads, Random Walk,
Apparent Size and Waiting Time distributions, Analysis Component, Consensus Pattern and Period Size, Redundancy).
Points used:

1. "A tandem repeat in DNA is two or more adjacent, approximate copies of a pattern of nucleotides."
2. Parameters `Match Mismatch Delta PM PI Minscore MaxPeriod`; "The recomended values for Match Mismatch and Delta
   are 2, 7, and 7"; "Probabilistic data is available for PM values of 80 and 75 and PI values of 10 and 20";
   "The alignment of a tandem repeat must meet or exceed this alignment score [Minscore] to be reported".
3. Sum of heads: "R(d,k,PM) = the total number of heads in head runs of length k or longer in an iid Bernoulli
   sequence of length d … well approximated by the normal distribution and its exact mean and variance can be
   calculated … the largest number, x, such that 95% of the time R ≥ x"; "We require that the smallest pattern
   for tuple size k have a sum-of-heads criterion of at least k+1".
4. Random walk: `Δd_max = floor(2.3·√(PI·d))`; apparent size estimated "by simulation"; waiting time recursion.
5. Analysis: candidate pattern `j+1..i` "aligned with the surrounding sequence using wraparound dynamic
   programming (WDP). If at least two copies of the pattern are aligned with the sequence, the tandem repeat is
   reported"; narrow band for patterns > 20; consensus "by majority rule", realignment; "Period size is defined
   as the most common matching distance between corresponding characters in the alignment and may not be
   identical to consensus size"; multiples limited "to, at most, three pattern sizes".
6. Table explanation: indices, period, "Number of copies aligned with the consensus pattern", consensus size,
   "Percent of matches between adjacent copies overall", "Percent of indels between adjacent copies overall",
   score, composition, "Entropy measure based on percent composition". Alignment explanation: statistics refer
   to "the matches, mismatches and indels overall between adjacent copies in the sequence, not between the
   sequence and the consensus pattern".
7. Redundancy: "restricted to the three best scoring period sizes … the same period size may be detected more
   than once with different scores and slightly different indices."
8. Version 4.10.0: MaxPeriod above 2000 is an error ("very large TRs are outside the scope of the TRF statistical
   models"); expected tables for `test_seqs` s1–s3.

### TRF 4.10.0 source — github.com/Benson-Genomics-Lab/TRF (git clone, commit 355c1f9, 2020-06-29; AGPL-3.0)

Compiled with `./configure && make` and run as `trf seq.fa 2 7 7 80 10 <minscore> <maxperiod> -h -d`. An
instrumented copy (scratch only: two `fprintf(stderr, …)` lines printing the trigger `(i, d, consensus size)` and
the raw `match / mismatch / indel`, score and copy number in `get_statistics`) provided the per-candidate
reference. Facts read from the code: `init_sm` (only identical A/C/G/T score Match — "changed to use Similarity
Matrix to avoid N matching itself"); `newwrap` (two-pass WDP, backward scan to the leftmost best row, forward
redo from one pattern length earlier, zero cells killed to −1000 beyond the candidate); traceback preference and
copy-number accumulation (`get_pair_alignment_with_copynumber`); non-weighted `get_consensus`
(`WEIGHTCONSENSUS 0`); `get_statistics` (adjacent-copy cursors, most-common distance, composition, entropy, integer
truncation of percentages); `newtupbo` (tuple sizes 4/5/7 for d ≤ 29 / 30–159 / ≥ 160, `Min_Distance_Window = 20`,
copy-number rule 1.9 / ramp / 1.8, distance-seen array, `SMALLDISTANCE = 20` full vs narrow band);
`GetTopPeriods` / `multiples_criteria_4` (three best periods, period-1 composition ≥ 80 %); `trfclean.h`
(`RemoveBySize`, index sort, `RemoveRedundancy` 90 % overlap + `IsRedundant`); `sumdata80` / `waitdata80` tables.

Licence note: TRF is AGPL-3.0 and this library is MIT, so TRF code (and its simulated tables) is not ported
line-by-line; the C# implementation is written from the published description, and the sum-of-heads cut-offs are
recomputed independently (below).

---

## Stage A findings (description)

| # | Old description (TestSpec / XML / algorithm doc) | Source says | Verdict |
|---|---|---|---|
| A1 | % matches / % indels = columns of the alignment of the window **against the consensus** | between **adjacent copies** (README Alignment Explanation item 6) | wrong → corrected |
| A2 | Copy number = observed bases / period | copies aligned with the consensus = aligned consensus columns / consensus size (TRF: 10.0 for the 29-bp deletion tract, old 9.667) | wrong → corrected |
| A3 | Consensus size = period; period = candidate size | period = most common matching distance, consensus size may differ | wrong → corrected |
| A4 | Global alignment of a window against tiled consensus | local wraparound DP | wrong → corrected |
| A5 | N treated as an ordinary symbol (all-N = perfect repeat) | N never matches | wrong → corrected |
| A6 | Overlap handling: drop anything contained in a higher-scoring repeat | TRF reports overlapping repeats of different periods; only redundant multiples / same-period duplicates (≥ 90 % overlap) are removed | wrong → corrected |
| A7 | "k-tuple seeding tables non-redistributable" | TRF is open source (AGPL-3.0); the sum-of-heads table is reproducible from the published definition | corrected |
| A8 | Scoring 2/7/7, Minscore 50, PM .80 / PI .10, majority consensus | as source | confirmed |

---

## Reference cross-checks (all numbers from this session)

1. **Sum-of-heads criterion** — exact mean/variance of `R(d,k,0.8)` by a run-length Markov chain, cut-off
   `max(k+1, ⌊μ − 1.65σ⌋)`: equal to TRF's `sumdata80` for **2000 / 2000** distances (floor(μ − 1.645σ) differs at
   412, the exact 5 % quantile at 128 — TRF uses z = 1.65). The README's illustrative example (PM .75, k 5, d 100 →
   26) is not reproduced by either the normal approximation (μ = 46.27, σ = 11.60 → 27) or the exact quantile (27);
   the PM = 80 table actually used by TRF is reproduced exactly.
2. **Analysis component** — for each candidate `(i, d)` that TRF analysed (instrumented binary), the C# private
   analysis at the same `(i, d)`: **1 524 / 1 524 identical** for patterns ≤ 20 (indices, score, consensus size,
   copy number, adjacent-copy match/mismatch/indel counts), 827 / 959 for patterns > 20 (TRF narrow band;
   600 random sequences with embedded imperfect repeats, seeds 31 and 32). A Python prototype of the same
   algorithm agreed with the binary on 1 460 / 1 460 small-pattern candidates before the C# port, and the C# port
   equals the prototype on 700 / 700 whole-sequence runs.
3. **Whole pipeline vs TRF `.dat`** (700 random sequences with embedded imperfect repeats, 1–4 arrays each,
   substitution 0–15 %, indel 0–10 %, N in every 7th sequence; `2 7 7 80 10 50 500`):

   | Embedded periods | TRF rows | C# rows | exact rows (indices, period, size, score, consensus) | region level (period ±1, ≥ 90 % mutual overlap) | C# rows confirmed by TRF |
   |---|---|---|---|---|---|
   | ≤ 20 (500 seqs) | 1 154 | 1 152 | 1 069 (92.6 %) | 1 108 (96.0 %) | 1 109 (96.3 %) |
   | ≤ 100 (200 seqs) | 348 | 333 | 280 (80.5 %) | 324 (93.1 %) | 324 (97.3 %) |

   On the 1 349 exact rows the other fields (copies to TRF's 0.1, % matches / % indels / composition to TRF's
   integers, entropy to 0.01 for N-free regions) agree on 1 347 (two equal-score traceback ties reached from a
   different trigger position). Entropy differs on 29 of 77 regions containing N (TRF divides the four counts by
   the length including N; the library uses the canonical normalised Shannon entropy).
4. **Crafted cases** (all identical to TRF): TRF README test_seqs s1/s2/s3; `CACACACACA`; `CAGCAGCAGTAGCAGCAG`
   (86.67 %); `CACACATACACA` (nothing: sum of heads 4 < 5); the 29-bp deletion tract (three overlapping rows at
   minscore 10: periods 3, 14, 11; two at 50); flanked + lowercase version; `NNNNNNNN` (nothing); N inside a CAG
   array (a mismatch); homopolymer with one G (period 1, 96.43 %); `(CA)×30` (only period 2); the 18-cell
   combinatorial grid of `RepeatsCombinatorialTests` (all 18 equal TRF, including the one cell TRF does not report).
5. **Bernoulli statistics** — `ComputeBernoulliStatistics` on `CACACACACA`, `CAGCAGCAGTAGCAGCAG`, the deletion tract:
   8/0/0, 13/2/0, 25/0/2 = TRF's statistics lines (old implementation gave frame-shifted counts for the deletion).
6. **Timing** (Release, this container): 100 kb random, maxPeriod 6 / 500 / 2000 → 0.06 / 0.34 / 0.59 s; 100 kb of
   mixed embedded repeats, maxPeriod 500 → 1.3 s, 545 rows (TRF 549 rows, 0.22 s); pathological 100 kb perfect
   `(CA)n`, maxPeriod 500 → 54 s (compiled TRF 74 s; every even distance re-aligns the whole array in both).
   The old exhaustive window scan was O(n²·P·L²) (each (start, period, window) a global alignment).

---

## Convention adopted

- Analysis exactly as TRF (see algorithm doc §4.1); full WDP for every pattern size (TRF: band for > 20).
- Detection: k-tuple trigger + sum-of-heads criterion + distance-seen suppression + copy rule + three best periods.
  Not reproduced (declared): apparent-size criterion (simulated cut-offs), best-period list (d > 250), narrow band.  **Resolved by WP7** (all three implemented; see §WP7 revision).
  Random-walk range distances: implemented in the parameter-set API (WP6 revision below); the legacy overloads keep
  the 2026-09-30 behaviour.
- Output: 0-based `Start`, exact percentages, ordered by (start, end, period); `minPeriod` applied after redundancy
  elimination; `maxPeriod ≤ 2000`, `minScore ≥ 1`; eager validation on both overloads.
- `ComputeBernoulliStatistics(tract, period)`: TRF analysis of the tract with the last `period` bases as the
  candidate; PM = matches / comparisons, PI = indels / comparisons (equal to TRF's statistics for that region).

---

## WP6 revision (2026-10-01, B04 completeness audit L8–L10)

### Sources opened (this session)

- TRF 4.10.0 README (git clone, commit 355c1f9): parameter list ("Match, Mismatch, and Delta … in the range of 3 to 7
  … recommended 2, 7, and 7"; "Probabilistic data is available for PM values of 80 and 75 and PI values of 10 and
  20"; Minscore; MaxPeriod 1..2000), options `-m` ("every location that occurred in a tandem repeat changed to the
  letter 'N'"), `-f` ("500 nucleotides on each side"), `-r`, `-l <n>` ("longest TR array expected … n million bp",
  default 2), `-ngs` ("Short 50 flanks are appended"); "Table Explanation" item 9 (entropy based on percent
  composition); "Alignment Explanation"; Waiting Time / Sum of Heads / Random Walk sections and Table 1 image
  (`images_for_readme/table1.gif`: k 4 / 5 / 7 for 1–29 / 30–159 / 160–500 — PM = .80 only).
- TRF 4.10.0 source, read to understand behaviour (not transcribed): `tr30dat.c` `get_statistics` (composition
  denominator = all non-gap symbols of the aligned repeat, entropy summed largest-p first, `log(p)/log(2)`),
  the PM-dependent initialisation (PM = 75: tuple sizes 3/4/5/7 up to distances 29/43/159/MAXDISTANCE; PM = 80:
  4/5/7 up to 29/159), `sumdata75` / `sumdata80` tables (used only as the comparison oracle), `d_range` =
  ⌊2.3·√(Pindel·d)⌋ with `Pindel = (float)PI/100`, `new_meet_criteria_3` (range summation), `newwrap` row limits
  (`maxwraplength`), `MakeMaskedFile`, `print_flanking_sequence`; `trfrun.h` MAXDISTANCE = max(200, min(max(MaxPeriod,
  500), ⌊0.6·n⌋)); `trf.c` option parsing (PM other than 75/80 → "No sum table file"; note: 4.10.0 parses the `-l`
  value from the MaxPeriod argument — `ParseUInt(av[8]…)` — a TRF bug; the library implements the documented meaning).

### Findings

- **L8:** TRF's entropy column uses p_b = count_b / (region length including N); the library's `Entropy` normalises
  over A/C/G/T, so N-containing regions differ (e.g. U1: 1.8425 vs TRF 1.83). New `EntropyTrf`.
- **L9:** weights, PM, PI, Minscore, MaxPeriod, `-l`, `-r` were hard-coded. The sum-of-heads cut-off for PM = 75 is
  reproduced by the same exact-moment derivation (tuple sizes 3/4/5/7) with floor max(k + 1, 5): **2000/2000** d
  equal to `sumdata75` (and still 2000/2000 for `sumdata80`; the floor 5 binds only for k = 3, d ≤ 20). The tuple
  ranges are TRF's configuration; each starts where the cut-off already reaches k + 1 (first such d at PM = .75:
  k = 4 → 27, 5 → 41, 7 → 91), consistent with the paper's rule. PI enters only through the random-walk range.
- **L10:** masked sequence, flanking sequence and the alignment rows were not exposed.

### Reference cross-checks (compiled TRF 4.10.0 in `scratchpad/trf/src/trf`; harness `scratchpad/wp6/`)

Test set: 700 sequences (200–1 500 bp random background, 1–3 embedded arrays of period 1–250, substitution 0–15 %,
indel 0–8 %, 25 % of sequences with N inside arrays / N runs; `gen6.py`, seed 606; 843 837 bp). TRF:
`trf set700.fa M X D PM PI S P -h -d -ngs`; C#: `FindApproximateTandemRepeats(seq, TandemRepeatsFinderParameters)`.
"Exact" = every .dat field equal (copies as TRF's float `%.1f`, truncated percentages, entropy `%.2f`, consensus);
"region" = some C# row covers ≥ 50 % of the TRF row.

| Parameters | TRF rows | exact | region | p ≤ 20 exact / region |
|---|---|---|---|---|
| 2 7 7 80 10 50 500 (recommended) | 1 305 | 87.8 % | 99.4 % | 98.3 % / 99.8 % |
| 2 5 7 80 10 50 2000 | 1 409 | 88.9 % | 99.8 % | 98.0 % / 99.9 % |
| 2 3 5 80 10 40 200 | 1 703 | 84.3 % | 99.1 % | 95.8 % / 99.9 % |
| 2 7 7 75 20 50 500 | 1 357 | 87.4 % | 99.6 % | 97.5 % / 99.9 % |
| 2 5 5 75 10 30 100 | 1 759 | 87.3 % | 99.9 % | 95.8 % / 99.9 % |
| 2 3 3 80 20 50 500 | 2 414 | 64.2 % | 93.8 % | 92.0 % / 99.0 % |
| 3 7 7 80 10 60 50 | 1 135 | 92.2 % | 99.7 % | 96.9 % / 100 % |
| 2 7 7 80 10 50 500 `-r` | 2 628 | 89.1 % | 99.7 % | 98.6 % / 99.9 % |
| 2 5 7 75 10 50 500 `-r` | 3 057 | 89.1 % | 99.8 % | 97.0 % / 99.9 % |

For comparison, the legacy overload (`maxPeriod` 500) on the same set: 83.1 % / 96.4 %; the new API without the
random-walk range: 83.6 % / 96.4 % (PM 75 / PI 20: 83.7 % / 98.2 % → 87.4 % / 99.6 % with it). The residual
differences are the declared ones (narrow-band WDP for consensus > 20 — consensus ≤ 20 rows are 98.3 % exact at the
recommended set — and the simulated waiting-time criterion, which matters most with the permissive 2 3 3 weights).
Fields of rows at the same locus agree except 1–6 copy-number / % indel values per set (different co-optimal paths).

- **Entropy (L8):** every same-locus row containing N has `EntropyTrf` printing TRF's value (229/229 recommended;
  269, 284, 243, 266, 271, 206, 444, 566 for the other sets — 0 mismatches); the legacy `Entropy` differs on
  152/229.
- **Masking (L10):** TRF `-m` vs `MaskApproximateTandemRepeats`: identical for every sequence whose reported loci are
  identical (601/601, 565/565, 582/582, 461/461 sequences for 2 7 7 80 10 50 500, 2 5 7 80 10 50 2000,
  2 7 7 75 20 50 500, 2 3 5 80 10 40 200; 0 mismatching positions); over all 700 sequences 0.33–0.60 % of positions
  differ, all from the locus differences above. Soft mask lower-cases exactly the hard-masked positions (700/700).
- **Flanks (L10):** 50-bp flanks equal to TRF `-ngs` on every same-locus row (1 147/1 147 … 2 723/2 723);
  500-bp flanks equal to TRF `-f` (alignment file, 200-sequence subset) 314/314.
- **Alignment rows (L10):** equal to the sequence / consensus lines of TRF's alignment file for 191/191 rows with
  consensus ≤ 20 and 116/123 above (TRF's narrow band picks another co-optimal path).
- **`-l` (L9):** TRF built with the `-l` value read in bp (only change): caps 60 / 120 / 250 bp → periods ≤ 20
  96.5 % / 97.4 % / 98.6 % exact (as without a cap); larger periods differ because TRF's narrow-band pass starts its
  forward rows at the backward optimum. Default (2 000 000) = uncapped for every input here.
- **Timing (Release):** 100 kb random, recommended set: 0.51 s (PM 75 / PI 20: 1.1 s).
- **Crafted rows locked in tests:** U1–U5 (`RepeatFinder_TrfParameters_Tests.cs`) — all TRF rows reproduced except
  two consensus > 20 rows (asserted on indices / period / score only; the U4 period-49 row is fully locked since WP7).


## WP7 revision (2026-10-01, B04 completeness audit L7 — TRF detection pipeline)

### Sources opened (this session)

- **TRF 4.10.0 README** (repository clone, commit 355c1f9): "Apparent Size Distribution" — *S = the distance between
  the first and last run of k heads in an iid Bernoulli sequence of length d … We estimate the distribution of S by
  simulation because we make it conditional on first meeting the sum-of-heads criterion … we determine the maximum
  number y such that 95 % of the time S > y … if PM = .75, k = 5 and d = 100, then the criterion is 56. In order to
  test the apparent-size criterion, we compute the distance between the first and last tuple on list D_d*;
  "Detection" (sum-of-heads and apparent-size tests on the distance lists, nearby distances via the random walk);
  "Narrow Band Alignment" — *we limit WDP calculations to a narrow diagonal band … for patterns larger than 20
  characters. In accordance with the random walk results, the band radius is Δd_max. The band is periodically
  recentered around a run of matches in the current best alignment*; "Multiple Reporting of Repeat at Different
  Pattern Sizes"; What's New 4.04 (*widened radius of narrowband alignment*), 4.07b (*changed alignment to go further
  when score drops to 0*). The paper itself (NAR 27:573) was not reachable from this sandbox (publisher / PMC hosts
  blocked by the egress proxy); the README reproduces its method sections.
- **TRF source, read for behaviour only** (AGPL — nothing transliterated): `newtupbo` (candidate loop, order of the
  tests, which width a pattern is aligned with), `new_meet_criteria_3` (apparent-size test as "first tuple of the
  window lies within W of the window start", the 35 % rule for range distances, range sums over *linked* distances),
  `link_Distance_window` / `no_matches_so_unlink_Distance` (a distance list becomes active when tested and inactive
  when found empty), `add_to_bestperiodlist` / `adjust_bestperiod_entry` / `search_for_range_in_bestperiodlist`
  (5 best periods per analysed region; span test i − 2d + 1 + W .. i), `multiples_criteria_4` / `GetTopPeriods`
  (3 of 5 tested), `narrowbandwrap` / `get_narrowband_pair_alignment_with_copynumber` (radii max(6, Δd_max) and
  min(2·max(6, Δd_max), ⌊size/3⌋), recentring after 3 diagonal row maxima, anchored forward start, tie rules,
  traceback through genuine zeros), `waitdata80` / `waitdata75` (used only as the oracle for the derived table).

### Derivation of the apparent-size criterion (exact, replaces TRF's simulation)

With f = position of the k-th head of the first k-run and E = the last head of the last k-run in a Bernoulli(PM)
sequence of length L, the sequence splits into a prefix ending at f (probability F(f) = p^k·[f = k] + q·p^k·N(f−k−1),
N(m) = P(no k-run in m tosses), contributing exactly k heads to R), a middle part f+1..E (a run-length Markov chain
started inside a counted run and ending inside one after S = E − f steps, accumulating the other heads R′ of R) and a
tail (a tail toss, then no k-run: T(m) = q·N(m−1), T(0) = 1). Hence
P(S = s, R ≥ x) = M(s, R′ ≥ x − k) · Σ_f F(f)·T(L − s − f), where M does not depend on L, so one chain run per tuple
size gives every d. y = max{y : P(S > y | R ≥ x) ≥ 0.95}; TRF's test (first tuple ≤ window start + W) is S ≥ y + 1,
i.e. W = max(d, 20) − y − 1. Results:

- README example (PM .75, k 5, d 100): **y = 56** (reproduced).
- C# table = independent NumPy prototype (`scratchpad/wp7/cond.py`) for **4000/4000** (d, PM).
- vs TRF 4.10.0 `waitdata80`: equal at **825/2000** d, |Δ| ≤ 1 at 1 795, Δ ∈ [−2, +3], mean +0.03 (k 4) /
  +0.28 (k 5) / +0.47 (k 7); `waitdata75`: equal 713/2000, |Δ| ≤ 1 at 1 640, Δ ∈ [−3, +4]. TRF's table is visibly
  Monte-Carlo (non-monotone in d, e.g. 62, 63, 63, 62 at d 160–163). Variants checked on a sample of d (`variants.py`): conditioning on R > x (PM 80 mean Δ −0.25 / +0.37 for k ≤ 5 / k 7,
  exact 90/159 vs 116/159 for the README definition), 95.5 % / 96 % quantiles (bias +0.8 … +5.7), L = d below 20 (no
  better); none fits clearly better, so the README's definition (S > y, R ≥ x) is kept.
- Unconditional variants (no sum-of-heads conditioning, or the plain waiting time of the first run) do not match
  the table at all (constant in d where TRF's values grow with d).

### Reference cross-checks (compiled TRF 4.10.0; harness `scratchpad/wp7/`: `run7.py`, `cmp7.py`, `rl7.py`,
`maskcmp7.py`, `bisect7.py`, `ablate.py`, C# driver `xc7`)

Same 700-sequence set as WP6 (843 837 bp). Exact = every .dat field; region = ≥ 50 % overlap.

| Parameters | TRF rows | WP6 exact / region | **WP7 exact / region** | with TRF's waitdata substituted |
|---|---|---|---|---|
| 2 7 7 80 10 50 500 | 1 305 | 87.8 / 99.4 | **99.8 / 100** (1 303) | 100 / 100 |
| 2 5 7 80 10 50 2000 | 1 409 | 88.9 / 99.8 | **99.9 / 100** (1 407) | 100 / 100 |
| 2 3 5 80 10 40 200 | 1 703 | 84.3 / 99.1 | **99.9 / 100** (1 701) | 100 / 100 |
| 2 7 7 75 20 50 500 | 1 357 | 87.4 / 99.6 | **99.9 / 100** (1 356) | 100 / 100 |
| 2 5 5 75 10 30 100 | 1 759 | 87.3 / 99.9 | **99.8 / 100** (1 755) | 100 / 100 |
| 3 7 7 80 10 60 50 | 1 135 | 92.2 / 99.7 | **100 / 100** | 100 / 100 |
| 2 3 3 80 20 50 500 | 2 414 | 64.2 / 93.8 | **99.8 / 100** (2 409) | 100 / 100 |
| 2 7 7 80 10 50 500 `-r` | 2 628 | 89.1 / 99.7 | **99.8 / 100** | 100 / 100 |
| 2 5 7 75 10 50 500 `-r` | 3 057 | 89.1 / 99.8 | **99.9 / 100** | 100 / 100 |
| `-l` 60 / 120 / 250 bp (TRF `-l`-in-bp build) | 2 766 / 2 458 / 1 545 | 68.3 / 50.3 / 57.0 | **100 / 100 / 99.8** | 100 / 100 / 100 |
| legacy overload, maxPeriod 500 | 1 305 | 83.1 / 96.4 | **99.8 / 100** | — |
| 1 Mb single sequence (1 400 concatenated records), recommended | 1 544 | — | **99.9 / 100** (1 543) | 100 / 100 |

- **Masks** (`-m`): identical for 700/700, 700/700, 700/700, 699/700 sequences (2 7 7 80 10 50 500, 2 5 7 80 10 50
  2000, 2 7 7 75 20 50 500, 2 3 5 80 10 40 200; WP6: 641, 631, 633, 536); 700/700 with TRF's table.
- **Alignment rows** vs TRF alignment file (200-sequence subset): consensus ≤ 20 192/192, **> 20 149/149** (WP6
  116/123); 500-bp flanks 341/341.
- **Ablation** (each component switched off in a scratch build; exact % on the seven sets in table order):
  no apparent-size test 92.4 / 94.6 / 93.0 / 91.5 / 93.6 / 95.4 / 92.0; full WDP instead of the band 94.2 / 93.5 /
  90.1 / 94.8 / 92.7 / 95.9 / 69.5; all distances active 99.8 / 99.9 / 99.9 / 99.9 / 99.7 / 100 / 99.8; no best-period
  list: unchanged on this set and on the 1 Mb sequence; sequences where it matters were found by a targeted search
  (two or three adjacent arrays, one with period > 250: 2 of the first 200 generated sequences, e.g. D9).
- **Residual, row by row** (`bisect7.py`: the TRF table is substituted only for d in a range, bisected to the single
  entry that restores TRF's output). All 25 non-identical (set, sequence) cases over the seven sets (16 TRF rows not reproduced exactly) are explained by exactly one
  table entry each, and in every case TRF's simulated W is one below the exact value (the exact criterion is
  marginally more permissive there):

  | Sets | Sequence | entry (d, PM) | exact W | TRF W |
  |---|---|---|---|---|
  | 2 7 7 80 10 50 500; 2 5 7 80 10 50 2000; 2 3 5 80 10 40 200; 2 3 3 80 20 50 500 | s166 | 145, 80 | 32 | 31 |
  | 2 7 7 80 10 50 500; 2 5 7 80 10 50 2000 | s342 | 43, 80 | 30 | 29 |
  | 2 7 7 80 10 50 500; 2 5 7 80 10 50 2000; 2 3 5 80 10 40 200; 2 3 3 80 20 50 500 | s599 | 119, 80 | 32 | 31 |
  | 2 3 5 80 10 40 200; 2 3 3 80 20 50 500 | s416 | 41, 80 | 30 | 29 |
  | 2 7 7 75 20 50 500; 2 5 5 75 10 30 100 | s411 | 96, 75 | 43 | 42 |
  | 2 7 7 75 20 50 500 | s514 | 157, 75 | 44 | 43 |
  | 2 5 5 75 10 30 100 | s350 | 24, 75 | 16 | 15 |
  | 2 5 5 75 10 30 100 | s457, s614 | 29, 75 | 16 | 15 |
  | 2 3 3 80 20 50 500 | s18, s253, s530 | 29, 80 | 20 | 19 |
  | 2 3 3 80 20 50 500 | s443 | 39, 80 | 29 | 28 |
  | 2 3 3 80 20 50 500 | s451 | 43, 80 | 30 | 29 |
  | 2 3 3 80 20 50 500 | s549 | 122, 80 | 32 | 31 |
  | 2 3 3 80 20 50 500 | s631 | 32, 80 | 28 | 27 |

- **Timing** (Release, 1 Mb single sequence, `scratchpad/wp7/r1m.fa`): recommended set 15.1 s (WP6) → **5.1 s**;
  permissive 2 3 3 80 20 50 500 55.5 s → **14.4 s** (compiled TRF: 1.8 s). The band replaces O(region · pattern) by
  O(region · band) for patterns > 20 and the apparent-size test removes most spurious candidates. One-off table cost:
  sum-of-heads + apparent-size tables for both PM ≈ 150 ms on first use (sum-of-heads now filled in one chain pass per
  tuple size: identical values, 2000/2000 vs `sumdata80/75`).
- **Tests locked to TRF** (`RepeatFinder_TrfDetection_Tests.cs`, D1–D9): selected from seeded random sequences
  (`search.py`, `bpl.py`, `trim.py` shortened the linking / best-period cases while the property held) such that the
  WP6 code disagrees with TRF and the disagreement disappears only with the named component (ablation builds
  `v_fullwdp`, `v_nobestlist`, `v_nolink`, reflection switch for the apparent-size table).

## WP14 revision (2026-10-01, B04 completeness audit 4a′ + 9 — caller-supplied apparent-size table, TRF output formats)

### Sources opened (this session)

- TRF 4.10.0 README (`scratchpad/trfsrc/README.md`): "Data file" (`-d`), `-h` (no HTML), `-ngs` ("more compact .dat
  output on multisequence files, … flanking sequence"), "Table Explanation", "Apparent Size Distribution".
- TRF 4.10.0 source, read for the output layout and field types only (AGPL-3.0; nothing copied):
  `trfrun.h` 179–300 (single-sequence `.dat` / `-ngs` writer: header once, `\n\nSequence: …\n\n\n\nParameters: …\n\n\n`,
  row format `%d %d %d %.1f %d %d %d %d %d %d %d %d %.2f %s`, 50-bp flanks with '.' at sequence ends), 330–600
  (multi-sequence writer, `.summary.html`), IL struct 60–80 (`float copies`, `float entropy`, int percentages);
  `trfclean.h` 99–140 (`TRFClean`: RemoveBySize → SortByIndex → RemoveRedundancy → **SortByCount**, i.e. rows in report
  order), 668–697 (`MakeFileName`), 699–835 (`OutputHTML`, `EO_MAX_TBL` 120, heading every 22 rows, `OutputHeading`);
  `tr30dat.c` 3240–3260 (`(int)(100*(float)match/x)`, `(int)(100*(double)ACGTcount/count)`, anchor label
  `%d--%d,%d,%3.1f,%d,%d` with the **double** `Copynumber` and `OUTPUTcount`), 3980 / 4239 (`OUTPUTcount` reset per
  sequence, incremented for every reported alignment before MaxPeriod filtering and redundancy elimination), 3299 /
  3563 / 3912 (`waitdata80/75[2001]`, `waiting_time_criteria = waitdata[min(2000, d)]`).
- The TRF waiting-time tables are not monotonic (PM 80: 598 decreases, PM 75: 677), so `ApparentSizeTable` is not
  required to be monotonic; validated range: 0 ≤ y(d) ≤ max(d, 20) − 1 (TRF's own entries satisfy it: y ≥ 1 / 4).

### Changes

- `TandemRepeatsFinderParameters.ApparentSizeTable` (int[2001], copied on assignment; null = exact table),
  `ExactApparentSizeTable(pm)`, `ApparentSizeTableFromWaitingTimes(w)` (y = max(d,20) − w − 1); `TrfModel` uses the
  supplied table for both apparent-size uses (criteria test and best-period-list window).
- Result fields `CopyMatches / CopyMismatches / CopyIndels` (adjacent-copy counts, TRF's integer percentages need
  them in single precision) and `OutputIndex` (= TRF `OUTPUTcount`).
- `FormatTrfDatFileHeader`, `FormatTrfDatLines(…, TrfDatLayout.Dat | Ngs)`, `FormatTrfHtmlTables`,
  `FormatTrfHtmlSummary`; C `printf("%.Nf")` rounding (exact binary value, ties to even — .NET "F" rounds ties away).
- MCP `find_approximate_tandem_repeats`: `format` (json | dat | ngs | html), `sequenceName`, `apparentSizeTable`
  (comma/space-separated string) + `apparentSizeTableKind` (apparent | trfWaitingTimes); `mask_approximate_tandem_repeats`:
  the two table parameters.

### Reference cross-checks (compiled TRF 4.10.0; harness `scratchpad/wp14/`: C# driver `xc14` using only the public
API, `cmp14.py`, `cmphtml.py`; results `res_dat.txt`)

Same 700 sequences as WP6/WP7. "TRF table" = TRF's `waitdata80/75` loaded at run time from the local TRF build
(`scratchpad/trftables.json`) and passed through `ApparentSizeTableFromWaitingTimes` → `ApparentSizeTable`.

| Parameters | `.dat` file, TRF table | `-ngs`, TRF table | same-locus `.dat` rows byte-identical, exact table | HTML pages (tables + summary), TRF table / exact table |
|---|---|---|---|---|
| 2 7 7 80 10 50 500 | identical (700/700 blocks) | identical (595/595) | 1 303/1 303 (blocks 697/700) | 596/596 / 591/596 |
| 2 5 7 80 10 50 2000 | identical | identical (619/619) | 1 407/1 407 | 620/620 / — |
| 2 3 5 80 10 40 200 | identical | identical (657/657) | 1 701/1 701 | 658/658 / — |
| 2 7 7 75 20 50 500 | identical | identical (600/600) | 1 356/1 356 | 601/601 / — |
| 2 5 5 75 10 30 100 | identical | identical (643/643) | 1 755/1 755 | 644/644 / 637/644 |
| 3 7 7 80 10 60 50 | identical | identical (550/550) | 1 135/1 135 (file identical) | 551/551 / — |
| 2 3 3 80 20 50 500 | identical | identical (678/678) | 2 409/2 409 | 679/679 / 668/679 |
| 2 7 7 80 10 50 500 `-r` | identical | identical | 2 624/2 624 | — |

- Single-sequence files (`one.fa` 516 bp; `many.fa` 150 repeats → two HTML pages; `none.fa` no repeats): `.dat`,
  `-ngs` and HTML byte-identical with both tables.
- The first HTML run differed on 7/596 pages, all in the anchor label's copy number (e.g. TRF `…,15.7,…` vs `15.6`):
  TRF prints the double-precision `Copynumber` there and the float `copies` in the cell; fixed, then 596/596.
- Supplying `ExactApparentSizeTable(pm)` equals `null` on 700/700 sequences for PM 80 and PM 75.
- The 25 WP7 residual (set, sequence) cases: all reach TRF parity with TRF's table through the public parameter
  (included in the byte-identical files above). Unit tests use only single bisected entries (s342: d 43, y 13; s350:
  d 24, w 15), synthetic all-zero / strictest tables and the exact table — TRF's tables are not embedded anywhere.
- ~~Not generated: TRF's alignment page (`.txt.html`)~~ — generated since WP17 (`FormatTrfAlignmentPages`, see §WP17).

## WP17 revision (2026-10-01, B04 completeness audit 3 item 4 — TRF alignment pages)

### Sources opened (this session)

- TRF 4.10.0 README ("Alignment explanation", `-f`, `-r`, output file names).
- TRF 4.10.0 source, read for the output layout and index conventions only (AGPL; no code copied):
  `trfrun.h` 726–760 / 799–810 / 875–879 (alignment-file heading, "Length", "ACGTcount", non-ACGT warning, "Done."),
  1086–1150 (FASTA reader: description truncated to 199 characters, letters upper-cased, non-letters dropped);
  `tr30dat.c` 1530–1600 / 1239–1290 (traceback: index of a gap column = the next sequence / pattern position),
  1674–1821 (`shift_pattern_indices`, `alt3_print_alignment`: 65 columns per row, a new row at each copy for patterns
  > 6, otherwise a blank between copies while 2·pattern columns remain, '*' line, 10-base context), 1822–1860
  (`print_alignment_headings`: anchor `%d--%d,%d,%3.1f,%d,%d`), 2760–2840 (flanks, consensus pattern 65 per line),
  2884–3215 (`get_statistics`: Matches / Mismatches / Indels, `%0.2f` of float fractions, distance table
  `\n %3d  %3d  %0.2f`, ACGTcount over the repeat), 3837–3890 (Pmatch / Pindel / tuple sizes / tuple distances with
  MAXDISTANCE), 4241 ("Found at i:%d original size:%d final size:%d" for every reported alignment);
  `trfclean.h` 520–660 (`CleanAlignments`: dropped alignments removed from their "Found at" line, so a blank line
  survives only after the last kept one; `BreakAlignments`: > 120 alignments → heading repeated, "File k of N",
  199-character `fgets` chunks, sections cut at lines starting with 'F' / 'D', closing "\nDone.").

### Changes

- `RepeatFinder.FormatTrfAlignmentPages(sequence, repeats, name, parameters, prefix, outputCount?, version)` →
  `<prefix>.<parameters>.N.txt.html` pages (same parameter style as `FormatTrfHtmlTables`; the table anchors resolve).
  Statistics are recomputed from the alignment rows by the same `CompareAdjacentCopies` routine (now also returning the
  matching-distance distribution).
- Result fields `DetectionPosition` (0-based i), `DetectionDistance` (candidate distance, "original size") and
  `OutputCount` (TRF's final `OUTPUTcount`); overload `FindApproximateTandemRepeats(string, parameters, out int
  outputCount, minPeriod)` for the empty case.
- MCP `find_approximate_tandem_repeats`: `format=html` also returns `alignmentPages`; items gain `outputCount`,
  `detectionPosition`, `detectionDistance` (additive).

### Reference cross-check (compiled TRF 4.10.0; harness `scratchpad/wp17/`: public-API driver `xc17`, `cmp17.py`;
results `res_main.txt`, `res_single.txt`, `res_singles.txt`)

| Run (700 sequences, TRF's table supplied) | alignment pages TRF / byte-identical |
|---|---|
| 2 3 3 80 20 50 500 (and `-f`) | 678 / 678 (678 / 678) |
| 2 3 5 80 10 40 200 (and `-f`) | 657 / 657 (657 / 657) |
| 2 5 5 75 10 30 100 (and `-f`) | 643 / 643 (643 / 643) |
| 2 5 7 80 10 50 2000 (and `-f`) | 619 / 619 (619 / 619) |
| 2 7 7 75 20 50 500 (and `-f`) | 600 / 600 (600 / 600) |
| 2 7 7 80 10 50 500 (and `-f`) | 595 / 595 (595 / 595) |
| 3 7 7 80 10 60 50 (and `-f`) | 550 / 550 (550 / 550) |
| 2 7 7 80 10 50 500 `-r` (and `-r -f`) | 595 / 595 (595 / 595) |

- Total 9 874 / 9 874 pages byte-identical, first run; table / summary pages of the same runs stay identical.
- Single-sequence files: `many.fa` (150 repeats → 2 pages, BreakAlignments path) and `none.fa`, 3 parameter sets with
  and without `-r -f`: 18 / 18. 60 set700 sequences as separate files (every third lower-cased, every fifth with a
  260-character description), MaxPeriod 3 and 500: 120 / 120 after passing TRF's 199-character description (the
  first run differed on the 24 long-name pages only, in the "Sequence:" line — TRF's reader truncates; the WP14 table
  pages differ the same way, so the caller passes the truncated name). These include 42 pages with no kept repeat but
  dropped alignments (extra blank line from `OutputCount`) and 22 without any alignment.
- Exact table instead of TRF's (2 7 7 80 10 50 500 / 2 7 7 75 20 50 500): 541 / 595 and 548 / 600 identical; 51 / 49 of
  the others differ only in "Found at i:" (the apparent-size criterion fires at a different scan position — exactly
  what the table governs), anchor counters or the trailing blank line; the remaining 3 / 3 are the F46 residual
  sequences whose rows differ.
- `BreakAlignments` D/F cut — reproduced and verified (WP18, F65). The split reads 199-character chunks and starts
  a new alignment at any chunk beginning with 'F' and ends the page at any chunk beginning with 'D'. TRF's reader
  keeps every letter (`trfrun.h` `LoadSequenceFromFileBenson`, A–Z / a–z upper-cased), and `-f` flank lines start at
  column 0 (`tr30dat.c` 2782–2788), so a flank line starting with IUPAC 'D' or a letter 'F' cuts the pages in files
  with more than 120 repeats. The chunk emulation (`SplitTrfAlignmentFile`) reproduces this. Reference check
  (`scratchpad/wp18/`, compiled TRF 4.10.0, `-f` and `-r -f`, 2 7 7 80 10 50 500 and 2 5 5 75 10 30 100):
  - 20 random sequences with D and F background letters (260 microsatellites each) × 2 sets × (`-f`, `-r -f`):
    342 / 342 alignment pages byte-identical.
  - 5 more random files and 3 LCG files (2 7 7 80 10 50 500 `-f`): 14 / 14 byte-identical.
  - The cuts are real. Page 1 holds 0–13 "Found at" instead of 120, and alignments after the last page are lost.
  - Test O19 locks one LCG case (both pages, TRF's table and the exact table give the same pages).
- `CleanAlignments` `strstr` anchor matching — unreachable, so no emulation is needed. Proof:
  - The kept reference is `first--last,period,copies,size,OUTPUTcount` (`tr30dat.c` 3243). The anchor of every
    section has the same form with its own count d.
  - OUTPUTcount is unique and increases in file order. The kept list is sorted by count (`SortByCount`,
    `trfclean.h` 119), so the sections tested against kept count c' before its own section all have d < c'.
  - Each string contains a single "--", so the fields must line up. Containment then needs d to begin with the
    digits of c' while d ≠ c', so d ≥ 10·c' > c'. That contradicts d < c'.
  - The "Fo" / "Do" tests in `CleanAlignments` need a lower-case 'o', which an upper-cased flank line never has.
