# Evidence Artifact: REP-APPROX-001

**Test Unit ID:** REP-APPROX-001
**Algorithm:** Approximate (imperfect / interrupted) tandem-repeat detection with the Tandem Repeats Finder model (`RepeatFinder.FindApproximateTandemRepeats`) and the TRF Bernoulli statistics (`RepeatFinder.ComputeBernoulliStatistics`)
**Date Collected:** 2026-09-30 (campaign 2026-09, batch B04)

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
  Not reproduced (declared): apparent-size criterion (simulated cut-offs), best-period list (d > 250), narrow band.
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
  two consensus > 20 rows (asserted on indices / period / score only).

