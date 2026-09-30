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
  Not reproduced (declared): apparent-size criterion (simulated cut-offs), random-walk range distances, best-period
  list (d > 250), narrow band.
- Output: 0-based `Start`, exact percentages, ordered by (start, end, period); `minPeriod` applied after redundancy
  elimination; `maxPeriod ≤ 2000`, `minScore ≥ 1`; eager validation on both overloads.
- `ComputeBernoulliStatistics(tract, period)`: TRF analysis of the tract with the last `period` bases as the
  candidate; PM = matches / comparisons, PI = indels / comparisons (equal to TRF's statistics for that region).
