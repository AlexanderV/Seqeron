# Evidence Artifact: REP-STR-001

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
