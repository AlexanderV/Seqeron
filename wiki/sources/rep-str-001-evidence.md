---
type: source
title: "Evidence: REP-STR-001 (Microsatellite/STR detection — approximate TRF model)"
tags: [validation, genomic-analysis]
doc_path: docs/Evidence/REP-STR-001-Evidence.md
sources:
  - docs/Evidence/REP-STR-001-Evidence.md
source_commit: 32e9f11dfac8dcd66d61e1eff994bb086758eb81
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: REP-STR-001

The validation-evidence artifact for test unit **REP-STR-001** — **Microsatellite / Short
Tandem Repeat (STR) detection**: a perfect detector (default `FindMicrosatellites`) plus the
opt-in **approximate / imperfect / interrupted** detector `FindApproximateTandemRepeats` and
`ComputeBernoulliStatistics`, modelled on Benson's **Tandem Repeats Finder** (1999). One
instance of the templated per-algorithm [[algorithm-validation-evidence|evidence artifact]]
pattern; the algorithm sits in the [[repetitive-element-detection]] repeats family (tandem
sub-problem, microsatellite/STR by unit length). See [[test-unit-registry]] for how units are
tracked.

This unit is the concrete implementation that **closes the "exact-copies-only" gap** the
repeats family previously documented as a Framework/Simplified limitation: the approximate
detector reports an interrupted tract as ONE repeat with quantified imperfection, where the
perfect detector fragments it.

## What this file records

- **Online sources** (accessed 2026-06-24):
  - **Benson, G. (1999) — Tandem Repeats Finder**, *Nucleic Acids Research* 27(2):573–580
    (rank 1, canonical TRF model) — the *approximate*-repeat definition ("two or more
    contiguous, approximate copies of a pattern"), the seven reported statistics (indices,
    period size, copy number, consensus size, percent matches, percent indels, alignment
    score), Smith-Waterman-style **wraparound dynamic programming (WDP)** alignment, and
    **consensus by majority rule** from the alignment.
  - **TRF definitions / desc pages** (`tandem.bu.edu`, rank 2) — verbatim statistic
    definitions; the **Bernoulli model** (aligning two adjacent copies = *n* independent
    coin-tosses), **PM = matching probability** (average percent identity between copies),
    **PI = indel probability**; statistics are between **adjacent copies**, not copy-vs-
    consensus; defaults **PM = .80, PI = .10**. The sum-of-heads `R(d,k,pM)` and random-walk
    `W(d,pI)` k-tuple *seeding* machinery is noted as the non-reproducible residual.
  - **TRF reference implementation** (Benson-Genomics-Lab/TRF, rank 3) — usage
    `trf File Match Mismatch Delta PM PI Minscore MaxPeriod`, recommended `(Match,Mismatch,
    Delta) = (2, 7, 7)`, `Minscore` default 50 (≈25 aligned chars at match weight 2).
  - **Wikipedia: Microsatellite** (rank 4) — unit "ten nucleotides or less" (library 1–6 bp),
    "repeated 5–50 times"; worked examples `TATATATATA` (di-), `GTCGTCGTCGTCGTC` (tri-).
- **Datasets (documented oracles)** — worked scoring with match `+2`, mismatch/indel `−7`:
  - Perfect dinucleotide `CACACACACA` → period 2, `CA`×5, score 20, 100% matches, 0% indels.
  - Interrupted trinucleotide `CAGCAGCAGTAGCAGCAG` (copy 4 `C→T`) → period 3, `CAG`×6, 17
    match/1 mismatch, score 27, 94.4% matches; perfect detector reports only `CAG`×3.
  - Interrupted dinucleotide `CACACATACACA` → period 2, `CA`×6, score 15, 91.6% matches.
  - Single-base deletion `CAGCAGCAGCAGCAGAGCAGCAGCAGCAG` (29 bp) → 29 match/1 indel, score
    **51 ≥ 50** (reported at default Minscore), 96.6% matches, 3.3% indels, 9.6 copies.
  - Bernoulli PM/PI between adjacent copies: perfect `CA` → PM 1.00, PI 0; `CAG…TAG` → PM
    13/15 ≈ 0.867; `CA…T…` → PM 0.80 (exactly on the default 0.80 threshold, inclusive);
    unrelated `ACACTGTG` → PM 0.00. E[matches] = PM·d is the only source-stated moment.
- **Recommended coverage** — MUST: perfect-dinucleotide control, interrupted-substitution
  one-repeat vs perfect-fragmentation, single-base-deletion percent-indels, Minscore gate;
  Bernoulli PM=1/PI=0, adjacent-copy PM distinct from consensus percent-matches, indel PI>0.
  SHOULD: scoring constants `(+2,−7,−7)`, PM threshold + exposed defaults. COULD: determinism
  / null / empty / invalid-parameter guards.

## Corner cases and assumptions

- **Approximate copies are intrinsic** — substitutions and indels within a tract are expected,
  quantified by percent-matches / percent-indels, not treated as errors.
- **Consensus vs period** — Benson allows `ConsensusSize ≠ Period`; the implemented subset
  builds the consensus by majority rule over period-aligned columns, so `ConsensusSize ==
  Period` by construction.
- **ASSUMPTION (deterministic exhaustive scan):** the subset enumerates every (start, period)
  window and scores by alignment instead of TRF's probabilistic k-tuple seeding and the sum-of-
  Bernoulli significance test — the honest residual in `LIMITATIONS.md` / algorithm doc §5.3.
  It changes *which windows are examined*, not the *statistics* of a reported repeat.
- **ASSUMPTION (percentage denominator = total alignment columns):** Benson gives no verbatim
  formula for percent-matches/indels; the fraction-of-columns convention reproduces the worked
  numbers. Match/mismatch/indel weights, Minscore, and statistic *names* are source-verbatim.
- **ASSUMPTION (adjacent-copy segmentation for Bernoulli PM/PI):** exact for substitution-only
  and perfect tracts; for indel-bearing tracts only the qualitative property (PI>0; match+
  mismatch+indel fractions partition the trials) is asserted, not a fragile exact PI. The
  `R(d,k,pM)` / `W(d,pI)` closed-form percentiles are NOT reproduced (non-redistributable TRF
  simulation tables — the genome-scale seeding residual).

No source contradictions. Wikipedia's encyclopedic definition, the Benson 1999 paper, the TRF
tool docs, and the reference implementation agree; the perfect vs approximate detectors are
complementary (the approximate one is opt-in).
