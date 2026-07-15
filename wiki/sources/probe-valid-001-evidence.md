---
type: source
title: "Evidence: PROBE-VALID-001 (gapped Smith–Waterman off-target probe-specificity scan)"
tags: [validation, primer, alignment]
doc_path: docs/Evidence/PROBE-VALID-001-Evidence.md
sources:
  - docs/Evidence/PROBE-VALID-001-Evidence.md
source_commit: 4de32c233ad726853dbae99f237ce61d34c3b01a
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PROBE-VALID-001

The validation-evidence artifact for test unit **PROBE-VALID-001** — a **hybridization-probe
validation / off-target specificity scan** built on a **gapped (Smith–Waterman) local
alignment** of the candidate probe against a pooled reference, replacing an earlier ungapped
Hamming scan. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for
unit tracking. The synthesized algorithm lives on the concept page
[[probe-offtarget-specificity-scan]]. It is the **specificity-checking sibling** of the
composition-rule unit [[taqman-probe-design-rules|PROBE-DESIGN-001]] in the MolTools
reagent-design family — a genuinely different algorithm (alignment-based off-target detection),
not a re-validation of the TaqMan composition rules.

## What this file records

- **Online sources (all authority rank 1, primary literature via encyclopedic/tutorial hosts):**
  - **Smith & Waterman (1981)** — the **local-alignment recurrence** with the defining
    **zero floor** (`H(i,j) = max{diag+s, gap_up, gap_left, 0}`): negative cells reset to 0,
    traceback starts at the max and stops at 0, so it returns the best-scoring *subsequence*
    (with indels via the gap-length maxima) rather than an end-to-end alignment. Motivates the
    gapped scan and the **local trimming** of mismatched tails.
  - **Altschul et al. (1990) BLAST** and **Altschul et al. (1997) Gapped BLAST/PSI-BLAST** —
    the **off-target / homology-search rationale** and the **gapped-vs-ungapped** improvement:
    a single gapped alignment detects insertions/deletions the original ungapped HSPs missed.
    The unit realizes the *gapped local alignment* property (indel-aware), **not** the
    seed-and-extend heuristic or a genome-scale index — a full SW scan is exact, not seeded.
  - **Karlin & Altschul (1990) / Altschul et al. (1990) statistics** (via NCBI's own tutorial
    and a CMU 03-711 lecture quoting the primaries) — the **E-value / bit-score / λ** machinery
    for the opt-in `ComputeLambdaNucleotide` / `ComputeKarlinAltschul` path (see formulas below).
  - **Kane et al. (2000)** 50-mer oligo microarray specificity (rank 1/5) — the empirical
    **>75% identity** cross-hybridization threshold: an off-target is called at ≥ 0.75 identity
    over the probe length (complemented by a **<14–15 contiguous complementary bases** caveat).
    0.75 is the caller-configurable default, not hard-coded.
- **Karlin–Altschul formulas (verbatim):** `E = K·m·n·e^(−λS)`; λ = unique positive root of
  `Σ p_i p_j e^(λ s_ij) = 1`; bit score `S' = (λS − ln K)/ln 2`; `E = m·n·2^(−S')`.
  Preconditions: expected per-pair score **< 0** and **at least one positive score** (so the
  positive root exists) — otherwise the theory breaks down.
- **Datasets (hand-derived, cross-checked by an independent Python re-implementation):**
  - **Indel-only off-target** — probe `ACGTACGTACGT` (12 nt); an off-target `ACGTACTGTACGT`
    (a `T` inserted after position 6) aligns `ACGTAC-GTACGT` with 12/12 identical columns and one
    gap → identity 1.0, `HasGaps = true`. The **ungapped Hamming scan (maxMismatches = 3) misses
    it** (every fixed 12-window has ≥ 6 mismatches); the gapped scan finds both the exact
    on-target at start 5 and this off-target.
  - **Indel + mismatch off-target** — off-target region `ACGTACTGTACTT`; SW **trims** the
    mismatched `TT` tail (zero floor) → 10 identical aligned columns → identity 10/12 = 0.8333;
    admitted at `minIdentity = 0.75`, rejected at 0.90.
  - **Karlin–Altschul worked example** — scheme +1/−3, uniform 0.25 freqs → **λ = 1.3740631…**
    (≈ published NCBI blastn 1.37; expected per-pair score −2.0 < 0, positive +1 exists). With
    K = 0.711, S = 30, m = 20, n = 1000 → **bit score S' = 59.9627**, **E = 1.7802e−14**
    (`K·m·n·e^(−λS)` equals `m·n·2^(−S')`). E decreases as S rises (E(31) < E(30)) and scales
    linearly in m·n (E(n=2000) = 2·E(n=1000)).
- **Corner cases / failure modes:** SW local trimming excludes a trailing mismatched tail
  (identity over the matched core); ungapped scans miss any indel-reachable off-target (a single
  indel frame-shifts every downstream base); lowering the identity cutoff admits more
  lower-similarity off-targets.

## Deviations and assumptions

- **ASSUMPTION: on-target = the first perfect ungapped full-coverage exact match** (identity 1.0,
  coverage 1.0, no gaps); additional perfect repeats and all imperfect/indel hits are off-targets.
  The literature defines specificity as intended-vs-non-intended signal but prescribes no
  algorithmic on/off label for pooled references; this is an API/labelling choice (the intended
  hybridization site is the exact complement), exposed transparently — not a sourced numeric
  constant.
- **Scoring:** BLAST DNA scoring (+2/−3, gap −2, `SequenceAligner.BlastDna`) for the SW scan;
  **published K supplied by the caller** (its closed form needs the Karlin–Altschul
  score-probability lattice machinery), λ solved from the score scheme.

No source contradictions — the SW recurrence, BLAST gapped/ungapped distinction, Karlin–Altschul
statistics, and Kane's 75% threshold are mutually consistent. Recommended coverage — **MUST:**
an indel-only off-target found by the gapped scan but missed by the ungapped Hamming scan; the
perfect on-target classified as on-target and excluded from the off-target count; exact
identity/coverage on the two hand-derived hits (1.0 and 0.8333); λ for +1/−3 equals 1.374 to
≤ 1e-6; bit-score/E-value identities and E's monotonic decrease in S / linearity in m·n.
**SHOULD:** the identity threshold gates the 0.8333 hit (in at 0.75, out at 0.90). **COULD:**
null/empty probe and null-references guards.
