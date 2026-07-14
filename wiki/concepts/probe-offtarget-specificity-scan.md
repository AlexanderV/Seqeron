---
type: concept
title: "Probe off-target specificity scan (gapped Smith–Waterman)"
tags: [primer, alignment, algorithm, validation]
sources:
  - docs/Evidence/PROBE-VALID-001-Evidence.md
  - docs/algorithms/MolTools/Probe_Validation.md
source_commit: f54e82403f6b4c465967bf24771be22113f31606
created: 2026-07-10
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: probe-valid-001-evidence
      evidence: "Test Unit ID: PROBE-VALID-001 ... Algorithm: Hybridization Probe Validation — Gapped (Smith–Waterman) Off-Target Scan"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:taqman-probe-design-rules
      source: probe-valid-001-evidence
      evidence: "Specificity-checking sibling of PROBE-DESIGN-001 in the MolTools reagent-design PROBE family; validates a candidate probe against a pooled reference rather than checking composition rules"
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:alignment-statistics
      source: probe-valid-001-evidence
      evidence: "Smith & Waterman (1981) local-alignment recurrence with the zero floor; the gapped scan reports per-hit identity/coverage over the locally-aligned columns"
      confidence: medium
      status: current
---

# Probe off-target specificity scan (gapped Smith–Waterman)

The **specificity / off-target check** for a hybridization probe (test unit
**PROBE-VALID-001**): align the candidate probe against a **pooled reference** with a
**gapped local alignment (Smith–Waterman)** and report every site whose identity clears a
threshold. This is the **specificity-checking sibling** of the composition-rule unit
[[taqman-probe-design-rules|PROBE-DESIGN-001]] — a genuinely **distinct algorithm** in the
MolTools reagent-design family, not a re-validation of the TaqMan rules. Literature-traced in
[[probe-valid-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern.

## The default validation record (`ValidateProbe`)

The per-algorithm spec (`docs/algorithms/MolTools/Probe_Validation.md`, unit PROBE-VALID-001)
documents the **default** `ValidateProbe(probe, references, maxMismatches=3,
selfComplementarityThreshold=0.3)` path — the ungapped scan whose behaviour is **unchanged**; the
gapped Smith–Waterman scan below is an **opt-in supplement**, not a rewrite of the default. The
default builds a small validation record:

- **Specificity score** from the total hit count `h` across all references:
  `h == 0 → 0.0`, `h == 1 → 1.0`, `h > 1 → 1.0 / h` (a cross-hybridization multiplicity penalty).
  Invariant `0 ≤ SpecificityScore ≤ 1`. Note this default `OffTargetHits` **pools** the intended
  on-target with off-targets — the on/off separation only exists in `ScanOffTargetsGapped`.
- **Self-complementarity** — fraction of aligned matches of the probe against its own reverse
  complement (`0 ≤ value ≤ 1`); recorded as an issue when it exceeds `selfComplementarityThreshold`
  (default 0.3).
- **Secondary structure** — a sequence-level **hairpin-potential** screen (always run), reported
  as the `HasSecondaryStructure` flag.
- **`IsValid`** is true when **no issues** were recorded, *or* under the fallback rule
  `OffTargetHits ≤ 1 && SelfComplementarity ≤ 0.4`.
- **Empty probe** → a structured *invalid* result (`SpecificityScore = 0.0`, `OffTargetHits = 0`,
  issue `"Empty probe sequence"`) rather than an exception; a **null** probe or reference
  collection throws `ArgumentNullException`. Inputs are uppercased before analysis.

A separate exact-hit helper, **`CheckSpecificity(probe, ISuffixTree)`**, counts exact occurrences
of the probe in a pre-built suffix tree — `O(m)` in the probe length, exact hits only (no
mismatch/indel tolerance), for uniqueness checks against an indexed genome.

## Why gapped local alignment (the core improvement)

The default `ValidateProbe` uses an **ungapped Hamming scan** (fixed-width windows,
`maxMismatches = 3`). Its limitation: an off-target reachable only through an **insertion or
deletion** is invisible, because a single indel frame-shifts every downstream position so every
fixed window mismatches heavily. The opt-in **Smith–Waterman** scan (Smith & Waterman
1981) closes this gap — the recurrence's gap-length maxima admit indels, and the **zero floor**
(negative cells reset to 0; traceback starts at the max score and stops at 0) makes it *local*,
returning the best-scoring sub-alignment rather than an end-to-end one. This is exactly the
**gapped-vs-ungapped BLAST** improvement (Altschul et al. 1990 → 1997): "gapped BLAST produces a
single alignment with gaps, detecting insertions and deletions the original version missed." The
unit realizes the **gapped-local-alignment** property (indel-aware, and *exact* — a full SW scan,
not the seed-and-extend heuristic or a genome-scale index).

### Consequence: local trimming

Because traceback stops at 0, a **trailing mismatched tail is excluded** from the reported hit —
identity is over the matched core. On the hand-derived indel+mismatch off-target
`ACGTACTGTACTT`, SW trims the `TT` tail and reports 10 identical aligned columns → identity
10/12 = **0.8333**, not a full-length figure.

## The off-target call (identity threshold)

An aligned site is called an **off-target** when its identity over the probe length clears a
cutoff. The default **0.75** comes from **Kane et al. (2000)**: on 50-mer oligo microarrays, a
non-target transcript **>75% similar** over the probe target "may show cross-hybridization"
(complemented by a **<14–15 contiguous complementary bases** caveat). 0.75 is a
**caller-configurable parameter**, not a hard-coded constant — the 0.8333 hit is admitted at
0.75 and rejected at 0.90.

**On/off labelling (assumption).** The literature defines specificity as intended-vs-non-intended
signal but prescribes no algorithmic label for a pooled reference. The unit classifies the
**first perfect ungapped full-coverage exact match** (identity 1.0, coverage 1.0, no gaps) as the
intended **on-target**; additional perfect repeats and all imperfect/indel hits are off-targets.
This is an API/labelling choice (the intended hybridization site is the exact complement), exposed
transparently — not a sourced numeric constant.

### Worked oracle

Probe `ACGTACGTACGT` (12 nt) vs a reference embedding an exact copy at start 5 and an
indel copy `ACGTACTGTACGT` (a `T` inserted after position 6) at start 27:

- The exact copy → on-target, identity 1.0, coverage 1.0, no gaps.
- The indel copy → gapped alignment `ACGTAC-GTACGT`, 12/12 identical columns, one gap, identity
  1.0, `HasGaps = true` — found by the gapped scan but **missed** by the ungapped Hamming scan
  (every fixed 12-window over the indel region has ≥ 6 mismatches).

Scoring: BLAST DNA (+2/−3, gap −2, `SequenceAligner.BlastDna`).

## Opt-in Karlin–Altschul significance statistics

An opt-in path (`ComputeLambdaNucleotide` / `ComputeKarlinAltschul`) attaches BLAST-style
significance to a raw score, per **Karlin & Altschul (1990)** / **Altschul et al. (1990)**:

- **E-value** `E = K·m·n·e^(−λS)` (m = query length, n = database length, S = raw score).
- **λ** is the unique positive root of `Σ p_i p_j e^(λ s_ij) = 1`.
- **Bit score** `S' = (λS − ln K)/ln 2`, giving the equivalent `E = m·n·2^(−S')`.
- **Preconditions:** the expected per-pair score must be **< 0** and at least one score must be
  **positive** (so the positive root exists); otherwise the statistical theory breaks down.
- **K** is caller-supplied (its closed form needs the Karlin–Altschul score-probability lattice);
  λ is solved from the scoring scheme.

**Cross-check:** for +1/−3 scoring with uniform 0.25 base frequencies, the solver gives
**λ = 1.3740631…**, matching NCBI blastn's published 1.37 (expected per-pair score
0.25·1 + 0.75·(−3) = −2.0 < 0, and +1 exists). With K = 0.711, S = 30, m = 20, n = 1000 →
bit score **59.9627**, **E = 1.7802e−14**; E decreases as S rises and scales linearly in m·n.

## Distinction from the TaqMan design rules

[[taqman-probe-design-rules]] (PROBE-DESIGN-001) checks **sequence-composition** constraints of a
5'-nuclease hydrolysis probe (no 5'-G, more Cs than Gs, GC/length windows, probe-Tm gate).
PROBE-VALID-001 does something orthogonal: it takes a candidate probe and **aligns it against a
reference** to find cross-hybridizing off-targets. Same PROBE family and reagent-design goal
(a probe that binds only its intended target), different machinery — composition rules vs
gapped local alignment plus BLAST/Karlin–Altschul statistics.

## Change history captured

The gapped SW scan + on/off-target separation were added (2026-06-24) as an **opt-in supplement**
to the default ungapped-Hamming `ValidateProbe` (whose behaviour is unchanged); the
Karlin–Altschul λ / bit-score / E-value statistics were added the same day as a further opt-in
extension.
