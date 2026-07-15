---
type: concept
title: "Position weight matrix (log-odds PWM/PSSM construction and threshold scanning)"
tags: [motif, algorithm]
sources:
  - docs/algorithms/Pattern_Matching/Position_Weight_Matrix.md
source_commit: 19070d6ba2f6b3d30d50a67db9183a714db89787
created: 2026-07-15
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: position-weight-matrix
      evidence: "Test Unit ID: PAT-PWM-001 ... Algorithm Group: Pattern Matching; Position Weight Matrix (MotifFinder.CreatePwm / MotifFinder.ScanWithPwm)."
      confidence: high
      status: current
    - predicate: alternative_to
      object: concept:consensus-from-alignment
      source: position-weight-matrix
      evidence: "Section 7.1 lists the consensus sequence ('one best character per column') as a related motif representation; the PWM's own Consensus field (INV-02) is exactly that per-column best base, but the PWM keeps the full weighted log-odds score per column rather than collapsing to a single character — same aligned-sequence input, richer model."
      confidence: high
      status: current
---

# Position weight matrix (PWM / PSSM)

A **position weight matrix** (PWM, a.k.a. **position-specific scoring matrix / PSSM**) models a
conserved DNA motif as **per-position, per-base log-odds scores**. It is the **weighted /
probabilistic** branch of the motif family: where the exact [[known-motif-search]] tests character
equality and [[iupac-degenerate-matching]] tests allowed-base-set membership, a PWM assigns each
`(base, column)` a **real-valued score** and ranks candidate windows by their **summed** score.
Seqeron exposes it as `MotifFinder.CreatePwm(...)` (build the matrix from aligned sequences) and
`MotifFinder.ScanWithPwm(...)` (slide it across a subject and emit above-threshold windows).
Validated under test unit **PAT-PWM-001**; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the per-unit artifact pattern.

## Construction — from an alignment to log-odds (the three matrices)

Building follows the standard PFM → PPM → PWM chain (Wikipedia "Position weight matrix"):

1. **PFM (counts).** `PFM[k,j] = Σᵢ 1(Xᵢⱼ = k)` — occurrences of base `k` in column `j` over the
   `N` aligned training sequences.
2. **PPM (smoothed frequency).** `PPM[k,j] = (PFM[k,j] + p) / (N + |Σ|·p)`, where `p` is the
   **pseudocount** (default `0.25`) and `|Σ| = 4` for DNA. Pseudocounts avoid `log 0` for bases
   unseen in a column (Nishida et al. 2008).
3. **PWM (log-odds).** `PWM[k,j] = log₂(PPM[k,j] / b_k)`, with the background **fixed at
   `b_k = 0.25`** for every base.

The window score is the sum of the per-column contributions:
`Score(S) = Σⱼ PWM[Sⱼ, j]`.

## Matrix layout and derived properties

The matrix is a fixed `double[4, Length]` with **Row 0 = A, 1 = C, 2 = G, 3 = T** (`Length` = motif
width = the aligned training length, INV-01). The `PositionWeightMatrix` carries three derived
fields:

| Field | Definition | Invariant |
|---|---|---|
| `Consensus` | highest-scoring base in each column | INV-02 — `GenerateConsensus()` picks the per-column maximum |
| `MaxScore` | Σ of per-column **maxima** | INV-03 — aggregate columnwise maxima |
| `MinScore` | Σ of per-column **minima** | INV-03 — aggregate columnwise minima |

`MaxScore`/`MinScore` bound any achievable window score, so the `threshold` can be read relative to
the matrix's own dynamic range.

## Scanning contract

`ScanWithPwm(sequence, pwm, threshold)` slides the PWM over every window and emits a `MotifMatch`
for each window scoring **`>= threshold`** (default `0.0`):

| Aspect | Behaviour |
|---|---|
| `Position` | 0-based window start |
| `MatchedSequence` | the scored subject substring |
| `Pattern` | set to `pwm.Consensus` (the reported hit is labelled with the motif's consensus) |
| `Score` | the summed log-odds value for that window |
| Match rule | `score >= threshold` (inclusive `>=`) |
| Invalid base | a window is dropped if any character maps to a negative base index (outside A/C/G/T) |
| Complexity | `CreatePwm` `O(N·L)`, `ScanWithPwm` `O(S·L)`, `O(1)` scan auxiliary |

## Preconditions and edge cases

| Case | Behaviour |
|---|---|
| `sequences` null | `ArgumentNullException` |
| empty training set | `ArgumentException` (≥1 aligned sequence required) |
| unequal training lengths | `ArgumentException` (columns require alignment) |
| non-ACGT training char | `ArgumentException` (strict validation, consistent with `CreateConsensusFromAlignment`) |
| single training sequence | valid PWM (`Count = 1` permitted) |
| null `sequence` / `pwm` in scan | `ArgumentNullException` |
| subject shorter than PWM | no matches (the scan loop never runs) |

Training sequences are upper-cased before counting.

## Scope, simplifications, and siblings

This is the **weighted-scoring** member of the MotifFinder motif family. Deliberately simplified
(from the spec's conformance section):

- **Uniform background fixed at `0.25`** — callers cannot model GC-biased / non-uniform backgrounds.
- **DNA-only, 4-row matrix** — no protein alphabets or arbitrary symbol sets.
- **No gap / profile-HMM state transitions** — insertions/deletions are out of scope; use a
  profile-HMM model where gap-aware scoring is required.
- Returns **raw log-odds scores**, not calibrated p-values or probabilities.

Siblings in the motif family (`MotifFinder`):

- **Best-base-per-column collapse** — [[consensus-from-alignment]] (`CreateConsensusFromAlignment`)
  and the ambiguity-coded [[iupac-degenerate-consensus]] (`GenerateConsensus`) reduce an alignment to
  a single string; the PWM keeps the full per-column weight vector. The PWM's own `Consensus` field is
  exactly the best-base-per-column reduction, which is why this unit is the `alternative_to` the plain
  plurality consensus — same aligned input, richer model.
- **Exact / degenerate matching** — [[known-motif-search]] (equality) and
  [[iupac-degenerate-matching]] (allowed-base-set membership) are the unweighted match/non-match
  branches; the PWM replaces the boolean per-position test with an additive score.
- A worked biological instance of PWM scanning is [[splice-acceptor-site-prediction]] (the
  `AcceptorPwm` 3′ splice-site scorer).

**References:** Wikipedia "Position weight matrix" (PFM→PPM→PWM chain), Nishida et al. 2008
(pseudocounts), Kel et al. 2003 (MATCH TFBS scanning), Stormo 2000 (DNA binding-site representation),
Rosalind CONS (profile/consensus). **No source contradictions**; the two API-shape assumptions —
consensus used as the `MotifMatch.Pattern` label and the inclusive `>=` threshold — are documented
implementation choices, not deviations from the cited theory.
