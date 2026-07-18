---
type: concept
title: "Shannon entropy profile (general-alphabet sliding-window entropy series)"
tags: [analysis, algorithm]
sources:
  - docs/algorithms/Statistics/Entropy_Profile.md
  - docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md
source_commit: dd6f9c3fb21684add3c59107cae5ea989bbd3315
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-entropy-profile-001-evidence
      evidence: "Test Unit ID: SEQ-ENTROPY-PROFILE-001; Algorithm Group: Statistics; Algorithm: Shannon Entropy Profile (sliding-window per-symbol Shannon entropy); Implementation Status: Production"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:shannon-entropy
      source: seq-entropy-profile-001-evidence
      evidence: "The profile is the windowed consumer of the scalar per-symbol Shannon entropy: each window is delegated to SequenceStatistics.CalculateShannonEntropy (the general-alphabet, all-letters entry point on shannon-entropy), yielding H = −Σ pᵢ log₂ pᵢ bits per window."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:windowed-sequence-complexity-profile
      source: seq-entropy-profile-001-evidence
      evidence: "Sibling sliding-window profiler in the same complexity/entropy family, but a distinct method: SequenceStatistics.CalculateEntropyProfile emits entropy-only over ALL letters (default W=50, step=1, IEnumerable<double>), whereas SequenceComplexity.CalculateWindowedComplexity emits a ComplexityPoint (DNA-canonical A/T/G/C Shannon + linguistic complexity, default W=64, step=10)."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:k-mer-statistics
      source: seq-entropy-profile-001-evidence
      evidence: "Per-window entropy is over the single-symbol (k=1) composition of the window; the k-mer/block Shannon k-entropy H_k = −Σ p(kmer) log₂ p(kmer) owned by k-mer-statistics (CalculateKmerEntropy, SEQ-COMPLEX-KMER-001) is the correlation-aware higher-order form the spec directs callers to for order-sensitive analysis."
      confidence: medium
      status: current
---

# Shannon entropy profile (general-alphabet sliding-window entropy series)

The **Shannon entropy profile** slides a fixed-width window along a sequence and emits, per window,
the per-symbol **Shannon entropy** `H = −Σ pᵢ log₂ pᵢ` in **bits** — a **series** of complexity values
that traces local composition along the sequence. Low values mark low-complexity / repetitive tracts
(a homopolymer window is `0` bits); high values mark compositionally diverse regions (a uniform
4-nucleotide window is the `2`-bit DNA maximum). It is the **windowed consumer** of the scalar
[[shannon-entropy]] measure — specifically its **general-alphabet** entry point
`SequenceStatistics.CalculateShannonEntropy` — driven per window by
`SequenceStatistics.CalculateEntropyProfile`. Validated as test unit **SEQ-ENTROPY-PROFILE-001**
(CLEAN, PASS/PASS 2026-06-16); [[test-unit-registry]] tracks the unit,
[[seq-entropy-profile-001-evidence]] holds the source trace and worked oracles, and
[[algorithm-validation-evidence]] describes the artifact pattern. The implementation is
[[research-grade-limitations|research-grade]].

## Where it sits — the general-alphabet *entropy-only* profiler

This is a **distinct entry point** from its neighbours, not a re-derivation:

- vs. the **scalar** [[shannon-entropy]] — that reduces a whole distribution to one number; this scans
  the *same* `H` across sliding windows and returns the **series**. It reuses the scalar's
  general-alphabet kernel (`SequenceStatistics.CalculateShannonEntropy`, counting **all letters**), not
  the DNA-canonical `SequenceComplexity.CalculateShannonEntropy` (A/T/G/C-only) path.
- vs. the **windowed complexity profile** [[windowed-sequence-complexity-profile]] — that is a
  **different method** (`SequenceComplexity.CalculateWindowedComplexity`) emitting a `ComplexityPoint`
  per window (DNA-canonical Shannon **and** linguistic complexity, defaults `W = 64`, `step = 10`).
  This unit is the **entropy channel alone**, over the **general alphabet**, with defaults `W = 50`,
  `step = 1`, returning a bare `IEnumerable<double>` — no linguistic complexity, no coordinates, no
  ambiguity-code filtering.
- vs. the **k-mer k-entropy** of [[k-mer-statistics]] — per-window entropy here is over the
  **single-symbol (k = 1)** composition of the window; higher-order (block / n-mer) entropy is *not
  implemented* here, and the spec directs callers to `SequenceComplexity.CalculateKmerEntropy` for
  correlation-aware analysis.

## Core model

For each window `x` the entropy of its **symbol-frequency distribution** is

```
H = −Σᵢ pᵢ · log₂ pᵢ          (bits; convention 0·log 0 = 0)
```

where `pᵢ = countᵢ / total` over the **letters** in the window (case-folded; non-letter characters are
ignored). Base-2 logs make the unit bits. For `k` equally-likely symbols `H` attains its maximum
`log₂ k`; the DNA 4-letter alphabet peaks at `log₂ 4 = 2` bits. The profile is the ordered list of one
`H` per window of width `W` slid along the sequence.

## Method contract

`SequenceStatistics.CalculateEntropyProfile(string sequence, int windowSize = 50, int stepSize = 1)`
in `Seqeron.Genomics.Analysis` (`SequenceStatistics.cs`) is the sliding-window driver; each window is
materialized as a substring and delegated to the per-window kernel
`SequenceStatistics.CalculateShannonEntropy(string)`.

| Name | Type | Default | Constraints |
|------|------|---------|-------------|
| `sequence` | `string` | required | letters counted (case-folded via upper-case); non-letters ignored |
| `windowSize` | `int` | `50` | a window is produced only when `W ≤ length` |
| `stepSize` | `int` | `1` | window advance in symbols, `≥ 1` |

**Output:** `IEnumerable<double>` — one Shannon entropy value (bits) per window, in 0-based **offset
order**. The enumerable is **lazy** (deferred, streaming): windows are produced on demand at offsets
`0, step, 2·step, …` while `offset ≤ length − W`.

**Alphabet handling (distinct from the DNA-canonical path):** counting is over every `char.IsLetter`
symbol, so degenerate / `N` symbols count as their **own** symbol, **protein** windows can exceed
`2` bits (`k > 4`), and there is **no T↔U normalization** — `U` and `T` are distinct symbols if both
appear. Null/empty `sequence` or `windowSize > length` yields an **empty** profile (no exception).

## Invariants

- **INV-01** Every profile value `H ≥ 0` (each `−pᵢ log₂ pᵢ ≥ 0`).
- **INV-02** Every profile value `H ≤ log₂ k` (`k` = distinct symbols in the window; `≤ 2` bits for DNA,
  larger for protein).
- **INV-03** A homopolymer window yields `H = 0` (single symbol, `p = 1`, `log₂ 1 = 0`).
- **INV-04** A window with all symbols equally frequent yields `H = log₂ k` (uniform maximum).
- **INV-05** Number of windows = `⌊(n − W)/step⌋ + 1` when `W ≤ n`, else `0` (offsets
  `0, step, 2·step, … ≤ n − W`).

## Edge cases

| Case | Result | Rationale |
|------|--------|-----------|
| null / empty sequence | empty profile | guarded input |
| `windowSize > length` | empty profile | no full window exists (INV-05) |
| `windowSize == length` | single value | exactly one window |
| homopolymer window | `0.0` bits | INV-03 |
| uniform DNA window | `2.0` bits | INV-04 / INV-02 |

## Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Full profile | `O(n × W)` | `O(σ)` | `n` = length, `W` = window width, `σ` = alphabet size; one `O(W)` frequency pass per window. A suffix tree does **not** apply — this is per-window frequency counting, not substring-occurrence search. |

## Worked oracle

`CalculateEntropyProfile("AAATGC", windowSize: 4, stepSize: 1)` → `[0.8112781244591328, 1.5, 2.0]`:

- `AAAT` → `−(3/4 log₂ 3/4 + 1/4 log₂ 1/4) = 0.8112781244591328`
- `AATG` → `−(1/2 log₂ 1/2 + 2·(1/4 log₂ 1/4)) = 1.5`
- `ATGC` → `log₂ 4 = 2.0`

## Limitations

Single-symbol (k = 1) entropy does **not** distinguish sequences with identical composition but
different order (`AATT` and `ATAT` both yield `1` bit) — higher-order correlations require the k-mer
entropy of [[k-mer-statistics]]. It is alphabet-sensitive (protein windows can exceed `2` bits since
`k > 4`), and no statistical-significance / background-model correction is applied.

## References

Shannon, C.E. (1948) "A Mathematical Theory of Communication", *Bell System Technical Journal*
27(3):379–423. Wikipedia contributors, *Entropy (information theory)*. *Entropy-Based Biological
Sequence Study*, IntechOpen (Eq. 3, `yᵢ = −Σⱼ pᵢⱼ log pᵢⱼ`). A
[[research-grade-limitations|research-grade]] windowed application of standard Shannon entropy.
