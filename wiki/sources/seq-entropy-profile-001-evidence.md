---
type: source
title: "Evidence: SEQ-ENTROPY-PROFILE-001 (sliding-window Shannon entropy profile)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md
sources:
  - docs/Evidence/SEQ-ENTROPY-PROFILE-001-Evidence.md
source_commit: 60364fa35e17472ed4b4847deceae1f24784348f
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-ENTROPY-PROFILE-001

The validation-evidence artifact for test unit **SEQ-ENTROPY-PROFILE-001** — a **sliding-window
per-symbol Shannon entropy profile**: for each fixed-width window slid along the sequence it emits
`H = −Σ pᵢ log₂ pᵢ` in **bits** over the window's **single-character (k = 1, mono-nucleotide)**
composition. This is the **same per-window Shannon measure** already documented on
[[windowed-sequence-complexity-profile]] as the Shannon component of the `ComplexityPoint` profile —
a **second entry point** that emits the entropy profile *alone* (no linguistic complexity), not a new
measure. The method itself is written up on that concept page; this file records the source trace and
worked oracles. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for how units
are tracked, and [[k-mer-statistics]] for the k-mer (k ≥ 1) generalization of the same Shannon
formula.

## What this file records

- **Online sources:**
  - **Shannon, C. E. (1948)** *A Mathematical Theory of Communication*, *Bell System Technical
    Journal* 27(3):379–423 (DOI 10.1002/j.1538-7305.1948.tb01338.x; primary, rank 1) — the
    originating definition of entropy `H = −Σ pᵢ log pᵢ` as the measure of choice/uncertainty of a
    discrete source.
  - **Entropy (information theory)** — Wikipedia (rank 4, citing Shannon 1948): `H(X) = −∑ₓ p(x)
    log_b p(x)`; **b = 2 ⇒ bits (shannons)**; **maximum `log_b(n)`** attained by the uniform
    distribution over `n` outcomes.
  - **Entropy-Based Biological Sequence Study** — IntechOpen chapter 75997 (rank 3–4, DNA
    application): DNA entropy `yᵢ = −Σⱼ pᵢⱼ log pᵢⱼ` over the four-nucleotide alphabet
    {A,C,G,T}; **max entropy = 2 bits (log₂4)** at equal nucleotide probability; the method slides a
    **counter of width `W`** over the sequence, per-window symbol frequencies becoming the `pᵢ`.
- **Datasets (hand-derived per-window oracles, `H = −Σ pᵢ log₂ pᵢ`):**

  | Window | Counts | H (bits) |
  |--------|--------|----------|
  | `AAAA` | A=4 | **0.0** (homopolymer, deterministic) |
  | `AATT` | A=2,T=2 | **1.0** |
  | `ATGC` | A=T=G=C=1 | **2.0** (= log₂4, uniform max) |
  | `AAAT` | A=3,T=1 | **0.8112781244591328** (3:1 skew) |
  | `AATG` / `GCAA` | 2,1,1 | **1.5** |
  | `AAATTC` | A=3,T=2,C=1 | **1.4591479170272448** |

- **Sliding-window profiles** (window offsets stepped by `stepSize`; only full windows):

  | Sequence | windowSize | stepSize | Expected profile (bits) |
  |----------|-----------|----------|--------------------------|
  | `AAATGC` | 4 | 1 | `AAAT`,`AATG`,`ATGC` → [0.8112781244591328, 1.5, 2.0] |
  | `AAATGCAA` | 4 | 2 | `AAAT`,`ATGC`,`GCAA` → [0.8112781244591328, 2.0, 1.5] |

- **Documented corner cases / failure modes:** zero-probability convention `0·log 0 ≡ 0` ⇒
  homopolymer windows give H = 0; a window with every symbol equally frequent attains H = log₂(k)
  for k distinct symbols present (2 bits over the full 4-letter alphabet); **`windowSize > length`
  ⇒ no full window exists** (empty profile); `windowSize == length` ⇒ a single value.
- **Recommended coverage:** per-window uniform (2), two-symbol equal (1), homopolymer (0), skewed
  (3:1 → 0.8112781…); full profile at stepSize = 1 and > 1 with exact per-window values and the
  number/order of windows; `windowSize > length → empty`, `== length → single value`; invariants
  `0 ≤ H ≤ log₂ k` and `H ≥ 0`; case-insensitivity (case-folded before counting).

## Deviations and assumptions

One **ASSUMPTION**, non-value-affecting: the implementation computes `pᵢ` from **per-symbol (k = 1)
mono-nucleotide** frequencies of the letters present in the window (case-folded, non-letters
ignored). The cited sources define `H` over *any* symbol probability distribution; the mono-symbol
alphabet choice (rather than k-mer/block entropy) is the implementation's modelling choice, consistent
with the four-letter DNA application (max 2 bits) in the IntechOpen chapter. It changes only the
alphabet over which `pᵢ` is taken, not the formula. (The **k-mer / block** generalization of the same
Shannon formula is the separate SEQ-COMPLEX-KMER-001 / `AnalyzeKmers.Entropy` k-entropy — see
[[k-mer-statistics]].)

## Contradictions

No source contradictions — Shannon 1948 (`H = −Σ p log p`, `0·log 0 = 0`), the Wikipedia
information-theory definition (base-2 ⇒ bits, uniform max `log_b n`), and the IntechOpen DNA
application (4-letter alphabet, max 2 bits, sliding window of width W) are mutually consistent; the
hand-derived per-window and sliding-profile oracles follow directly from those definitions.
</content>
</invoke>
