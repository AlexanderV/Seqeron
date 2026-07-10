---
type: source
title: "Evidence: SEQ-COMPLEX-DUST-001 (DUST triplet-frequency low-complexity score)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md
sources:
  - docs/Evidence/SEQ-COMPLEX-DUST-001-Evidence.md
source_commit: dfe172057b03d328680fbf256469d7675c2604a4
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-COMPLEX-DUST-001

The validation-evidence artifact for test unit **SEQ-COMPLEX-DUST-001** — the **DUST**
low-complexity score, the triplet-frequency member of the complexity/entropy family and the
DNA-side low-complexity masker. One instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the method itself is written up on
the concept page [[dust-low-complexity-score]]. This file records the source trace and worked
oracles. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Morgulis, Gertz, Schäffer & Agarwala (2006)** *J Comput Biol* 13(5):1028–1040 (primary,
    rank 1; PMID 16796549, DOI 10.1089/cmb.2006.13.1028) — the symmetric DUST paper. DUST "is a
    heuristic algorithm that employs a scoring function based on counting nucleotide **triplet
    frequencies in 64-base windows**"; the new implementation reuses the same scoring function and
    only changes the masking rule (now **symmetric and context-insensitive**).
  - **Li, H. (2025)** *Finding low-complexity DNA sequences with longdust* (arXiv:2509.07357,
    rank 1; author of minimap2/samtools) — restates the SDUST score **verbatim**: for a string `x`
    of length `L`, `score = ∑_t c_x(t)(c_x(t)−1)/2 / (L−2)`, where `c_x(t)` is the count of triplet
    `t` and the sum runs over all triplets. Normalization is by `(L−2)` = number of triplets
    `ℓ = L−2` (k = 3 hardcoded). Scoring function written `S_S(c_x) = 1/ℓ(x) · ∑_t c_x(t)(c_x(t)−1)/2 − T`
    (threshold `T` subtracted from the raw score). Defaults: window `w = 64`; threshold =
    complexity level 20, score **2.0**. **Direction: a HIGH score indicates LOW complexity** — the
    high-scoring (low-complexity) regions are the ones masked.
  - **lh3/sdust `sdust.c`** (rank 3, Heng Li reference C implementation) — incremental accumulation
    `++*L, *rw += cw[t]++, *rv += cv[t]++;` adds the *current* count `cw[t]` to the running score
    `rw` **before** incrementing, so summing pre-increment counts `0+1+…+(c−1)` over each triplet's
    occurrences yields exactly `c(c−1)/2` and `rw = ∑_t c(c−1)/2`. Threshold comparison
    `if (rw * 10 > L * T)` ⇔ `rw/L > T/10`; with `T = 20` this is `score > 2.0`. Defaults `int W = 64;`
    and `int T = 20;`.
- **Datasets (hand-derived oracles, k = 3, divisor = number of triplets = L−2):**

  | Input | L | Σ c(c−1)/2 | L−2 | Score |
  |-------|---|-----------|-----|-------|
  | `ATGC` | 4 | 0 (ATG=1, TGC=1) | 2 | **0.0** |
  | `ACGTACGT` | 8 | 2 (ACG=2, CGT=2) | 6 | **0.333…** |
  | `AAAAAA` | 6 | 6 (AAA=4) | 4 | **1.5** |
  | `ACACACAC` | 8 | 6 (ACA=3, CAC=3) | 6 | **1.0** |
  | `AAAAAAAAAA` | 10 | 28 (AAA=8) | 8 | **3.5** |

  Threshold reference: default window **64**, default mask-if-score-> threshold **2.0**, triplet
  word size **k = 3**.
- **Documented corner cases:** below triplet length (`L < 3`) → score `∑c(c−1)/2 / (L−2)` undefined
  (no triplets, `L−2 ≤ 0`); all-distinct triplets → every `c(c−1)/2 = 0` → score **0** (maximum
  complexity); maximally repetitive homopolymer of length `L` → one triplet repeated `L−2` times →
  maximal score `(L−2)(L−3)/2 / (L−2) = (L−3)/2`.
- **Recommended coverage:** the five hand-derived oracle scores exactly; DnaSequence/string
  overloads agree; case-insensitivity (`ToUpperInvariant`); null DnaSequence → ArgumentNullException,
  null/empty string → 0; input shorter than wordSize → 0.

## Deviations and assumptions

Two **ASSUMPTIONs**, both flagged in the artifact (the paper and reference impl hardcode k = 3):

1. **General word size `wordSize`** — the repository method exposes a `wordSize` parameter; for
   `wordSize = w` the normalization generalizes to the number of words `L − w + 1` (= `L − 2` when
   `w = 3`). Consistent with the formula, but **only k = 3 is source-backed**; tests assert exact
   source-derived values only for k = 3.
2. **Input shorter than one word (`L < wordSize`)** — neither source defines a score when no word
   exists; the implementation returns **0** (no repeats ⇒ minimal complexity). A defined-output
   convention, not a source value.

## Contradictions

No source contradictions — Morgulis et al. 2006, Li 2025, and lh3/sdust agree on the
`∑_t c(c−1)/2 / (L−2)` triplet-frequency score, the `w = 64` / `T = 20` (score 2.0) defaults, and
the HIGH-score ⇒ LOW-complexity direction. The lh3/sdust incremental accumulation is proven
algebraically equal to the closed-form `c(c−1)/2` sum.
</content>
</invoke>
