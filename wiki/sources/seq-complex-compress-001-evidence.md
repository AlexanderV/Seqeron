---
type: source
title: "Evidence: SEQ-COMPLEX-COMPRESS-001 (Lempel–Ziv compression-based sequence complexity)"
tags: [validation, analysis]
doc_path: docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md
sources:
  - docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md
source_commit: c2d2b19e7359c655c98c0b2b7fc08aadd63ff843
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: SEQ-COMPLEX-COMPRESS-001

The validation-evidence artifact for test unit **SEQ-COMPLEX-COMPRESS-001** — sequence complexity
by **Lempel–Ziv (LZ76) complexity** `c(S)`, the compression-based member of the complexity/entropy
family. One instance of the templated per-algorithm [[algorithm-validation-evidence|evidence
artifact]] pattern; the method itself is written up on the concept page
[[sequence-complexity-compression-lempel-ziv]]. This file records the source trace and worked
oracles. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Lempel & Ziv (1976)** *IEEE Trans. Inf. Theory* 22(1):75–81 (primary, rank 1; paper text
    paywalled) — complexity ties to "the gradual buildup of new patterns along the sequence", i.e.
    the number of distinct phrases produced by a left-to-right parse.
  - **Wikipedia "Lempel–Ziv complexity"** (rank 4, cites the 1976 primary) — verbatim definition
    ("number of different sub-strings encountered as the binary sequence is viewed as a stream from
    left to right"), the delimiter parsing rule, and the standard `O(n)` pointer-scan pseudocode
    returning component count `C` (with the trailing `if v ≠ 1 then C := C+1`).
  - **Naereen/Lempel-Ziv_Complexity** (rank 3 reference impl, MIT) — the `set`-of-seen-substrings
    factorization and four verbatim binary doctests with exact return values.
  - **entropy / AntroPy `lziv_complexity`** (rank 3, cites Lempel-Ziv 1976 & Zhang 2009) — same
    "number of different substrings" definition plus the normalization `LZ_n = LZ / (n / log_b n)`
    (b = number of unique characters), motivated by raw LZ's length dependence.
  - **Asymptotic upper bound** (rank 4 search synthesis) — `b(n) = n / log_σ n` for uniform
    symbols; `γ = c(n)/b(n) → 1` for a maximally complex (random) sequence.
- **Datasets (documented oracles):**
  - *Naereen binary doctests (raw c):* `1001111011000010` → **8**; `1010101010101010` → **7**;
    `1001111011000010000010` → **9**; `100111101100001000001010` → **10** (component lists given).
  - *Normalized LZ (derived from the entropy formula):* `1001111011000010` (n=16, b=2, c=8) →
    `8/(16/log₂16) = 8/4` = **2.0**; `"0"×16` (b=1→2 clamped, c=5) → `5/4` = **1.25**.
  - *Traced raw values:* `"AAAA"` → c=2 (`A/AA`); `"ACGT"` → c=4 (`A/C/G/T`); `""` → c=0;
    `"0"×16` → c=5 (`0/00/000/0000/00000`, general homopolymer `c = ⌊(√(8n+1)−1)/2⌋`).
- **Documented corner cases:** empty input → 0; homopolymer is low but **not** 1; trailing
  incomplete component counted by the Wikipedia pseudocode but not by the Naereen set variant;
  raw LZ grows with `n` (normalize for cross-length comparison); single distinct symbol (b<2)
  makes `log_b n` undefined.
- **Recommended coverage:** the four Naereen raw doctest values exactly; normalized 2.0 for
  `1001111011000010`; homopolymer `"0"×16` raw 5 + b<2 normalized 1.25; `"ACGT"` → 4;
  monotonicity (more repetitive ⇒ ≤ complexity at equal length); DNA b=4 normalization uses log₄.

## Deviations and assumptions

Two **ASSUMPTIONs**, both flagged in the artifact, neither invented:

1. **Trailing-component convention** — the exhaustive-history parse can end on a partial substring
   that never became "new". Wikipedia's pseudocode adds 1 for it; the Naereen set-based reference
   does not. The unit adopts the **Naereen contract** because it has exact reproducible doctest
   values; the choice affects at most the last component by 1.
2. **Normalization log base = alphabet size (b ≥ 2)** — entropy/antropy use `b` = number of unique
   symbols present (≤4 for DNA). For **b < 2** (single-symbol input) `log_b n` is undefined; the
   entropy/antropy reference **clamps the base to 2** and returns the *normalized* value, not the
   raw count (`"0"×16` → 1.25). This corrects an earlier "return raw count" reading (implementation
   and test M8 corrected to 1.25 on 2026-06-16).

## Contradictions

No source contradictions among the definition sources — Lempel-Ziv 1976, Wikipedia, Naereen, and
entropy/antropy agree on the "number of distinct substrings" LZ76 measure and the
`n / log_b n` normalization. The only internal correction was the b<2 clamp reading (raw → 1.25),
resolved by reading the antropy source directly.
