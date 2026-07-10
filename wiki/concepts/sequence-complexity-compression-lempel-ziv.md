---
type: concept
title: "Sequence complexity via compression (Lempel–Ziv)"
tags: [analysis, algorithm]
sources:
  - docs/Evidence/SEQ-COMPLEX-COMPRESS-001-Evidence.md
source_commit: c2d2b19e7359c655c98c0b2b7fc08aadd63ff843
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: seq-complex-compress-001-evidence
      evidence: "Test Unit ID: SEQ-COMPLEX-COMPRESS-001 ... Algorithm: Lempel–Ziv complexity (compression-based sequence complexity)"
      confidence: high
      status: current
---

# Sequence complexity via compression (Lempel–Ziv)

Measuring how **compressible / patterned** a sequence is by the **Lempel–Ziv (LZ76)
complexity** `c(S)` — the number of distinct phrases produced when the string is parsed
left-to-right, adding a new component each time the scanned suffix stops matching anything
already seen in its own history. It is the compression-based member of the
**sequence complexity / entropy family** (Shannon entropy, DUST, linguistic, k-mer, compression).
Validated as test unit **SEQ-COMPLEX-COMPRESS-001**;
the source trace and worked oracles are in [[seq-complex-compress-001-evidence]], and
[[test-unit-registry]] tracks the unit. See [[algorithm-validation-evidence]] for the artifact pattern.

## Where it sits in the complexity family

This is a **distinct complexity measure**, not a repeat detector or an entropy of composition:

- vs. **Shannon entropy of composition** — the protein SEG detector
  [[protein-low-complexity-seg]] scores a window by `H = −Σ pᵢ·log₂ pᵢ` (bits/residue over the
  *composition* of a fixed window); LZ complexity instead counts **new substrings along the whole
  sequence**, so it is sensitive to *order/pattern buildup*, not just symbol frequencies. A
  perfectly periodic string (`1010…`) has flat composition entropy per window yet a small,
  slowly-growing LZ count.
- vs. **k-mer Shannon k-entropy** — [[k-mer-statistics]] reduces a k-mer count profile to a
  fixed-`k` Shannon entropy; LZ uses **variable-length** phrases discovered adaptively rather than
  a fixed word length.
- vs. **repeat detection** — [[repetitive-element-detection]] finds and types explicit repeated
  substrings; LZ collapses repetitiveness into a single scalar (repetitive ⇒ low `c`). LZ
  complexity is *low* exactly where those repeat/low-complexity tracts are, so the two are
  complementary views of the same signal.
- vs. **DUST triplet score** — [[dust-low-complexity-score]] is the other `SEQ-COMPLEX-*` family
  member and the DNA low-complexity **masker**: `∑_t c(c−1)/2 / (L−2)` over a fixed `k = 3` triplet
  count distribution, where a *high* score means *low* complexity (opposite numeric direction to LZ)
  and a published mask threshold (2.0, window 64). LZ uses adaptive variable-length phrases rather
  than fixed triplet counts.

## Definition and parsing rule (LZ76 exhaustive history)

`c(S)` = "the number of different sub-strings encountered as the sequence is viewed as a stream
from left to right" (Lempel & Ziv 1976). The parse walks a growing history: at each position take
the shortest new substring that is **not** already a substring of the scanned history; when found,
it becomes a new component and the cursor advances past it. Wikipedia gives the standard `O(n)`
pointer-scan pseudocode returning the component count `C`; Naereen's reference implementation
realises the same factorization with a `set` of seen substrings (grow the candidate while it is
already in the set, else add it and reset). Both are recorded verbatim in
[[seq-complex-compress-001-evidence]].

**Raw-count oracle (Naereen binary doctests):**

| Input | Components (LZ76) | c |
|-------|-------------------|---|
| `1001111011000010` | `1/0/01/11/10/110/00/010` | 8 |
| `1010101010101010` | `1/0/10/101/01/010/1010` | 7 |
| `1001111011000010000010` | `…/00/010/000` | 9 |
| `100111101100001000001010` | `…/000/0101` | 10 |

DNA traces: `ACGT` → `A/C/G/T` → c = 4 (every symbol its own component — the parsing-rule
boundary); `AAAA` → `A/AA` → c = 2; empty → c = 0.

## Normalization for cross-length comparison

Raw `c` grows with length `n` (a longer sequence yields a higher count), so cross-length
comparison uses the entropy/antropy normalization `LZ_n = c / (n / log_b n)`, where `b` = number
of **distinct symbols actually present** (alphabet size; ≤ 4 for DNA). The denominator
`b(n) = n / log_b n` is the asymptotic upper bound for a maximally complex (random) sequence, so
the normalized value `γ = c / b(n) → 1` as randomness increases. Worked: `1001111011000010`
(n = 16, b = 2, c = 8) → `8 / (16 / log₂16) = 8/4 = 2.0`.

## Corner cases and the two documented conventions

- **Homopolymer (low but not minimal complexity):** a run of one symbol adds one longer component
  per step, so `"0"×n` gives `c = ⌊(√(8n+1)−1)/2⌋` — e.g. `"0"×16` → `0/00/000/0000/00000` → **c = 5**
  (much lower than a random string of the same length, but *not* 1).
- **Single distinct symbol (b < 2):** `log_b n` is undefined for b = 1. The entropy/antropy
  reference **clamps the log base to 2** and still returns the *normalized* value, not the raw
  count: `"0"×16` → `5 / (16 / log₂16) = 1.25`. (An earlier "return raw count" reading was
  corrected to 1.25 during validation — see the [[seq-complex-compress-001-evidence|Evidence]]
  correction note.)
- **Trailing-component convention (ASSUMPTION):** the parse can end on a partial substring that
  never became "new". The Wikipedia pseudocode adds 1 for it (`if v ≠ 1 then C := C+1`); the
  Naereen set-based reference does not. Seqeron adopts the **Naereen contract** because its
  doctests give exact reproducible expected values, and the choice affects at most the last
  component by 1.
- **Monotonicity invariant:** a more repetitive sequence has `c` ≤ that of a less repetitive
  sequence of equal length (the "productivity buildup" property) — a good property-based test.

## References

Lempel A. & Ziv J. (1976) *IEEE Trans. Inf. Theory* 22(1):75–81; normalization from Vallat
`entropy`/AntroPy `lziv_complexity` (citing Zhang et al. 2009, *J. Math. Chem.* 46(4):1203–1212);
reference implementation Naereen/Lempel-Ziv_Complexity. Full citations in
[[seq-complex-compress-001-evidence]] (do not duplicate here). A
[[research-grade-limitations|research-grade]] implementation of the standard LZ76 measure.
