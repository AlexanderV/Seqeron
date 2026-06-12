# Validation Report: SEQ-ENTROPY-001 — Shannon Entropy

- **Validated:** 2026-06-12   **Area:** Composition
- **Canonical method(s):** `SequenceComplexity.CalculateShannonEntropy(DnaSequence)`,
  `SequenceComplexity.CalculateShannonEntropy(string)`,
  `SequenceComplexity.CalculateKmerEntropy(DnaSequence, k)`
  (wrapper smoke: `SequenceStatistics.CalculateShannonEntropy(string)` — separate impl, all-letters alphabet)
- **Stage A verdict:** PASS-WITH-NOTES (one rounded value typo in spec's cross-verification table)
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia — Entropy (information theory)** (https://en.wikipedia.org/wiki/Entropy_(information_theory)):
  confirmed verbatim:
  - Formula: `H(X) = −Σ p(x) log_b p(x)`.
  - Log base determines units; base 2 → **bits**.
  - 0·log0 convention: "the value of the corresponding summand 0 log_b(0) is taken to be 0,
    which is consistent with the limit lim p→0+ p log(p) = 0."
  - Maximum: "The maximal entropy of an event with n different outcomes is log_b(n);
    it is attained by the uniform probability distribution."
- **Shannon (1948), A Mathematical Theory of Communication** — original definition of H and the
  0·log0 = 0 limit convention. (cited in spec)
- DNA specialisation (n = 4): H_max = log2(4) = **2 bits** — consistent with Wikipedia "Sequence logo"
  information content R_i = log2(4) − H_i.

### Formula check
Spec formula `H = −Σ p_i log2 p_i` (bits, base 2) matches the authoritative source **exactly**,
including log base and the 0·log0 = 0 convention. Invariants INV-ENT-001..005 are all genuine
mathematical properties (range 0..2 for DNA, max at uniform = log2(n), homopolymer = 0,
k-mer bound log2(4^k) = 2k, empty sum = 0 / expansibility axiom).

### Edge-case semantics check
Empty, single symbol, homopolymer→0, uniform→max, non-DNA chars excluded, null guard,
invalid k guard — all defined and sourced in the spec.

### Independent cross-check (hand / Python, exact)
| Input | Formula | Computed |
|-------|---------|----------|
| `ATGCATGCATGCATGC` (uniform 4) | −4·(0.25·log2 0.25) | **2.0** |
| `AAAAAAA` (homopolymer) | −(1·log2 1) | **0.0** |
| `ATATATAT` (50/50) | −2·(0.5·log2 0.5) | **1.0** |
| `ATGATGATG` (3 uniform) | −3·(⅓·log2 ⅓) | **1.5849625…** = log2(3) |
| k-mer `ATCG` k=2 (AT,TC,CG) | −3·(⅓·log2 ⅓) | **1.5849625…** = log2(3) |
| k-mer `ATGCATGCATGCATGC` k=2 (AT4,TG4,GC4,CA3 /15) | −3·(4/15·log2 4/15) − (3/15·log2 3/15) | **1.9898981…** |

### Findings / divergences
- **Doc typo (minor, fixed):** the cross-verification table previously listed the
  `ATGCATGCATGCATGC` k=2 value as `≈ 1.98082`. The formula shown is correct; its actual value is
  **≈ 1.98990**. Corrected the displayed value in `tests/TestSpecs/SEQ-ENTROPY-001.md`.
  The corresponding test (S04) asserts the **formula** directly, so it was never affected.

---

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs:78-161`
  - `CalculateShannonEntropyCore` (l.93): fixed alphabet {A,T,G,C}, counts only those chars,
    `total==0 → 0`, `count>0` guard implements 0·log0=0, uses `Math.Log2` → **bits**.
  - `CalculateKmerEntropyCore` (l.136): overlapping windows, `length<k → 0`, same Log2/guard logic.
  - Guards: `ArgumentNullException.ThrowIfNull` (l.80,130); `k<1 → ArgumentOutOfRangeException` (l.131).
  - String overload uppercases via `ToUpperInvariant()` (l.90); DnaSequence overload feeds
    `sequence.Sequence` (already canonical uppercase) — both reach the same Core.

### Formula realised correctly?
Yes. `entropy -= p * Math.Log2(p)` over `p = count/total` is exactly `−Σ p log2 p`.
Log base = 2 (bits) and 0·log0 handled by the `count>0` / `total==0` guards — both match Stage A.

### Cross-verification table recomputed vs code (tests run)
All values above are asserted by the canonical tests; the 42 Entropy-filtered tests pass,
including exact `2.0`, `0.0`, `1.0`, `log2(3)`, and the formula-based k=2 dinucleotide case.

### Variant/delegate consistency
String vs DnaSequence overload bitwise-equal (S02 test). `SequenceStatistics` wrapper is a
distinct implementation (counts all letters via `char.IsLetter`, not a fixed DNA alphabet) and is
correctly scoped "smoke only" in the spec; its 3 smoke tests (uniform→~2, homopolymer→0, empty→0) pass.

### Numerical robustness
`p ∈ (0,1]` always; `Math.Log2` exact at powers of two; no overflow/div-by-zero on stated ranges.

### Test quality audit
Assertions check exact sourced values (not "no throw"): `Is.EqualTo(2.0)` / `0.0` / `1.0` exact,
`log2(3)` within 1e-10, k-mer formula within 1e-10. Edge cases (empty, null, k<1, length<k,
homopolymer, non-DNA) all covered. Deterministic.

### Findings / defects
None in code or tests.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — formula/log base/conventions all match authoritative sources;
  one rounded display value in the spec table was wrong (1.98082 → 1.98990) and has been corrected.
- **Stage B: PASS** — code faithfully realises `H = −Σ p log2 p` in bits with correct 0·log0 and
  edge handling; all worked examples reproduce.
- **End state: CLEAN.** Build succeeds (0 warn/err); Entropy filter 42/42 pass;
  full suite **4461 passed, 0 failed**.
- Files changed: `tests/TestSpecs/SEQ-ENTROPY-001.md` (doc value correction only; no source/test code changed).
