# Validation Report: SEQ-ENTROPY-PROFILE-001 — Shannon Entropy Profile

- **Validated:** 2026-06-16   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateEntropyProfile(string, int windowSize=50, int stepSize=1)`; per-window kernel `SequenceStatistics.CalculateShannonEntropy(string)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Wikipedia — Entropy (information theory)** (WebFetch, 2026-06-16). Confirmed verbatim:
  - Formula: `H(X) := −∑ₓ p(x) log_b p(x)`.
  - Log base for bits: "Base 2 gives the unit of bits (or 'shannons')" → binary log.
  - Maximum: "The maximal entropy of an event with n different outcomes is log_b(n): it is attained by the uniform probability distribution." → log₂(n) for bits.
  - Zero-probability convention: "the value of the corresponding summand 0 log_b(0) is taken to be 0 … consistent with the limit lim_{p→0⁺} p log(p) = 0."
- **SciPy `scipy.stats.entropy(base=2)`** (reference implementation, run this session) and an independent from-scratch Python implementation of `H = −Σ pᵢ log₂ pᵢ` — used to cross-check every expected number, NOT the repo code.

### Formula check
The unit computes `H = −Σ pᵢ log₂ pᵢ` (bits) over per-symbol (mono-nucleotide) frequencies of each
sliding window of width W advanced by `stepSize`. This matches the Wikipedia/Shannon-1948 definition
exactly: symbols, base-2 log → bits, uniform → log₂ k maximum, 0·log0 ≡ 0. The mono-symbol alphabet
choice is a recorded modelling assumption that does not alter the formula (only the distribution over
which pᵢ is taken); it is consistent with the 4-letter / 2-bit DNA application.

### Edge-case semantics check
- Homopolymer window → H = 0 (zero-prob convention). Sourced & confirmed.
- Uniform k-symbol window → H = log₂ k (max-entropy property). Sourced & confirmed (k=4 → 2.0; k=8 → 3.0).
- W > length → no full window → empty profile. Consistent with the sliding-window definition.
- W == length → exactly one window. Consistent.
- null/empty → empty profile (guarded). Reasonable, sourced as guarded input.

### Independent cross-check (numbers — all recomputed externally this session)
| Window / config | External value (Python `−Σ pᵢ log₂ pᵢ` + SciPy base 2) |
|---|---|
| `AAAA` (A=4) | 0.0 |
| `AATT` (A2T2) | 1.0 |
| `ATGC` (1,1,1,1) | 2.0 |
| `AAAT` (A3T1) | 0.8112781244591328 (SciPy: 0.8112781244591328) |
| `AATG` (A2T1G1) | 1.5 |
| `GCAA` (G1C1A2) | 1.5 |
| `AAATTC` (A3T2C1) | 1.4591479170272448 (SciPy: …446 — last-ULP fp diff, within 1e-10) |
| profile `AAATGC`,4,1 | [0.8112781244591328, 1.5, 2.0] |
| profile `AAATGCAA`,4,2 | [0.8112781244591328, 2.0, 1.5] |
| window counts `AAATGCAA` w4 s1/s2/s3 | 5 / 3 / 2 |
| `AAA`,4,1 (W>n) | [] |
| `AATT`,4,1 (W==n) | [1.0] |
| `ACGN` / `ATUG` / `A-C-G-T` (4 distinct) | 2.0 each |
| `ACDEFGHI` (8 distinct) | 3.0 (= log₂ 8) |

Every value the spec/tests assert traces to an external reference computed this session.

### Findings / divergences
None. The description (Evidence doc, TestSpec, algorithm doc) is mathematically correct and the cited
sources genuinely say what is claimed. INV-01…INV-05 are all true properties.

## Stage B — Implementation

### Code path reviewed
- `SequenceStatistics.CalculateShannonEntropy` — `SequenceStatistics.cs:732`. Case-folds via
  `ToUpperInvariant`, counts only `char.IsLetter`, computes `−Σ (count/total)·Math.Log2(count/total)`,
  returns 0 for null/empty or total==0. Realises the formula exactly (base-2, zero-prob safe).
- `SequenceStatistics.CalculateEntropyProfile` — `SequenceStatistics.cs:957`. Guards null/empty and
  `windowSize > length` (→ empty); iterates offsets `0, step, 2·step, … ≤ n−W`; yields
  `CalculateShannonEntropy(window)` per window. Window enumeration matches INV-05 exactly.

### Formula realised correctly?
Yes. Verified by reading the code and by recomputing every cross-check value with an independent
reference (Python + SciPy) — the code agrees with the external reference to within 1e-10 on all rows.

### Cross-verification table recomputed vs code
The 16 test assertions (after additions) all pass against the actual code, and every expected value
equals the externally-computed reference above. No code-echoed expectations: each MUST/SHOULD value
is the externally-derived Shannon number, not whatever the implementation happened to return.

### Variant/delegate consistency
The profile driver delegates to the same kernel used directly by the M1–M5 kernel tests, so the
sliding profile and the per-window kernel are consistent by construction (confirmed by M6/M7 matching
the kernel values).

### Test quality audit (HARD gate)
- **Sourced, not echoed:** every expected value is the externally-computed `−Σ pᵢ log₂ pᵢ` result; a
  deliberately-wrong implementation (e.g. base-e, or T↔U merging) would fail M1–M7, C2–C4.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` on all known values; invariant tests
  (S4/S5) are in addition to, not in place of, exact-value tests. No skips/weakened asserts/widened
  tolerances introduced.
- **Coverage gaps found & fixed:** the original 14 tests omitted documented behaviors from §3.3/§6.2:
  (1) `N` counted as its own symbol, (2) no T↔U normalization (T/U distinct), (3) non-letters ignored,
  (4) the general INV-2 with k>4 (protein window exceeding 2 bits). Added **C2 (`ACGN`→2.0),
  C3 (`ATUG`→2.0), C4 (`A-C-G-T`→2.0), and an 8-symbol case (`ACDEFGHI`→3.0)** with values sourced
  from the Shannon formula and confirmed externally this session. These exercise the kernel's
  alphabet/filtering branches that were previously untested.
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6617` (was 6613; +4). Changed test file
  builds warning-free (the 4 build warnings are pre-existing NUnit2007 warnings in unrelated files).

### Findings / defects
- **Stage-B test gap (fixed this session):** missing coverage of documented `N`/T-vs-U/non-letter and
  k>4 behaviors. Fixed by adding 4 sourced exact-value tests. No code defect — the implementation was
  already correct for all of these; the tests merely lock the documented behavior to external values.

## Verdict & follow-ups
- **Stage A:** PASS — formula, log base, maximum, and zero-prob convention all confirmed against
  Wikipedia/Shannon-1948; all numeric expectations independently reproduced (SciPy + hand).
- **Stage B:** PASS — code faithfully realises the validated formula and window enumeration; all
  cross-check values match the external reference.
- **End-state:** CLEAN — no code defect; the one test-coverage gap was completely fixed in-session,
  full suite green (Failed 0).
- **Test-quality gate:** PASS.
