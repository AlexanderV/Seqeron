# Validation Report: SEQ-COMPLEX-DUST-001 — DUST Score (triplet-frequency low-complexity score)

- **Validated:** 2026-06-16   **Area:** Complexity
- **Canonical method(s):** `SequenceComplexity.CalculateDustScore(DnaSequence, int)`, `SequenceComplexity.CalculateDustScore(string, int)` (`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs:346,361,368`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (with extracted numbers)

1. **Li & Li (2025), "Finding low-complexity DNA sequences with longdust", arXiv:2509.07357** — full PDF retrieved this session (WebFetch, cached, read pages 1–6).
   - §2.1 Notations (verbatim): "ℓ(x) = Σ_t c_x(t) = |x| − k + 1 is the total number of k-mers in x."
   - §2.5 The SDUST scoring (verbatim equation): `S_S(c⃗_x) = (1/ℓ(x)) · Σ_t [ c_x(t)(c_x(t)−1)/2 ] − T`, attributed to SDUST (Morgulis et al., 2006).
   - §2.5: "SDUST hardcodes k = 3 and uses w = 64 by default." Threshold default discussion: longdust default T = 0.6; SDUST's own threshold is a separate level.
   - Direction (§2.2/§2.5): repeated k-mers raise Σ c(c−1)/2, so a HIGH score means LOW complexity.
2. **lh3/sdust reference C implementation (`sdust.c`, master)** — retrieved this session (WebFetch of raw source).
   - Score accumulation (verbatim): `*rw += cw[t]++;` — adds the *current* count before incrementing, so summing 0+1+…+(c−1) over a triplet's occurrences yields exactly c(c−1)/2; thus `rw = Σ_t c(c−1)/2`.
   - Triplet count (verbatim): `++*L;` — L is incremented once per triplet added, so L = number of triplets in the window.
   - Threshold (verbatim): `if (rw * 10 > L * T)` ⟺ `rw/L > T/10`; with default `T = 20` ⇒ score > 2.0. Defaults (verbatim): `W = 64, T = 20`.
3. **Morgulis, Gertz, Schäffer, Agarwala (2006), J Comput Biol 13(5):1028–40** — primary paper. Direct PDF (MSU mirror) was unreachable from the sandbox (curl/WebFetch timed out), but its score function is restated verbatim and attributed by source (1) above (Li 2025 §2.5) and realised verbatim by source (2) (lh3/sdust). WebSearch summary independently confirmed: "scoring function based on counting nucleotide triplet frequencies in 64-base windows", k=3, w=64, higher score ⇒ lower complexity.

### Formula check

Repo formula (DUST_Score.md §2.2, code): `S(x) = ( Σ_t c_t(c_t−1)/2 ) / (L − 2)` for triplets, generalized in code to divide by the word count `L − wordSize + 1`.

Cross-checked against source (1) §2.5 with ℓ(x) = |x| − k + 1: for k = 3, ℓ = L − 3 + 1 = **L − 2**. The repo divisor `L − wordSize + 1` equals `L − k + 1 = ℓ(x)` exactly — it IS the number-of-k-mers normalization, not an approximation. **Match confirmed.** (The repo's `MaskLowComplexity` uses threshold 2.0 = SDUST level T=20, consistent with source (2).)

The `−T` term in the SDUST scoring is the *masking threshold* applied separately (in `MaskLowComplexity`), not part of the per-sequence complexity score this unit returns; the unit returns the un-thresholded `(1/ℓ)Σ c(c−1)/2`, which is the correct decomposition.

### Edge-case semantics check

- L < k (no k-mer exists): sources do not define a score; repo convention returns 0 (documented ASSUMPTION). Defensible.
- All-distinct k-mers ⇒ every c(c−1)/2 = 0 ⇒ score 0 (max complexity). Sourced (INV-2).
- Homopolymer length L ⇒ one triplet repeated L−2 times ⇒ (L−2)(L−3)/2 / (L−2) = (L−3)/2. Sourced derivation (INV-5).

### Independent cross-check (hand computation from source (1) formula, k=3, divisor = L−2)

| Input | L | Triplet counts | Σ c(c−1)/2 | L−2 | Score |
|-------|---|----------------|-----------|-----|-------|
| `AAAAAA` | 6 | AAA=4 | 6 | 4 | 1.5 |
| `ACGTACGT` | 8 | ACG=2,CGT=2,GTA=1,TAC=1 | 2 | 6 | 0.3333… |
| `ATGC` | 4 | ATG=1,TGC=1 | 0 | 2 | 0.0 |
| `ACACACAC` | 8 | ACA=3,CAC=3 | 6 | 6 | 1.0 |
| `AAAAAAAAAA` | 10 | AAA=8 | 28 | 8 | 3.5 |
| `AATAATAA` | 8 | AAT=2,ATA=2,TAA=2 | 3 | 6 | 0.5 |

All values trace to the source-(1) formula, not to code output. INV-1..INV-5 are genuine mathematical properties.

### Findings / divergences (Stage A)

- The spec's prior "known failure mode" of dividing by `(L−wordSize)` i.e. `L−3` (one less than the k-mer count) would be wrong; both sources fix the divisor at `ℓ = L−k+1 = L−2`. The current description/code use the correct `L−2`.
- Only k=3 is source-backed; the exposed `wordSize` parameter is a documented, accepted extrapolation. No Stage-A defect.

**Stage A verdict: PASS** — formula, normalization (1/ℓ, ℓ = L−k+1), direction, defaults (w=64, threshold 2.0/level 20), and edge cases all confirmed against two retrieved authoritative sources and hand computation.

## Stage B — Implementation

### Code path reviewed

`SequenceComplexity.cs:368` `CalculateDustScoreCore`: counts overlapping words over `wordCount = L − wordSize + 1` positions (`:375`), sums `count*(count-1)/2.0` (`:391`), returns `sum / wordCount` (`:396`). Validation: null DnaSequence throws (`:348`), wordSize<1 throws (`:349,363`), null/empty string ⇒ 0 (`:364`), L<wordSize ⇒ 0 (`:370`).

### Formula realised correctly?

Yes. `sum / wordCount` = `(Σ_t c_t(c_t−1)/2) / (L − wordSize + 1)` = source (1)'s `(1/ℓ)Σ c(c−1)/2` with ℓ = L−k+1. Exact, no precision loss (integer counts in `double`, divisor > 0 on the reachable path since `L ≥ wordSize ⇒ wordCount ≥ 1`).

### Cross-verification table recomputed vs code

Ran the canonical test fixture (18 cases) — all pass with the hand-computed exact values above (M1 1.5, M2 1/3, M3 0.0, M4 1.0, M5 3.5, M7 0.5), each asserted `Within(1e-10)`.

### Variant/delegate consistency

`DnaSequence` and `string` overloads share `CalculateDustScoreCore`; M6 asserts they agree (1.5). `MaskLowComplexity` reuses the same core with k=3 and threshold 2.0 (level 20) per source (2). Consistent.

### Test quality audit (HARD gate)

- **Sourced, not code-echoed:** every MUST expected value is derived from source (1)'s formula by hand; the file header explicitly states "A wrong divisor (e.g. L−3) would FAIL these tests." A deliberately-wrong divisor (L−3) gives 2.0 for `AAAAAA` (vs asserted 1.5) — the tests would catch it. ✔
- **No green-washing:** exact-value asserts with tight `1e-10` tolerance, not ranges/Greater/Contains, for all formula cases. INV-1 (≥0) is asserted *in addition to* the exact value (M7), not in place of it. ✔
- **Coverage:** both public overloads; all five Stage-A formula cases; all-distinct (INV-2); two homopolymers (INV-5); dinucleotide repeat; overload agreement; case-insensitivity; null DnaSequence→throw; null/empty string→0; L<wordSize→0; wordSize=0→throw on both overloads. Every Stage-A branch and documented edge/error case is exercised. ✔
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6598`; `dotnet build` 0 errors. (The 4 NUnit2007 warnings are pre-existing in `ApproximateMatcher_EditDistance_Tests.cs`, unrelated and untouched.) ✔

Minor (non-defect) note: there is no exact-value test for a non-default `wordSize` (e.g. k=2), because no external source defines a value for k≠3 (ASSUMPTION 1); asserting one would violate the "sourced expectations" rule. The default k=3 path is fully covered. The MCP smoke test (`Seqeron.Mcp.Sequence.Tests/ComplexityDustScoreTests.cs`) is a separate binding unit using relational asserts and a different empty-input contract (throws on `""`); out of scope here, no change made.

### Findings / defects (Stage B)

None. Implementation faithfully realises the validated formula; tests assert sourced exact values and cover all documented branches.

## Verdict & follow-ups

- **Stage A:** PASS. **Stage B:** PASS. **End-state: CLEAN.**
- **Test-quality gate:** PASS (sourced exact values, no green-washing, full branch/edge coverage, honest green 6598/0).
- No code or test changes were required this session.
- Follow-up (optional, non-blocking): MaskLowComplexity uses a fixed-window threshold scan rather than the SDUST perfect-interval rule — already documented as an intentional simplification, not a defect for this unit (which validates the score, not the masker).
