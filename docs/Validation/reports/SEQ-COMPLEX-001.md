# Validation Report: SEQ-COMPLEX-001 — Sequence Complexity Metrics

- **Validated:** 2026-06-24   **Area:** Sequence Composition
- **Canonical method(s):** `SequenceComplexity.CalculateLinguisticComplexity` (+ Shannon entropy,
  k-mer entropy, windowed complexity, low-complexity regions, DUST score, masking, compression ratio)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs`
Tests: `SequenceComplexityTests.cs`, `SequenceComplexity_CalculateDustScore_Tests.cs`,
`SequenceComplexity_CalculateKmerEntropy_Tests.cs`,
`SequenceComplexity_CalculateWindowedComplexity_Tests.cs`,
`SequenceComplexity_EstimateCompressionRatio_Tests.cs` (126 tests).

> Note on scope: since the 2026-06-12 archived report, the DUST and compression metrics were
> re-implemented and re-validated under their own sub-units (`SEQ-COMPLEX-DUST-001`,
> `SEQ-COMPLEX-COMPRESS-001`, `SEQ-COMPLEX-KMER-001`, `SEQ-COMPLEX-WINDOW-001`). The original
> `SEQ-COMPLEX-001` TestSpec markdown lagged behind on those two; this session re-validated the
> code against authoritative sources and **fixed the stale description** (no code change).

---

## Stage A — Description

### Sources opened (this session, external)
- **Wikipedia, "Linguistic sequence complexity"** — confirms vocabulary usage
  `U_i = actual / max`, `max = min(4^i, N−i+1)`, and presents the **product** aggregate
  `C = ∏ U_i`. Worked examples `ACGGGAAGCTGATTCCA` (U₂=14/16, U₃=15/15, U₄=14/14) and
  `ACACACACACACACACA` (U₁=1/2, U₂=2/16, U₃=2/15).
- **Troyanskaya et al. (2002), Bioinformatics 18(5):679–688** — confirms the **summation**
  aggregate `LC = A(s)/M(s) = Σ_i A_i / Σ_i M_i` computed via suffix trees. This is the form
  the spec selects and the code implements.
- **Li (2025), "Finding low-complexity DNA sequences with longdust", Bioinformatics
  42(3):btag112, §2.5** — authoritative restatement of the original Morgulis 2006 SDUST score.
  Exact formula (read from the PDF, line for line):
  `S_S(c⃗_x) = (1/ℓ(x)) · Σ_t c_x(t)(c_x(t)−1)/2 − T`, with `ℓ(x) = |x|−k+1` = **number of
  k-mers**. The normalization denominator is **ℓ(x) (the number of triplets), not (ℓ(x)−1)**.
- **Shannon (1948) / Wikipedia "Entropy"** — `H = −Σ p_i log₂ p_i`, max for 4-letter alphabet = 2.

### Formula check
| Measure | Spec formula | Authoritative source | Status |
|---|---|---|---|
| Linguistic complexity | `ΣA / ΣM`, `M_i = min(4^i, N−i+1)` | Troyanskaya 2002 (summation) | ✅ matches |
| Shannon entropy | `H = −Σ p_i log₂ p_i` | Shannon 1948 | ✅ matches |
| k-mer entropy | Shannon H over k-mer frequencies | Shannon 1948 / Li 2025 §2.6 | ✅ matches |
| DUST score | (was `Σ/(w−1)`) → **`Σ/ℓ(x)`** | Li 2025 §2.5 (Morgulis 2006) | ⚠️ **spec was wrong; fixed** |
| Compression | (was unique-substring 14/27) → **normalized LZ-76** | Lempel–Ziv 1976 / Zhang 2009 | ⚠️ **spec was stale; fixed** |

### Findings / divergences (Stage A)
1. **Product-vs-summation LC** (pre-existing, benign): two published LC definitions exist.
   The spec correctly selects and cites Troyanskaya's **summation** form, which the code
   implements. Recorded as a note, not a defect.
2. **DUST normalization in the spec was wrong.** The original `SEQ-COMPLEX-001.md` §3.4/§5.4
   stated `DUST = Σ/(w−1)` (denominator `N−3`), giving 8.0 / 6/13 / 2.5. The authoritative
   Li (2025) §2.5 restatement of Morgulis divides by `ℓ(x) = N−2` (the **number of triplets**),
   giving 7.5 / 6/14 / 2.0. **The code and the SEQ-COMPLEX-DUST-001 unit already use the
   correct `/ℓ(x)` divisor.** I corrected the stale TestSpec (§3.4, §4.6, §5.4) to the sourced
   `/ℓ(x)` convention.
3. **Compression-ratio section was superseded.** §4.8/§5.5 still described the old
   unique-substring heuristic (14/27, 5/112). `EstimateCompressionRatio` is now normalized
   Lempel–Ziv complexity, validated under SEQ-COMPLEX-COMPRESS-001. I updated §4.8/§5.5 to
   point at the LZ metric and its sourced worked values (2.0, 1.125, 1.25).

### Independent cross-check (hand / Python recomputation — all exact)
| Input | Measure | Recomputed | Code / current test expects |
|---|---|---|---|
| `ATGCTAGCATGCAATG` | LC (mw10) | 91/103 = 0.883495 | 91/103 ✅ |
| `AAAAAAAAAAAAAAAA` | LC (mw10) | 10/103 = 0.097087 | 10/103 ✅ |
| `A` / `ATGC` | LC | 1.0 / 1.0 | 1.0 / 1.0 ✅ |
| `ACGGGAAGCTGATTCCA` | LC (mw4) | 47/49 | 47/49 ✅ |
| `ACACACACACACACACA` | LC (mw10) | 20/112 = 5/28 | 5/28 ✅ |
| `ATGCATGCATGCATGC` | Shannon | 2.0 | 2.0 ✅ |
| `ATATATAT` / `ATGATGATG` | Shannon | 1.0 / log₂3 | 1.0 / log₂3 ✅ |
| `ATCG` | k-mer H (k2) | log₂3 = 1.5849625 | log₂3 ✅ |
| `ATGCATGCATGCATGC` | k-mer H (k2) | 1.9898981 | exact formula ✅ |
| `AAAAAAAAAAAAAAAAAA` | DUST | 120/16 = 7.5 | 7.5 ✅ (matches `/ℓ(x)`) |
| `ATGCTAGCATGCTAGC` | DUST | 6/14 = 3/7 | 6/14 ✅ |
| `AAAAAAA` | DUST | 10/5 = 2.0 | 2.0 ✅ |
| `1001111011000010` | norm. LZ | 8/(16/log₂16)=2.0 | 2.0 ✅ |
| `ACGTACGTACGTACGT` | norm. LZ | 9/(16/log₄16)=1.125 | 1.125 ✅ |

**Stage A verdict: PASS-WITH-NOTES** — every formula now matches its authoritative source
exactly. Two stale TestSpec sections (DUST divisor, compression metric) were corrected to the
already-implemented, sourced behaviour; the LC product-vs-summation distinction is a documented
benign note.

---

## Stage B — Implementation

### Code path reviewed
- **LC**: `CalculateLinguisticComplexityCore` (lines 39–66) — for word length 1..min(maxWord,N):
  distinct substrings via HashSet, `maxPossible = min(4^w, N−w+1)`; returns `Σobs/Σmax`.
  Realises Troyanskaya summation exactly.
- **Shannon**: 93–120 — frequencies over fixed `{A,T,G,C}` (non-ATGC excluded from both
  numerator and denominator), `−Σ p log₂ p`.
- **k-mer entropy**: 164–189 — H over overlapping k-mer counts; `len < k` → 0.
- **DUST**: `CalculateDustScoreCore` (368–401) — triplet counts, `Σ c(c−1)/2` (count promoted
  to `double` to avoid int overflow), normalized by `wordCount = N−wordSize+1 = ℓ(x)`.
  **Divisor = ℓ(x) ⇒ matches Li 2025 §2.5 / Morgulis.** `len < wordSize` → 0.
- **Lempel–Ziv** (`CalculateLempelZivComplexityCore` 528–552, normalized 554–576) — exhaustive-
  history parse; base clamped to 2 when <2 symbols. `EstimateCompressionRatio` delegates here.
- Windowed / low-complexity-region / masking reuse the Shannon and DUST cores.

### Formula realised correctly?
Yes. An independent Python reimplementation mirrors the C# control flow and produces identical
numbers to both the code and the test assertions for every row of the cross-check table.

### Cross-verification recomputed vs code
All 126 SequenceComplexity tests assert exact externally-sourced values (`Within(1e-10)` or
exact equality). DUST asserts 7.5 / 6/14 / 2.0 / 1.5 / 3.5 / 0.5 — all reproduced by the
`/ℓ(x)` divisor. LZ asserts 2.0 / 1.125 / 1.25 — all reproduced.

### Variant/delegate consistency
String overloads upper-case then call the same `*Core` as the `DnaSequence` overloads (parity
tests pass). `EstimateCompressionRatio(string|DnaSequence)` delegates to
`CalculateNormalizedLempelZivComplexity` (delegation tests pass). Windowed/region/mask reuse the
shared cores, so variants agree by construction.

### Numerical robustness
DUST promotes `count` to `double` before `count·(count−1)` (documented int-overflow guard for
long homopolymers). Div-by-zero guarded (`possibleTotal>0`, `wordCount>0` via length guard,
`total==0` early returns, `logBaseN<=0` degenerate guard). No overflow on stated ranges.

### Test quality audit
Assertions check exact sourced values (Troyanskaya / Wikipedia / Shannon / Li-2025-restated
Morgulis / Naereen LZ doctests), not tautologies. Guard clauses (null, k<1, maxWord<1,
windowSize<1, stepSize<1, wordSize<1) tested. Range invariants (LC∈[0,1], H∈[0,2], k-mer
H∈[0,log₂4^k]) tested. Deterministic.

### Findings / defects (Stage B)
None. The DUST divisor that the old TestSpec markdown described as `/(w−1)` is, in the code,
correctly `/ℓ(x)` — i.e. the **code was already right** and only the description was stale.
No code change required.

---

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES**, **Stage B: PASS**, **End state: CLEAN.**
- Code changed: **no**. Description changed: **yes** — `tests/TestSpecs/SEQ-COMPLEX-001.md`
  §3.4/§4.6/§4.8/§5.4/§5.5 corrected from the superseded DUST `/(w−1)` divisor and the obsolete
  unique-substring compression heuristic to the authoritative `/ℓ(x)` DUST normalization
  (Li 2025 §2.5) and the normalized Lempel–Ziv compression metric, matching the
  already-validated SEQ-COMPLEX-DUST-001 / SEQ-COMPLEX-COMPRESS-001 sub-units.
- Full unfiltered suite after edits: **18208 passed, 0 failed**.
