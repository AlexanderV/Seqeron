# Validation Report: SEQ-COMPLEX-COMPRESS-001 — Compression-based sequence complexity (Lempel–Ziv 1976)

- **Validated:** 2026-06-16   **Area:** Complexity
- **Canonical method(s):** `SequenceComplexity.CalculateLempelZivComplexity(string|DnaSequence)`,
  `CalculateNormalizedLempelZivComplexity(string|DnaSequence)`,
  `EstimateCompressionRatio(string|DnaSequence)` (delegate → normalized LZ)
- **Stage A verdict:** FAIL → corrected (2026-09: the "LZ76" description was the Naereen LZ78 incremental parse)
- **Stage B verdict:** FAIL → fixed (review campaign 2026-09, batch B04: F1) — supersedes the 2026-06 verdicts and values below

## Review campaign 2026-09 revision (B04) — supersedes the 2026-06 report below

- **End-state: FIXED** (batch report `docs/Validation/review-2026-09/B04.md`, F1).
- **F1 (defect):** `CalculateLempelZivComplexity` used the Naereen "set of seen phrases" parse — Ziv–Lempel **1978** incremental parsing — while documented as LZ76. Every LZ value in the 2026-06 report below (doctests 8/7/9/10, `"0"×16` → 5 / 1.25, `1001111011000010` → 2.0, `ACGT×4` → 9 / 1.125) is the LZ78 count and is **wrong**. Correct LZ76 (exhaustive history): `1001111011000010` → 6 (1/0/01/1110/1100/0010), normalized 1.5; `0001101001000101` → 6 (Lempel & Ziv 1976); `010011101101100` → 6 (Estévez-Rams et al.); `1010…`×16 → 3; `"0"×16` → 2, normalized 0.5 (b < 2 clamp to 2 kept); `ACGT×4` → 5, normalized 0.625; `HELLO WORLD!`×4 → 11 / 0.38596001132145313; A..Z → 26 / 1.0.
- **Sources opened:** antropy 0.2.2 `entropy.py` (`_lz_complexity`, `lziv_complexity` doctests, PyPI wheel); NeuroKit2 `complexity_lempelziv.py`; Naereen `lempel_ziv_complexity` 0.2.2 (identified as LZ78); Lempel & Ziv 1976 example, Estévez-Rams et al. arXiv:1311.0546, Wikipedia LZ complexity (snippets).
- **Reference cross-check:** brute-force exhaustive-history definition = antropy on 20 000 random strings (0 mismatches); 6 DNA strings + 20 kb DNA (c = 2756, normalized 0.984423382950957) = antropy to 1e-12; independent Kaspar–Schuster port agrees.
- **Fix:** LZ76 from the Longest-Previous-Factor array (Crochemore & Ilie 2008; shared suffix-array + Kasai LCP helpers), O(n log² n). Normalization c/(n/log_b n) (Zhang 2009) unchanged. `EstimateCompressionRatio`, MCP `compression_ratio`, `complexity_compression_ratio` inherit the fix; raw + normalized LZ76 also exposed as MCP `lempel_ziv_complexity` (F49).
- **Tests:** `SequenceComplexity_EstimateCompressionRatio_Tests.cs` re-sourced (M1–M9, LZ1976, Estévez-Rams, antropy text and DNA datasets, 20 kb); fuzz reference = brute-force LZ76; MCP `CompressionRatioTests` → 1.5 / 0.625. TestSpec `tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md`, Evidence, `docs/algorithms/Complexity/Lempel_Ziv_Complexity.md` corrected.
- **Current code path:** `SequenceComplexity.cs:1712/1581` (raw), :1599/1608 (normalized), :1621/1631 (`EstimateCompressionRatio`), core :1636, `ComputeLongestPreviousFactor` :1667.

*Historical 2026-06-16 report (superseded — its LZ values are LZ78, see above):*

## Stage A — Description

### Sources opened this session (retrieved, not trusted from the label)

1. **Naereen/Lempel-Ziv_Complexity** `lempel_ziv_complexity.py` — WebFetched raw source.
   Docstring: "the number of different substrings encountered as the stream is viewed from
   begining to the end." Algorithm body matches the implementation's core byte-for-byte.
   Doctests: `1001111011000010`→**8**, `1010101010101010`→**7**, `1001111011000010000010`→**9**,
   `100111101100001000001010`→**10**.
2. **Wikipedia — Lempel–Ziv complexity** — WebFetched. Definition "number of different
   sub-strings encountered as the binary sequence is viewed as a stream"; O(n) pseudocode.
3. **antropy `lziv_complexity`** (`antropy/src/antropy/entropy.py`) — WebFetched raw source.
   Normalization block (verbatim):
   `base = sum(np.bincount(s) > 0); base = 2 if base < 2 else base; return _lz_complexity(s) / (n / log(n, base))`.
4. **entropy `lziv_complexity`** (older sibling) — WebFetched raw source. Identical
   normalization block (clamps base to 2). Both cite Lempel–Ziv 1976 and Zhang et al. 2009.
5. Lempel & Ziv (1976) IEEE TIT 22(1):75–81 — primary, paywalled; the parsing convention and
   worked numbers are taken from sources #1–#4 which cite it.

### Formula check

- **Raw LZ (c):** exhaustive-history left-to-right set-based parse → number of distinct
  components. Matches Naereen #1 exactly. ✅
- **Normalized LZ:** `LZn = c / (n / log_b n)`, b = number of distinct symbols. Matches antropy
  #3 / entropy #4 exactly for b ≥ 2. ✅
- **b < 2 (single-symbol) case:** the reference **clamps base to 2** and returns
  `c/(n/log_2 n)`. The TestSpec/Evidence claimed "returns the raw count" and attributed it to
  source #4 — **this is wrong**; source #4 does the opposite. Stage-A divergence (corrected).

### Edge-case semantics

- Empty → 0 components (Naereen `len(set)`=0). ✅
- Homopolymer `"0"×16` → `0/00/000/0000/00000` = 5 (traced). ✅
- Single base → 1 component. ✅
- b<2 normalization → clamp base to 2 (corrected from "raw count").

### Independent cross-check (numbers I computed this session)

Independent Python reimplementation of the Naereen set parser (`/tmp/lz_verify.py`):

| Input | My c | Spec c | Match |
|-------|------|--------|-------|
| `1001111011000010` | 8 | 8 | ✅ |
| `1010101010101010` | 7 | 7 | ✅ |
| `1001111011000010000010` | 9 | 9 | ✅ |
| `100111101100001000001010` | 10 | 10 | ✅ |
| `0000000000000000` | 5 | 5 | ✅ |
| `ACGT` | 4 | 4 | ✅ |
| `AAAA` | 2 | 2 | ✅ |
| `ACGTACGTACGTACGT` | 9 | 9 (for C2) | ✅ |

Normalized (my computation): `1001111011000010`→**2.0**; `ACGTACGTACGTACGT`→**1.125**;
`"0"×16` antropy-style clamp→**1.25** (NOT 5.0).

### Stage A findings

- **PASS-WITH-NOTES.** Raw-LZ description and all b≥2 normalization values are correct and
  externally confirmed. One divergence: the b<2 normalization convention in the
  TestSpec/Evidence ("return raw count", value 5.0) contradicts its own cited reference, which
  clamps the base to 2 (value 1.25). Description corrected in TestSpec §1.3/§4 (M8)/§6 (A2) and
  Evidence (normalized table + corner cases).

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceComplexity.cs:601-573` *(2026-06 line numbers and parser — superseded; current: public overloads :1571–1631, LZ76 core `CalculateLempelZivComplexityCore` :1636, `ComputeLongestPreviousFactor` :1667, normalization :1795)*.
- `CalculateLempelZivComplexityCore` (522): set-based exhaustive-history parser — identical to
  the Naereen reference. ✅
- `CalculateNormalizedLempelZivComplexityCore` (548): computed `c/(n/log_b n)` for b≥2 (✅) but
  **returned the raw count for b<2** — the defect.
- `EstimateCompressionRatio` (507/517): thin delegate to normalized LZ. ✅

### Defect found & fixed

**b<2 normalization diverged from the cited reference.** Code returned `c` (raw) when fewer than
two distinct symbols were present; the antropy/entropy reference clamps the base to 2 and returns
the normalized value. Fix: replaced `if (b < MinAlphabetForNormalization) return c;` with
`if (b < MinAlphabetForNormalization) b = MinAlphabetForNormalization;` (clamp), keeping the
`n==1 ⇒ log_b(1)=0` div-by-zero guard. For `"0"×16` the result is now **1.25** (= 5/(16/log₂16)),
matching the reference.

### Cross-verification table recomputed vs the fixed code

| Method | Input | Expected (sourced) | Code (after fix) |
|--------|-------|--------------------|------------------|
| raw LZ | doctests 1–4 | 8/7/9/10 (Naereen) | 8/7/9/10 ✅ |
| raw LZ | `"0"×16` / `ACGT` / `AAAA` / `A` / `""` | 5/4/2/1/0 (traced) | 5/4/2/1/0 ✅ |
| normalized | `1001111011000010` | 2.0 (antropy) | 2.0 ✅ |
| normalized | `"0"×16` (b<2 clamp) | **1.25** (antropy code) | 1.25 ✅ |
| normalized | `ACGTACGTACGTACGT` (b=4) | 1.125 (antropy) | 1.125 ✅ |
| EstimateCompressionRatio | `1001111011000010` / DnaSeq ACGT×4 | 2.0 / 1.125 | 2.0 / 1.125 ✅ |

### Variant/delegate consistency

String and DnaSequence overloads agree (raw + normalized); `EstimateCompressionRatio` equals
`CalculateNormalizedLempelZivComplexity` for both overloads — now exercised by tests.

### Test quality audit (HARD gate)

- **Sourced expectations, not code echoes:** raw values come from the Naereen doctests; normalized
  from the antropy code; the b<2 value was corrected from the unsourced 5.0 to the sourced 1.25. A
  deliberately-wrong parser (entropy-convention 6 for doctest 1) fails M1.
- **No green-washing:** M8 was *not* adjusted to match code output — the **code** was fixed and the
  test locked to the sourced 1.25. No assertion weakened, no tolerance widened, no skip.
- **Cover all logic:** added 5 tests closing untested overloads/branches — normalized DnaSequence
  parity (1.125), normalized null-DnaSequence throws, normalized empty→0, normalized n==1 degenerate
  guard (raw count 1.0), `EstimateCompressionRatio(DnaSequence)` delegation (1.125). Fixture 15→20.
- **Honest green:** full unfiltered `Seqeron.Genomics.Tests` = **Failed: 0, Passed: 6603**; build
  **0 errors**; changed files warning-free (the 4 pre-existing NUnit2007 warnings are in an unrelated
  file). **Gate: PASS.**
- **MCP note:** `Seqeron.Mcp.Sequence.Tests/ComplexityCompressionRatioTests.cs` targets net9.0 and
  cannot be built by the installed net8 SDK (pre-existing environmental limit). It is a binding/
  ordering smoke test, unchanged, and uses b=4 inputs unaffected by the b<2 fix.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (raw-LZ + b≥2 normalization correct; b<2 convention mis-sourced →
  description corrected).
- **Stage B:** FAIL → FIXED (b<2 clamp-to-2 corrected to match the cited reference; tests locked to
  sourced values; coverage gaps closed).
- **End-state:** ✅ CLEAN — defect completely fixed in-session; build + full suite green.
- **Logged:** FINDINGS_REGISTER A29.
