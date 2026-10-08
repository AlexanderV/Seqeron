# Evidence Artifact: SEQ-COMPLEX-COMPRESS-001

**Test Unit ID:** SEQ-COMPLEX-COMPRESS-001
**Algorithm:** Lempel–Ziv complexity (compression-based sequence complexity)
**Date Collected:** 2026-06-14 (corrected 2026-09-28, review-2026-09 B04)

---

## ⚠ Correction 2026-09-28 — LZ78 parse mislabelled as LZ76

The 2026-06 version of this artifact took the Naereen `lempel_ziv_complexity` set-based parse
(values 8/7/9/10, homopolymer 5) as the Lempel–Ziv (1976) measure. That parse is the Ziv–Lempel
**1978** incremental parsing (each phrase = a previously seen phrase + one symbol), not the LZ76
exhaustive history. Sources actually opened on 2026-09-28:

- **antropy 0.2.2** (PyPI wheel, `antropy/entropy.py`, `_lz_complexity` + `lziv_complexity`): the
  Kaspar–Schuster (1987) scan, citing Lempel & Ziv 1976 and Zhang 2009. Doctests: `1001111011000010`
  → **6** (`1 / 0 / 01 / 1110 / 1100 / 0010`), normalized **1.5**; `HELLO WORLD! ×4` → 11 /
  0.38596001132145313; `A..Z` → 26 / 1.0. Normalization `c/(n/log(n, base))`, `base = 2 if base < 2`.
- **lempel_ziv_complexity 0.2.2** (Naereen, PyPI wheel): set-based parse, `1001111011000010` → 8 — LZ78-style.
- **NeuroKit2** `complexity_lempelziv.py` (raw.githubusercontent.com): same Kaspar–Schuster scan,
  citing Lempel & Ziv 1976 and Kaspar & Schuster 1987.
- **WebSearch snippets**: Lempel & Ziv 1976 classic example `0001101001000101` has complexity 6
  (`0·001·10·100·1000·101`); Estévez-Rams et al. arXiv:1311.0546: "an exhaustive history of the
  sequence is one where any factor is not a substring of the string preceding it, while the factor
  with the last symbol removed is … E(u)=0.1.00.11.101.101100, C(u)=6" for u = 010011101101100.
  Wikipedia (snippet) gives `1001111011000010` = 1/0/01/1110/1100/0010, complexity 6.
- **Numerical cross-check:** a brute-force Python implementation of the exhaustive-history
  definition agreed with antropy `_lz_complexity` on 20 000 random strings (n = 1–60, alphabets
  {A}, {A,T}, {A,C,G}, {A,C,G,T}): 0 mismatches.

Corrected reference values (antropy 0.2.2 = brute-force definition):

| Input | c (LZ76) | normalized | (old LZ78 value) |
|---|---|---|---|
| `1001111011000010` | 6 | 1.5 | 8 / 2.0 |
| `0001101001000101` | 6 | 1.5 | — |
| `010011101101100` | 6 | 1.5627562382434075 | — |
| `1010101010101010` | 3 | 0.75 | 7 |
| `1001111011000010000010` | 7 | 1.4189100604755036 | 9 |
| `100111101100001000001010` | 8 | 1.5283208335737188 | 10 |
| `0`×16 | 2 | 0.5 | 5 / 1.25 |
| `AAAA` | 2 | 1.0 | 2 |
| `ACGT` | 4 | 1.0 | 4 |
| `ACGTACGTACGTACGT` | 5 | 0.625 | 9 / 1.125 |
| `ATGCGATCGATTAGCCGATAGCTAGGCTAACGTTAGCATG` | 16 | 1.0643856189774725 | — |
| `CAG`×20 | 4 | 0.2484555351907228 | — |
| 20 000-base LCG DNA (see test R2) | 2756 | 0.984423382950957 | — |

For n = 1 antropy raises ZeroDivisionError (`log(1, b) = 0`); this library returns the raw count 1.
Everything below this section is the original 2026-06 record; where it conflicts with this
section (Naereen values, "set-based contract", homopolymer 5, ASSUMPTION 1), this section wins.

---

## Online Sources

### Lempel, A. & Ziv, J. (1976) — "On the Complexity of Finite Sequences" (primary, citation index)

**URL:** https://www.mindat.org/reference.php?id=14583094 (bibliographic record); DOI: https://doi.org/10.1109/TIT.1976.1055501
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed, IEEE Trans. Inf. Theory 22(1):75–81)

**Key Extracted Points:**

1. **Origin / definition basis:** Retrieved via WebSearch query *"Lempel Ziv 1976 On the Complexity of Finite Sequences complexity measure c(S) definition"*. The result records the publication (IEEE Trans. on Information Theory, 22(1), 75–81, 1976, doi:10.1109/tit.1976.1055501) and states the measure links complexity "to the gradual buildup of new patterns along the given sequence," i.e. **the number of distinct phrases produced as the sequence is parsed left-to-right**. (The full paper text behind a paywall; the precise parsing rule and worked numbers are taken from the rank-3/4 sources below, all of which cite this paper.)

### Wikipedia — "Lempel–Ziv complexity" (cites the 1976 primary)

**URL:** https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the primary Lempel & Ziv 1976)

**Key Extracted Points:**

1. **Definition (verbatim):** the complexity is "the number of different sub-strings (or sub-words) encountered as the binary sequence is viewed as a stream (from left to right)." "The Lempel–Ziv complexity corresponds to the number of iterations needed to finish this procedure."
2. **Parsing rule (verbatim):** "We have to move the delimiter (starting in position 1) the furthest possible to the right, so that the sub-word between position 1 and the delimiter position be a word of the sequence that starts before the position 1 of the delimiter. As soon as the delimiter is set on a position where this condition is not met, we stop, move the delimiter to this position, and start again by marking this position as a new initial position." (= shortest new substring not already encountered in the scanned history → start a new component.)
3. **Algorithm (verbatim pseudocode):** an O(n) left-to-right scan returning component count `C`:
   ```
   i := 0; C := 1; u := 1; v := 1; vmax := v
   while u + v <= n do
      if S[i + v] ≠ S[u + v] then v := v + 1
      else
         vmax := max(v, vmax); i := i + 1
         if i ≠ u then C := C + 1; u := u + vmax; v := 1; i := 0; vmax := v
         else v := 1
      end if
   end while
   if v ≠ 1 then C := C + 1
   ```

### Naereen / Lempel-Ziv_Complexity — reference implementation (Python, MIT)

**URL:** https://raw.githubusercontent.com/Naereen/Lempel-Ziv_Complexity/master/src/lempel_ziv_complexity.py
**Accessed:** 2026-06-14
**Authority rank:** 3 (well-maintained open-source reference implementation; links to the Wikipedia/primary)

**Key Extracted Points:**

1. **Definition (verbatim docstring):** "It is defined as the number of different substrings encountered as the stream is viewed from begining to the end."
2. **Reference algorithm (verbatim):** parse with a `set` of seen substrings; at index `ind` take `sub = sequence[ind:ind+inc]`; if `sub` already in the set grow `inc += 1`, else add `sub`, advance `ind += inc`, reset `inc = 1`; return `len(set)`. (Exhaustive-history LZ76 factorization.)
3. **Worked examples (verbatim doctests with exact return values):**
   - `lempel_ziv_complexity('1001111011000010')` → **8** ; components `1 / 0 / 01 / 11 / 10 / 110 / 00 / 010`
   - `lempel_ziv_complexity('1010101010101010')` → **7** ; components `1 / 0 / 10 / 101 / 01 / 010 / 1010`
   - `lempel_ziv_complexity('1001111011000010000010')` → **9** ; components `1 / 0 / 01 / 11 / 10 / 110 / 00 / 010 / 000`
   - `lempel_ziv_complexity('100111101100001000001010')` → **10** ; components `… / 0101`

### AntroPy / `antropy.lziv_complexity` — reference implementation (cites Lempel-Ziv 1976 & Zhang 2009)

**URL:** https://raphaelvallat.com/antropy/build/html/generated/antropy.lziv_complexity.html (HTTP 404 at fetch time); definition/normalization corroborated by the sibling `entropy.lziv_complexity` doc below.
**Accessed:** 2026-06-14
**Authority rank:** 3

**Key Extracted Points:**

1. Surfaced via WebSearch as the canonical Python LZ-complexity reference; confirms the "number of different substrings encountered" definition and a `normalize` option dividing by `n / log_b(n)`.

### `entropy.lziv_complexity` documentation (cites Lempel-Ziv 1976 & Zhang et al. 2009)

**URL:** https://raphaelvallat.com/entropy/build/html/generated/entropy.lziv_complexity.html
**Accessed:** 2026-06-14
**Authority rank:** 3

**Key Extracted Points:**

1. **Definition (verbatim):** "the number of different substrings encountered as the sequence is viewed from beginning to the end."
2. **Normalization formula (verbatim):** with `normalize=True`, `LZn = LZ / (n / log_b(n))` where `LZ` = raw complexity count, `n` = sequence length, `b` = number of unique characters (alphabet size). Rationale: "raw LZ is heavily influenced by sequence length (longer sequence will result in higher LZ)."
3. **References given:** Lempel, A., & Ziv, J. (1976), IEEE Trans. Inf. Theory 22(1):75–81; Zhang, Y. et al. (2009), J. Math. Chem. 46(4):1203–1212.

### Asymptotic upper bound b(n) (WebSearch corroboration)

**URL:** WebSearch query *"normalized Lempel-Ziv complexity formula C(n) divide n / log_alpha(n) upper bound b(n) Lempel Ziv 1976"*
**Accessed:** 2026-06-14
**Authority rank:** 4 (search synthesis of primary-citing sources: arXiv nlin/0608049, arXiv 1311.0546, Zhang 2009)

**Key Extracted Points:**

1. **Upper bound (verbatim from results):** "for uniformly distributed symbols … b(n) = n/log n, and c(n) is normalized to b(n) resulting in γ = c(n)/b(n)." Also "C(u) < N/((1−ε_N) log_σ N) … leading to an asymptotic value C(u) < N/log N." Here the log base σ is the alphabet size; γ = c(n)/b(n) → 1 for a maximally complex (random) sequence.

---

## Documented Corner Cases and Failure Modes

### From Naereen reference implementation

1. **Empty input:** `len(set)` over an empty scan = 0 components → complexity 0.
2. **Single symbol / homopolymer:** a run of one symbol grows one extra-long component each step, so the set is `{0, 00, 000, …}`; for `n` identical symbols the component count is the largest `k` with `1+2+…+k ≤ n`, i.e. `c = ⌊(√(8n+1)−1)/2⌋`. E.g. `"0"×16` → components `{0,00,000,0000,00000}` → **c = 5** (verified by tracing the reference parser). `"AAAA"` → `{A, AA}` → **c = 2**. This is much lower than a random string of the same length (low complexity), but is NOT 1.
3. **Trailing incomplete component:** when the final substring extends to the end without ever becoming "new" on its own, the standard pseudocode adds 1 (`if v ≠ 1 then C := C+1`); the Naereen set-based variant does NOT count a trailing partial substring that never completed (it only counts substrings actually added to the set). The Naereen doctest values are the contract used here.

### From entropy/antropy

1. **Length dependence:** raw LZ grows with `n`; normalization by `n/log_b(n)` is required for cross-length comparison.
2. **Single distinct symbol (b=1):** `log_b(n)` with b=1 is undefined (log base 1). The entropy/antropy reference handles this by clamping the base to 2 (`base = 2 if base < 2 else base`) and returning the normalized value `c/(n/log_2 n)` — NOT the raw count. (Corrected 2026-06-16 from the earlier "return raw count" reading.)

---

## Test Datasets

### Dataset: Naereen reference doctests (binary alphabet) — *superseded 2026-09-28: LZ78 counts, not LZ76; use the correction table at the top*

**Source:** Naereen/Lempel-Ziv_Complexity `src/lempel_ziv_complexity.py` (retrieved 2026-06-14)

| Input string | Components (LZ76 exhaustive history) | Raw complexity c |
|--------------|--------------------------------------|------------------|
| `1001111011000010` | `1/0/01/11/10/110/00/010` | 8 |
| `1010101010101010` | `1/0/10/101/01/010/1010` | 7 |
| `1001111011000010000010` | `1/0/01/11/10/110/00/010/000` | 9 |
| `100111101100001000001010` | `…/000/0101` | 10 |

### Dataset: Normalized LZ (derived from entropy formula `LZn = c / (n / log_b n)`) — *superseded 2026-09-28: built on the LZ78 counts; correct values in the table at the top*

**Source:** entropy/antropy `lziv_complexity` normalization (retrieved 2026-06-14); derivation shown.

| Input string | n | b (alphabet) | c | log_b(n) | b(n)=n/log_b(n) | normalized = c/b(n) |
|--------------|---|--------------|---|----------|-----------------|---------------------|
| `1001111011000010` | 16 | 2 | 8 | log₂16 = 4 | 16/4 = 4 | 8/4 = **2.0** |
| `0000000000000000` | 16 | 1→2 (clamped) | 5 | log₂16 = 4 | 16/4 = 4 | 5/4 = **1.25** |

**Correction (2026-06-16):** the entropy/antropy reference does NOT return the raw count for b<2.
Its source (`antropy/src/antropy/entropy.py`, function `lziv_complexity`) reads
`base = 2 if base < 2 else base; return _lz_complexity(s) / (n / log(n, base))` — i.e. it clamps the
log base to 2. For `"0"×16` this gives `5/(16/log₂16) = 5/4 = 1.25`, not 5.0. The implementation and
test M8 were corrected to 1.25 during validation of SEQ-COMPLEX-COMPRESS-001.

Additional traced raw values (reference parser, retrieved 2026-06-14): `"AAAA"` → c=2 (`A/AA`); `"ACGT"` → c=4 (`A/C/G/T`); `""` → c=0.

---

## Assumptions

1. **ASSUMPTION: trailing-component convention** — The exhaustive-history factorization can end with a final partial substring that never became "new" before the sequence ended. The Wikipedia pseudocode adds 1 for it (`if v≠1`); the Naereen set-based reference does not. We adopt the **Naereen set-based contract** because it has explicit, reproducible doctest expected values that we can encode exactly; the chosen convention is consistent across all four Naereen worked examples. This affects at most the last component count by 1 and is documented, not invented.
2. **ASSUMPTION: normalization log base = alphabet size (distinct symbols actually present, b≥2)** — entropy/antropy use `b` = "number of unique characters" in the sequence. For DNA this is ≤4. When b<2 (single-symbol input) `log_b n` is undefined; we return the raw count in that degenerate case. Source-backed by the entropy doc formula; the b<2 fallback is the documented degenerate handling.

---

## Recommendations for Test Coverage

*(2026-06 recommendations; the expected values cited below are the superseded LZ78 ones — current tests lock the LZ76 values of the correction table, see `tests/TestSpecs/SEQ-COMPLEX-COMPRESS-001.md`.)*

1. **MUST Test:** raw LZ76 complexity equals the four Naereen doctest values exactly (8, 7, 9, 10). — Evidence: Naereen `lempel_ziv_complexity.py` doctests.
2. **MUST Test:** normalized LZ for `1001111011000010` = 2.0 via `c/(n/log_b n)`. — Evidence: entropy normalization formula + derivation.
3. **MUST Test:** homopolymer `"0"×16` yields raw complexity 5 (components `0/00/000/0000/00000`) and the b<2 normalization fallback returns the raw count (5). — Evidence: traced reference parser; entropy b<2 undefined-log rule.
4. **SHOULD Test:** maximally complex short string `"ACGT"` (all distinct symbols) → each symbol its own component → c=4. — Rationale: boundary of the parsing rule.
5. **SHOULD Test:** monotonic-style invariant — a more repetitive sequence has ≤ complexity of a less repetitive one of equal length. — Rationale: stated property of the measure (productivity buildup).
6. **COULD Test:** DNA alphabet sequence (b=4) normalization uses log base 4. — Rationale: real bio-sequence application (Zhang 2009).

---

## References

1. Lempel, A., & Ziv, J. (1976). On the Complexity of Finite Sequences. IEEE Transactions on Information Theory, 22(1), 75–81. https://doi.org/10.1109/TIT.1976.1055501
2. Wikipedia contributors. Lempel–Ziv complexity. https://en.wikipedia.org/wiki/Lempel%E2%80%93Ziv_complexity (accessed 2026-06-14)
3. Besson, L. (Naereen). Lempel-Ziv_Complexity (Python reference implementation). https://github.com/Naereen/Lempel-Ziv_Complexity/blob/master/src/lempel_ziv_complexity.py (accessed 2026-06-14)
4. Vallat, R. entropy / AntroPy — `lziv_complexity`. https://raphaelvallat.com/entropy/build/html/generated/entropy.lziv_complexity.html (accessed 2026-06-14)
5. Zhang, Y., Hao, J., Zhou, C., & Chang, K. (2009). Normalized Lempel-Ziv complexity and its application in bio-sequence analysis. Journal of Mathematical Chemistry, 46(4), 1203–1212. https://doi.org/10.1007/s10910-008-9512-2 (abstract paywalled; normalization formula taken from ref. 4 which cites it)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Correction — LZ78 set parse replaced by LZ76 exhaustive history (antropy / Kaspar–Schuster); reference values recomputed (see top section).
