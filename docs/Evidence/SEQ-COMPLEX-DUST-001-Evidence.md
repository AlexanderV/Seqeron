# Evidence Artifact: SEQ-COMPLEX-DUST-001

**Test Unit ID:** SEQ-COMPLEX-DUST-001
**Algorithm:** DUST Score (triplet-frequency low-complexity score of Morgulis et al. 2006 SDUST/DUST)
**Date Collected:** 2026-06-14 (revised 2026-09-28 and 2026-09-30, review campaign B04)

---

## 2026-09-30 Revision (B04 completeness audit WP4) — triplets only; longdust; linker; N input

**Sources opened (2026-09-30):**
- lh3/sdust `sdust.c` (git clone of github.com/lh3/sdust, compiled with `make`): `#define SD_WLEN 3`;
  line 68 `TODO: is this right for SD_WLEN!=3?`; `sdust_core` comment "N or the end of sequence; N effectively
  breaks input into pieces of independent sequences" — but the `else` branch resets only `l` and `t`, not the
  triplet window `w`, `L`, `rw`, `rv`, `cw`, `cv`.
- NCBI C++ Toolkit `src/algo/dustmask/symdust.cpp` / `include/algo/dustmask/symdust.hpp`
  (`triplet_type`; `DEFAULT_LINKER = 1`; constructor clamps linker to 1–32; `save_masked_regions`:
  `if (s + linker_ >= b1.first) res.back().second = max(s, b1.second)`), `src/app/dustmask/dust_mask_app.cpp`
  (`GetDustMasks_SkipNs`: only runs of N ≥ window are cut out). Binary: Debian `ncbi-blast+` 2.12.0 `dustmasker`.
- lh3/longdust 1.4-r97 (git clone, commit 9491215, MIT): `README.md` ("[SDUST] hardcodes k=3";
  `S_L(x) = Σ log c_x(t)! − f(ℓ(x)/4^k)`), `tex/longdust.tex` (Li H, Li B: `Q = Σ log c! − f(ℓ)`,
  `S = Q − T·ℓ`, T = 0.6; `f(ℓ;q) ≈ Σ_t e^{−ℓq_t} Σ_n log n!·(ℓq_t)^n/n!`), `longdust.c` (`ld_opt_init`,
  `ld_cal_f`, `ld_cal_f2`, `ld_dust1`, `ld_dust2`, `ld_extend`, `ld_if_backward`, `ld_backward`, `ld_forward`),
  compiled with `make` (and an AddressSanitizer build).

**F34 — `wordSize ≠ 3` rejected; longdust added.** Every source defines DUST on triplets only (above), so the
former extrapolation (divisor L − wordSize) is removed; the sourced k-mer generalisation is longdust.
Cross-check of `FindLongdustRegions` vs the compiled `longdust` binary: 1 500 random repeat-rich inputs
(5 bp–12 kb, 4 507 434 bp, 410 containing N/IUPAC; k ∈ {3..8}, w ∈ {50, 100, 300, 1000, 2000, 5000},
T ∈ {0.3, 0.45, 0.6, 0.8, 1.0}, `-e` 0/10/50/200, `-b` 2/3, `-f`, `-a`, `-g` 0.3/0.41/0.6) → 12 818 reference
intervals, **0 mismatches**; 121 homopolymer edge cases around w + k − 1 plus one 2 Mb sequence
(`longdust_asan`, 0 ASan errors) → **0 mismatches**. `CalculateLongdustScore` vs a C helper linking
`longdust.c` (`ld_cal_f`/`ld_cal_f2` + the same Σ log c accumulation): 3 000 inputs (k 1–10, GC off/0.3/0.41)
→ **3 000 bit-identical**; vs an independent Python lgamma/Poisson evaluation (1 057 uniform cases)
→ max relative difference 8.7·10⁻⁸.

**F35 — dustmasker linker + soft mask.** C# `FindLowComplexityIntervals(…, linker)` / `MaskLowComplexity(…,
softMask: true)` vs `dustmasker -window W -level T -linker L -outfmt interval|fasta` on 1 500 ACGT inputs
(10–500 bp; W ∈ {8, 16, 30, 64}; level ∈ {2, 10, 15, 20, 30}; linker ∈ {1, 2, 3, 5, 10, 20, 32}; the linker
changed the output in 261 cases): **0 interval mismatches, 0 soft-mask mismatches**. Second run (1 000
inputs): symdust's rule applied to dustmasker's own linker-1 output reproduces dustmasker(L) in 1 000/1 000;
C#(L) = symdust rule applied to sdust in 1 000/1 000; dustmasker's core differs from sdust in 10 cases,
all W = 8 / level 30 (symdust `thresholds_` has only W − 3 entries), so C#(L) = dustmasker(L) in 990/1 000.
Worked row (dustmasker 2.12.0, `-window 64 -level 20`, 96-bp A×16 / (CA)×8 / A×13 sequence): linker 1–17 →
`10-25, 43-58, 81-93`; 18–19 → `10-58, 81-93`; 32 → `10-93` (closed).

**F36 — string overload with N.** Upstream sdust on
`ACGTNNAAAAAAAAAAAANACGTACACACACACACACANNGGGCCCTAGGTCA` (53 bp, `-w 64 -t 20`) prints `8 20 / 21 34 /
35 72` — an interval ending past the sequence end, caused by the window leaking across the N. The C# string
path resets the window at every non-ACGT symbol, i.e. it equals sdust run on each maximal ACGT piece
(`[6,18) [23,38)`). 3 000 random inputs (1 037 with N/IUPAC/lower case; W 3–100, T 1–30): C# vs per-piece
sdust **0 mismatches**; ACGT-only inputs vs whole-input sdust **0 mismatches**; upstream whole-input sdust
differs from its own per-piece output on 274 of the N inputs. dustmasker is not the reference on N input
(it treats IUPAC codes as bases and only cuts N runs ≥ window).

---

## 2026-09-28 Revision (review campaign B04) — normaliser corrected to ℓ − 1; SDUST masking

The 2026-06 reading below of the longdust restatement (`/(L−2)`, i.e. divide by the number of
triplets ℓ) is **superseded**. Both reference implementations normalise by **ℓ − 1**:

- **lh3/sdust `sdust.c`** (opened: https://raw.githubusercontent.com/lh3/sdust/master/sdust.c, 2026-09-28),
  `find_perfect`: `new_r = r, new_l = kdq_size(w) - i - 1; if (new_r * 10 > T * new_l)` — the interval
  from triplet i to the window end contains `kdq_size(w) − i` = ℓ triplets, so `new_l = ℓ − 1`.
  (`rw * 10 > L * T` in `sdust_core` is only a pre-filter deciding whether to search.)
- **NCBI C++ Toolkit dustmasker** (opened: `src/algo/dustmask/symdust.cpp`,
  `include/algo/dustmask/symdust.hpp` on raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public, 2026-09-28):
  `thresholds_[0] = 1; thresholds_[i] = i*level_`, test `score*10 > thresholds_[count]` — the
  threshold scales with ℓ − 1.
- **Morgulis et al. 2006** (WebSearch snippet of the JCB PDF, 2026-09-28): "ℓ = n − 2 is the number of
  triplets"; score normalised by (ℓ − 1). Full text not openable (publisher/mirror blocked).
- **Numerical confirmation** — sdust compiled with gcc (`make`, sources above), default `-w 64 -t 20`:

| Input | ℓ | Σ c(c−1)/2 | Σ/(ℓ−1) | Σ/ℓ | sdust output |
|---|---|---|---|---|---|
| `ACGTTGCAGTCATGCGATC`+A×7+`TGCATCGGATCCTAGGCTAA` (A-run 7) | 5 | 10 | 2.5 | 2.0 | masked [19,26) |
| same with A-run 6 | 4 | 6 | 2.0 | 1.5 | none |
| `GGTTCATGGCTT`+(AC)×6+`GGTTCATGGCTT` | 10 | 20 | 2.22 | 2.0 | masked [12,24) |
| `GGTTCATGGCTT`+`ACACACACACA`+`GGTTCATGGCTT` | 9 | 16 | 2.0 | 1.78 | none |
| `AAAAAAA` (standalone) | 5 | 10 | 2.5 | 2.0 | masked [0,7) |

Only the ℓ − 1 normaliser reproduces the reference (score > 2.0 ⇔ masked). The longdust README
(opened: https://raw.githubusercontent.com/lh3/longdust/master/README.md) writes `S_D(x) = Σ c(c−1)/2 / ℓ(x)`;
this informal restatement does not match the SDUST code and is not used for SDUST.

**SDUST masking cross-check:** the C# port of `sdust_core` (`MaskLowComplexity`) was compared with the
compiled binary on 3,000 random repeat-rich sequences (3–400 bp; W ∈ {3,5,8,16,30,64,100};
T ∈ {10,12,15,20,25,30}) and one 1-Mb sequence: 0 mismatches.

## Online Sources

### Morgulis, Gertz, Schäffer & Agarwala (2006). A fast and symmetric DUST implementation to mask low-complexity DNA sequences. J Comput Biol 13(5):1028–1040.

**URL:** https://pubmed.ncbi.nlm.nih.gov/16796549/
**Accessed:** 2026-06-14 (fetched via WebFetch of the PubMed abstract page)
**Authority rank:** 1 (peer-reviewed primary paper defining the algorithm)

**Key Extracted Points:**

1. **Citation:** WebFetch of the PubMed page returned: Title "A fast and symmetric DUST implementation to mask low-complexity DNA sequences"; Authors Aleksandr Morgulis, E Michael Gertz, Alejandro A Schäffer, Richa Agarwala; Journal of Computational Biology; 13(5):1028-40; 2006; DOI 10.1089/cmb.2006.13.1028.
2. **Method:** WebSearch result block (researchr/PubMed summary) states DUST "is a heuristic algorithm that employs a scoring function based on counting nucleotide triplet frequencies in 64-base windows" and that the new implementation "uses the same function to assign a complexity score to a sequence" while changing only the masking rule (now symmetric and context-insensitive).

### Li, H. (2025). Finding low-complexity DNA sequences with longdust (arXiv:2509.07357).

**URL:** https://arxiv.org/pdf/2509.07357
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (primary preprint by the author of minimap2/samtools; restates the SDUST score verbatim)

**Key Extracted Points:**

1. **Score formula (verbatim from WebFetch):** the SDUST complexity score for a string `x` of length `L` is `∑_t c_x(t)(c_x(t)-1)/2 / (L-2)`, where `c_x(t)` is the count of triplet `t` in `x` and the sum runs over all triplets.
2. **Normalization:** division is by `(L-2)`, i.e. the number of triplets ℓ = L−2 (k = 3 hardcoded).
3. **Thresholding form:** the internal scoring function is written `S_S(c_x) = 1/ℓ(x) ∑_t [c_x(t)(c_x(t)-1)/2] − T` (WebSearch result block), confirming both the `1/ℓ` normalization and the threshold `T` subtracted from the raw score.
4. **Default parameters:** default window size `w = 64`; threshold corresponds to complexity level 20, score 2.0 (WebFetch).
5. **Direction:** repeated triplets raise `∑ c(c-1)/2`, so a HIGH score indicates LOW complexity; low-complexity (high-scoring) regions are the ones masked (WebSearch result block).

### lh3/sdust — reference C implementation (`sdust.c`, master).

**URL:** https://raw.githubusercontent.com/lh3/sdust/master/sdust.c
**Accessed:** 2026-06-14 (fetched via WebFetch of the raw source file)
**Authority rank:** 3 (reference implementation by Heng Li, reimplementing the symmetric DUST algorithm)

**Key Extracted Points:**

1. **Incremental score accumulation (verbatim line):** `++*L, *rw += cw[t]++, *rv += cv[t]++;` — on adding triplet `t`, the *current* count `cw[t]` is added to the running score `rw` before the count is incremented. Summing the pre-increment counts 0+1+2+…+(c−1) over all occurrences of a triplet yields exactly `c(c−1)/2`, so `rw = ∑_t c(c−1)/2` and `L` counts the triplets added.
2. **Threshold comparison (verbatim line):** `if (rw * 10 > L * T)` — equivalent to `rw / L > T / 10`; with default `T = 20` this is `score > 2.0`, confirming the threshold value 2.0.
3. **Default parameters (verbatim):** `int W = 64;` and `int T = 20;` — default window width 64 and threshold level 20 (score 2.0).

---

## Documented Corner Cases and Failure Modes

### From Li (2025) / Morgulis et al. (2006)

1. **Below triplet length (L < 3):** the score `∑ c(c-1)/2 / (L-2)` is undefined when L−2 ≤ 0 (no triplets exist); the operation has no defined complexity for such inputs.
2. **All-distinct triplets:** when every triplet occurs exactly once, every `c(c-1)/2 = 0`, so the score is 0 (maximum complexity).
3. **Maximally repetitive (single triplet):** a homopolymer of length L has one triplet repeated L−2 times, giving the maximal score `(L−2)(L−3)/2 / (L−2) = (L−3)/2`.

---

## Test Datasets

> **Note (2026-09-28):** the hand-derived table below used the superseded ℓ divisor. Correct values
> (divisor ℓ − 1): `ATGC` 0; `ACGTACGT` 2/5 = 0.4; `AAAAAA` 6/3 = 2.0; `ACACACAC` 6/5 = 1.2;
> `AAAAAAAAAA` 28/7 = 4.0; `AATAATAA` 3/5 = 0.6; `AAAAAAA` 10/4 = 2.5; (AC)×6 20/9; `AAA` (ℓ = 1) ⇒ 0.

### Dataset: Hand-derived worked examples (k = 3, divisor = number of triplets = L−2)

**Source:** Derived directly from the Li (2025) formula `∑_t c_t(c_t−1)/2 / (L−2)` and the lh3/sdust accumulation.

| Input | L | Triplets (count) | Σ c(c−1)/2 | L−2 | Score |
|-------|---|------------------|------------|-----|-------|
| `ATGC` | 4 | ATG=1, TGC=1 | 0 | 2 | 0.0 |
| `ACGTACGT` | 8 | ACG=2,CGT=2,GTA=1,TAC=1 | 1+1=2 | 6 | 0.3333333333… |
| `AAAAAA` | 6 | AAA=4 | 4·3/2=6 | 4 | 1.5 |
| `ACACACAC` | 8 | ACA=3,CAC=3 | 3+3=6 | 6 | 1.0 |
| `AAAAAAAAAA` | 10 | AAA=8 | 8·7/2=28 | 8 | 3.5 |

### Dataset: Threshold reference

**Source:** lh3/sdust `if (rw * 10 > L * T)` with `T = 20`; Li (2025) score 2.0 / level 20.

| Parameter | Value |
|-----------|-------|
| Default window size | 64 |
| Default threshold (mask if score >) | 2.0 |
| Triplet/word size k | 3 |

---

## Assumptions

1. **ASSUMPTION: General word size `wordSize`** — The paper and reference implementation hardcode k = 3 (triplets). The repository method exposes a `wordSize` parameter; for `wordSize = w` the normalization generalizes to (number of words − 1) = `L − w` (= ℓ − 1 = `L − 3` when w = 3). This generalization is consistent with the formula but only k = 3 is source-backed; tests assert exact source-derived values only for k = 3.
2. **ASSUMPTION: Fewer than two words (ℓ ≤ 1)** — The normaliser ℓ − 1 is 0 and no word pair exists; the implementation returns 0; this is a defined-output convention, not a source value.

---

## Recommendations for Test Coverage

1. **MUST Test:** Homopolymer `AAAAAA` ⇒ score 1.5 (k=3) — Evidence: Li (2025) formula, hand-derived.
2. **MUST Test:** `ACGTACGT` ⇒ score 0.333… (k=3) — Evidence: Li (2025) formula, hand-derived.
3. **MUST Test:** All-distinct-triplet input `ATGC` ⇒ score 0.0 — Evidence: Li (2025) formula; all c(c−1)/2 = 0.
4. **MUST Test:** Repetitive dinucleotide `ACACACAC` ⇒ score 1.0 — Evidence: Li (2025) formula, hand-derived.
5. **MUST Test:** Longer homopolymer `AAAAAAAAAA` ⇒ score 3.5 — Evidence: Li (2025) formula, hand-derived.
6. **MUST Test:** DnaSequence and string overloads agree (same input ⇒ same score) — Evidence: both wrap one core.
7. **SHOULD Test:** Case-insensitivity (string overload upper-cases) — Rationale: DnaSequence normalizes to upper-case; string overload documents `ToUpperInvariant`.
8. **SHOULD Test:** Null DnaSequence ⇒ ArgumentNullException; null/empty string ⇒ 0 — Rationale: documented validation, matches sibling methods.
9. **COULD Test:** Input shorter than wordSize ⇒ 0 — Rationale: defined-output convention (ASSUMPTION 2).

---

## References

1. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. (2006). A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13(5):1028–1040. https://doi.org/10.1089/cmb.2006.13.1028 (abstract retrieved via https://pubmed.ncbi.nlm.nih.gov/16796549/)
2. Li H. (2025). Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
3. Li H. sdust — Symmetric DUST for finding low-complexity regions in DNA sequences (reference C implementation). https://raw.githubusercontent.com/lh3/sdust/master/sdust.c (accessed 2026-06-14)
4. NCBI C++ Toolkit dustmasker `symdust.cpp`/`symdust.hpp` — https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/dustmask/symdust.cpp (accessed 2026-09-28)

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Normaliser corrected to ℓ − 1 (sdust/dustmasker source + binary); SDUST masking evidence added.
- **2026-09-30**: F34 triplets-only + longdust, F35 dustmasker linker/soft mask, F36 N-splitting string overload (B04 WP4); cross-check numbers above.
