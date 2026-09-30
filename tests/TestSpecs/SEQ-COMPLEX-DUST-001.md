# Test Specification: SEQ-COMPLEX-DUST-001

**Test Unit ID:** SEQ-COMPLEX-DUST-001
**Area:** Complexity
**Algorithm:** DUST Score (triplet-frequency low-complexity score) + symmetric DUST (SDUST) masking (+ dustmasker linker / soft mask) + longdust k-mer generalisation
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-09-30

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Morgulis et al. (2006). A fast and symmetric DUST implementation. J Comput Biol 13(5):1028–1040. | 1 | https://doi.org/10.1089/cmb.2006.13.1028 (abstract: https://pubmed.ncbi.nlm.nih.gov/16796549/) | 2026-06-14 |
| 2 | Li, H. & Li, B. (2025). Finding low-complexity DNA sequences with longdust. arXiv:2509.07357 (manuscript `tex/longdust.tex` in the repo). | 1 | https://arxiv.org/abs/2509.07357 | 2026-09-30 |
| 3 | lh3/sdust — reference C implementation (`sdust.c`), compiled and run. | 3 | https://raw.githubusercontent.com/lh3/sdust/master/sdust.c | 2026-09-28 |
| 4 | NCBI dustmasker `symdust.cpp`/`symdust.hpp`, `dust_mask_app.cpp`; BLAST+ 2.12.0 `dustmasker` binary run. | 3 | https://raw.githubusercontent.com/ncbi/ncbi-cxx-toolkit-public/master/src/algo/dustmask/symdust.cpp | 2026-09-30 |
| 5 | lh3/longdust 1.4-r97 (`longdust.c`, MIT) — compiled and run. | 3 | https://github.com/lh3/longdust | 2026-09-30 |

### 1.2 Key Evidence Points

1. Score `S(x) = (Σ_t c_t(c_t−1)/2) / (ℓ−1)`, ℓ = L−2 triplets (k = 3). — sdust `find_perfect` (`new_r*10 > T*new_l`, `new_l = ℓ−1`), dustmasker (`thresholds_[i] = i·level`), Morgulis (2006); confirmed with the compiled sdust binary (Evidence, 2026-09-28 revision). The 2026-06 `/(L−2)` reading of Li (2025) is superseded.
2. The reference implementation accumulates `r += c[t]++`, i.e. `r = Σ_t c_t(c_t−1)/2`. — lh3/sdust.
3. Default window size = 64; threshold = 2.0 (level 20, since `T=20` and `rw/L > T/10`). — Li (2025) & lh3/sdust.
4. Higher score ⇒ lower complexity (repeated triplets raise the sum); high-scoring regions are masked. — Li (2025) / Morgulis et al. (2006).
5. DUST is defined for triplets only: symdust `triplet_type`, sdust `SD_WLEN 3`, longdust README "It [SDUST] hardcodes k=3"; sdust.c `TODO: is this right for SD_WLEN!=3?`. Other word sizes are rejected (F34).
6. dustmasker linker: merge when `prev.last + linker ≥ next.first` (closed coords), linker ∈ 1–32, default 1 (= sdust adjacency merge). — symdust `save_masked_regions`.
7. sdust: "N effectively breaks input into pieces of independent sequences" (sdust.c comment); upstream code leaks the window across N (e.g. interval 35–72 on a 53-bp input) — each ACGT run is scanned independently here (F36).
8. longdust `S_L(x) = Σ log c(t)! − f(ℓ/4^k)`, low complexity iff `S_L − T·ℓ > 0`; defaults k 7, w 5000, T 0.6, X-drop 50, b 3. — Li & Li (2025), longdust.c.

### 1.3 Documented Corner Cases

- All-distinct triplets ⇒ every `c(c−1)/2 = 0` ⇒ score 0 (maximum complexity).
- Homopolymer of length L ≥ 4 ⇒ score `(L−2)/2`.
- ℓ ≤ 1 (L ≤ 3) ⇒ no triplet pair exists; score undefined in the sources (implementation returns 0 by convention).

### 1.4 Known Failure Modes / Pitfalls

1. Dividing by the number of words ℓ instead of ℓ − 1 deflates the score (a 7-A run would score 2.0 and not be masked, contrary to sdust/dustmasker).
2. Treating high score as high complexity inverts the masking decision. — Li (2025).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `CalculateDustScore(DnaSequence, int)` | SequenceComplexity | Canonical | Triplet score over DnaSequence; default wordSize 3. |
| `CalculateDustScore(string, int)` | SequenceComplexity | Canonical | String overload; upper-cases; null/empty ⇒ 0. |
| `MaskLowComplexity(DnaSequence, int, double, char)` | SequenceComplexity | Canonical | SDUST perfect-interval masking; port of lh3/sdust `sdust_core`. |
| `MaskLowComplexity(DnaSequence, int, double, char, int linker, bool softMask)` | SequenceComplexity | Canonical | + dustmasker linker / soft mask. |
| `MaskLowComplexity(string, …, linker, softMask)` | SequenceComplexity | Canonical | Raw string with N/IUPAC (independent ACGT runs). |
| `FindLowComplexityIntervals(string / DnaSequence, int, double, int)` | SequenceComplexity | Canonical | SDUST intervals `[start, end)`. |
| `CalculateLongdustScore(string, int, double?)` | SequenceComplexity | Canonical | longdust `S_L(x)`. |
| `FindLongdustRegions(string, …)` | SequenceComplexity | Canonical | Port of longdust `ld_dust1/ld_dust2`. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score ≥ 0 for all inputs (sum of binomial coefficients ≥ 0; divisor > 0). | Yes | Formula |
| INV-2 | All-distinct words ⇒ score = 0. | Yes | Li (2025) (each c(c−1)/2 = 0) |
| INV-3 | Score = (Σ_t c_t(c_t−1)/2) / (ℓ−1) for wordSize 3. | Yes | sdust / dustmasker source |
| INV-4 | DnaSequence and string overloads return identical scores for the same (upper-case) sequence. | Yes | Both wrap one core |
| INV-5 | Homopolymer of length L ≥ 4 (k=3) ⇒ score = (L−2)/2. | Yes | Derived |
| INV-6 | `MaskLowComplexity` output equals lh3/sdust intervals (masked with maskChar); same length as input. | Yes | sdust binary |
| INV-7 | SDUST masking is reversal-symmetric. | Yes | Morgulis (2006) |
| INV-8 | linker = 1 output ≡ default (sdust) output; larger linker only merges (masked set grows monotonically with linker). | Yes | symdust rule |
| INV-9 | For input with non-ACGT symbols, intervals = sdust run on each maximal ACGT run, shifted; no interval covers a non-ACGT position. | Yes | sdust contract |
| INV-10 | Soft mask: `ToUpper(soft) == ToUpper(input)`; lower-case positions = hard-masked positions. | Yes | dustmasker -outfmt fasta |
| INV-11 | longdust intervals equal compiled longdust; forward ⊆ both-strand union. | Yes | longdust binary |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Homopolymer AAAAAA | ℓ=4, AAA=4, Σ=6, /3 | 2.0 | sdust formula |
| M2 | ACGTACGT | ACG=2,CGT=2,GTA=1,TAC=1, Σ=2, /5 | 0.4 | sdust formula |
| M3 | All-distinct ATGC | ATG=1,TGC=1, Σ=0 | 0.0 | INV-2 |
| M4 | ACACACAC | ACA=3,CAC=3, Σ=6, /5 | 1.2 | sdust formula |
| M5 | Homopolymer AAAAAAAAAA | ℓ=8, Σ=28, /7 | 4.0 | INV-5 |
| M6 | Overload agreement | DnaSequence vs string, AAAAAA | equal (2.0) | INV-4 |
| M7 | Score ≥ 0 (INV-1) | AATAATAA | 0.6 | formula |
| M8 | sdust threshold boundary | A×7, A×6, (AC)×6, ACACACACACA | 2.5, 2.0, 20/9, 2.0 | sdust binary (masked iff > 2.0) |
| M9 | Single triplet | AAA (ℓ=1) | 0.0 | convention |
| M10 | SDUST masking = sdust | 15 reference rows (W 5–64, T 10–30) | exact masked strings | sdust binary |
| M11 | Reversal symmetry | mask(rev(x)) = rev(mask(x)) | equal | INV-7 |
| M12 | Mask validation | W < 3, threshold < 0 / NaN, null, linker ∉ 1–32 | ArgumentOutOfRange / ArgumentNull | contract |
| M13 | wordSize ≠ 3 rejected | 0, −1, 1, 2, 4, 7, 8 | ArgumentOutOfRange (MCP: ArgumentOutOfRange / ArgumentException) | F34 sources |
| M14 | Linker = dustmasker | 96-bp A/(CA)/A sequence, linker 1, 17, 18, 19, 32 | [10,26)[43,59)[81,94) ×2; [10,59)[81,94) ×2; [10,94) | dustmasker 2.12.0 binary |
| M15 | Soft mask = dustmasker -outfmt fasta | same sequence, linker 1 and 32; lower-case input | exact strings | dustmasker binary |
| M16 | N input = sdust per ACGT run | 5 rows (N, R, n; W 16/64) | exact intervals + masked strings | sdust binary per piece |
| M17 | longdust regions | VNTR 408 bp, STR+N 307 bp; defaults, -k5 -w100, -f, -k4 -w200 -t1.0 | (120,288); (80,138)(202,227); (117,291); …; (120,288) | longdust binary |
| M18 | longdust score | 6 rows incl. GC 0.41 | 17-digit values from longdust's `ld_cal_f/ld_cal_f2` | longdust source (compiled helper) |
| M19 | longdust validation | k 0/15, w 0/65535, T 0/NaN, X-drop −1, b 1, GC 0/1, null | exceptions | longdust asserts/options |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null DnaSequence | `CalculateDustScore((DnaSequence)null!)` | ArgumentNullException | Matches sibling methods |
| S2 | Null/empty string | `CalculateDustScore("")`, `null` | 0.0 | Documented validation |
| S3 | Case-insensitivity | "aaaaaa" via string overload | 2.0 (= "AAAAAA") | ToUpperInvariant |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Shorter than wordSize | "AT", wordSize 3 | 0.0 | ASSUMPTION 2 convention |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No `SequenceComplexity_CalculateDustScore_Tests.cs` existed prior to this unit. An implementation of `CalculateDustScore` existed in `SequenceComplexity.cs` but with a wrong divisor (`total - 1` instead of the triplet count `total = L−2`), and no dedicated tests.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | No prior test |
| M2 | ❌ Missing | No prior test |
| M3 | ❌ Missing | No prior test |
| M4 | ❌ Missing | No prior test |
| M5 | ❌ Missing | No prior test |
| M6 | ❌ Missing | No prior test |
| M7 | ❌ Missing | No prior test |
| S1 | ❌ Missing | No prior test |
| S2 | ❌ Missing | No prior test |
| S3 | ❌ Missing | No prior test |
| C1 | ❌ Missing | No prior test |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceComplexity_CalculateDustScore_Tests.cs` — all cases above.
- **Remove:** nothing (no prior DUST tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceComplexity_CalculateDustScore_Tests.cs` | Canonical | 11 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented overload-agreement test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented non-negativity + exact test | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented exception test | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented null/empty test | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented case-insensitivity test | ✅ Done |
| 11 | C1 | ❌ Missing | Implemented short-input test | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact 1.5 |
| M2 | ✅ | Exact 1/3 |
| M3 | ✅ | Exact 0.0 |
| M4 | ✅ | Exact 1.0 |
| M5 | ✅ | Exact 3.5 |
| M6 | ✅ | Overload equality |
| M7 | ✅ | ≥ 0 + exact |
| S1 | ✅ | ArgumentNullException |
| S2 | ✅ | 0.0 |
| S3 | ✅ | 1.5 (case-insensitive) |
| C1 | ✅ | 0.0 |

In-scope cases: 11. ✅ count: 11.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | ~~General `wordSize` extrapolation~~ — removed (F34): `wordSize ≠ 3` throws | M13 |
| 2 | ℓ ≤ 1 ⇒ 0 (defined-output convention) | C1, M9 |

---

## 7. Open Questions / Decisions

1. The source hardcodes k = 3; the repository exposes a `wordSize` parameter. Decision (2026-09-30, F34): keep the parameter for source compatibility but reject every value other than 3; the sourced k-mer generalisation is longdust (`CalculateLongdustScore` / `FindLongdustRegions`).
2. **2026-09-28 (review campaign B04):** the normaliser was corrected from ℓ to ℓ − 1 (sdust `find_perfect`, dustmasker `thresholds_`, confirmed with the compiled sdust binary), and `MaskLowComplexity` was replaced by the real SDUST perfect-interval algorithm (identical to sdust on 3,000 random cases). Tests M1–M7/S3 were re-locked to the reference-confirmed values; §5 audit rows above describe the 2026-06 state.
3. **2026-09-30 (B04 completeness audit WP4):** F34 `wordSize ≠ 3` rejected + longdust added (1 500 inputs / 12 818 intervals, 0 mismatches vs compiled longdust; `S_L` bit-identical on 3 000); F35 dustmasker `linker` + soft mask (1 500 inputs, 0 interval / 0 soft-mask mismatches vs dustmasker 2.12.0); F36 `string` overload + `FindLowComplexityIntervals` with N splitting (3 000 inputs, 1 037 with non-ACGT, 0 mismatches vs sdust per ACGT run). Tests M13–M19.
