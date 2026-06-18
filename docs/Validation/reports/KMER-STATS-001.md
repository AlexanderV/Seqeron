# Validation Report: KMER-STATS-001 — K-mer Statistics

- **Validated:** 2026-06-16   **Area:** K-mer
- **Canonical method(s):** `KmerAnalyzer.AnalyzeKmers(string sequence, int k)` → `KmerStatistics(TotalKmers, UniqueKmers, MaxCount, MinCount, AverageCount, Entropy)` (with helpers `CountKmers`, `CalculateKmerEntropy`, `GetKmerFrequencies`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect)

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)

1. **Wikipedia — K-mer** (https://en.wikipedia.org/wiki/K-mer) — WebFetched.
   Confirms verbatim: a sequence of length L has **L − k + 1** k-mers; k-mer universe = **nᵏ**.
   Worked example **GTAGAGCTGT**: k=1 → G,T,A,C; k=2 → GT,TA,AG,GA,AG,GC,CT,TG,GT (9, with GT & AG doubled);
   k=3 → 8 distinct trimers. **AGAT** → three 2-mers AG, GA, AT.
2. **BioInfoLogics — k-mer counting, part I** (Clavijo 2018) — WebFetched.
   Confirms **L − k + 1**; **ATCGATCAC** k=3 count table ATC=2, TCG=1, CGA=1, GAT=1, TCA=1, CAC=1 →
   7 total, **6 distinct**, 5 unique; and the **distinct ("counted once") vs unique ("appear only once")** distinction.
3. **Manca / Bonnici / Franco (2021), "Spectral concepts in genome informational analysis", arXiv:2106.15351** —
   PDF downloaded and text-extracted this session. Verbatim definition:
   > *"the k-entropy E_k(G) is defined as the Shannon entropy of the probability distribution assigning to any k-mer its frequency:*
   > *E_k(G) = Σ_{α∈D_k(G)} p(α) log₂ p(α)   where p(α) = mult_G(α) / (|G| − k + 1)."*
   Log base **2** (bits); p(α) is the multiplicity over the **L − k + 1** windows. (Primary, authority rank 1.)

### Formula check

| Quantity | Validated formula | Source |
|----------|-------------------|--------|
| TotalKmers | L − k + 1 = Σ_α mult(α) | Wikipedia; BioInfoLogics |
| UniqueKmers | distinct k-mer count | Wikipedia tables; BioInfoLogics |
| Max/Min | extremes of multiplicity table | Wikipedia table (GTAGAGCTGT k=1 → max 4, min 1) |
| AverageCount | Total/Distinct (mean multiplicity) | arithmetic; derived from totals |
| Entropy | −Σ p(α) log₂ p(α), p(α)=mult/(L−k+1), bits | Manca 2021 (primary), corroborated by arXiv:2511.05300 |

### Edge-case semantics (all sourced)

- empty / null sequence, k > L ⇒ L−k+1 ≤ 0 ⇒ no k-mers ⇒ all-zero statistics.
- k ≤ 0 ⇒ invalid (ArgumentOutOfRangeException via CountKmers).
- single distinct k-mer (homopolymer) ⇒ entropy = −1·log₂1 = **0** (minimum diversity).
- all windows distinct ⇒ entropy = **log₂(L−k+1)** (maximum diversity).

### Independent cross-check (numbers, hand-computed in Python this session — not from the repo)

| Seq | k | Total | Distinct | Max | Min | Avg | Entropy (bits) |
|-----|---|-------|----------|-----|-----|-----|----------------|
| GTAGAGCTGT | 1 | 10 | 4 | 4 | 1 | 2.5 | 1.846439344671 |
| GTAGAGCTGT | 2 | 9 | 7 | 2 | 1 | 9/7→1.29 | 2.725480556998 |
| GTAGAGCTGT | 3 | 8 | 8 | 1 | 1 | 1.0 | 3.000000000000 (=log₂8) |
| ATCGATCAC | 3 | 7 | 6 | 2 | 1 | 7/6→1.17 | 2.521640636343 |
| AGAT | 2 | 3 | 3 | 1 | 1 | 1.0 | 1.584962500721 (=log₂3) |
| AAAA | 2 | 3 | 1 | 3 | 3 | 3.0 | 0 |

All match the TestSpec / Evidence expected values exactly, and the count tables match the
Wikipedia / BioInfoLogics tables.

### Findings / divergences

- **Sign convention (noted, not a defect):** the Manca paper prints `E_k = Σ p log₂ p` **without** the
  leading minus sign, which as written would make E_k ≤ 0. Standard Shannon entropy negates the sum
  (`−Σ p log₂ p`) to yield a non-negative quantity; the same paper's equipartition-maximum proposition
  (k-hapax genome has the *maximum* E_k) only holds for the non-negative form. The implementation
  correctly uses `−Σ p log₂ p`. This is a documented sign-convention nuance, not a description error.

**Stage A verdict: PASS** — every formula, definition, edge case, and worked number traces to an
authoritative source retrieved this session.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/KmerAnalyzer.cs`:
- `AnalyzeKmers` (L528–560): builds the count table via `CountKmers`; all-zero record on empty table
  (L532–542); `TotalKmers = values.Sum()`, `UniqueKmers = counts.Count`, `MaxCount = values.Max()`,
  `MinCount = values.Min()`, `AverageCount = Math.Round(values.Average(), 2)`, `Entropy = CalculateKmerEntropy`.
- `CountKmers` (L20–42): upper-cases input; `k≤0` throws; `k>L` or null/empty ⇒ empty dict; otherwise
  L−k+1 overlapping windows.
- `CalculateKmerEntropy` (L331–346) over `GetKmerFrequencies` (L177–189): p(α) = count/Σcount =
  count/(L−k+1); `entropy −= freq·Math.Log2(freq)` with the `freq>0` guard implementing 0·log0 = 0.

### Formula realised correctly?

Yes. `values.Sum()` = Σ mult = L−k+1, `counts.Count` = distinct, `values.Average()` = Total/Distinct,
and the entropy term is exactly `−Σ p log₂ p` with p = mult/(L−k+1) and the 0·log0 = 0 convention.
Log base is `Math.Log2` (bits). Matches the validated description.

### Cross-verification table recomputed vs code

The 15-test fixture was run (`--filter ...AnalyzeKmers`): **15 passed, 0 failed**. Each asserted value
equals the independently hand-computed value above (entropies within 1e-10; integers/records exact).

### Variant/delegate consistency

`AnalyzeKmers` reuses `CountKmers` and `CalculateKmerEntropy`; M6/M7 cross-check Total = Σcounts and
Unique = distinct count directly against `CountKmers`, confirming agreement.

### Numerical robustness

`int` totals are bounded by sequence length (no overflow on stated ranges); div-by-zero avoided by the
empty-table early return; `Math.Log2(freq)` guarded by `freq>0`. AverageCount rounded to 2 dp for
display (documented; exact ratio also verifiable). No precision concerns.

### Test quality audit (HARD gate)

- **Sourced expectations, not code echoes:** entropy constants (1.846439344671, 2.725480556998, 3.0,
  2.521640636343, 1.584962500721, 0) and all totals/distinct/max/min are asserted as exact sourced
  values — a deliberately-wrong implementation fails the absolute assertion, not an echo. ✅
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` for entropy and exact record equality for
  edge cases; S2/S3 are genuine INV-4/INV-5 invariant bounds (supplementary to the exact entropies
  already locked by M1–M5, not a substitute for them). No skip/ignore, no widened tolerance, no
  weakened assertion. ✅
- **Cover all the logic:** worked examples (M1–M5), homopolymer→entropy 0 (S1), all-distinct→log₂D
  (M3,M5), empty (M8), null (C2), k>L (M9), k≤0 throws (M10), case-insensitivity (C1), INV-1/2/3
  cross-checks (M6,M7). All Stage-A branches and documented edge/error cases exercised. ✅
- **Honest green:** full **unfiltered** `dotnet test` = **6607 passed, 0 failed, 1 skipped** (the skip is
  the unrelated `MFE_Benchmark_AllScenarios`, not this unit); `dotnet build` = **0 errors** (4 pre-existing
  warnings in unrelated files; no file was changed in this unit). ✅

### Findings / defects

None. No code, test, or spec change was required.

## Verdict & follow-ups

- **Stage A: PASS · Stage B: PASS · End-state: ✅ CLEAN.**
- Test-quality gate: **PASS** (exact sourced values, full coverage, honest unfiltered green 6607/0).
- One BY-DESIGN note logged (Manca formula sign convention) in the findings register; no action needed.
- No follow-ups.
