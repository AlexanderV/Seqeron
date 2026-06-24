# Validation Report: POP-DIV-001 — Diversity Statistics (π, Watterson's θ, Tajima's D, heterozygosity)

- **Validated:** 2026-06-24   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateNucleotideDiversity`, `CalculateWattersonTheta`, `CalculateTajimasD`, `CalculateDiversityStatistics` (+ private `CalculateObservedHeterozygosity`/`CalculateExpectedHeterozygosity`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_Diversity_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened (live web, 2026-06-24)
- Wikipedia "Tajima's D" — confirms the main formula `D = (k̂ − S/a₁) / √(e₁S + e₂S(S−1))` and **every** variance constant a₁, a₂, b₁, b₂, c₁, c₂, e₁, e₂ character-for-character, plus the n=5 worked example (S=4, a₁=2.08, pairwise diffs {3,2,2,3,1,3,2,2,1,1}, k̂=2.0, S/a₁=1.92, d=0.08).
- Wikipedia "Watterson estimator" — confirms θ̂_w = K/aₙ with aₙ = Σ_{i=1}^{n−1} 1/i (the (n−1)th harmonic number; summation bound is n−1, NOT n).
- TestSpec POP-DIV-001.md, Evidence POP-DIV-001-Evidence.md (cite Nei & Li 1979, Watterson 1975, Tajima 1989, Nei 1978).

### Formula check
- **Nucleotide diversity:** π = Σ_{i<j} d_ij / (C(n,2)·L), d_ij = # differing sites; **per-site** normalization (÷L). Matches Nei & Li / Wikipedia.
- **Watterson:** θ_w = S/(aₙ·L), aₙ summed to n−1. Confirmed exact, correct harmonic bound.
- **Tajima's D:** every constant matches Wikipedia/Tajima(1989) exactly. k̂ is the **unnormalized** average pairwise differences (per-pair, not per-site).
- **Heterozygosity:** H_exp = (1−Σp²) per site (basic gene diversity, Wikipedia Zygosity); H_obs = Nei (1978) unbiased n/(n−1)·(1−Σp²) per site. Identity H_obs = n/(n−1)·H_exp holds. For haploid data H_obs ≡ π (per-site).

### Edge-case semantics (all sourced)
n<2 → π=0; S=0 → θ=0,D=0; n<2 or L≤0 → θ=0; n<3 → D=0 (Tajima requires n≥3); variance≤0 → D=0 guard. All defined and standard.

### Independent cross-check — hand computation (Python, this session)
Tajima Wikipedia example n=5, S=4, k̂=2.0:

| Quantity | Hand value |
|----------|-----------|
| a₁ | 2.083333 |
| a₂ | 1.423611 |
| b₁ / b₂ | 0.500000 / 0.366667 |
| c₁ / c₂ | 0.020000 / 0.022667 |
| e₁ / e₂ | 0.009600 / 0.003933 |
| S/a₁ | 1.920000 |
| d | 0.080000 |
| V | 0.085590 |
| **D** | **0.273450** |
| π = k̂/L | 0.1000 |

Watterson WT-M01 (S=10, n=10, L=1000): a₉ = 2.82897 → θ = 10/(2.82897·1000) = 0.0035349 ≈ 0.00353. Confirmed.

**Findings:** None. Description is mathematically exact.

## Stage B — Implementation

### Code path reviewed
- π — `CalculateNucleotideDiversity` (L205–234): all C(n,2) pairs, count differing sites, `totalDiff/(comparisons*length)`. Per-site. Exact.
- θ — `CalculateWattersonTheta` (L239–252): `a1` loop `i=1; i<sampleSize` (= Σ_{i=1}^{n−1}); returns `S/(a1*L)`. Correct harmonic bound; n<2 or L≤0 → 0.
- D — `CalculateTajimasD` (L266–304): a1,a2 to n−1; b1,b2,c1,c2,e1,e2 (L288–293) match validated formulas symbol-for-symbol; variance = e1·S + e2·S(S−1); guard `variance<=0 → 0`; `d = k̂ − S/a1`. Exact. Guards S=0 or n<3 → 0.
- `CalculateDiversityStatistics` (L309–350): counts S correctly; feeds **unnormalized** k̂ = π·L to Tajima (correct); n<2 → all zeros.
- Heterozygosity (L358–415): H_obs = Nei n/(n−1)·(1−Σp²) per site; H_exp = (1−Σp²) per site. Match.

### Cross-verification recomputed vs code
- DS-M01 (seqs ACGTACGT/ACGTATGT/ACGTACGA, n=3, L=8) recomputed in Python: S=2, π=1/6, θ_W=1/6, D=0, H_exp=1/9, H_obs=1/6 — all match the test's exact assertions and the code output.
- Tajima example D=0.273450 → tests TD-C01/TD-C02 assert 0.273±0.005: PASS.
- Watterson WT-M01 θ=0.0035349 → asserts 0.00353±0.0005: PASS.

### Variant/delegate consistency
`CalculateDiversityStatistics` reuses the canonical `CalculateNucleotideDiversity`/`CalculateWattersonTheta`/`CalculateTajimasD` directly; no divergent duplicate logic. `DiversityStatistics` record fields populated correctly (SampleSize = input count; `SegregratingSites` — known intentional typo, kept for API compat per spec note 2).

### Test quality audit
29 tests, all assert exact sourced values (Wikipedia D≈0.273, π=0.1, θ≈0.00353, exact fractions 1/6, 1/9, 2/12, S/L), plus edge cases (empty/single/identical/all-different) and range invariants. Not tautological. Deterministic. Full filtered run: 29 passed / 0 failed.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS, Stage B: PASS — State CLEAN.** All formulas verified against live Wikipedia (Tajima's D, Watterson) and independent Python hand computation; per-site vs per-sequence normalization is correct throughout (π and θ per-site ÷L; Tajima fed unnormalized k̂). No code changes required. No defects.
