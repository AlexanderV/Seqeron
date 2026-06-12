# Validation Report: POP-DIV-001 — Diversity Statistics (π, Watterson's θ, Tajima's D)

- **Validated:** 2026-06-12   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateNucleotideDiversity`, `CalculateWattersonTheta`, `CalculateTajimasD`, `CalculateDiversityStatistics` (+ private `CalculateObserved/ExpectedHeterozygosity`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_Diversity_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened
- Wikipedia "Tajima's D" (fetched 2026-06-12) — confirms the main formula and every variance constant.
- Wikipedia "Watterson estimator" (fetched 2026-06-12) — confirms θ̂ = K/aₙ and aₙ = Σ_{i=1}^{n−1} 1/i (the (n−1)th harmonic number).
- TestSpec POP-DIV-001.md, Evidence POP-DIV-001-Evidence.md (cites Nei & Li 1979, Watterson 1975, Tajima 1989).

### Formula check (verified EXACTLY against Wikipedia)
- **Nucleotide diversity:** π = Σ_{i<j} d_ij / (C(n,2)·L), d_ij = # differing sites. Per-site normalization (÷L). Confirmed.
- **Watterson:** θ_w = S/aₙ, aₙ = Σ_{i=1}^{n−1} 1/i. Per-site: S/(aₙ·L). Harmonic sum runs to n−1 (NOT n). Confirmed.
- **Tajima's D:** D = (k̂ − S/a₁) / √(e₁S + e₂S(S−1)), with:
  - a₁ = Σ_{i=1}^{n−1} 1/i ; a₂ = Σ_{i=1}^{n−1} 1/i²
  - b₁ = (n+1)/(3(n−1)) ; b₂ = 2(n²+n+3)/(9n(n−1))
  - c₁ = b₁ − 1/a₁ ; c₂ = b₂ − (n+2)/(a₁·n) + a₂/a₁²
  - e₁ = c₁/a₁ ; e₂ = c₂/(a₁² + a₂)
  
  Every constant matches the Wikipedia / Tajima (1989) definition character-for-character. No off-by-one in either harmonic sum.

### Edge-case semantics (all sourced)
n<2 → π=0; S=0 → θ=0,D=0; n<2 or L≤0 → θ=0; n<3 → D=0 (Tajima requires n≥3); variance≤0 → D=0 guard. All defined and standard.

### Independent cross-check — hand computation (Wikipedia worked example n=5, S=4, k̂=2.0)
| Quantity | Hand value | Source |
|----------|-----------|--------|
| a₁ | 25/12 = 2.083333 | Wiki "a₁=2.08" |
| a₂ | 205/144 = 1.423611 | — |
| b₁ | 0.500000 | (n+1)/(3(n−1)) |
| b₂ | 0.366667 | 2(n²+n+3)/(9n(n−1)) |
| c₁ | 0.020000 | b₁−1/a₁ |
| c₂ | 0.022667 | — |
| e₁ | 0.009600 | c₁/a₁ |
| e₂ | 0.003933 | c₂/(a₁²+a₂) |
| S/a₁ | 1.920000 | Wiki "1.92" |
| d | 0.080000 | Wiki "d = 2 − 1.92 = .08" |
| V | 0.085590 | e₁S + e₂S(S−1) |
| **D** | **0.273450** | matches spec "D ≈ 0.273" |

Watterson cross-check (WT-M01): S=10, n=10, L=1000 → a₉=2.82897 → θ = 10/(2.82897·1000) = 0.0035349 ≈ 0.00353. Confirmed.

**Findings:** None. Description is mathematically exact.

## Stage B — Implementation

### Code realises the formula (file:line)
- π — `CalculateNucleotideDiversity` (L179–208): all C(n,2) pairs, count differing sites, `totalDiff / (comparisons * length)`. Exact.
- θ — `CalculateWattersonTheta` (L213–226): a1 loop `i=1; i<sampleSize` (= Σ_{i=1}^{n−1}), returns `S/(a1*L)`. Exact, correct harmonic bound.
- D — `CalculateTajimasD` (L240–278): a1,a2 loop to n−1; b1,b2,c1,c2,e1,e2 (L262–267) match validated formulas symbol-for-symbol; variance = e1·S + e2·S(S−1); guard `variance<=0 → 0`; `d = k̂ − S/a1`. Exact.
- `CalculateDiversityStatistics` (L283–324): counts S correctly, computes k̂ = π·L and feeds unnormalized k̂ to Tajima's D (correct — Tajima uses average pairwise diffs, not per-site). H_obs = Nei (1978) n/(n−1)·(1−Σp²); H_exp = (1−Σp²); identity H_obs = n/(n−1)·H_exp holds.

### Cross-verification recomputed vs code
- Hand D=0.273450 → test TD-C01/TD-C02 assert 0.273 ±0.005: PASS.
- θ_W(n=5,S=4,L=20) = 0.0960 → TD-C02 asserts 0.096: PASS.
- DS-M01 hand values (S=2, π=1/6, θ=1/6, D=0, H_exp=1/9, H_obs=1/6) independently recomputed in Python — all match the test assertions exactly.

### Variant/delegate consistency
Single canonical path; `CalculateDiversityStatistics` reuses the public π/θ/D methods, no divergent duplicate.

### Test quality audit
29 spec tests (99 total in the Diversity file). Assertions check exact sourced values (1/6, 25/12, 0.273, 0.096, Nei identity) — not tautologies. All Stage-A edge cases (empty, single, n<3, S=0, L≤0, identical) covered. Deterministic.

**Findings / defects:** None.

## Verdict & follow-ups
- Stage A PASS, Stage B PASS. No code changes required.
- Tests: `--filter ~Diversity` 99/99 passed; full suite 4484/4484 passed (baseline held).
- **State: CLEAN.**
