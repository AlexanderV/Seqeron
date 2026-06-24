# Validation Report: POP-LD-001 — Linkage Disequilibrium (D, D', r²)

- **Validated:** 2026-06-24   **Area:** Population Genetics (PopGen)
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateLD(var1Id, var2Id, genotypes, distance)`, `PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs` (lines 729–846)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Linkage disequilibrium"** (fetched 2026-06-24) — confirmed verbatim:
  - `D = p_AB − p_A·p_B`
  - Lewontin's `D' = D / D_max` with the **sign-dependent** normalizer:
    - `D < 0` → `D_max = min{ p_A·p_B , (1−p_A)(1−p_B) }`
    - `D > 0` → `D_max = min{ p_A(1−p_B) , p_B(1−p_A) }`
    - range −1 ≤ D' ≤ 1
  - `r² = D² / (p_A(1−p_A)·p_B(1−p_B))`, range [0,1]
  - "LD for diploid frequencies": diploid dosage covariance = D/2 (with the indicator scaling), and after variance normalization the **diploid correlation equals the haploid r_AB** ("Surprisingly, this is identical to the haploid LD correlation rAB").
- Cross-referenced **Lewontin (1964)** (D' normalization), **Hill & Robertson (1968)** (r²), **Gabriel et al. (2002)** (block thresholding) via the Wikipedia citations.

### Diploid-dosage convention (key nuance vs the "haplotype-frequency" framing)
The unit's API does **not** take haplotype frequencies p_AB/p_A/p_B; it takes **diploid genotype dosages (0/1/2)** and estimates LD phase-free from the genotype covariance. This is a sourced design (Wikipedia "LD for diploid frequencies"). With the 0/1/2 dosage encoding, under random union of gametes `Cov(X₁,X₂) = 2D`, so `D = Cov/2`. Wikipedia's "Cov = D/2" uses 0/½/1 indicator scaling — the factor-2 difference is purely the variable scaling and cancels in the normalized r². The code's `d = cov/2` is correct for 0/1/2 dosages.

### Edge-case semantics (sourced)
- Empty data → D=0, r²=0, D'=0 (no association).
- Monomorphic locus → variance = 0 → r² undefined; convention returns 0 (guarded, not NaN).
- Complete LD → |D'| = 1; linkage equilibrium → D=0 → D'=0, r²=0.

### Independent cross-check (hand computation)
Abstract haplotype-freq example `p_A=p_B=0.5, p_AB=0.4`:
- D = 0.4 − 0.25 = **0.15**; D>0 → D_max = min(0.25, 0.25) = **0.25** → D' = **0.6**; r² = 0.15²/0.0625 = **0.36**. ✓

Concrete genotype dataset (Perfect LD) `(0,0)(0,0)(1,1)(1,1)(2,2)(2,2)`:
- mean₁=mean₂=1; cov = var₁ = var₂ = 2/3 → r² = (2/3)²/((2/3)²) = **1.0**.
- d = (2/3)/2 = 1/3; p₁=p₂=0.5 → D_max = min(0.25,0.25)=0.25 → D' = (1/3)/0.25 = 1.333 → clamp **1.0**. ✓

No-LD 3×3 balanced design → cov = 0 → r² = 0, D' = 0. ✓

### Findings / divergences
None. TestSpec + Evidence match the authoritative sources exactly, including the sign-dependent D_max and the r² denominator `p_A·p_a·p_B·p_b`.

## Stage B — Implementation

### Code path reviewed
`CalculateLD` (PopulationGeneticsAnalyzer.cs:729–781), `FindHaplotypeBlocks` (787–846).

### Formula realised correctly? (evidence)
- Allele freqs `p1 = Σgeno1/(2n)`, `p2 = Σgeno2/(2n)` (B-allele frequency from 0/1/2 dosage). ✓
- Population `cov`, `var1`, `var2` (÷n). `rSquared = cov²/(var1·var2)` only when `var1>0 && var2>0`, else 0; clamped [0,1] → squared Pearson correlation of dosages = haplotype r². ✓
- `d = cov/2` → correct D for 0/1/2 encoding (Cov = 2D). ✓
- `dMax = d ≥ 0 ? min(p1·q2, q1·p2) : min(p1·p2, q1·q2)` → **sign-dependent D_max exactly matching Wikipedia/Lewontin**. ✓
- `dPrime = dMax > 1e-10 ? |d|/dMax : 0`, then `min(dPrime, 1.0)` → div-by-zero guarded, clamped to [0,1]. ✓
- Monomorphic locus → variance 0 → r²=0 via guard; D'=0 via dMax guard → no NaN/Infinity. ✓

### Cross-verification table recomputed vs code
| Case | Inputs | Expected | Code result |
|------|--------|----------|-------------|
| Perfect LD | (0,0)(0,0)(1,1)(1,1)(2,2)(2,2) | r²=1.0, D'=1.0 | cov=var=2/3 → r²=1.0; d=1/3, dMax=0.25 → 1.333→clamp 1.0 ✓ |
| No LD (3×3) | balanced | r²=0, D'=0 | cov=0 → r²=0, d=0 → D'=0 ✓ |
| Monomorphic locus 1/2 | one locus constant | r²=0, D'=0, no NaN | var=0 → guard → 0; p=0/1 → dMax guard → 0 ✓ |
| Empty | [] | r²=0, D'=0 | early return ✓ |
| Abstract haplo-freq | p_A=p_B=.5, p_AB=.4 | D'=0.6, r²=0.36 | matches formula by hand (Stage A) |

### Variant/delegate consistency
`FindHaplotypeBlocks` delegates pairwise LD to `CalculateLD` (adjacent-pair r² ≥ threshold, default 0.7; simplified Gabriel 2002). Variants ordered by position; blocks require ≥2 variants; Start ≤ End. Consistent with canonical.

### Test quality audit
31 LD tests assert exact sourced values (r²=1.0 / D'=1.0 within 1e-10 for perfect LD; 0.0 for no-LD), explicit NaN rejection on monomorphic loci, parametrized [0,1] range invariants, and block formation/ordering/position/threshold cases. Real assertions, deterministic, cover all Stage-A edge cases.

### Findings / defects
None. The three protocol-flagged pitfalls are handled correctly:
- D_max **is** sign-dependent (correct D').
- r² denominator is the full `Var₁·Var₂` (= p_A·p_a·p_B·p_b after normalization) — no missing factor.
- Monomorphic locus guarded (`var>0`, `dMax>1e-10`) — no division by zero / NaN.

## Verdict & follow-ups
- **Stage A: PASS** — formulas, sign-dependent D_max, r² denominator match Wikipedia/Lewontin/Hill&Robertson; both abstract and concrete examples reproduced by hand.
- **Stage B: PASS** — code faithfully realises the validated formulas with correct guards; cross-verification reproduced; tests assert exact sourced values.
- **State: CLEAN** — no defect. No code changes. LD filter: 31 passed / 0 failed.
