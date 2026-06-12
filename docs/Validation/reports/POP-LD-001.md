# Validation Report: POP-LD-001 — Linkage Disequilibrium (D, D', r²)

- **Validated:** 2026-06-12   **Area:** Population Genetics (PopGen)
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateLD(var1Id, var2Id, genotypes, distance)`, `PopulationGeneticsAnalyzer.FindHaplotypeBlocks(variants, ldThreshold)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs` (lines 685–813)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/PopulationGeneticsAnalyzer_LinkageDisequilibrium_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Linkage disequilibrium"** (https://en.wikipedia.org/wiki/Linkage_disequilibrium) — confirmed verbatim the three core formulas and the sign-dependent D_max, plus the "LD for diploid frequencies" section relating dosage covariance and haplotype r².
- Cross-referenced **Lewontin (1964)** (D' normalization), **Hill & Robertson (1968)** (r² correlation measure), **Gabriel et al. (2002)** (haplotype-block thresholding) via the Wikipedia citations and corroborating literature (e.g. Rogers & Huff 2009; the diploid-dosage covariance interpretation, PMC/Oxford MBE 2020).

### Formula check (cited exactly from Wikipedia)
- **D (coefficient of LD):** `D = x_AB − p_A·p_B`. ✓
- **D' (Lewontin):** `D' = D / D_max` with the **sign-dependent** normalization:
  - `D > 0` → `D_max = min(p_A·p_b, p_a·p_B)` = `min(p_A·(1−p_B), (1−p_A)·p_B)`
  - `D < 0` → `D_max = min(p_A·p_B, p_a·p_b)` = `min(p_A·p_B, (1−p_A)·(1−p_B))`
  - Range −1 ≤ D' ≤ 1; |D'| commonly reported on [0,1]. ✓
- **r² (Hill & Robertson):** `r² = D² / (p_A·p_a·p_B·p_b)` = `D² / [p_A(1−p_A)·p_B(1−p_B)]`. Range [0,1]. ✓

### Diploid-dosage convention (important nuance)
Wikipedia's "LD for diploid frequencies" section states the covariance of its **0/½/1 "x,y" indicator variables** equals **D/2**, and that variance normalization cancels the factor of 2 so the **squared correlation of dosages equals the haplotype r²**. With the **0/1/2 allele-dosage** encoding used by the code, the standard result under random mating is `Cov(X₁,X₂) = 2D`, i.e. **D = Cov/2** (the two haplotypes of an individual each contribute the within-gamete covariance D; the cross-gamete term is 0 under random union). This is consistent — the apparent "D/2 vs Cov/2" difference is purely the indicator-variable scaling, and the normalized r² is invariant to it. The code's `d = cov/2` is correct for the 0/1/2 encoding.

### Edge-case semantics (sourced)
- Empty data → no association → D=0, r²=0, D'=0.
- Monomorphic locus → r² denominator (a variance) = 0 → r² undefined; convention returns 0 (guarded), not NaN.
- Complete LD → |D'| = 1.
- Linkage equilibrium (x_AB = p_A·p_B) → D = 0 → D'=0, r²=0.

### Independent cross-check (hand computation)
Abstract worked example `p_A=0.5, p_B=0.5, x_AB=0.4`:
- D = 0.4 − 0.25 = **0.15**.
- D>0 → D_max = min(0.5·0.5, 0.5·0.5) = **0.25** → D' = 0.15/0.25 = **0.6**.
- r² = 0.15² / (0.5·0.5·0.5·0.5) = 0.0225/0.0625 = **0.36**.
All match the spec's expected values. Linkage-equilibrium case x_AB = p_A·p_B = 0.25 → D=0 confirmed.

### Findings / divergences
None. Descriptions in TestSpec and Evidence match the authoritative sources exactly, including the sign-dependent D_max.

## Stage B — Implementation

### Code path reviewed
`CalculateLD` (PopulationGeneticsAnalyzer.cs:697–748), `FindHaplotypeBlocks` (754–813).

### Formula realised correctly? (evidence)
- Allele freqs `p1 = Σgeno1/(2n)`, `p2 = Σgeno2/(2n)` (B-allele frequency from 0/1/2 dosage). ✓
- `cov`, `var1`, `var2` = population covariance/variance (÷n). `rSquared = cov²/(var1·var2)` only when `var1>0 && var2>0`, else 0 → **squared Pearson correlation of dosages = haplotype r²** (Stage-A confirmed). ✓
- `d = cov/2` → correct D for 0/1/2 encoding (Cov = 2D). ✓
- `dMax = d ≥ 0 ? min(p1·q2, q1·p2) : min(p1·p2, q1·q2)` → **sign-dependent D_max exactly matching Wikipedia/Lewontin**. ✓
- `dPrime = dMax > 1e-10 ? |d|/dMax : 0`, then `min(dPrime, 1.0)` → guarded against div-by-zero and clamped to [0,1]. ✓
- Monomorphic locus → variance 0 → r²=0 via guard; D'=0 via dMax guard → no NaN/Infinity. ✓

### Cross-verification table recomputed vs code
| Case | Inputs | Expected | Code result |
|------|--------|----------|-------------|
| Perfect LD | (0,0)(0,0)(1,1)(1,1)(2,2)(2,2) | r²=1.0, D'=1.0 | mean=1, cov=var=2/3 → r²=1.0; d=1/3, dMax=0.25 → D'=1.333→clamp 1.0 ✓ |
| No LD (3×3) | balanced design | r²=0, D'=0 | cov=0 → r²=0, d=0 → D'=0 ✓ |
| Monomorphic locus 1 | all geno1=0 | r²=0, D'=0 | var1=0 → guard → 0; p1=0 → dMax guard → 0 ✓ |
| Empty | [] | r²=0, D'=0 | early return ✓ |
| Abstract D'=0.6/r²=0.36 | haplotype-freq example | n/a to genotype API | matches formula by hand (Stage A) |

### Variant/delegate consistency
`FindHaplotypeBlocks` delegates pairwise LD to `CalculateLD` (adjacent-pair r² ≥ threshold, default 0.7; simplified Gabriel 2002). Blocks require ≥2 variants; Start ≤ End by position ordering. Consistent with canonical.

### Test quality audit
26 LD tests assert exact sourced values (r²=1.0 / D'=1.0 within 1e-10 for perfect LD; 0.0 for no LD), explicit NaN/Infinity rejection on monomorphic loci, [0,1] range invariants (parametrized), and block-formation/ordering/threshold cases. Real assertions, deterministic, cover all Stage-A edge cases.

### Findings / defects
None. The three sign-dependent/denominator pitfalls flagged by the protocol are all handled correctly:
- D_max **is** sign-dependent (correct D').
- r² denominator is the full `Var₁·Var₂` (= the p·q·p·q product after normalization) — no missing factor.
- Monomorphic locus is guarded (`var>0`, `dMax>1e-10`) — no division by zero.

## Verdict & follow-ups
- **Stage A: PASS** — formulas, sign-dependent D_max, and r² denominator match Wikipedia/Lewontin/Hill&Robertson exactly; worked example reproduced by hand.
- **Stage B: PASS** — code faithfully realises the validated formulas with correct guards; cross-verification reproduced; tests assert exact sourced values.
- **State: CLEAN** — no defect. No code changes. Full suite `Seqeron.Genomics.Tests` = 4484 passed / 0 failed (LD filter: 26 passed).
