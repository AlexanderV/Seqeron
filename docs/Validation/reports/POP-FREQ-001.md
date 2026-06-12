# Validation Report: POP-FREQ-001 — Allele Frequencies (allele freq, MAF, genotype freq)

- **Validated:** 2026-06-12   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(int, int, int)`,
  `CalculateMAF(IEnumerable<int>)`, `FilterByMAF(IEnumerable<Variant>, double, double)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs:112-170`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia "Allele frequency" (en.wikipedia.org/wiki/Allele_frequency) | p = f(AA) + ½f(AB), q = f(BB) + ½f(AB); equivalent to (2·hom + het)/(2N); p + q = 1. Diploid example: 6 AA, 3 AB, 1 BB → A copies = 6×2+3 = 15, B copies = 1×2+3 = 5, total = 20 → **p = 15/20 = 0.75, q = 5/20 = 0.25**. |
| Wikipedia "Minor allele frequency" (en.wikipedia.org/wiki/Minor_allele_frequency) | MAF = frequency of the **second most common** allele; therefore MAF ≤ 0.5; conceptually MAF = min of the two allele frequencies, i.e. min(p, q). |
| Wikipedia "Genotype frequency" (en.wikipedia.org/wiki/Genotype_frequency) | Genotype frequency = (individuals with that genotype) / (total individuals). Four-o'clock example: 49 AA, 42 Aa, 9 aa → freq(a) = (42 + 2×9)/200 = 60/200 = 0.30, **p = 0.70, q = 0.30**. Allele and genotype frequencies sum to 1. |
| Gillespie (2004) *Population Genetics: A Concise Guide*, ISBN 978-0-8018-8008-7 | Standard diploid allele-counting convention (homozygote contributes 2 copies, heterozygote 1) — consistent with the above. |

### Formula check (exact, as cited)

- **Allele frequency (biallelic, diploid):** with genotype counts (AA, Aa, aa) and N = AA+Aa+aa individuals,
  - p = freq(A) = (2·AA + Aa) / (2N)
  - q = freq(a) = (2·aa + Aa) / (2N)
  - p + q = 1
  Each allele copy counted: homozygote → 2 copies, heterozygote → 1 copy each. Denominator is **2N (diploid)**.
- **MAF** = min(p, q) (the less common allele; ≤ 0.5). Equivalently min(altFreq, 1−altFreq).
- **Genotype frequency** = genotype count / N (observed); sums to 1.

### Edge-case semantics (sourced / derived)

- Monomorphic locus → MAF = 0 (one allele fixed).
- All heterozygous → p = q = 0.5.
- N = 0 → no defined frequency; spec/Evidence prescribe graceful (0, 0).
- Negative genotype counts are not biologically meaningful → rejected (`ArgumentOutOfRangeException`).

### Independent cross-check (exact numbers)

- **Protocol worked example:** AA=30, Aa=50, aa=20, N=100 →
  p = (2·30 + 50)/200 = 110/200 = **0.55**; q = (2·20 + 50)/200 = 90/200 = **0.45**; MAF = min(0.55, 0.45) = **0.45**;
  genotype freqs 30/100 = 0.30, 50/100 = 0.50, 20/100 = 0.20 (sum = 1). Matches spec AF-M12 (NDSU molecular example).
- **Wikipedia diploid (6-3-1)** → p = 0.75, q = 0.25 (AF-M10). ✓
- **Wikipedia flower (49-42-9)** → p = 0.70, q = 0.30 (AF-M02). ✓
- **NDSU blood type (1787-3039-1303)**, 2N = 12258 → f(M) = 6613/12258 ≈ 0.5394, f(N) = 5645/12258 ≈ 0.4606 (AF-M11). ✓

### Findings / divergences

None. All formulas, the 2N diploid denominator, the allele-copy counting, MAF = min(p,q), and the
sum-to-1 invariants are confirmed against authoritative sources.

---

## Stage B — Implementation

### Code path reviewed

`PopulationGeneticsAnalyzer.cs`:
- `CalculateAlleleFrequencies` — lines 112-133
- `CalculateMAF` — lines 138-151
- `FilterByMAF` — lines 156-170

### Formula realised correctly? (evidence)

- **CalculateAlleleFrequencies:** `totalAlleles = 2*(homMaj+het+homMin)` (= 2N, diploid ✓);
  `majorAlleles = 2*homMaj + het`, `minorAlleles = 2*homMin + het` (homozygote 2 copies, het 1 copy ✓);
  returns `(major/total, minor/total)`. Negative counts each throw `ArgumentOutOfRangeException` (lines 117-122).
  `totalAlleles == 0` guard returns `(0, 0)` (line 126). By construction major+minor = (2·hom_maj+het+2·hom_min+het)/2N = 2N/2N = 1.
- **CalculateMAF:** `altFreq = sum(genotypes)/(2*count)` with 0/1/2 encoding, then `Math.Min(altFreq, 1-altFreq)` = min(p,q) ✓; empty → 0 (line 142).
- **FilterByMAF:** `maf = Math.Min(AF, 1-AF)`; include iff `minMAF ≤ maf ≤ maxMAF`; `yield return` (lazy, order-preserving) ✓.

### Cross-verification table recomputed vs code

| Input | Formula result | Code result | Test |
|-------|----------------|-------------|------|
| AA=30, Aa=50, aa=20 | p=0.55, q=0.45, MAF=0.45 | 0.55 / 0.45 | AF-M12 ✓ |
| AA=6, Aa=3, aa=1 | p=0.75, q=0.25 | 0.75 / 0.25 | AF-M10 ✓ |
| AA=49, Aa=42, aa=9 | p=0.70, q=0.30 | 0.70 / 0.30 | AF-M02 ✓ |
| MM=1787, MN=3039, NN=1303 | 0.5394 / 0.4606 | 6613/12258, 5645/12258 | AF-M11 ✓ |
| genotypes [0×5,1×4,2×1] | altFreq=6/20=0.3, MAF=0.3 | 0.3 | MAF-M01 ✓ |
| genotypes [0,1×4,2×5] | altFreq=14/20=0.7, MAF=0.3 | 0.3 | MAF-M02 ✓ |
| all 0 / all 2 | MAF=0 | 0 | MAF-M05/M06 ✓ |
| [0,1,1,2] | altFreq=0.5, MAF=0.5 | 0.5 | MAF-M08 ✓ |
| N=0 | (0,0) | (0,0) | AF-M07 ✓ |

### Variant/delegate consistency

`TestHardyWeinberg` (line 413) and `CalculateLD` (lines 714-715) independently use the same
(2·hom + het)/2N convention — consistent. No `*Fast`/duplicate allele-frequency variant exists.

### Test quality audit

`PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs` (53 cases) asserts exact sourced values
(`Within(1e-10)` or exact fractions such as `6613.0/12258.0`), covers all Stage-A edge cases
(monomorphic, all-het, N=0, single sample, large population, negative-count throws, boundary MAF).
`Properties/PopulationGeneticsProperties.cs` (23 cases) covers AF-C01 (sum=1), AF-C02 (range [0,1]),
MAF-C01 (MAF ∈ [0, 0.5]) as randomized property tests. Assertions are real, deterministic, non-tautological.

### Findings / defects

None. No wrong denominator (correctly 2N, not N), heterozygote correctly counted as 1 copy each,
MAF correctly uses `min` (not `max`).

---

## Verdict & follow-ups

- **Stage A: PASS** — description matches Wikipedia (Allele/Minor allele/Genotype frequency) and Gillespie (2004) exactly.
- **Stage B: PASS** — code realises the validated formulas; all worked examples recomputed and match.
- **State: CLEAN** — no defect found; no code changes required.
- **Tests:** `--filter ~AlleleFrequency` → 53 passed, 0 failed; `~PopulationGeneticsProperties` → 23 passed, 0 failed;
  full `Seqeron.Genomics.Tests` suite → **4484 passed, 0 failed** (baseline preserved).
