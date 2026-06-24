# Validation Report: POP-FREQ-001 — Allele Frequencies (allele freq, MAF, genotype freq)

- **Validated:** 2026-06-24   **Area:** Population Genetics
- **Canonical method(s):** `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies(int, int, int)`,
  `CalculateMAF(IEnumerable<int>)`, `FilterByMAF(IEnumerable<Variant>, double, double)`
  (`src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs:138-196`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened (this session) & what they confirm

| Source | Confirms |
|--------|----------|
| Wikipedia "Allele frequency" (en.wikipedia.org/wiki/Allele_frequency) — fetched 2026-06-24 | p = f(AA) + ½f(AB), q = f(BB) + ½f(AB); denominator is **2N** (chromosome copies); diploid worked example **6 AA, 3 AB, 1 BB → A copies = 6×2+3 = 15, B copies = 1×2+3 = 5, total = 20 → p = 0.75, q = 0.25**; p + q = 1. |
| Wikipedia "Minor allele frequency" (en.wikipedia.org/wiki/Minor_allele_frequency) — fetched 2026-06-24 | MAF = frequency of the **second most common** allele → MAF ≤ 0.5 by definition; equivalently the smaller of the two allele frequencies, i.e. min(p, 1−p). |
| Wikipedia "Genotype frequency" (per Evidence + prior report) | Four-o'clock example 49 AA, 42 Aa, 9 aa → freq(a) = (42 + 2×9)/200 = 0.30, p = 0.70. |
| Gillespie (2004) *Population Genetics: A Concise Guide*, ISBN 978-0-8018-8008-7 | Standard diploid allele-counting convention (homozygote = 2 copies, heterozygote = 1 each). |

### Formula check (exact, as cited)

- **Allele frequency (biallelic diploid)**, counts (AA, Aa, aa), N = AA+Aa+aa:
  - p = (2·AA + Aa) / (2N), q = (2·aa + Aa) / (2N), p + q = 1. Denominator **2N (diploid)** — confirmed.
- **MAF** = min(altFreq, 1 − altFreq) = min(p, q) ≤ 0.5.
- **Genotype frequency** = genotype count / N; sums to 1.

### Edge-case semantics (sourced / derived)

- Monomorphic locus → MAF = 0 (one allele fixed).
- All heterozygous → p = q = 0.5.
- N = 0 → no defined frequency; spec/Evidence prescribe graceful (0, 0).
- Negative genotype counts not biologically meaningful → rejected (`ArgumentOutOfRangeException`).

### Independent cross-check (exact numbers, computed by hand this session)

- **Wikipedia diploid (6-3-1):** A=15, B=5, 2N=20 → p = 15/20 = **0.75**, q = 5/20 = **0.25**. ✓ (AF-M10; matches fetched Wikipedia value)
- **Wikipedia flower (49-42-9):** (2·49+42)/200 = 140/200 = **0.70**; (2·9+42)/200 = 60/200 = **0.30**. ✓ (AF-M02)
- **NDSU molecular (30-50-20):** (2·30+50)/200 = 110/200 = **0.55**; (2·20+50)/200 = 90/200 = **0.45**; MAF = 0.45. ✓ (AF-M12)
- **NDSU blood type (1787-3039-1303):** 2N = 12258 → f(M) = 6613/12258 ≈ **0.5394**, f(N) = 5645/12258 ≈ **0.4606**. ✓ (AF-M11)
- **MAF** genotypes [0×5,1×4,2×1]: altFreq = 6/20 = 0.3 → MAF = **0.3**. ✓; genotypes with altFreq 0.7 → MAF = 1−0.7 = **0.3** ✓; all-0 / all-2 → MAF = **0** ✓.

### Findings / divergences

None. The 2N diploid denominator, allele-copy counting, MAF = min(p,q), and the sum-to-1 /
range invariants are all confirmed against authoritative sources.

---

## Stage B — Implementation

### Code path reviewed

`PopulationGeneticsAnalyzer.cs`:
- `CalculateAlleleFrequencies` — lines 138-159
- `CalculateMAF` — lines 164-177
- `FilterByMAF` — lines 182-196

Source unchanged since the last code touch (`b240ab1c`, which added the negative-count validation); current code matches the prior validated state.

### Formula realised correctly? (evidence)

- **CalculateAlleleFrequencies:** `totalAlleles = 2*(homMaj+het+homMin)` (= 2N ✓);
  `majorAlleles = 2*homMaj + het`, `minorAlleles = 2*homMin + het` (homozygote 2 copies, het 1 each ✓);
  returns `(major/total, minor/total)`. Each negative count throws `ArgumentOutOfRangeException` (lines 143-148).
  `totalAlleles == 0` guard returns `(0, 0)` (line 152-153). By construction major+minor = 2N/2N = 1.
- **CalculateMAF:** `altFreq = sum(genotypes)/(2*count)` with 0/1/2 encoding; `Math.Min(altFreq, 1-altFreq)` = min(p,q) ✓; empty → 0 (line 168-169).
- **FilterByMAF:** `maf = Math.Min(AF, 1-AF)`; include iff `minMAF ≤ maf ≤ maxMAF`; `yield return` (lazy, order-preserving) ✓.

### Cross-verification table recomputed vs code (via tests)

| Input | Formula result | Test |
|-------|----------------|------|
| AA=6, Aa=3, aa=1 | p=0.75, q=0.25 | AF-M10 ✓ |
| AA=49, Aa=42, aa=9 | p=0.70, q=0.30 | AF-M02 ✓ |
| AA=30, Aa=50, aa=20 | p=0.55, q=0.45 | AF-M12 ✓ |
| MM=1787, MN=3039, NN=1303 | 6613/12258, 5645/12258 | AF-M11 ✓ |
| genotypes [0×5,1×4,2×1] | MAF=0.3 | MAF-M01 ✓ |
| genotypes altFreq=0.7 | MAF=0.3 | MAF-M02 ✓ |
| all 0 / all 2 | MAF=0 | MAF-M05/M06 ✓ |
| [0,1,1,2] | MAF=0.5 | MAF-M08 ✓ |
| N=0 | (0,0) | AF-M07 ✓ |

### Variant/delegate consistency

`TestHardyWeinberg` and `CalculateLD` independently use the same (2·hom + het)/2N convention. No `*Fast`/duplicate allele-frequency variant exists.

### Test quality audit

`PopulationGeneticsAnalyzer_AlleleFrequency_Tests.cs` asserts exact sourced values (`Within(1e-10)` or
exact fractions such as `6613.0/12258.0`), covers all Stage-A edge cases (monomorphic, all-het, N=0,
single sample, large population, negative-count throws, boundary MAF). `Properties/PopulationGeneticsProperties.cs`
covers AF-C01 (sum=1), AF-C02 (range [0,1]), MAF-C01 (MAF ∈ [0, 0.5]) as randomized property tests.
Assertions are real, deterministic, non-tautological.

### Findings / defects

None. Correct 2N denominator (not N), heterozygote correctly counted as 1 copy each, MAF correctly uses `min` (not `max`).

---

## Verdict & follow-ups

- **Stage A: PASS** — description matches Wikipedia (Allele/Minor allele/Genotype frequency) and Gillespie (2004) exactly; re-fetched live this session.
- **Stage B: PASS** — code realises the validated formulas; all worked examples recomputed and match.
- **State: CLEAN** — no defect found; no code changes required.
- **Tests:** `--filter ~AlleleFrequency | ~PopulationGeneticsProperties` → **88 passed, 0 failed**. No code changed, so full-suite run not required.
