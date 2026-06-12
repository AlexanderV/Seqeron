# Validation Report: CHROM-CENT-001 — Centromere position / chromosome classification by centromere

- **Validated:** 2026-06-12   **Area:** Chromosome
- **Canonical method(s):** `ChromosomeAnalyzer.AnalyzeCentromere(chromosomeName, sequence, windowSize, minAlphaSatelliteContent)` → private `DetermineCentromereType(chromosomeLength, centStart, centEnd)`
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Centromere_Tests.cs`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia "Centromere"** (fetched 2026-06-12) — reproduces the Levan, Fredga & Sandberg (1964) table verbatim:

  | Centromere position | Arm ratio (q/p) | Sign | Classification |
  |---|---|---|---|
  | Medial *sensu stricto* | 1.0–1.6 | M | Metacentric |
  | Medial region | 1.7 | m | Metacentric |
  | Submedial | 3.0 | sm | Submetacentric |
  | Subterminal | 3.1–6.9 | st | Subtelocentric |
  | Terminal region | 7.0 | t | Acrocentric |
  | Terminal *sensu stricto* | ∞ | T | Telocentric |

  Confirms: p = short arm ("petit"), q = long arm, **arm ratio r = q/p (long/short)**.
- **Independent secondary source** (ResearchGate / Pensoft summaries of Levan 1964) — confirms the **interval** reading used by implementations: m = 1.0–1.7, sm = 1.7–3.0, st = 3.0–7.0, a (t) > 7.0, telocentric = single arm only.
- **Primary reference:** Levan A, Fredga K, Sandberg AA (1964). "Nomenclature for centromeric position on chromosomes." *Hereditas* 52(2):201–220.

### Formula check
- **Arm ratio:** r = q/p (long arm / short arm), r ≥ 1. ✔
- **Centromeric index:** CI = (short arm p) / (total length) × 100, range 0–50. ✔
  CI is consistent with r via CI = 100/(1+r):
  - r = 1.0 → CI = 50 (perfectly metacentric)
  - r = 1.7 → CI = 37.04 (≈ the 37.5 metacentric/submetacentric border; literature rounds the border ratio to 1.7)
  - r = 3.0 → CI = 25.0 (submetacentric/subtelocentric border) ✔
  - r = 7.0 → CI = 12.5 (subtelocentric/acrocentric border) ✔
  - r = ∞ (p = 0) → CI = 0 (telocentric) ✔
  The prompt's idealized CI bands (37.5 / 25 / 12.5) match these; the only divergence is the 37.5↔1.667 vs Levan's rounded 1.7, which is standard in the literature.

### Boundary mapping (spec & sources agree)
- **r = 1.7 → Metacentric** (Levan sign "m")
- **r = 3.0 → Submetacentric** (sign "sm")
- **r = 7.0 → Acrocentric** (sign "t")
- **p = 0 (r = ∞) → Telocentric** (sign "T")

### Edge-case semantics
- p = q → r = 1.0 → Metacentric, CI = 50 (defined). ✔
- p = 0 → Telocentric (defined). ✔
- p > q never occurs because arms are assigned by min/max (swap-safe), so r ≥ 1 always. ✔

### Independent cross-check (hand computations)
- p = 40, q = 60, total = 100 → r = 1.5, CI = 40 → **Metacentric** ✔
- p = 10, q = 90, total = 100 → r = 9.0, CI = 10 → **Acrocentric** (Levan "t" region; true Telocentric only at p = 0) ✔
- IEEE-754 boundary check (verified numerically): 17/10 == literal 1.7 and ≤ 1.7 true; 3000/1000 ≤ 3.0 true; 7000/1000 < 7.0 false. ✔

### Findings / divergences (Stage A notes)
- **Note 1 (nomenclature, not a defect):** the r ≥ 7.0 / CI 0–12.5 class is labelled **Acrocentric** (Levan sign "t"), and **Telocentric** is reserved for p = 0 (sign "T", r = ∞). Some textbooks loosely call the r > 7 class "telocentric"; the spec and code follow Levan's own sign convention, which matches the Wikipedia table. Acceptable.
- **Note 2 (interval vs point table):** Levan's published table lists representative point ratios (1.7, 3.0, 7.0); the half-open interval scheme used by the spec/code (m ≤1.7, sm (1.7,3.0], st (3.0,7.0), a ≥7.0) is the standard operationalization and is corroborated by independent sources.

## Stage B — Implementation

### Code path reviewed
`ChromosomeAnalyzer.cs`:
- `AnalyzeCentromere` (line 360) — detection by sliding window (repeat content × low GC variability), then `DetermineCentromereType`.
- `DetermineCentromereType` (line 494):
  ```csharp
  int centMid = (centStart + centEnd) / 2;
  int pArm = Math.Min(centMid, chromosomeLength - centMid);   // short arm
  int qArm = Math.Max(centMid, chromosomeLength - centMid);   // long arm
  if (pArm == 0) return "Telocentric";
  double armRatio = (double)qArm / pArm;                       // r = q/p, swap-safe
  return armRatio switch {
      <= 1.7 => "Metacentric",
      <= 3.0 => "Submetacentric",
      <  7.0 => "Subtelocentric",
      _      => "Acrocentric"
  };
  ```

### Formula realised correctly?
- r = q/p with p = min arm, q = max arm → exactly the validated Levan arm ratio, swap-safe. ✔
- Thresholds map to the validated boundaries:
  - r = 1.7 → `<= 1.7` → Metacentric ✔
  - r = 3.0 → `<= 3.0` → Submetacentric ✔
  - r = 7.0 → not `< 7.0` → Acrocentric ✔
  - p = 0 → Telocentric ✔
  Boundary inclusivity matches the spec ("1.7 → Metacentric, 3.0 → Submetacentric, 7.0 → Acrocentric") exactly. No CI-uses-q-instead-of-p error; arm ratio is correctly long/short.

### Cross-verification table recomputed vs code
| p | q | total | r = q/p | Code class | Expected (Levan) | Match |
|---|---|-------|---------|-----------|------------------|-------|
| 50 | 50 | 100 | 1.0 | Metacentric | Metacentric | ✔ |
| 40 | 60 | 100 | 1.5 | Metacentric | Metacentric | ✔ |
| 10 | 17 | 27  | 1.7 | Metacentric | Metacentric | ✔ |
| 10 | 30 | 40  | 3.0 | Submetacentric | Submetacentric | ✔ |
| 10 | 50 | 60  | 5.0 | Subtelocentric | Subtelocentric | ✔ |
| 10 | 70 | 80  | 7.0 | Acrocentric | Acrocentric | ✔ |
| 10 | 90 | 100 | 9.0 | Acrocentric | Acrocentric | ✔ |
| 0  | 100| 100 | ∞ | Telocentric | Telocentric | ✔ |

### Variant/delegate consistency
- `IsAcrocentric` is set as `centType == "Acrocentric"` (line 426) — consistent with type (M8). ✔
- Note: a *separate* utility pair `CalculateArmRatio` / `ClassifyChromosomeByArmRatio` (lines 925–952) exists but is **not** part of CHROM-CENT-001 (different canonical method, different p/q-based scheme, tested in `ChromosomeAnalyzerTests.cs`). Out of scope here; the Levan classification under test lives solely in `DetermineCentromereType`.

### Numerical robustness
- Integer division for `centMid` and Min/Max guard before the double divide; `pArm == 0` short-circuits before division → no div-by-zero. ✔
- IEEE-754 boundary equality verified (17/10 ≤ 1.7, 3000/1000 ≤ 3.0, 7000/1000 ≮ 7.0). ✔

### Test quality audit
- 25 tests in `ChromosomeAnalyzer_Centromere_Tests.cs`, all passing. They assert exact class strings (Metacentric/Submetacentric/Subtelocentric/Acrocentric), the `IsAcrocentric` invariant, the valid type-set, and structural invariants (Start ≤ End, Length = End − Start), plus edge cases (empty/null/short/uniform/non-repetitive) and case-insensitivity.
- Gap (minor): tests reach the boundaries via the stochastic sliding-window detector (arm ratios ~1.0, ~2.0, ~5.0, ~21), so the **exact** boundary values r = 1.7 / 3.0 / 7.0 are not directly exercised. These were instead verified here by source-checked hand computation and IEEE-754 trace; the classification code is a pure deterministic switch on the ratio, so the boundary behavior is fully determined and correct. Not a defect.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES — Levan (1964) nomenclature, arm-ratio r = q/p, CI = p/total×100 ∈ [0,50], and boundary mapping (1.7→M, 3.0→sm, 7.0→t/Acrocentric, ∞→Telocentric) all confirmed against Wikipedia + independent sources. Notes: r>7 class is Levan "t"=Acrocentric (Telocentric only at p=0); interval operationalization of the point table is standard.
- **Stage B:** PASS — `DetermineCentromereType` realises r = q/p (swap-safe, no q/p inversion) with thresholds and boundary inclusivity matching the validated Levan boundaries; worked examples recompute correctly; edge cases handled.
- **State:** CLEAN — no defect found; no code change required. Full suite 4484/4484 passing; Centromere filter 25/25 passing.
- **No defects logged.**
