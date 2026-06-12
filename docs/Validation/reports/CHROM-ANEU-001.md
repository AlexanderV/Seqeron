# Validation Report: CHROM-ANEU-001 — Aneuploidy / Copy-Number Variation Detection

- **Validated:** 2026-06-12   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth, binSize)`, `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(copyNumberStates, minFraction)` (supporting type `CopyNumberState`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (lines 121–127, 826–916)
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Aneuploidy_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Aneuploidy"** (fetched 2026-06-12): aneuploidy = "the presence of an abnormal number of chromosomes in a cell" (e.g. 45 or 47 instead of 46), explicitly excluding whole extra haploid sets (that is euploid polyploidy). Copy-number terminology confirmed verbatim:
  - Nullisomy = 0 copies; Monosomy = 1; Disomy = 2 (normal diploid); Trisomy = 3; Tetrasomy = 4; Pentasomy = 5.
  - Named examples confirmed: **Trisomy 21 (Down syndrome)** = extra chr21; **Monosomy X (Turner syndrome)** = single X, karyotype 45,X.
- **Wikipedia, "Copy-number variation"** (fetched 2026-06-12): confirms NGS read-depth methods are standard for CNV detection. The standard read-depth model (cross-checked against bioinformatics convention): copy-number ratio = observed coverage / expected diploid baseline; log2 ratio = log2(ratio); deletion → ratio < 1 (log2 < 0); duplication → ratio > 1 (log2 > 0).
- **Griffiths et al. (2000), An Introduction to Genetic Analysis** (cited in Evidence doc): same aneuploidy nomenclature (mono-/tri-/tetra-somy) for abnormal chromosome number.

### Formula check
The spec/Evidence define the read-depth approach:
```
logRatio    = log2(observedDepth / medianDepth)
copyNumber  = round(2^logRatio × 2) = round((observedDepth / medianDepth) × 2)
confidence  = 1 − min(1, |copyNumber/2 − 2^logRatio|)
```
This is the canonical read-depth CNV model with a diploid baseline (×2). log base = 2 (correct). Disomy baseline of 2 copies at ratio 1.0 is correct. Threshold for calling = rounding to nearest integer copy number, which is the explicit design choice the spec documents (no fuzzy log2 ±0.3 gate; classification is integer-copy rounding). This is internally consistent and matches Wikipedia's integer-copy terminology.

### Edge-case semantics check
- Empty input / median ≤ 0 → empty result (avoids ÷0): defined and sourced (Implementation).
- Nullisomy: depth 0 → log2(0) = −∞, 2^(−∞) = 0 → CN 0; confidence = 1 − |0 − 0| = 1.0. Defined.
- Clamp to [0, 10] documented.
- Whole-chromosome classification requires a dominant CN over ≥ minFraction (default 0.8) of bins; CN=2 is never flagged (definition of aneuploidy = abnormal number).
- Documented limitation: sex chromosomes are not special-cased (a normal male monosomic X/Y would be flagged) — explicitly accepted in spec §8 and Evidence §7.

### Independent cross-check (hand computation, medianDepth = 30)
| Depth | Ratio | log2(ratio) | round(ratio×2) → CN | Call | Confidence |
|-------|-------|-------------|---------------------|------|-----------|
| 0.0  | 0.0 | −∞     | 0 | Nullisomy (deletion) | 1.0 |
| 15.0 | 0.5 | −1.000 | 1 | Monosomy (deletion)  | 1.0 |
| 30.0 | 1.0 | 0.000  | 2 | Disomy / Normal      | 1.0 |
| 45.0 | 1.5 | 0.585  | 3 | Trisomy (duplication)| 1.0 |
| 60.0 | 2.0 | 1.000  | 4 | Tetrasomy (duplication)| 1.0 |

All worked values match the spec's "Standard Depth Values" table and the Evidence "Depth to Copy Number Mapping" table. log2(1.5) ≈ 0.585 and log2(0.5) = −1 confirmed.

### Findings / divergences
None. The description is biologically and mathematically correct and faithfully sourced.

## Stage B — Implementation

### Code path reviewed
- `DetectAneuploidy` — `ChromosomeAnalyzer.cs:831–875`. Groups by chromosome, bins by `Position / binSize`, averages depth per bin, computes `logRatio = Math.Log2(meanDepth / medianDepth)`, `copyNumber = Round(Pow(2, logRatio) * 2)` clamped to [0,10], `confidence = 1 − min(1, |copyNumber/2 − Pow(2, logRatio)|)`. Output ordered by bin key.
- `IdentifyWholeChromosomeAneuploidy` — `ChromosomeAnalyzer.cs:880–916`. Per chromosome, takes the dominant CN; if its fraction ≥ minFraction and CN ≠ 2, emits a record with the classification name (0→Nullisomy, 1→Monosomy, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, else "Copy number = N").

### Formula realised correctly?
Yes. `Math.Pow(2, Math.Log2(x))` = x, so `copyNumber = Round((mean/median) × 2)` — exactly the validated formula with the correct log base (2) and the correct diploid ×2 baseline. Confidence formula matches the spec exactly. The classification switch maps each integer CN to the Wikipedia term; CN=2 is correctly excluded from aneuploidy calls.

### Cross-verification table recomputed vs code
The five worked rows above were re-derived against the code path and match the test assertions in `DetectAneuploidy_*` (Normal/Trisomy/Monosomy/Nullisomy/Tetrasomy) and `IdentifyWholeChromosomeAneuploidy_*`. LogRatio assertions check exact `Math.Log2` values (incl. `double.NegativeInfinity` for nullisomy), and confidence = 1.0 at every CN boundary — all pass.

### Variant/delegate consistency
The static `GetAneuploidyTerm` used by `AnalyzeKaryotype` (line 186) and the inline switch in `IdentifyWholeChromosomeAneuploidy` (line 902) use the same standard nomenclature. The karyotype path adds "Disomy" for 2 (not used by the aneuploidy classifier, which excludes 2 by definition) and "Polysomy (N copies)" vs "Copy number = N" for CN>5 — a cosmetic label difference between the two methods, both unambiguous; not a defect.

### Test quality audit
Tests assert exact sourced values (CN, exact `Math.Log2` log-ratios, confidence 1.0 at boundaries), cover all Stage-A edge cases (empty, zero/negative median, clamp to 10, very low → 0, single point, binning/averaging, multi-chromosome grouping, ordering, minFraction at/below/custom threshold, CN>5 formatting). No tautological "no-throw" assertions. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — terminology and read-depth/log2 model independently confirmed against Wikipedia (Aneuploidy, Copy-number variation) and Griffiths et al.; worked example reproduced.
- **Stage B: PASS** — code computes the validated formula exactly; tests lock the sourced values; build + full suite green.
- **State: CLEAN.** No defects. `dotnet build` succeeds (0 warnings/errors); Aneuploidy filter = 35 passed; full `Seqeron.Genomics.Tests` = 4484 passed, 0 failed.
- Accepted limitation (pre-documented, not a defect): sex chromosomes are not baseline-adjusted, so a normal male X/Y could read as monosomic.
