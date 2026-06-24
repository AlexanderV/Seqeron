# Validation Report: CHROM-ANEU-001 — Aneuploidy Detection (read-depth copy-number)

- **Validated:** 2026-06-24   **Area:** Chromosome Analysis
- **Canonical method(s):** `ChromosomeAnalyzer.DetectAneuploidy(depthData, medianDepth, binSize)`, `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy(copyNumberStates, minFraction)` (supporting type `CopyNumberState`)
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs:832–917`
- **Test file:** `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Aneuploidy_Tests.cs`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm
- **Wikipedia, "Aneuploidy"** (fetched 2026-06-24): aneuploidy = "the presence of an abnormal number of chromosomes in a cell" — explicitly *excluding* whole extra/missing haploid sets (those are euploid polyploidy). Copy-number terminology confirmed verbatim: Nullisomy = 0, Monosomy = 1, Disomy = 2 (normal diploid), Trisomy = 3, Tetrasomy = 4, Pentasomy = 5. Clinical examples confirmed: **Down = trisomy 21**, **Turner = monosomy X (45,X)**, **Edwards = trisomy 18**, **Patau = trisomy 13**.
- **Read-depth CNV model** (standard bioinformatics convention, also documented in Evidence via Wikipedia "Copy number variation"): copy-number ratio = observed coverage / expected diploid baseline; log2 ratio = log2(ratio); deletion → ratio < 1 (log2 < 0); duplication → ratio > 1 (log2 > 0).

### Formula check
The spec/Evidence define the read-depth approach (Evidence §2.1, spec §3.1/§9):
```
logRatio   = log2(observedDepth / medianDepth)
copyNumber = round(2^logRatio × 2) = round((observedDepth / medianDepth) × 2)
confidence = 1 − min(1, |copyNumber/2 − 2^logRatio|)
```
Canonical read-depth CNV model with a diploid (×2) baseline, log base 2. Disomy at ratio 1.0 → CN 2 is correct. The "call" is integer-copy rounding (no fuzzy log2 ±0.3 gate); this is the explicit, internally consistent design choice and matches Wikipedia's integer-copy terminology.

### Boundary vs CHROM-KARYO-001
KARYO-001 (`AnalyzeKaryotype`) produces the full karyotype string and uses a static `GetAneuploidyTerm` for labeling; ANEU-001 is the dedicated *aneuploidy call* from read depth (`DetectAneuploidy`) plus whole-chromosome classification (`IdentifyWholeChromosomeAneuploidy`). The two share nomenclature but the ANEU call is the read-depth→CN path here. No overlap conflict.

### Edge-case semantics check
- Empty input / median ≤ 0 → empty (avoids ÷0): defined, sourced (Implementation).
- Nullisomy: depth 0 → log2(0) = −∞, 2^(−∞) = 0 → CN 0; confidence = 1 − |0 − 0| = 1.0.
- Clamp to [0, 10] documented.
- Whole-chromosome classification: dominant CN over ≥ minFraction (default 0.8) of bins; CN=2 never flagged (definition of aneuploidy).
- Documented limitation (spec §8, Evidence §7): sex chromosomes not baseline-adjusted → a normal male monosomic X/Y would be flagged. Accepted, not a defect.

### Independent cross-check (hand computation, medianDepth = 30)
| Depth | Ratio | log2(ratio) | round(ratio×2) → CN | Call | Confidence |
|-------|-------|-------------|---------------------|------|-----------|
| 0.0  | 0.0 | −∞     | 0 | Nullisomy   | 1.0 |
| 15.0 | 0.5 | −1.000 | 1 | Monosomy    | 1.0 |
| 30.0 | 1.0 | 0.000  | 2 | Disomy/Normal | 1.0 |
| 45.0 | 1.5 | 0.585  | 3 | Trisomy     | 1.0 |
| 60.0 | 2.0 | 1.000  | 4 | Tetrasomy   | 1.0 |
All values match spec §3.1 and Evidence §2.1. log2(1.5) ≈ 0.585, log2(0.5) = −1 confirmed.

### Findings / divergences
None. Description is biologically and mathematically correct and faithfully sourced.

## Stage B — Implementation

### Code path reviewed
- `DetectAneuploidy` — `ChromosomeAnalyzer.cs:832–876`. Guards `data.Count==0 || medianDepth<=0` → yield break. Groups by chromosome, bins by `Position / binSize` (ordered), averages depth per bin, `logRatio = Math.Log2(meanDepth/medianDepth)`, `copyNumber = round(Pow(2,logRatio)*2)` clamped to [0,10], `confidence = 1 − min(1, |copyNumber/2 − Pow(2,logRatio)|)`. Bin span = `[key*binSize, (key+1)*binSize − 1]` so Start < End.
- `IdentifyWholeChromosomeAneuploidy` — `ChromosomeAnalyzer.cs:881–917`. Per chromosome takes dominant CN by fraction; if fraction ≥ minFraction and CN ≠ 2 emits a record with switch-mapped name (0→Nullisomy, 1→Monosomy, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, else "Copy number = N").

### Formula realised correctly?
Yes. `Math.Pow(2, Math.Log2(x)) == x`, so CN = `round((mean/median)×2)` — exactly the validated formula, correct log base (2), correct ×2 diploid baseline. Confidence formula matches the spec exactly. CN=2 correctly excluded from aneuploidy calls.

### Cross-verification table recomputed vs code
The five worked rows above were re-derived against the code path and match the `DetectAneuploidy_*` assertions (Normal/Trisomy/Monosomy/Nullisomy/Tetrasomy). LogRatio assertions check exact `Math.Log2` values incl. `double.NegativeInfinity`; confidence = 1.0 at every CN boundary. Clamp tests (10×→CN 10; 0.01×→CN 0) reproduced.

### Variant/delegate consistency
The static `GetAneuploidyTerm` used by `AnalyzeKaryotype` and the inline switch in `IdentifyWholeChromosomeAneuploidy` use the same nomenclature; cosmetic label difference for CN>5 ("Polysomy (N)" vs "Copy number = N") — not a defect.

### Test quality audit
31 tests assert exact sourced values (CN, exact `Math.Log2` log-ratios, confidence 1.0 at boundaries) and cover all Stage-A edge cases (empty, zero/negative median, clamp to 10, very low → 0, single point, binning/averaging, multi-chromosome grouping, ordering, minFraction below/at/custom threshold, CN>5 formatting). No tautological "no-throw" assertions. Deterministic.

### Findings / defects
None.

## Verdict & follow-ups
- **Stage A: PASS** — terminology, clinical examples, and read-depth/log2 model independently confirmed against Wikipedia ("Aneuploidy"); worked example reproduced by hand.
- **Stage B: PASS** — code computes the validated formula exactly; tests lock the sourced values.
- **State: CLEAN.** No defects. Build: 0 warnings/0 errors. Aneuploidy filter = 31 passed, 0 failed.
- Accepted pre-documented limitation (not a defect): sex chromosomes are not baseline-adjusted, so a normal male X/Y could read as monosomic.
