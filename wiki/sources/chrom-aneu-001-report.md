---
type: source
title: "Validation report: CHROM-ANEU-001 (aneuploidy detection — read-depth copy number, ChromosomeAnalyzer.DetectAneuploidy / IdentifyWholeChromosomeAneuploidy)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-ANEU-001.md
sources:
  - docs/Validation/reports/CHROM-ANEU-001.md
source_commit: 15b20caaeff245cda0742f1f97189ee5419ef422
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-ANEU-001

The two-stage **validation write-up** for test unit **CHROM-ANEU-001** (aneuploidy detection —
estimating per-bin copy number from sequencing read depth and classifying whole-chromosome ploidy
against the normal disomic CN=2 baseline), validated 2026-06-24. This is the *report* artifact that
feeds one row of the [[validation-ledger]]; it records the validator's independent **verdict** on both
the algorithm description (Stage A) and the shipped code (Stage B), and the wider campaign is
[[validation-and-testing]]. The depth→copy-number model, the classification rule, the confidence
formula, the documented clinical oracles and the two limitations are synthesized in the concept
[[aneuploidy-detection]] (the chromosome copy-number/ploidy family anchor); [[test-unit-registry]]
defines the unit. Distinct from [[chrom-aneu-001-evidence]] — the pre-implementation evidence artifact
sourced from `docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No code defect; no test change required. Build
**0 warnings / 0 errors**; the aneuploidy filter ran **31 passed, 0 failed**.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs:832–917`:

- `DetectAneuploidy(depthData, medianDepth, binSize)` (`:832–876`) — guards `data.Count==0 ||
  medianDepth<=0` → `yield break`; groups by chromosome, bins by `Position / binSize` (ordered),
  averages depth per bin, then `logRatio = Math.Log2(meanDepth/medianDepth)`,
  `copyNumber = round(2^logRatio × 2)` clamped to `[0, 10]`,
  `confidence = 1 − min(1, |copyNumber/2 − 2^logRatio|)`. Bin span `[key*binSize, (key+1)*binSize − 1]`
  so `Start < End`.
- `IdentifyWholeChromosomeAneuploidy(copyNumberStates, minFraction)` (`:881–917`) — per chromosome
  takes the dominant CN by bin fraction; if fraction ≥ `minFraction` (default 0.8) **and CN ≠ 2**, emits
  a switch-mapped record (0→Nullisomy, 1→Monosomy, 3→Trisomy, 4→Tetrasomy, 5→Pentasomy, else
  "Copy number = N").
- Supporting type `CopyNumberState`.
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_Aneuploidy_Tests.cs`.

## Stage A — description (algorithm faithfulness)

Confirmed against **Wikipedia "Aneuploidy"** (fetched 2026-06-24): aneuploidy = "an abnormal number of
chromosomes in a cell", explicitly *excluding* whole extra/missing haploid sets (euploid polyploidy).
Copy-number ladder confirmed verbatim (Nullisomy 0 / Monosomy 1 / Disomy 2 / Trisomy 3 / Tetrasomy 4 /
Pentasomy 5) along with the clinical oracles Down = trisomy 21, Turner = monosomy X (45,X), Edwards =
trisomy 18, Patau = trisomy 13. The read-depth CNV model (copy-number ratio = observed / expected
diploid baseline; log2 ratio; deletion ratio < 1, duplication > 1) is the standard convention.

The validated formula is the canonical read-depth CNV model with a diploid (×2) baseline, log base 2 —
disomy at ratio 1.0 → CN 2. The call is **integer-copy rounding** (no fuzzy log2 ±0.3 gate), an
explicit, internally consistent design choice matching the integer-copy terminology.

**Independent hand cross-check (medianDepth = 30):** depth 0.0 → ratio 0.0, log2 −∞, CN 0 Nullisomy;
15.0 → 0.5, −1.000, CN 1 Monosomy; 30.0 → 1.0, 0.000, CN 2 Disomy/Normal; 45.0 → 1.5, 0.585, CN 3
Trisomy; 60.0 → 2.0, 1.000, CN 4 Tetrasomy — all confidence 1.0 (expected = observed at every
integer-CN ratio). Every row matches the spec and Evidence. **Findings: none.**

## Boundary vs [[karyotype-analysis]] / CHROM-KARYO-001

CHROM-KARYO-001 (`AnalyzeKaryotype`) produces the full karyotype string and uses a static
`GetAneuploidyTerm` for labeling; CHROM-ANEU-001 is the dedicated **aneuploidy call from read depth**
(`DetectAneuploidy`) plus whole-chromosome classification (`IdentifyWholeChromosomeAneuploidy`). The two
share nomenclature but ANEU owns the read-depth→CN path — **no overlap conflict**.

## Stage B — implementation

`Math.Pow(2, Math.Log2(x)) == x`, so `copyNumber = round((mean/median)×2)` is exactly the validated
formula with the correct log base and ×2 diploid baseline; the confidence formula matches the spec
exactly; **CN=2 is correctly excluded** from aneuploidy calls. The five worked rows were re-derived
against the code and match the `DetectAneuploidy_*` assertions (Normal/Trisomy/Monosomy/Nullisomy/
Tetrasomy); log-ratio assertions check exact `Math.Log2` values including `double.NegativeInfinity`, and
the clamp tests (10× → CN 10, 0.01× → CN 0) reproduce.

**Variant/delegate consistency (not a defect):** the static `GetAneuploidyTerm` (used by
`AnalyzeKaryotype`) and the inline switch in `IdentifyWholeChromosomeAneuploidy` share the nomenclature;
only a cosmetic label difference for CN>5 ("Polysomy (N)" vs "Copy number = N").

**Test-quality audit:** 31 tests assert exact sourced values (CN, exact `Math.Log2` log-ratios,
confidence 1.0 at boundaries) and cover all Stage-A edge cases — empty, zero/negative median, clamp to
10, very low → 0, single point, binning/averaging, multi-chromosome grouping, ordering, `minFraction`
below/at/custom threshold, CN>5 formatting. No tautological "no-throw" assertions; deterministic.

## Findings

- **No code defect and no test change.** Terminology, clinical examples and the read-depth/log2 model
  were independently confirmed against Wikipedia "Aneuploidy" and the worked example reproduced by hand;
  the code computes the validated formula exactly and the tests lock the sourced values.
- **Accepted pre-documented limitation (not a defect):** sex chromosomes are **not** baseline-adjusted,
  so a normal male (single X / Y) could read as monosomic — see [[aneuploidy-detection]] for the full
  limitation writeup (sex-chromosome and partial-aneuploidy boundaries).
