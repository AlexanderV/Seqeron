# Validation Report: ONCO-MRD-001 — Minimal/Molecular Residual Disease (MRD) Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectMRD`, `OncologyAnalyzer.TrackVariantsOverTime`, `OncologyAnalyzer.IsVariantDetected` (reuses `CtDnaDetectionProbability`, ONCO-CTDNA-001)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Chen et al. (2021), "Commercial ctDNA Assays for Minimal Residual Disease"** (review PDF, learn.colontown.org) — independent review, NOT the repo's cited PMC9265001. Confirms verbatim: Signatera "selects up to 16 somatic SNVs"; a plasma sample is ctDNA-positive/MRD-positive when "at least 2 of the 16 variants are detected"; ~1 million reads/sample (ultra-deep, ~100,000×).
- **WebSearch — Signatera positivity rule** (Natera + multiple peer-reviewed sources): "Samples with at least two tumor-specific variants are defined as ctDNA-positive … those with fewer than 2 variants are considered negative." Top 16 somatic variants selected from tumour WES with matched normal/buffy-coat to remove germline + CHIP.
- **Reinert et al. (2019), JAMA Oncol 5(8):1124–1131** (PubMed 31070691) — primary tumour-informed multiplex-PCR MRD study (assay later commercialised as Signatera); serial ctDNA detects relapse ahead of imaging.
- **Wan et al. (2020), Sci Transl Med 12(548):eaaz8084 — INVAR** (search of GenomeWeb/EurekAlert/Science abstract): IMAF = "background-subtracted, depth-weighted mean allele fraction across patient-specific tumor-mutated loci"; IMAF as low as 0.0011 at ~97% specificity.

### Formula check
| Claim | Source | Status |
|-------|--------|--------|
| MRD-positive ⟺ ≥2 of tracked variants detected (τ=2 default) | Chen 2021; multiple independent | Confirmed verbatim |
| Panel size up to 16 SNVs | Chen 2021; Natera | Confirmed |
| IMAF = depth-weighted mean AF over loci | Wan 2020 | Confirmed (with note below) |
| Panel Poisson p = 1 − e^(−n·f·m) | Standard Poisson detection model; Avanzini 2020 | Confirmed |

### Edge-case semantics
- Exactly 1 detected ⇒ negative; 0 detected ⇒ negative — confirmed ("fewer than 2 … negative").
- Empty panel ⇒ invalid input — reasonable (nothing to interrogate).
- All total reads = 0 ⇒ IMAF = 0, p = 0 — consistent with definition.

### Independent cross-check (numbers)
- **M6 IMAF** = (3+1+0)/(200+150+180) = 4/530 = **0.007547169811320755** (Python recompute — matches).
- **M7 Poisson** p = 1 − e^(−16) = **0.9999998874648253** (Python recompute, IEEE-754 double — matches test).

### Findings / divergences (Stage A → PASS-WITH-NOTES)
1. **IMAF simplification (documented).** Wan's true IMAF is *background-subtracted* and *weighted by tumour AF / fragment-size likelihood*. The implementation computes a plain read-pooled Σalt/Σtotal depth-weighted mean. This is correctly labelled as an intentional simplification in `MRD_Detection.md §5.3` and the Evidence assumption register. It remains a depth-weighted mean AF; acceptable as a documented divergence, not a defect.
2. **Evidence-doc numeric typo (fixed this session).** The Evidence dataset table listed the m=16 Poisson value as `0.9999998874648379`; the correct IEEE-754 double (and the test's value) is `0.9999998874648253`. Corrected in `docs/Evidence/ONCO-MRD-001-Evidence.md`. The test was already correct (it asserts `1.0 - Math.Exp(-16.0)`).

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`
- `DetectMRD` L5734–5790; `TrackVariantsOverTime` L5805–5830; `IsVariantDetected` L5703–5712; `CtDnaDetectionProbability` L5423–5447.

### Formula realised correctly?
- Single pass: counts `trackedCount`, `detectedCount` (via `IsVariantDetected`: alt ≥ r_min), `altReadSum`, `totalReadSum`.
- IMAF = `altReadSum / totalReadSum` (0 when Σtotal=0). Negative read counts clamped to 0 via `Math.Max`.
- Poisson delegated to `CtDnaDetectionProbability(n, imaf, m)` ⇒ p = 1 − e^(−n·imaf·m). No duplicated formula.
- Status = Positive ⟺ `detectedCount ≥ positivityThreshold`. Default τ=2 (`DefaultMrdPositivityThreshold`). Matches sourced rule.
- Longitudinal: per-timepoint `DetectMRD`, preserves order, records earliest positive index (−1 if none). Correct.

### Cross-verification table recomputed vs code
| Case | Expected (source) | Code result | Match |
|------|-------------------|-------------|-------|
| 2 of 3 detected | Positive, D=2 | Positive, D=2 | ✅ |
| 1 of 3 | Negative, D=1 | Negative, D=1 | ✅ |
| 0 of 3 | Negative, D=0, IMAF=0 | same | ✅ |
| 16 markers, 2 det | Positive, tracked=16 | same | ✅ |
| IMAF (3,200)(1,150)(0,180) | 4/530 | 0.0075471698… | ✅ |
| Poisson n=1000 f=0.001 m=16 | 1−e^−16 = 0.9999998874648253 | same | ✅ |
| n=0 default | p=0 | p=0 | ✅ |
| Longitudinal [0,1,2,3] | FirstPositive=2 | 2 | ✅ |

### Variant/delegate consistency
`IsVariantDetected` is the per-locus primitive reused by `DetectMRD`/`TrackVariantsOverTime`; default r_min=1 consistent across all three. Poisson reuse of `CtDnaDetectionProbability` is exact.

### Test quality audit (HARD gate)
Pre-existing suite (16 tests) used exact sourced values throughout (no Greater/AtLeast/range green-washing, deterministic). **Gaps found and fixed this session — 7 tests added (now 23):**
- `IsVariantDetected` (public canonical method) was untested directly → added 3 tests: default-cutoff boundary (alt=1 detected, 0 not), custom-cutoff boundary (alt=3 detected, 2 not), invalid `minSupportingReads` throws.
- `DetectMRD` validation gaps: `minSupportingReads < 1` and `genomeEquivalents < 0` (both `ArgumentOutOfRangeException`) were untested → added C5, C6.
- Documented edge cases untested: default `genomeEquivalents=0 ⇒ p=0` (M7b) and "all total reads = 0 ⇒ IMAF=0, p=0" (M3b) → added.

All assertions trace to the description/source values (alt≥r_min rule from Wan 2020; n=0⇒λ=0⇒p=0 from the Poisson model), not to code echoes. Build warning-free for the changed file.

### Findings / defects
- **No code defect.** All formulas and edge cases realise the validated description.
- IMAF simplification is documented (Stage A note 1), not a code bug.

## Verdict & follow-ups
- **Stage A: PASS-WITH-NOTES** — biology/maths confirmed against independent sources; two documentation notes (IMAF simplification already disclosed; Evidence numeric typo fixed).
- **Stage B: PASS** — implementation faithful; test coverage extended to all public methods + documented edge/error cases.
- **End-state: CLEAN** — Evidence typo corrected; 7 coverage-gap tests added; `dotnet build` 0 errors / no new warnings; full unfiltered suite **6674 passed, 0 failed** (1 pre-existing skipped benchmark).
- **Test-quality gate: PASS** (after adding the 7 tests; previously LIMITED on public-method/edge coverage).
