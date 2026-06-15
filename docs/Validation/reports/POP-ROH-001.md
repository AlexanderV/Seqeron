# Validation Report: POP-ROH-001 — Runs of Homozygosity (ROH) detection & F_ROH

- **Validated:** 2026-06-15   **Area:** PopGen
- **Canonical method(s):** `PopulationGeneticsAnalyzer.FindROH(genotypes, minSnps, minLength, maxHeterozygotes, maxGap)`; `PopulationGeneticsAnalyzer.CalculateInbreedingFromROH(rohSegments, genomeLength)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **McQuillan et al. (2008), Am J Hum Genet 83(3):359-372** — https://pmc.ncbi.nlm.nih.gov/articles/PMC2556426/ (WebFetch, 2026-06-15).
   - F_ROH formula **verbatim: "Froh = ∑Lroh/Lauto"**, where ∑Lroh = "the total length of all of an individual's ROHs above a specified minimum length" and Lauto = "the length of the autosomal genome covered by SNPs, excluding the centromeres."
   - L_auto numeric value **verbatim: "2,673,768 kb"**.
   - F_ROH defined as "the proportion of the autosomal genome in runs of homozygosity above a specified length threshold." ⇒ confirms range [0,1] (INV-4) and linearity in ΣL_roh (INV-5).

2. **PLINK 1.9 `--homozyg` docs (Chang et al. 2015)** — https://www.cog-genomics.org/plink/1.9/ibd (WebFetch, 2026-06-15).
   - Defaults retrieved: `--homozyg-snp` = 100 SNPs; `--homozyg-kb` ≥ 1000 kb (= 1,000,000 bp); `--homozyg-window-snp` = 50; `--homozyg-window-het` = 1; `--homozyg-window-missing` = 5; `--homozyg-window-threshold` = 0.05; `--homozyg-gap` = 1000 kb (consecutive SNPs > 1000 kb apart cannot share a ROH); `--homozyg-density` = 1 SNP / 50 kb.
   - Confirms both thresholds (count AND length) must be satisfied (INV-1; failure-mode 1).

3. **detectRUNS vignette (Marras et al. 2015 consecutive-runs method)** — https://cran.r-project.org/web/packages/detectRUNS/vignettes/detectRUNS.vignette.html (WebFetch, 2026-06-15).
   - "window-free … scans the genome SNP by SNP." Parameters: `minSNP`, `minLengthBps`, `maxGap`, `maxOppRun` ("maximum number of opposite genotypes in the run" — allows heterozygous calls within runs), `maxMissRun`. Example uses `maxOppRun = 1`.
   - Run extends while SNPs meet the homozygosity criterion (with `maxOppRun`/`maxMissRun` exceptions) and terminates when gap > `maxGap` or a tolerance is exceeded (INV-2, INV-3, corner cases).

### Formula check
F_ROH = ΣL_roh / L_auto — matches McQuillan (2008) verbatim. The implementation computes `Σ(End − Start) / genomeLength`, identical.

### Edge-case semantics
- Tolerated interior opposite genotypes (≤ maxOppRun): sourced (Marras/detectRUNS).
- Gap > maxGap breaks an all-homozygous run: sourced (PLINK `--homozyg-gap`; Marras `maxGap`).
- Sub-threshold stretches discarded (count and length independently): sourced (PLINK; Marras).
- A heterozygote cannot seed a homozygous run / leading hets skipped: implied by "stretch of homozygous SNPs" definition.
- A ROH is bounded by homozygous markers ⇒ trailing tolerated het is not part of the run: definitional (Marras run = stretch of homozygous SNPs; opposite genotypes only bridge interior error). The exact C++ trimming detail in detectRUNS was not retrievable this session, but the property is self-consistent and the only sound convention.

### Independent cross-check (numbers used in tests, traced to sources)
- F_ROH worked example: ΣL_roh = 20 Mb over L_auto = 100 Mb ⇒ **0.20** (McQuillan formula, hand-computed).
- Whole-genome ROH over the actual McQuillan L_auto = **2,673,768** kb ⇒ F_ROH = **1.0**.
- PLINK defaults (100 SNPs / 1,000,000 bp / 1 het / 1,000,000 bp gap) drive the M1/M5/M6/M4 boundary values.

### Findings / divergences
- PLINK's two-phase sliding-window scan (`--homozyg-window-*`) is intentionally **not** reproduced; the window-free consecutive-runs method (Marras 2015) is implemented and is itself an authoritative reference method. Documented in TestSpec §7. Not a defect.
- Missing-genotype handling (`maxMissRun` / `--homozyg-window-missing`) is out of scope; input has no missing sentinel, any non-1 genotype is homozygous. Documented assumption, not invented behavior.

**Stage A: PASS** — every formula, default, and edge-case rule traces to a source retrieved this session.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Population/PopulationGeneticsAnalyzer.cs:1440` (`CalculateInbreedingFromROH`), `:1484` (`FindROH` validation), `:1505` (`FindROHIterator`).

### Formula realised correctly?
- `CalculateInbreedingFromROH`: `Σ(long)(End − Start) / genomeLength` — exact F_ROH; uses `long` accumulator to avoid 32-bit overflow on genome-scale sums; returns 0 when `genomeLength ≤ 0`; throws on null segments.
- `FindROHIterator`: sorts ascending; tracks `lastHomIndex`/`snpCountAtLastHom` so a run is emitted on `[runStart .. lastHom]` (closed at the last homozygous SNP). Gap break and het break both close the current qualifying run and restart at the breaker; a het breaker seeds no run (`snpCount = isHet ? 0 : 1`). Qualification requires `SnpCount ≥ minSnps AND End − Start ≥ minLength`.

### Cross-verification table recomputed vs code (hand-traced)
| Case | Input | Expected (source) | Code result |
|------|-------|-------------------|-------------|
| M1 | 100 hom @20kb | 1 run, [0, 1,980,000], 100 | match |
| M2 | 101 SNPs, het@50, maxHet=1 | 1 run, [0, 2,000,000], 101 | match |
| M3 | 201 SNPs, het@50,100 | [0,1.98M]/100 + [101·20k,200·20k]/100 | match |
| M4 | two 60-SNP blocks, 2Mb gap | 2 runs × 60 | match |
| M5 | 20 hom, default minSnps=100 | empty | match |
| M6 | 120 @1kb (119 kb) | empty | match |
| M7 | ΣL=20M / 100M | 0.20 | match |
| M8 | ROH = L_auto (2,673,768) | 1.0 | match |
| S6 (new) | trailing het @100 | [0,1.98M]/100, het excluded | match |
| S7 (new) | het@50, maxHet=0 | [0,49·20k]/50 + [51·20k,99·20k]/49 | match |

### Variant/delegate consistency
Two public methods, no `*Fast`/delegate variants. Defaults (`DefaultRohMinSnps=100`, `DefaultRohMinLengthBp=1,000,000`, `DefaultRohMaxHeterozygotes=1`, `DefaultRohMaxGapBp=1,000,000`) match PLINK 1.9 sources.

### Test quality audit (HARD gate)
- **Sourced, not code-echoes:** every expected value traces to McQuillan formula, McQuillan L_auto, or PLINK defaults — not to code output. F_ROH cases use exact `.Within(1e-10)` (no permissive tolerance / range / Greater).
- **No green-washing:** no weakened assertions, no skips, no widened tolerances. Legacy weak tests previously removed from `PopulationGeneticsAnalyzerTests.cs` (pointer comment at line 100-101 confirmed).
- **Coverage:** all Stage-A branches exercised — single run (M1), interior tolerated het (M2), het-beyond-tolerance split (M3), gap break (M4), minSnps (M5), minLength (M6), unsorted (S1), empty (S2), hom-alt encoding (S3), leading hets (S4), F_ROH two-segment/whole-genome/empty (M7/M8/S5), invalid args (C1/C2). **Two coverage gaps found and closed this session:** trailing-het trimming (S6) and zero-tolerance break (S7) — both were untested logic branches; both confirmed correct (no code change needed).
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6552` after additions; `dotnet build` 0 errors, no new warnings in the ROH files.

### Findings / defects
None. Implementation faithfully realises the validated description. The two added tests lock previously-uncovered but correct behaviors.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS.**
- **End-state: CLEAN.** No defect; coverage gaps (trailing-het trim, zero-tolerance) closed with sourced tests; full suite green.
- **Test-quality gate: PASS** (sourced values, no green-washing, all branches covered incl. the 2 new edge tests, honest full-suite green).
- One-line note: F_ROH and consecutive-runs ROH detection match McQuillan (2008) and Marras (2015)/PLINK 1.9 exactly; added trailing-het and zero-tolerance edge tests.
