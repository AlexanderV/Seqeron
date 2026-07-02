# seqeron-epigenetics — tool map (all 13)

Server: **annotation**. One backing class: `EpigeneticsAnalyzer.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/annotation/<tool>.md`.

> Coordinates are **0-based**; CpG-island `end` is **exclusive**. Methylation levels & β-values are
> fractions **0..1**; histone-mark & accessibility signals are **0..1** (`predict_chromatin_state`
> throws on out-of-range). Always confirm exact I/O in the tool doc.

## CpG detection / composition

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_cpg_sites` | `EpigeneticsAnalyzer.FindCpGSites` | 0-based start index of every `CG` dinucleotide (case-insensitive). | `find_cpg_sites.md` |
| `cpg_observed_expected` | `EpigeneticsAnalyzer.CalculateCpGObservedExpected` | CpG O/E = observed / ((C·G)/len); `0` when len<2 or no C/G. | `cpg_observed_expected.md` |
| `find_cpg_islands` | `EpigeneticsAnalyzer.FindCpGIslands` | Gardiner-Garden & Frommer (1987): sliding `minLength` windows with GC ≥ `minGc` **and** O/E ≥ `minCpGRatio`, merged; merged run must still meet both. Defaults 200 / 0.5 / 0.6. | `find_cpg_islands.md` |

## Methylation contexts & calling

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_methylation_sites` | `EpigeneticsAnalyzer.FindMethylationSites` | Classify each C into CpG / CHG / CHH (H=A/C/T) with up-to-3-base context. Sequence-only → `methylationLevel`/`coverage` = 0. | `find_methylation_sites.md` |
| `simulate_bisulfite_conversion` | `EpigeneticsAnalyzer.SimulateBisulfiteConversion` | Frommer 1992: unprotected C→T (case preserved); cytosines in `methylatedPositions` (0-based) survive; other bases unchanged; strand-specific; same length out. | `simulate_bisulfite_conversion.md` |
| `methylation_from_bisulfite` | `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite` | Bismark (Krueger & Andrews 2011): at each reference CpG, read `C`=methylated, `T`=unmethylated, else ignored; `level`=meC/(meC+unmeC), `coverage`=valid C/T calls; zero-coverage CpGs omitted. | `methylation_from_bisulfite.md` |
| `methylation_profile` | `EpigeneticsAnalyzer.GenerateMethylationProfile` | Schultz 2012 weighted level `Σ(level·cov)/Σcov` per context (global/CpG/CHG/CHH); `totalCpGSites`, `methylatedCpGSites` (level≥0.5), `methylationByPosition` sorted. | `methylation_profile.md` |

## Differential / comparative methylation

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_dmrs` | `EpigeneticsAnalyzer.FindDMRs` | methylKit tiling (Akalin 2012): fixed `windowSize` windows with ≥ `minCpGCount` CpGs and `|mean(s2−s1)| > minDifference` (strict); two-sided Fisher p on pooled counts; `Hypermethylated` (s2>s1) / `Hypomethylated`. Defaults 1000 / 0.25 / 3. **`sample1`=control, `sample2`=treatment.** | `find_dmrs.md` |
| `predict_imprinted_genes` | `EpigeneticsAnalyzer.PredictImprintedGenes` | Report a gene when `|maternal−paternal| ≥ minDifference` (default 0.4); `parentalOrigin` Maternal/Paternal, `imprintingScore` = `|diff|/(mat+pat+0.01)` capped 1, `hasDMR` when `|diff|>0.5`. | `predict_imprinted_genes.md` |

## Epigenetic age

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `epigenetic_age` | `EpigeneticsAnalyzer.CalculateEpigeneticAge` | Horvath 2013 clock: linear predictor `intercept + Σ coef·β` over shared CpGs → years via inverse calibration (adult.age=20; exp branch for x<0, linear for x≥0). Caller supplies the coefficient table (clock-agnostic wrapper). | `epigenetic_age.md` |

## Chromatin state / accessibility

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `predict_chromatin_state` | `EpigeneticsAnalyzer.PredictChromatinState` | ChromHMM-style (Ernst & Kellis 2012) over six marks binarised at 0.5, priority: BivalentPromoter → BivalentEnhancer → ActivePromoter → Active/WeakEnhancer → Transcribed → Repressed → Heterochromatin → LowSignal. Signals **0..1** or throws. | `predict_chromatin_state.md` |
| `annotate_histone_modifications` | `EpigeneticsAnalyzer.AnnotateHistoneModifications` | Per-interval `{start,end,mark,signal}` → canonical Roadmap state for one mark (H3K4me3→ActivePromoter, H3K4me1→WeakEnhancer, H3K27ac→ActiveEnhancer, H3K36me3→Transcribed, H3K27me3→Repressed, H3K9me3→Heterochromatin, H3K9ac→ActivePromoter); below 0.5 → LowSignal. | `annotate_histone_modifications.md` |
| `find_accessible_regions` | `EpigeneticsAnalyzer.FindAccessibleRegions` | Merge position-sorted signal ≥ `threshold`, split at gaps > `maxGap`, keep spans ≥ `minWidth`; `accessibilityScore`=max signal, `peakType` Strong(>0.8)/Moderate(>0.5)/Weak; `nearbyGenes` empty here. Defaults 0.5 / 100 / 50. | `find_accessible_regions.md` |

## Envelope

- **No guarded unit.** No `EpigeneticsAnalyzer.*` entry in
  [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) (verified 2026-07) —
  nothing throws a `MinimumMode` guard, so there is no STOP rule. Results are model/heuristic output
  under the stated thresholds; keep the alpha / not-for-clinical caveat (see [`../SKILL.md`](../SKILL.md)).

## Not a tool here (route elsewhere)

- **Chromosome-scale heterochromatin blocks / G-bands** → `find_heterochromatin_regions` (chromosome
  server) via [`bio-chromosome`](../../bio-chromosome/SKILL.md).
- **Nucleosome positioning / standalone enhancer caller:** **no tool** — enhancer *states* come out of
  `predict_chromatin_state`; do not fabricate one.
- **ORFs/genes/promoters, variants, motifs, repeats, k-mers, splicing, miRNA, RNA-seq** →
  [`bio-annotation`](../../bio-annotation/SKILL.md).
