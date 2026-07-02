---
name: seqeron-epigenetics
description: >-
  Analyze DNA methylation and chromatin state with Seqeron (MCP tools OR the C#
  API). Use to find CpG sites / CpG islands (Gardiner-Garden & Frommer) and CpG
  observed/expected ratio, classify CpG/CHG/CHH methylation contexts, call
  per-CpG methylation from bisulfite sequencing reads (Bismark-style) or
  forward-simulate bisulfite conversion, aggregate a weighted methylation
  profile, find differentially methylated regions (DMRs) between two samples
  (methylKit tiling windows + Fisher exact), estimate epigenetic age (Horvath
  DNA-methylation clock), predict chromatin state / annotate histone
  modifications (ChromHMM-style H3K4me3/H3K4me1/H3K27ac/H3K36me3/H3K27me3/H3K9me3),
  find accessible (ATAC-seq-like) chromatin regions, and predict imprinted genes
  from allele-specific methylation. Triggers: "find CpG islands", "CpG O/E
  ratio", "call methylation from bisulfite", "bisulfite-convert this sequence",
  "methylation profile", "find DMRs between groups", "hyper/hypomethylated
  regions", "epigenetic age / DNAm age / methylation clock", "predict chromatin
  state", "annotate histone marks", "accessible chromatin peaks", "imprinted
  genes". Input is DNA (methylation β-values / read piles for the calling tools).
  Server: annotation (all EpigeneticsAnalyzer.*).
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-epigenetics — methylation, bisulfite, DMRs, chromatin, epigenetic age

Routing + orchestration skill for the epigenetics / DNA-methylation family on the **annotation**
server (13 tools, one backing class: `EpigeneticsAnalyzer.*`). It picks the right tool for a
CpG / methylation / bisulfite / DMR / chromatin / epigenetic-age question and gives a **dual-mode**
recipe (MCP tool name **and** the equivalent `Seqeron.Genomics` C# `Method ID`).

- **Rigor is delegated.** Tool-only computation, provenance, envelope, cross-check, units, 0-based
  coords, and the alpha / not-for-clinical caveat are owned by
  **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server annotation`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/annotation/*.md`;
  algorithm invariants/formulas in `docs/algorithms/Epigenetics/*.md`.

## Scope boundary — this family vs neighbours

- **This skill OWNS** the DNA-methylation / chromatin family: CpG sites/islands/O-E, methylation
  contexts (CpG/CHG/CHH), bisulfite calling + simulation, methylation profiles, DMRs, epigenetic
  age, chromatin-state / histone-mark annotation, accessible regions, imprinted genes.
  **[`bio-annotation`](../bio-annotation/SKILL.md)** name-drops this family shallowly — route here.
- **Everything else annotation** (ORFs/genes/promoters, variant calling & effect, motifs/PWMs,
  repeats/low-complexity, k-mers, splicing, miRNA, RNA-seq) stays in
  **[`bio-annotation`](../bio-annotation/SKILL.md)**.
- **Heterochromatin blocks / G-bands at chromosome scale** → `find_heterochromatin_regions`
  (chromosome server) via **[`bio-chromosome`](../bio-chromosome/SKILL.md)** — distinct from the
  per-locus `predict_chromatin_state`/`Heterochromatin` label here.
- **No nucleosome-positioning and no enhancer-caller tool exists.** Enhancer *states*
  (Active/Weak/Bivalent enhancer) fall out of `predict_chromatin_state`; there is no
  standalone nucleosome/enhancer tool — do not fabricate one.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Positions of every CpG dinucleotide | `find_cpg_sites` / `EpigeneticsAnalyzer.FindCpGSites` |
| CpG observed/expected ratio for one sequence | `cpg_observed_expected` / `EpigeneticsAnalyzer.CalculateCpGObservedExpected` |
| CpG islands (GC ≥ minGc **and** O/E ≥ minCpGRatio, merged) | `find_cpg_islands` / `EpigeneticsAnalyzer.FindCpGIslands` |
| Classify C into CpG / CHG / CHH context (sequence-only) | `find_methylation_sites` / `EpigeneticsAnalyzer.FindMethylationSites` |
| Forward-simulate bisulfite conversion of a strand | `simulate_bisulfite_conversion` / `EpigeneticsAnalyzer.SimulateBisulfiteConversion` |
| Call per-CpG methylation level from bisulfite reads | `methylation_from_bisulfite` / `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite` |
| Aggregate sites → weighted global/CpG/CHG/CHH profile | `methylation_profile` / `EpigeneticsAnalyzer.GenerateMethylationProfile` |
| DMRs between two samples (hyper/hypomethylated) | `find_dmrs` / `EpigeneticsAnalyzer.FindDMRs` |
| Epigenetic (DNAm) age from clock CpGs | `epigenetic_age` / `EpigeneticsAnalyzer.CalculateEpigeneticAge` |
| Chromatin state from six histone marks | `predict_chromatin_state` / `EpigeneticsAnalyzer.PredictChromatinState` |
| Per-interval state from a single histone mark | `annotate_histone_modifications` / `EpigeneticsAnalyzer.AnnotateHistoneModifications` |
| Accessible (ATAC-seq-like) chromatin peaks | `find_accessible_regions` / `EpigeneticsAnalyzer.FindAccessibleRegions` |
| Imprinted genes from allele-specific methylation | `predict_imprinted_genes` / `EpigeneticsAnalyzer.PredictImprintedGenes` |

> Coordinates 0-based; island `end` is exclusive. Methylation levels & β-values are fractions
> 0..1; histone/accessibility signals 0..1 (out-of-range signal throws for `predict_chromatin_state`).
> Confirm exact I/O per tool doc.

## Canonical dual-mode pipelines

### (a) CpG islands + independent O/E cross-check
1. **[MCP]** `find_cpg_islands`(sequence, minLength?=200, minGc?=0.5, minCpGRatio?=0.6) → `islands[]` (0-based start, exclusive end, `gcContent`, `cpGRatio`).
2. **[MCP]** for a reported island, re-derive its ratio independently: `cpg_observed_expected`(sequence=island_substring) → `ratio` should match `cpGRatio`.
- **[C# API]** `EpigeneticsAnalyzer.FindCpGIslands(...)` → `.CalculateCpGObservedExpected(...)`.
```
Provenance
1) find_cpg_islands(sequence,minLength=200,minGc=0.5,minCpGRatio=0.6) → islands[]
2) cpg_observed_expected(island_substr) → ratio ≈ island.cpGRatio (independent recompute)
Caveat: alpha — validate before use.
```

### (b) Bisulfite reads → per-CpG methylation → profile
1. **[MCP]** `methylation_from_bisulfite`(referenceSequence, bisulfiteReads=[{readSequence,startPosition}]) → `sites[]` (`position`, `type=CpG`, `methylationLevel`, `coverage`); uncovered CpGs omitted.
2. **[MCP]** `methylation_profile`(sites) → weighted `globalMethylation`/`cpGMethylation`/`cHGMethylation`/`cHHMethylation`, `totalCpGSites`, `methylatedCpGSites`, `methylationByPosition`.
- Optional round-trip check: `simulate_bisulfite_conversion`(referenceSequence, methylatedPositions) forward-models the reads you then call — a self-consistency oracle.
- **[C# API]** `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(...)` → `.GenerateMethylationProfile(...)`; `.SimulateBisulfiteConversion(...)`.
```
Provenance
1) methylation_from_bisulfite(ref,reads) → per-CpG sites (level=meC/(meC+unmeC), coverage=C/T calls)
2) methylation_profile(sites) → weighted per-context levels (Schultz 2012)
Cross-check: simulate_bisulfite_conversion(ref,methylatedPositions) reproduces the read pile.
Caveat: alpha.
```

### (c) DMRs between two samples
1. Build each sample's `sites[]` (from pipeline (b), or caller-supplied).
2. **[MCP]** `find_dmrs`(sample1=control, sample2=treatment, windowSize?=1000, minDifference?=0.25, minCpGCount?=3) → `regions[]` (`meanDifference`=mean(s2−s1), Fisher `pValue`, `cpGCount`, `annotation`=Hyper/Hypomethylated).
- **[C# API]** `EpigeneticsAnalyzer.FindDMRs(...)`. Direction is `sample2 − sample1` — **pass control first**. p is two-sided Fisher on pooled counts, not multiple-testing corrected. (Provenance block → `reference/pipelines.md` §3.)

### (d) Chromatin state from histone marks (+ accessibility)
1. **[MCP]** `predict_chromatin_state`(h3k4me3,h3k4me1,h3k27ac,h3k36me3,h3k27me3,h3k9me3) → `state` (ChromHMM-priority: Bivalent* → Active* → Transcribed → Repressed → Heterochromatin → LowSignal). Signals must be 0..1.
2. **[MCP]** per-interval single-mark labels: `annotate_histone_modifications`(modifications=[{start,end,mark,signal}]) → `annotations[]` with `predictedState`.
3. **[MCP]** accessible peaks: `find_accessible_regions`(accessibilitySignal=[{position,signal}], threshold?=0.5, minWidth?=100, maxGap?=50) → `regions[]` (`accessibilityScore`, `peakType` Strong/Moderate/Weak).
- **[C# API]** `.PredictChromatinState(...)` · `.AnnotateHistoneModifications(...)` · `.FindAccessibleRegions(...)`.

### (e) Epigenetic age (Horvath clock)
1. **[MCP]** `epigenetic_age`(methylationAtClockCpGs={cpgId→β}, coefficients={cpgId→coef}, intercept?=0.0) → `age` (years, Horvath inverse calibration, adult.age=20). CpGs missing from either side contribute nothing.
- **[C# API]** `EpigeneticsAnalyzer.CalculateEpigeneticAge(...)`. You supply the published clock coefficient table (wrapper is clock-agnostic) — the result is only as good as that table. (Provenance block → `reference/pipelines.md` §6.)

## Envelope

No epigenetics unit is a guarded unit in [`docs/Validation/LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (no `EpigeneticsAnalyzer.*` `MinimumMode` guard; verified 2026-07). So there is **no STOP rule** specific to this family — but the algorithms are model/heuristic-based (Gardiner-Garden & Frommer thresholds, ChromHMM present/absent binarization at 0.5, methylKit tiling, Horvath calibration). Report results as *model output under the stated parameters/thresholds*, not ground truth, and keep the alpha / not-for-clinical caveat from `bio-rigor`. Rule-out results (empty `islands`/`regions`/`imprinted`) are threshold-dependent — record the parameters used.

## End-to-end grounded example

**Task.** From a reference and a pile of bisulfite reads for two conditions, characterise methylation
and find DMRs, cross-checking each layer.
1. `find_cpg_sites`(referenceSequence) → CpG positions (the set DMR windows are built over). (`.FindCpGSites`)
2. `methylation_from_bisulfite`(referenceSequence, controlReads) and again with `treatmentReads` → two `sites[]`. (`.CalculateMethylationFromBisulfite`)
3. `methylation_profile`(each sites[]) → per-condition weighted CpG methylation (sanity summary). (`.GenerateMethylationProfile`)
4. `find_dmrs`(control_sites, treatment_sites, windowSize=1000, minDifference=0.25, minCpGCount=3) → `regions[]`. (`.FindDMRs`)
5. Cross-check: each DMR's `cpGCount` and window bounds must fall within the CpG positions from step 1, and `meanDifference` sign must agree with the two profile levels from step 3.
```
Provenance
1) find_cpg_sites(ref) → CpG positions
2) methylation_from_bisulfite(ref, controlReads|treatmentReads) → sites per condition
3) methylation_profile(sites) → weighted per-context methylation (cross-check magnitudes)
4) find_dmrs(control,treatment,windowSize=1000,minDifference=0.25,minCpGCount=3) → Hyper/Hypomethylated regions + Fisher p
5) Cross-check: DMR windows ⊆ step-1 CpGs; meanDifference sign agrees with step-3 profiles
Caveat: alpha; thresholds/params are model choices — Fisher p not multiple-testing corrected. Not for clinical use.
```

## Reference

- **This family's tool map (all 13 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`CpG_Site_Detection.md`](../../../docs/algorithms/Epigenetics/CpG_Site_Detection.md) ·
  [`Methylation_Analysis.md`](../../../docs/algorithms/Epigenetics/Methylation_Analysis.md) ·
  [`Bisulfite_Sequencing_Analysis.md`](../../../docs/algorithms/Epigenetics/Bisulfite_Sequencing_Analysis.md) ·
  [`Differentially_Methylated_Regions.md`](../../../docs/algorithms/Epigenetics/Differentially_Methylated_Regions.md) ·
  [`Chromatin_State_Prediction.md`](../../../docs/algorithms/Epigenetics/Chromatin_State_Prediction.md) ·
  [`Epigenetic_Age_Estimation.md`](../../../docs/algorithms/Epigenetics/Epigenetic_Age_Estimation.md)
- **Operating envelope:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (no epigenetics guard).
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (overlap: it name-drops methylation/epigenetics; all other annotation tools live there).
