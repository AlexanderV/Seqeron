# seqeron-epigenetics — fuller pipelines & parameter guidance

Dual-mode recipes for the epigenetics / DNA-methylation family (annotation server,
`EpigeneticsAnalyzer.*`). Rigor (tool-only, provenance, cross-check, alpha caveat) is delegated to
[`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/annotation/<tool>.md`. Tool index:
[`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `find_cpg_islands`: `minLength` = 200 nt · `minGc` = 0.5 · `minCpGRatio` = 0.6 (Gardiner-Garden & Frommer).
- `find_dmrs`: `windowSize` = 1000 bp · `minDifference` = 0.25 (strict `>`, absolute) · `minCpGCount` = 3.
- `predict_imprinted_genes`: `minDifference` = 0.4; `hasDMR` when `|diff| > 0.5`.
- `find_accessible_regions`: `threshold` = 0.5 · `minWidth` = 100 bp · `maxGap` = 50.
- `predict_chromatin_state` / `annotate_histone_modifications`: presence threshold **0.5** (ChromHMM binarization).
- `epigenetic_age`: `intercept` default 0.0; caller supplies the coefficient table (adult.age = 20).
- Coordinates 0-based; island `end` exclusive. Levels/β/signals are fractions 0..1.

## 1. CpG landscape of a sequence (sites → O/E → islands)

**Goal.** Characterise the CpG content and locate islands, with an independent recompute.

1. **[MCP]** `find_cpg_sites`(sequence) → `positions[]` (0-based CG starts).
2. **[MCP]** `cpg_observed_expected`(sequence) → whole-sequence `ratio` (context for whether islands are plausible).
3. **[MCP]** `find_cpg_islands`(sequence, minLength=200, minGc=0.5, minCpGRatio=0.6) → `islands[]`.
4. Cross-check: for each island, `cpg_observed_expected`(island_substring) should reproduce its `cpGRatio`, and the count of `find_cpg_sites` positions inside the island should be consistent with a high O/E.

- **[C# API]** `EpigeneticsAnalyzer.FindCpGSites(...)` · `.CalculateCpGObservedExpected(...)` · `.FindCpGIslands(...)`.

```
Provenance
1) find_cpg_sites(seq) → CpG positions
2) cpg_observed_expected(seq) → whole-seq O/E
3) find_cpg_islands(seq,minLength=200,minGc=0.5,minCpGRatio=0.6) → islands[]
4) cross-check: cpg_observed_expected(island_substr) ≈ island.cpGRatio
Caveat: alpha; island calls are threshold-dependent — record minLength/minGc/minCpGRatio.
```

## 2. Bisulfite reads → methylation calls → profile (+ round-trip oracle)

**Goal.** Call per-CpG methylation from aligned bisulfite reads and summarise it.

1. **[MCP]** `methylation_from_bisulfite`(referenceSequence, bisulfiteReads=[{readSequence,startPosition}])
   → `sites[]` (`position`, `type=CpG`, `context`, `methylationLevel`=meC/(meC+unmeC), `coverage`=valid C/T calls). Zero-coverage CpGs omitted; read bases past the reference ignored.
2. **[MCP]** `methylation_profile`(sites) → weighted `globalMethylation`/`cpGMethylation`/`cHGMethylation`/`cHHMethylation` (Schultz 2012 `Σ(level·cov)/Σcov`), `totalCpGSites`, `methylatedCpGSites` (level≥0.5), `methylationByPosition`.
3. Optional self-consistency oracle: `simulate_bisulfite_conversion`(referenceSequence, methylatedPositions) forward-models the read strand you expect — the calls in step 1 should agree with the protected/unprotected pattern you simulated.

- **[C# API]** `EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(...)` → `.GenerateMethylationProfile(...)`; oracle `.SimulateBisulfiteConversion(...)`.

```
Provenance
1) methylation_from_bisulfite(ref,reads) → per-CpG sites (Bismark C=meth, T=unmeth)
2) methylation_profile(sites) → weighted per-context methylation
3) simulate_bisulfite_conversion(ref,methylatedPositions) → expected converted strand (oracle for step 1)
Caveat: alpha; sequence-only find_methylation_sites has level=coverage=0 — measured levels need reads.
```

## 3. Differentially methylated regions between two conditions

**Goal.** Locate hyper/hypomethylated regions between control and treatment.

1. Build `sites[]` for each condition (pipeline 2, or caller-supplied `{position,type,context,methylationLevel,coverage}`).
2. **[MCP]** `find_dmrs`(sample1=control, sample2=treatment, windowSize=1000, minDifference=0.25, minCpGCount=3)
   → `regions[]` (`start`/`end` = first/last window position, `meanDifference`=mean(s2−s1), Fisher `pValue`, `cpGCount`, `annotation` Hyper/Hypomethylated).
3. Direction check: `annotation`=Hypermethylated ⇔ `meanDifference` > 0 ⇔ treatment > control. **Order matters — pass control as `sample1`.**

- **[C# API]** `EpigeneticsAnalyzer.FindDMRs(...)`.

```
Provenance
1) sites per condition (methylation_from_bisulfite / caller-supplied)
2) find_dmrs(control,treatment,windowSize=1000,minDifference=0.25,minCpGCount=3) → regions[]
3) check: Hypermethylated ⇔ meanDifference>0 (sample2−sample1)
Caveat: alpha; two-sided Fisher p on pooled counts, NOT multiple-testing corrected; window/threshold are model choices.
```

## 4. Allele-specific methylation → imprinted genes

**Goal.** Flag candidate imprinted loci from maternal/paternal methylation.

1. **[MCP]** `predict_imprinted_genes`(genes=[{geneId,start,end,maternalMethylation,paternalMethylation}], minDifference=0.4)
   → `imprinted[]` (`parentalOrigin` Maternal when maternal higher else Paternal, `imprintingScore`=`|diff|/(mat+pat+0.01)` capped 1, `hasDMR` when `|diff|>0.5`).

- **[C# API]** `EpigeneticsAnalyzer.PredictImprintedGenes(...)`.

```
Provenance
predict_imprinted_genes(genes,minDifference=0.4) → imprinted[] (parentalOrigin, imprintingScore, hasDMR)
Caveat: alpha; empty result is threshold-dependent (record minDifference).
```

## 5. Chromatin state & accessibility

**Goal.** Assign chromatin state from histone marks and locate accessible peaks.

1. **[MCP]** `predict_chromatin_state`(h3k4me3,h3k4me1,h3k27ac,h3k36me3,h3k27me3,h3k9me3) → `state`.
   Priority (each mark binarised at 0.5): BivalentPromoter (K4me3+K27me3) → BivalentEnhancer
   (K4me1+K27me3) → ActivePromoter (K4me3) → ActiveEnhancer (K4me1+K27ac) / WeakEnhancer (K4me1) →
   Transcribed (K36me3) → Repressed (K27me3) → Heterochromatin (K9me3) → LowSignal. **Signals 0..1 or it throws.**
2. **[MCP]** `annotate_histone_modifications`(modifications=[{start,end,mark,signal}]) → `annotations[]`
   with per-interval `predictedState` from that single mark (below 0.5 → LowSignal).
3. **[MCP]** `find_accessible_regions`(accessibilitySignal=[{position,signal}], threshold=0.5, minWidth=100, maxGap=50)
   → `regions[]` (`accessibilityScore`=max signal, `peakType` Strong>0.8/Moderate>0.5/Weak, `nearbyGenes` empty).

- **[C# API]** `.PredictChromatinState(...)` · `.AnnotateHistoneModifications(...)` · `.FindAccessibleRegions(...)`.

```
Provenance
1) predict_chromatin_state(6 marks 0..1) → state (ChromHMM priority)
2) annotate_histone_modifications(intervals) → per-interval single-mark state
3) find_accessible_regions(signal,threshold=0.5,minWidth=100,maxGap=50) → peaks (Strong/Moderate/Weak)
Caveat: alpha; 0.5 binarization is a fixed model choice; enhancer/nucleosome standalone tools do not exist.
```

## 6. Epigenetic age (Horvath clock)

**Goal.** Estimate DNAm age from methylation at clock CpGs.

1. **[MCP]** `epigenetic_age`(methylationAtClockCpGs={cpgId→β 0..1}, coefficients={cpgId→coef}, intercept=0.0)
   → `age` (years). Linear predictor over CpGs shared by both maps; missing CpGs contribute nothing.
   Inverse calibration (adult.age=20): `(21)·exp(x)−1` for x<0, `(21)·x+20` for x≥0.

- **[C# API]** `EpigeneticsAnalyzer.CalculateEpigeneticAge(...)`.

```
Provenance
epigenetic_age(β-values, published Horvath coefficients, intercept) → DNAm age (years)
Caveat: alpha; result is only as good as the supplied clock coefficient table (wrapper is clock-agnostic).
```

## Scope reminders

- **Sequence-only vs measured methylation.** `find_methylation_sites` gives contexts only
  (`methylationLevel`/`coverage` = 0); measured levels require `methylation_from_bisulfite` or
  caller-supplied sites.
- **No nucleosome-positioning / standalone enhancer tool** — enhancer states come from
  `predict_chromatin_state`. **Chromosome-scale heterochromatin/G-bands** →
  `find_heterochromatin_regions` on the chromosome server ([`bio-chromosome`](../../bio-chromosome/SKILL.md)).
- **Everything else annotation** (ORFs/genes/promoters, variants, motifs, repeats, k-mers, splicing,
  miRNA, RNA-seq) → [`bio-annotation`](../../bio-annotation/SKILL.md).
```
