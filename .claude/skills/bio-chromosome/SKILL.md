---
name: bio-chromosome
description: >-
  Chromosome- and assembly-scale genome analysis with Seqeron (MCP tools OR the
  C# API). Use for karyotype & ploidy / aneuploidy (trisomy, monosomy) from
  depth or chromosome descriptors, centromere / telomere / heterochromatin and
  G-band prediction, arm-ratio classification (metacentric/acrocentric),
  synteny blocks & large structural rearrangements (inversion/translocation/
  deletion/duplication), assembly QC (N50/L50, auN, Nx curve, contiguity, gaps,
  contig extraction, completeness/BUSCO-like, assembly comparison), and
  chromosome-scale composition/GC-skew profiling. Triggers: "analyze this
  karyotype", "detect aneuploidy / ploidy", "find the centromere / telomeres",
  "N50 of this assembly", "assembly stats", "find synteny blocks", "detect
  rearrangements", "GC-skew of this chromosome", "predict the replication
  origin", "sliding-window GC / composition". Server: chromosome (+ analysis for
  GC-skew/origin).
allowed-tools: Read, Bash, Grep, Glob
---

# bio-chromosome — karyotype, centromere/telomere, assembly QC, synteny/SV, composition

Routing + orchestration skill for the **Chromosome** server (32 tools:
`ChromosomeAnalyzer.*` + `GenomeAssemblyAnalyzer.*`). It picks the right tool for a
chromosome- or assembly-scale question and gives a **dual-mode** recipe (MCP tool calls
**and** the equivalent `Seqeron.Genomics` C# `Method ID`).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units /
  0-based (or here often 1-based inclusive) coordinates, and the alpha / not-for-clinical
  caveat are owned by **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server chromosome`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/chromosome/*.md`;
  algorithm invariants in `docs/algorithms/{Chromosome_Analysis,Assembly,Extended_GC_Skew_Analysis}/*.md`.

## Scope boundary — this server vs neighbours

- **GC-skew, cumulative/windowed GC-skew, replication-origin prediction live on the *analysis*
  server**, not chromosome: `gc_skew` / `cumulative_gc_skew` / `windowed_gc_skew` /
  `predict_replication_origin` (`GcSkewCalculator.*`) — route via **[`bio-annotation`](../bio-annotation/SKILL.md)**.
  Composition *here* = G-bands (`predict_g_bands`), local GC/complexity (`local_quality`), suspicious windows (`find_suspicious_regions`).
- **Variant-level SV / CNV** (single breakpoints, per-locus CNV calls, VEP/ACMG effects) →
  **[`bio-annotation`](../bio-annotation/SKILL.md)**. *This* skill = CHROMOSOME-SCALE structure:
  synteny-based rearrangements and depth-based whole-chromosome aneuploidy/ploidy.
  Gene-order-permutation rearrangements, reversal distance, ANI/orthologs/RBH and
  dot-plots → [`seqeron-comparative-genomics`](../seqeron-comparative-genomics/SKILL.md).

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Karyotype from chromosome descriptors → aneuploidy flags | `analyze_karyotype` / `ChromosomeAnalyzer.AnalyzeKaryotype` |
| **Whole-chromosome** aneuploidy from read depth | `identify_whole_chromosome_aneuploidy` / `…IdentifyWholeChromosomeAneuploidy` |
| Per-bin copy number / aneuploidy from depth | `detect_aneuploidy` / `ChromosomeAnalyzer.DetectAneuploidy` |
| Overall ploidy from normalized depth | `detect_ploidy` / `ChromosomeAnalyzer.DetectPloidy` |
| Locate centromere (alpha-satellite scan) | `analyze_centromere` / `ChromosomeAnalyzer.AnalyzeCentromere` |
| Arm ratio from centromere position | `arm_ratio` / `ChromosomeAnalyzer.CalculateArmRatio` |
| Classify centromere position (meta/acro…) | `classify_chromosome_by_arm_ratio` / `…ClassifyChromosomeByArmRatio` |
| Telomere tract presence / critical shortening | `analyze_telomeres` / `ChromosomeAnalyzer.AnalyzeTelomeres` |
| Telomere length from qPCR T/S ratio | `estimate_telomere_length_from_ts_ratio` / `…EstimateTelomereLengthFromTSRatio` |
| Cell divisions from telomere shortening | `estimate_cell_divisions_from_telomere_length` / `…EstimateCellDivisionsFromTelomereLength` |
| Heterochromatin regions (repeat-content windows) | `find_heterochromatin_regions` / `…FindHeterochromatinRegions` |
| GC-based cytogenetic G-bands | `predict_g_bands` / `ChromosomeAnalyzer.PredictGBands` |
| Synteny blocks from ortholog pairs | `find_synteny_blocks` / `ChromosomeAnalyzer.FindSyntenyBlocks` |
| Rearrangements (inv/transloc/del/dup) from blocks | `detect_rearrangements` / `ChromosomeAnalyzer.DetectRearrangements` |
| Synteny blocks between two **assemblies** (k-mer) | `find_syntenic_blocks_assemblies` / `GenomeAssemblyAnalyzer.FindSyntenicBlocks` |
| Assembly stats (N50/L50, GC, gaps) | `assembly_statistics` / `GenomeAssemblyAnalyzer.CalculateStatistics` |
| Nx/Lx at a threshold · full Nx curve · auN | `nx_statistics` · `nx_curve` · `au_n` (`GenomeAssemblyAnalyzer.CalculateNx / CalculateNxCurve / CalculateAuN`) |
| Length distribution of contigs/scaffolds | `length_distribution` / `…CalculateLengthDistribution` |
| Find / distribution of N-gaps | `find_gaps` · `gap_distribution` (`…FindGaps` / `…AnalyzeGapDistribution`) |
| Split scaffolds → contigs | `extract_contigs` · `analyze_scaffolds` (`…ExtractContigs` / `…AnalyzeScaffolds`) |
| Completeness: marker-gene (BUSCO-like) · k-mer spectrum | `assess_completeness` · `estimate_completeness_from_kmers` |
| Local per-window quality (GC/N/complexity) | `local_quality` / `GenomeAssemblyAnalyzer.CalculateLocalQuality` |
| Suspicious windows (GC/complexity outliers) | `find_suspicious_regions` / `…FindSuspiciousRegions` |
| Tandem repeats · repetitive regions · repeat content | `find_tandem_repeats` · `find_repetitive_regions` · `repeat_content` |
| Compare two assemblies (shared k-mers) | `compare_assemblies` / `GenomeAssemblyAnalyzer.CompareAssemblies` |

## Canonical dual-mode pipelines

### (a) Chromosome sequence → centromere → arm ratio → classify
1. **[MCP]** `analyze_centromere`(chromosomeName, sequence, windowSize?=100000, minAlphaSatelliteContent?=0.3) → centromere `start`/`end` (null if none), `type`.
2. **[MCP]** `arm_ratio`(centromerePosition, chromosomeLength) → `p`, `q`, `armRatio`.
3. **[MCP]** `classify_chromosome_by_arm_ratio`(armRatio) → Metacentric / Submetacentric / Subtelocentric / Acrocentric / Telocentric.
- **[C# API]** `ChromosomeAnalyzer.AnalyzeCentromere(...)` → `.CalculateArmRatio(pos,len)` → `.ClassifyChromosomeByArmRatio(r)`.
```
Provenance
1) analyze_centromere(name,seq,windowSize=100000,minAlphaSatelliteContent=0.3) → start,end,type
2) arm_ratio(centromerePosition=start,chromosomeLength=len) → p,q,armRatio
3) classify_chromosome_by_arm_ratio(armRatio) → type
Cross-check: analyze_centromere's own `type` (Levan) should agree with step 3's classification.
Envelope: none of these guarded. Caveat: alpha — validate before decision use.
```

### (b) Read depth → whole-chromosome aneuploidy / ploidy
1. **[MCP]** `detect_aneuploidy`(depthData=[{chromosome,position,depth}], medianDepth, binSize?=1000000) → per-bin `copyNumber`, `logRatio`, `confidence`.
2. **[MCP]** `identify_whole_chromosome_aneuploidy`(depthData, medianDepth, binSize?, minFraction?) → per-chromosome dominant CN ≠ 2 as ISCN aneuploidy.
3. **[MCP]** (genome-wide) `detect_ploidy`(normalizedDepths, expectedDiploidDepth) → overall `ploidy`.
- **[C# API]** `ChromosomeAnalyzer.DetectAneuploidy(...)` → `.IdentifyWholeChromosomeAneuploidy(...)` → `.DetectPloidy(...)`.
- Descriptor-only alternative (no depth): `analyze_karyotype`(chromosomes=[{name,length,isSexChromosome}], expectedPloidyLevel?=2) → `abnormalities`.

### (c) Ortholog pairs → synteny blocks → large structural rearrangements
1. **[MCP]** `find_synteny_blocks`(orthologPairs=[{chr1,start1,end1,gene1,chr2,start2,end2,gene2}], minGenes?=3, maxGap?=10) → collinear blocks + strand.
2. **[MCP]** `detect_rearrangements`(syntenyBlocks) → inversion / translocation / deletion / duplication events + breakpoints.
- **[C# API]** `ChromosomeAnalyzer.FindSyntenyBlocks(...)` → `.DetectRearrangements(blocks)`.
- Assembly-vs-assembly (no genes): `find_syntenic_blocks_assemblies`(assembly1, assembly2, minBlockSize?) → k-mer-anchored blocks with orientation flags.
- **Variant-level breakpoints/CNV** → cross-link **[`bio-annotation`](../bio-annotation/SKILL.md)**.

### (d) Assembly QC — contiguity, gaps, completeness
1. **[MCP]** `assembly_statistics`(sequences=[{id,sequence}]) → N50/L50, N90/L90, GC, gap counts, largest/mean/median.
2. **[MCP]** `au_n`(lengths) → threshold-free contiguity; `nx_curve`(sortedLengths, totalLength, thresholds?) → full N10…N90 curve.
3. **[MCP]** `find_gaps`(sequences, minGapLength?=1) → `[start,end]` inclusive per N-run; `gap_distribution`(gaps) → summary.
4. **[MCP]** completeness: `assess_completeness`(markerGenes, assembly) *(BUSCO-like)* or `estimate_completeness_from_kmers`(kmerSpectrum).
- **[C# API]** `GenomeAssemblyAnalyzer.CalculateStatistics(...)` · `.CalculateAuN(...)` · `.CalculateNxCurve(...)` · `.FindGaps(...)` · `.AnalyzeGapDistribution(...)` · `.AssessCompleteness(...)`.
- **Note:** `nx_statistics`/`nx_curve` expect **descending-sorted lengths + precomputed totalLength**; feed them from `length_distribution`/your own sort, not raw sequences.

### (e) Chromosome-scale composition profile
1. **[MCP]** `predict_g_bands`(sequence, bandSize?, darkBandGcThreshold?, lightBandGcThreshold?) → gpos100/gpos50/gneg bands by windowed GC.
2. **[MCP]** `local_quality`(sequences, windowSize?) → per-window GC / N-count / linguistic complexity (step = windowSize/2).
3. **[MCP]** `find_heterochromatin_regions`(sequence, windowSize?, minRepeatContent?) → merged repeat-rich regions, classed Telomeric/Centromeric/Constitutive.
- **[C# API]** `ChromosomeAnalyzer.PredictGBands(...)` · `GenomeAssemblyAnalyzer.CalculateLocalQuality(...)` · `ChromosomeAnalyzer.FindHeterochromatinRegions(...)`.
- **GC-skew / cumulative skew / replication origin** are on the **analysis** server → use **[`bio-annotation`](../bio-annotation/SKILL.md)**: `gc_skew`, `cumulative_gc_skew`, `windowed_gc_skew`, `predict_replication_origin` (`GcSkewCalculator.*`).

### (f) Telomere biology
1. **[MCP]** `analyze_telomeres`(chromosomeName, sequence, telomereRepeat?=TTAGGG, searchLength?, minTelomereLength?, criticalLength?) → present? / length / purity / critically-short.
2. **[MCP]** `estimate_telomere_length_from_ts_ratio`(tsRatio, referenceRatio, referenceLength) → length from qPCR.
3. **[MCP]** `estimate_cell_divisions_from_telomere_length`(birthLength, currentLength, lossPerDivision) → replicative age proxy.
- **[C# API]** `ChromosomeAnalyzer.AnalyzeTelomeres(...)` · `.EstimateTelomereLengthFromTSRatio(...)` · `.EstimateCellDivisionsFromTelomereLength(...)`.

## Envelope — STOP rule (guarded unit in scope)

- **CHROM-CENT-001** (`AssignSuprachromosomalFamily`, output `Sf1OrSf2Dimeric`) is **Permissive-only**
  (blocked in Strict & Moderate; default is Moderate). Resolving SF1 vs SF2 needs an SF-resolved
  consensus-monomer reference the library does not ship. If a task asks to *assign suprachromosomal
  family / distinguish SF1 vs SF2*, **STOP** and report the envelope — do not force a result. The
  general centromere tools above (`analyze_centromere`, etc.) are **not** guarded.

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** Given a draft assembly (a few `{id,sequence}` scaffolds), QC it and interpret one chromosome:
(1) contiguity + gaps, (2) split its largest scaffold into contigs, (3) locate that chromosome's
centromere and classify it, (4) profile heterochromatin, (5) corroborate contiguity threshold-free.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `assembly_statistics`(sequences) → `n50`, `l50`, `gcContent`, `gapCount`. (`GenomeAssemblyAnalyzer.CalculateStatistics`)
2. `find_gaps`(sequences, minGapLength=10) then `gap_distribution`(gaps) → gap classes. (`…FindGaps` / `…AnalyzeGapDistribution`)
3. `extract_contigs`(scaffolds, minContigLength=500) → `{id}_contig{n}` sequences over gap-free runs. (`…ExtractContigs`)
4. `analyze_centromere`(chromosomeName, sequence) → `start`/`end`/`type`; then `arm_ratio`(start, length) → `armRatio`; then `classify_chromosome_by_arm_ratio`(armRatio). (`ChromosomeAnalyzer.AnalyzeCentromere` → `.CalculateArmRatio` → `.ClassifyChromosomeByArmRatio`)
5. `find_heterochromatin_regions`(sequence) → Centromeric region should overlap the centromere from step 4 (independent repeat-content path). (`…FindHeterochromatinRegions`)
6. Cross-check contiguity: `au_n`(lengths) — threshold-free auN should track N50 direction. (`GenomeAssemblyAnalyzer.CalculateAuN`)

```
Provenance
1) assembly_statistics(sequences) → n50,l50,gcContent,gapCount
2) find_gaps(sequences,minGapLength=10)+gap_distribution → gap classes ([start,end] inclusive)
3) extract_contigs(scaffolds,minContigLength=500) → contigs
4) analyze_centromere→arm_ratio→classify_chromosome_by_arm_ratio → centromere type
5) find_heterochromatin_regions → Centromeric region overlaps step-4 centromere (cross-check)
6) au_n(lengths) → threshold-free contiguity consistent with N50
Envelope: none of the above guarded (CHROM-CENT-001 applies only to SF1/SF2 assignment — not invoked here).
Caveat: alpha software; not for clinical use — independently validate before relying on any call.
```

## Reference

- **Full domain tool index (all 32, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`).
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (all 32 by group, one-liners + Method ID + doc):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`Karyotype_Analysis.md`](../../../docs/algorithms/Chromosome_Analysis/Karyotype_Analysis.md) ·
  [`Centromere_Analysis.md`](../../../docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md) ·
  [`Telomere_Analysis.md`](../../../docs/algorithms/Chromosome_Analysis/Telomere_Analysis.md) ·
  [`Synteny_Analysis.md`](../../../docs/algorithms/Chromosome_Analysis/Synteny_Analysis.md) ·
  [`Aneuploidy_Detection.md`](../../../docs/algorithms/Chromosome_Analysis/Aneuploidy_Detection.md) ·
  [`Assembly_Statistics.md`](../../../docs/algorithms/Assembly/Assembly_Statistics.md) ·
  [`Comprehensive_GC_Analysis.md`](../../../docs/algorithms/Extended_GC_Skew_Analysis/Comprehensive_GC_Analysis.md)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (GC-skew/replication-origin **and** variant-level SV/CNV).
