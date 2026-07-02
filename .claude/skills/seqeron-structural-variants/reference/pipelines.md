# seqeron-structural-variants — fuller pipelines & parameter guidance

Dual-mode recipes for the germline read-evidence SV/CNV family (annotation server,
`StructuralVariantAnalyzer.*`). Rigor (tool-only, provenance, cross-check, alpha caveat) is
delegated to [`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/annotation/<tool>.md`.
Tool index: [`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `find_discordant_pairs`: `expectedInsertSize` = 400, `insertSizeStdDev` = 50, `maxInsertSize` = 10000
  (discordant if interchromosomal, non-FR, `|insert − mean| > 3·sd`, or `> maxInsertSize`).
- `cluster_discordant_pairs`: `clusterDistance` = 500, `minSupport` = 3.
- `find_split_reads`: `minClipLength` = 20 · `cluster_split_reads`: `clusterDistance` = 10, `minSupport` = 2.
- `assemble_breakpoint_sequence`: `minOverlap` = 10 (reserved; longest-clip heuristic today).
- `find_microhomology`: `maxLength` = 20 (suffix-of-left ∩ prefix-of-right).
- `segment_copy_number`: `changeThreshold` = 0.3, `minProbes` = 5 · `identify_cnvs`: `normalCopyNumber` = 2, `minLength` = 10000.
- `filter_svs`: `minQuality` = 20, `minSupport` = 2, `minLength` = 50, `maxLength` = 1e8 · `merge_overlapping_svs`: `overlapFraction` = 0.5.
- `genotype_sv`: altFraction thresholds 0.1 / 0.3–0.7 / 0.9; quality capped 99; `totalReads=0` → `./.`.

## 1. Paired-end SV calling (discordant signatures → clustered SVs)

**Goal.** Turn a paired-end library into typed SV candidates.

1. **[MCP]** `find_discordant_pairs`(readPairs, expectedInsertSize=400, insertSizeStdDev=50, maxInsertSize=10000)
   → discordant `pairs[]` (each `isDiscordant=true`): interchromosomal → translocation; large span → deletion;
   small span → insertion; same-strand (non-FR) → inversion/duplication.
2. **[MCP]** `cluster_discordant_pairs`(discordantPairs, clusterDistance=500, minSupport=3)
   → SV `variants[]` (id, chromosome, start/end, type, length, quality, supportingReads). SV type from signature.
3. Cross-check: a real event should also leave a split-read breakpoint (pipeline 2) inside the SV span.

- **[C# API]** `StructuralVariantAnalyzer.FindDiscordantPairs(...)` → `.ClusterDiscordantPairs(...)`.

```
Provenance
1) find_discordant_pairs(readPairs,expectedInsertSize=400,insertSizeStdDev=50,maxInsertSize=10000) → discordant pairs
2) cluster_discordant_pairs(pairs,clusterDistance=500,minSupport=3) → SV variants (type,span,support,quality)
Cross-check: SV span should contain an independent split-read breakpoint (pipeline 2).
Envelope: NO guarded unit. Model floor: Simplified signature-then-cluster (SV-DETECT-001); not exhaustive assembly.
Caveat: alpha — validate before use.
```

## 2. Split-read breakpoints → junction sequence + microhomology

**Goal.** Localise a breakpoint to single-base resolution and characterise the junction.

1. **[MCP]** `find_split_reads`(alignments, minClipLength=20) → split `reads[]` (primaryPosition,
   supplementaryPosition, clipLength, clippedSequence). Only soft-clips (`S`) ≥ `minClipLength` count.
2. **[MCP]** `cluster_split_reads`(splitReads, clusterDistance=10, minSupport=2) → `breakpoints[]`
   (position = mean of primary+supplementary; quality = `min(support·15, 100)`).
3. **[MCP]** `assemble_breakpoint_sequence`(splitReads, minOverlap=10) → junction `sequence`
   (= clipped sequence of the largest-`clipLength` read; `null` if no reads). Heuristic, not true OLC.
4. **[MCP]** `find_microhomology`(leftFlank, rightFlank, maxLength=20) → longest shared suffix/prefix
   at the junction → `microhomologyLength` + upper-cased `sequence` (MMBIR/MMEJ signature; 0 if none).

- **[C# API]** `.FindSplitReads` → `.ClusterSplitReads` → `.AssembleBreakpointSequence` → `.FindMicrohomology`.

```
Provenance
1) find_split_reads(alignments,minClipLength=20) → split reads
2) cluster_split_reads(splitReads,clusterDistance=10,minSupport=2) → breakpoints (single-base)
3) assemble_breakpoint_sequence(splitReads,minOverlap=10) → junction sequence (longest-clip heuristic)
4) find_microhomology(leftFlank,rightFlank,maxLength=20) → microhomology length + sequence
Envelope: NO guarded unit. Model floor: Simplified (SV-BREAKPOINT-001); assembly is longest-clip, not OLC.
Caveat: alpha.
```

## 3. Read-depth CNV calling (probes → segments → CNVs)

**Goal.** Call focal deletion/duplication CNVs from copy-number probe data.

1. **[MCP]** `segment_copy_number`(probes, changeThreshold=0.3, minProbes=5) → `segments[]`
   (span, meanLogR, `copyNumber = round(2·2^meanLogR)` clamped 0–10, mean BAF, probeCount).
   CBS-*like*: new segment at chr change or when `|logR − segMean| > changeThreshold` (after ≥ `minProbes`).
2. **[MCP]** `identify_cnvs`(segments, normalCopyNumber=2, minLength=10000) → CNV `variants[]`
   (Deletion below baseline / Duplication above; id `CNV1…`; quality = `|logRatio|·50`; supportingReads = probeCount).

- **[C# API]** `StructuralVariantAnalyzer.SegmentCopyNumber(...)` → `.IdentifyCNVs(...)`.
- **Model floor / scope:** one **non-allelic** copy number per segment for a **germline diploid** —
  no tumor purity, no GC-bias correction, no major/minor allelic split. For tumor allele-specific /
  ASCAT copy number → [`seqeron-oncology`](../../seqeron-oncology/SKILL.md); for arm-scale / aneuploidy
  → [`bio-chromosome`](../../bio-chromosome/SKILL.md).

```
Provenance
1) segment_copy_number(probes,changeThreshold=0.3,minProbes=5) → segments (span,logRatio,copyNumber,BAF,probeCount)
2) identify_cnvs(segments,normalCopyNumber=2,minLength=10000) → Deletion/Duplication CNVs
Envelope: NO guarded unit. Model floor: Simplified read-depth (SV-CNV-001); CBS-like, no purity/GC model, non-allelic.
Caveat: alpha.
```

## 4. SV post-processing (filter → merge → genotype → annotate)

**Goal.** Turn a raw candidate set into a filtered, deduplicated, genotyped, annotated call set.

1. **[MCP]** `filter_svs`(variants, minQuality=20, minSupport=2, minLength=50, maxLength=100000000)
   → keep SVs passing **all four** thresholds.
2. **[MCP]** `merge_overlapping_svs`(variants, overlapFraction=0.5) → `merged`: consecutive same-chr,
   same-type SVs with overlap ÷ smaller-length ≥ `overlapFraction` (union span, first id, summed
   support, max quality). Only adjacent-in-sort entries combine — pre-sort matters.
3. **[MCP]** `genotype_sv`(sv, refReads, altReads, totalReads) → `genotype` + `quality`
   (altFraction < 0.1 → 0/0; > 0.9 → 1/1; 0.3–0.7 → 0/1; else 0/1; `totalReads=0` → `./.` q0).
4. **[MCP]** `annotate_svs`(variants, genes) → per-SV `affectedGenes`, `affectedExons` (`GENE:exonN`),
   `functionalImpact` (HIGH/MODERATE/MODIFIER/LOW), `isPathogenic` (HIGH∪MODERATE), `populationFrequency` (always 0).

- **[C# API]** `.FilterSVs` → `.MergeOverlappingSVs` → `.GenotypeSV` → `.AnnotateSVs`.

```
Provenance
1) filter_svs(variants,minQuality=20,minSupport=2,minLength=50,maxLength=1e8) → passing SVs
2) merge_overlapping_svs(variants,overlapFraction=0.5) → merged (same-chr, same-type, adjacent)
3) genotype_sv(sv,refReads,altReads,totalReads) → genotype + quality (cap 99)
4) annotate_svs(variants,genes) → affectedGenes/exons, functionalImpact, isPathogenic
Envelope: NO guarded unit. annotate_svs populationFrequency is always 0 (no population DB attached).
Caveat: alpha; not for clinical use.
```

## Scope reminders

- **Tumor allele-specific / clonal copy number (ASCAT, major/minor CN, LOH, HRD)** →
  [`seqeron-oncology`](../../seqeron-oncology/SKILL.md). This family is germline, single-copy-number.
- **Chromosome-arm / whole-chromosome amplification-deletion, aneuploidy, karyotype** →
  [`bio-chromosome`](../../bio-chromosome/SKILL.md). CNVs here are focal, not arm-scale ploidy.
- **SNP / indel base-level calling, VEP-like effect, ACMG classification** →
  [`bio-annotation`](../../bio-annotation/SKILL.md). This family is large rearrangements from read signatures.
