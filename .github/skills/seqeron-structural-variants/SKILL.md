---
name: seqeron-structural-variants
description: >-
  Detect germline structural variants (SVs) and copy-number variants (CNVs) from
  sequencing read evidence with Seqeron (MCP tools OR the C# API). Use to flag
  discordant read pairs and split reads (breakpoint signatures), cluster them
  into SV candidates / breakpoints, assemble the breakpoint-junction sequence and
  find microhomology at it, segment copy-number probes and call CNVs
  (deletion / duplication), and genotype / filter / merge / gene-annotate SVs.
  Triggers: "find structural variants", "find discordant pairs / split reads",
  "cluster reads into SV candidates", "assemble the breakpoint", "microhomology
  at this junction", "call CNVs / segment copy number", "genotype this SV",
  "filter / merge / annotate SVs", "deletion / duplication / inversion /
  translocation from reads". (Tumor allele-specific copy number / ASCAT →
  seqeron-oncology; chromosome-arm-scale amplification / aneuploidy →
  bio-chromosome; SNP / indel calling → bio-annotation.) Server: annotation
  (all StructuralVariantAnalyzer.*).
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-structural-variants — germline read-evidence SV/CNV workflow

Routing + orchestration skill for the germline, **alignment-evidence** SV/CNV family on the
**annotation** server (12 tools, one backing class: `StructuralVariantAnalyzer.*`). It picks the
right tool for a discordant-pair / split-read / breakpoint / copy-number / genotype-filter-merge
question and gives a **dual-mode** recipe (MCP tool name **and** the equivalent `Seqeron.Genomics`
C# `Method ID`).

- **Rigor is delegated.** Tool-only computation, provenance, envelope, cross-check, units/0-based
  coords, and the alpha / not-for-clinical caveat are owned by
  **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server annotation`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/annotation/*.md`;
  algorithm invariants in `docs/algorithms/StructuralVar/*.md`.

## Scope boundary — this family vs neighbours

This skill OWNS **germline SV/CNV from read/alignment evidence**: paired-end discordant signatures +
split reads → clustered breakpoints → junction assembly / microhomology; read-depth-probe copy-number
segmentation → deletion/duplication CNV calls; and SV genotype / filter / merge / gene-annotation.
**[`bio-annotation`](../bio-annotation/SKILL.md)** name-drops this family shallowly (line "SVs, CNVs,
breakpoints…") — route here.

Three boundaries to draw:

- **Tumor allele-specific copy number / clonal layer → [`seqeron-oncology`](../seqeron-oncology/SKILL.md).**
  ASCAT-like purity/ploidy, allele-specific (major/minor) copy number, LOH/HRD. Here `segment_copy_number`
  emits a *single* copy number per segment (no purity model, no allelic split) for a germline diploid.
- **Chromosome-arm-scale amplification / deletion / aneuploidy → [`bio-chromosome`](../bio-chromosome/SKILL.md).**
  Whole-arm / whole-chromosome gain-loss, trisomy/monosomy, karyotype/ploidy, band-scale rearrangements.
  Here CNVs are focal segment-level events, not arm-scale ploidy calls.
- **SNP / indel calling (ref-vs-query substitutions & small indels) → [`bio-annotation`](../bio-annotation/SKILL.md).**
  This skill handles large rearrangements (≥ ~50 bp) from read *signatures*, not base-level variant calling.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Flag discordant read pairs (SV signatures) | `find_discordant_pairs` / `StructuralVariantAnalyzer.FindDiscordantPairs` |
| Cluster discordant pairs → SV candidates | `cluster_discordant_pairs` / `StructuralVariantAnalyzer.ClusterDiscordantPairs` |
| Find split reads (soft-clip CIGAR) | `find_split_reads` / `StructuralVariantAnalyzer.FindSplitReads` |
| Cluster split reads → breakpoints | `cluster_split_reads` / `StructuralVariantAnalyzer.ClusterSplitReads` |
| Assemble breakpoint-junction sequence | `assemble_breakpoint_sequence` / `StructuralVariantAnalyzer.AssembleBreakpointSequence` |
| Microhomology at a junction | `find_microhomology` / `StructuralVariantAnalyzer.FindMicrohomology` |
| Segment copy-number probes → log-ratio plateaus | `segment_copy_number` / `StructuralVariantAnalyzer.SegmentCopyNumber` |
| Segments → deletion/duplication CNVs | `identify_cnvs` / `StructuralVariantAnalyzer.IdentifyCNVs` |
| Genotype an SV (0/0, 0/1, 1/1, ./.) | `genotype_sv` / `StructuralVariantAnalyzer.GenotypeSV` |
| Filter SVs (quality / support / length) | `filter_svs` / `StructuralVariantAnalyzer.FilterSVs` |
| Merge overlapping same-type SVs | `merge_overlapping_svs` / `StructuralVariantAnalyzer.MergeOverlappingSVs` |
| Annotate SVs with genes/exons + impact | `annotate_svs` / `StructuralVariantAnalyzer.AnnotateSVs` |

## Canonical dual-mode pipelines

### (a) Paired-end SV calling: reads → discordant pairs → clustered SVs
1. **[MCP]** `find_discordant_pairs`(readPairs, expectedInsertSize=400, insertSizeStdDev=50, maxInsertSize=10000) → discordant `pairs[]`.
2. **[MCP]** `cluster_discordant_pairs`(discordantPairs, clusterDistance=500, minSupport=3) → SV `variants[]` (type inferred from signature).
- **[C# API]** `StructuralVariantAnalyzer.FindDiscordantPairs(...)` → `.ClusterDiscordantPairs(...)`.

### (b) Split-read breakpoints → junction sequence + microhomology
1. **[MCP]** `find_split_reads`(alignments, minClipLength=20) → split `reads[]`.
2. **[MCP]** `cluster_split_reads`(splitReads, clusterDistance=10, minSupport=2) → `breakpoints[]` (single-base resolution).
3. **[MCP]** `assemble_breakpoint_sequence`(splitReads, minOverlap=10) → junction `sequence` (longest-clip heuristic).
4. **[MCP]** `find_microhomology`(leftFlank, rightFlank, maxLength=20) → shared microhomology (MMBIR/MMEJ signature).
- **[C# API]** `.FindSplitReads` → `.ClusterSplitReads` → `.AssembleBreakpointSequence` → `.FindMicrohomology`.

### (c) Read-depth CNV calling: probes → segments → CNVs
1. **[MCP]** `segment_copy_number`(probes, changeThreshold=0.3, minProbes=5) → `segments[]` (CBS-like; copyNumber = round(2·2^meanLogR), clamped 0–10).
2. **[MCP]** `identify_cnvs`(segments, normalCopyNumber=2, minLength=10000) → CNV `variants[]` (Deletion below / Duplication above baseline).
- **[C# API]** `StructuralVariantAnalyzer.SegmentCopyNumber(...)` → `.IdentifyCNVs(...)`.
- Model floor: single (non-allelic) copy number per segment, no purity/GC model — allele-specific/tumor CN → seqeron-oncology.

### (d) Post-processing: filter → merge → genotype → annotate
1. **[MCP]** `filter_svs`(variants, minQuality=20, minSupport=2, minLength=50, maxLength=1e8) → passing `variants` (all 4 conditions).
2. **[MCP]** `merge_overlapping_svs`(variants, overlapFraction=0.5) → `merged` (adjacent, same-chr, same-type only).
3. **[MCP]** `genotype_sv`(sv, refReads, altReads, totalReads) → `genotype` + `quality` (altFraction thresholds; ./. if totalReads=0).
4. **[MCP]** `annotate_svs`(variants, genes) → per-SV `affectedGenes`, `affectedExons`, `functionalImpact`, `isPathogenic`.
- **[C# API]** `.FilterSVs` → `.MergeOverlappingSVs` → `.GenotypeSV` → `.AnnotateSVs`.

## Envelope — no guarded unit (verified)

**This family has NO guarded unit.** `docs/Validation/LIMITATIONS.md` has **no** SV/CNV entry
(the algorithm docs' Test-Unit IDs `SV-DETECT-001`, `SV-BREAKPOINT-001`, `SV-CNV-001` do **not**
appear there) — so nothing throws a `SeqeronLimitationException` and there is **no per-call STOP
rule to fabricate**. Report the **model floors** instead: all three algorithms are
Implementation-Status **Simplified** — signature-then-cluster heuristics, not exhaustive assembly
or probabilistic segmentation; `segment_copy_number` is CBS-*like*, not statistical change-point;
`assemble_breakpoint_sequence` is a longest-clip heuristic (no true overlap-layout); no purity/GC
model. The general `bio-rigor` envelope + alpha / not-for-clinical caveat still apply, and the STOP
rule re-engages if a downstream chain reaches a genuinely guarded unit.

## End-to-end grounded example

**Task.** From a paired-end library, call a candidate deletion, refine its breakpoint with split
reads, then filter and genotype it.

1. `find_discordant_pairs`(readPairs) → discordant pairs (large-span → deletion signature). (`.FindDiscordantPairs`)
2. `cluster_discordant_pairs`(pairs, minSupport=3) → one Deletion SV with support + span. (`.ClusterDiscordantPairs`)
3. `find_split_reads`(alignments) → soft-clip split reads at the same locus; `cluster_split_reads` → single-base breakpoint that should sit inside the step-2 SV span (independent evidence). (`.FindSplitReads` → `.ClusterSplitReads`)
4. `filter_svs`(variants, minQuality=20, minSupport=2, minLength=50) → the SV passes. (`.FilterSVs`)
5. `genotype_sv`(sv, refReads=10, altReads=10, totalReads=20) → `0/1`, quality 40. (`.GenotypeSV`)
```
Provenance
1) find_discordant_pairs(readPairs,expectedInsertSize=400,insertSizeStdDev=50) → discordant pairs
2) cluster_discordant_pairs(pairs,clusterDistance=500,minSupport=3) → Deletion SV (span,support)
3) find_split_reads(alignments,minClipLength=20) → cluster_split_reads(...,minSupport=2) → breakpoint ⊂ step-2 span (cross-check)
4) filter_svs(variants,minQuality=20,minSupport=2,minLength=50) → SV passes
5) genotype_sv(sv,refReads=10,altReads=10,totalReads=20) → 0/1, quality 40
Envelope: NO guarded unit (no SV/CNV row in LIMITATIONS.md). Model floor: Simplified signature-then-cluster heuristics.
Caveat: alpha software; not for clinical use — independently validate before relying on any call.
```

## Reference

- **This family's tool map (all 12 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`SV_Detection.md`](../../../docs/algorithms/StructuralVar/SV_Detection.md) ·
  [`Breakpoint_Detection.md`](../../../docs/algorithms/StructuralVar/Breakpoint_Detection.md) ·
  [`Copy_Number_Variation.md`](../../../docs/algorithms/StructuralVar/Copy_Number_Variation.md)
- **Operating envelope:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (no SV/CNV entry — none guarded)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) ·
  [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup)
- **Neighbours (scope boundaries):** [`bio-annotation`](../bio-annotation/SKILL.md) (SNP/indel core + shallow SV name-drop) ·
  [`seqeron-oncology`](../seqeron-oncology/SKILL.md) (tumor allele-specific / clonal CN) ·
  [`bio-chromosome`](../bio-chromosome/SKILL.md) (arm-scale amplification/aneuploidy)
