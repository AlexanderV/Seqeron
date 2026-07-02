---
name: seqeron-comparative-genomics
description: >-
  Compare whole genomes / gene sets with Seqeron (MCP tools OR the C# API). Use
  to compute Average Nucleotide Identity (ANI / genome identity, ANIb species
  boundary ~0.95), detect orthologs and reciprocal best hits (RBH) between two
  genomes, find syntenic blocks (collinear runs of orthologous genes),
  detect whole-genome rearrangements (inversions / transpositions from the
  signed gene-order permutation), estimate reversal (inversion) distance, find
  conserved gene clusters (common intervals), generate a dot-plot of shared
  k-mers, and run the end-to-end genome-comparison pipeline (core/dispensable
  partition). Triggers: "compute ANI between these genomes", "how similar are
  these two genomes", "find orthologs / RBH between genome A and B", "syntenic
  blocks between these genomes", "detect rearrangements between gene orders",
  "reversal distance", "dot-plot of these two sequences", "conserved gene
  clusters", "compare these two genomes". Server: analysis (all
  ComparativeGenomics.*). Chromosome-scale synteny/SV → bio-chromosome.
allowed-tools: Read, Bash, Grep, Glob
---

# seqeron-comparative-genomics — ANI, orthologs/RBH, synteny, rearrangements, dot-plots

Routing + orchestration skill for the genome-comparison family on the **analysis** server
(10 tools, one backing class: `ComparativeGenomics.*`). It picks the right tool for an
ANI / ortholog / RBH / synteny / rearrangement / reversal-distance / dot-plot question and
gives a **dual-mode** recipe (MCP tool name **and** the equivalent `Seqeron.Genomics` C#
`Method ID`). Carved out of `bio-annotation`, which name-drops this family shallowly.

- **Rigor is delegated.** Tool-only computation, provenance, envelope, cross-check, units /
  0-based coords, and the alpha / not-for-clinical caveat are owned by
  **[`bio-rigor`](../bio-rigor/SKILL.md)** — it applies by default; do not restate.
- **Don't know a tool name?** Use **[`seqeron-discovery`](../seqeron-discovery/SKILL.md)**
  (`python3 scripts/skills/find-tool.py <kw> --server analysis`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/analysis/*.md`;
  algorithm invariants/citations in `docs/algorithms/Comparative_Genomics/*.md`.

## Scope boundary — this family vs neighbours (IMPORTANT: shared tool names)

Two tool names live on **both** the Analysis and Chromosome servers. This skill owns the
**sequence / gene-order (genome-comparison) level**; **[`bio-chromosome`](../bio-chromosome/SKILL.md)**
owns the **chromosome / assembly scale**. Disambiguate by `Method ID` and inputs:

| Question | This skill (Analysis, `ComparativeGenomics.*`) | bio-chromosome (Chromosome) |
|---|---|---|
| Synteny | `find_syntenic_blocks` / `ComparativeGenomics.FindSyntenicBlocks` — chains an **ortholog map** (`gene id → gene id`) over `Gene` records | `find_synteny_blocks` / `ChromosomeAnalyzer.FindSyntenyBlocks` — groups **ortholog pairs by chromosome**, gap in **Mb**; also `find_syntenic_blocks_assemblies` / `GenomeAssemblyAnalyzer.FindSyntenicBlocks` (k-mer, assemblies) |
| Rearrangements | `detect_rearrangements` / `ComparativeGenomics.DetectRearrangements` — **breakpoints of the signed gene-order permutation** (Inversion / Transposition), needs genes + `orthologMap` | `detect_rearrangements` / `ChromosomeAnalyzer.DetectRearrangements` — inv/transloc/**del/dup** from **synteny blocks** |

- **Route to bio-chromosome** when the input is chromosomes / assemblies / synteny blocks, when
  you need deletion/duplication/translocation calls, or when gaps are in Mb.
- **Stay here** when the input is two gene sets / an ortholog map / permutations / raw sequences.
- **Pan-genome** (core/accessory across *many* genomes, Heaps' law) → **[`bio-metagenomics`](../bio-metagenomics/SKILL.md)**
  (`PanGenomeAnalyzer.*`). `compare_genomes` here does a **two-genome** core/dispensable partition only.
- **bio-annotation** name-drops ANI/orthologs/synteny — route here for the real recipes.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Genome-wide identity between two sequences (ANIb, fraction; ~0.95 = species) | `calculate_ani` / `ComparativeGenomics.CalculateANI` |
| Ortholog pairs (best-hit, 5-mer Jaccard) between two gene sets | `find_orthologs` / `ComparativeGenomics.FindOrthologs` |
| Reciprocal best hits (stricter orthology) | `find_reciprocal_best_hits` / `ComparativeGenomics.FindReciprocalBestHits` |
| Syntenic blocks (collinear ortholog runs) from an ortholog map | `find_syntenic_blocks` / `ComparativeGenomics.FindSyntenicBlocks` |
| Rearrangements = breakpoints of signed gene-order permutation | `detect_rearrangements` / `ComparativeGenomics.DetectRearrangements` |
| Lower-bound reversal (inversion) distance of two permutations | `reversal_distance` / `ComparativeGenomics.CalculateReversalDistance` |
| Conserved gene clusters (common intervals) across ≥2 genomes | `find_conserved_clusters` / `ComparativeGenomics.FindConservedClusters` |
| Dot-plot: matching k-mer coordinates of two sequences | `generate_dot_plot` / `ComparativeGenomics.GenerateDotPlot` |
| End-to-end two-genome comparison (RBH+synteny+rearrangements+core/dispensable) | `compare_genomes` / `ComparativeGenomics.CompareGenomes` |

## Canonical dual-mode pipelines

### (a) Genome identity — ANI + dot-plot cross-check
1. **[MCP]** `calculate_ani`(genome1Sequence, genome2Sequence, fragmentSize?=1000, minFragmentIdentity?=0.7) → `ani` (fraction in [0,1]; ~0.95 = 70% DDH species boundary; 0 if nothing qualifies).
2. **[MCP]** corroborate visually: `generate_dot_plot`(sequence1, sequence2, wordSize?=10, stepSize?=1) → `points`; a strong main diagonal is consistent with high ANI.
- **[C# API]** `ComparativeGenomics.CalculateANI(...)` · `.GenerateDotPlot(...)`.
```
Provenance
1) calculate_ani(g1,g2,fragmentSize=1000,minFragmentIdentity=0.7) → ani (ANIb, Goris 2007)
2) generate_dot_plot(g1,g2,wordSize=10,stepSize=1) → points (diagonal consistent w/ ani)
Caveat: alpha — not for clinical use; validate before relying on any call.
```

### (b) Orthology — RBH (strict) with best-hit corroboration
1. **[MCP]** `find_reciprocal_best_hits`(genome1Genes, genome2Genes, minIdentity?=0.3) → RBH pairs `{gene1Id, gene2Id, identity, coverage, alignmentLength}`.
2. **[MCP]** compare against the looser one-directional set: `find_orthologs`(genome1Genes, genome2Genes, minIdentity?=0.3, minCoverage?=0.5). RBH ⊆ best-hits — RBH should be a subset.
- **[C# API]** `ComparativeGenomics.FindReciprocalBestHits(...)` · `.FindOrthologs(...)`. Each `Gene` = `{id, genomeId, start, end, strand, sequence}`; `sequence` must be populated.
```
Provenance
1) find_reciprocal_best_hits(g1Genes,g2Genes,minIdentity=0.3) → RBH pairs (Moreno-Hagelsieb 2008)
2) find_orthologs(g1Genes,g2Genes,minIdentity=0.3,minCoverage=0.5) → best-hit pairs (RBH ⊆ this)
Caveat: alpha — validate before use.
```

### (c) Synteny → rearrangements → reversal distance (gene-order level)
1. Build an ortholog map (gene1 id → gene2 id) from pipeline (b)'s RBH pairs.
2. **[MCP]** `find_syntenic_blocks`(genome1Genes, genome2Genes, orthologMap, minBlockSize?=3, maxGap?=5) → blocks `{genome1Id,start1,end1,genome2Id,start2,end2,isInverted,geneCount,identity}`.
3. **[MCP]** `detect_rearrangements`(genome1Genes, genome2Genes, orthologMap) → breakpoint events `{type=Inversion|Transposition, genomeId, position, ...}`.
4. **[MCP]** quantify: `reversal_distance`(permutation1, permutation2) → `distance = ⌈breakpoints/2⌉` (lower bound, Bafna & Pevzner 1998).
- **[C# API]** `.FindSyntenicBlocks(...)` → `.DetectRearrangements(...)` → `.CalculateReversalDistance(...)`.
- **Route to [`bio-chromosome`](../bio-chromosome/SKILL.md)** instead if you need del/dup/translocation calls or chromosome-scale blocks (see scope table).
```
Provenance
1) orthologMap from RBH (pipeline b)
2) find_syntenic_blocks(g1Genes,g2Genes,orthologMap,minBlockSize=3,maxGap=5) → blocks
3) detect_rearrangements(g1Genes,g2Genes,orthologMap) → Inversion/Transposition breakpoints
4) reversal_distance(perm1,perm2) → ⌈breakpoints/2⌉ (lower bound)
Caveat: alpha — validate before use.
```

### (d) One-shot two-genome comparison
1. **[MCP]** `compare_genomes`(genome1Genes, genome2Genes, minOrthologIdentity?=0.3, minSyntenicBlockSize?=3) → `syntenicBlocks`, `orthologs` (RBH), `rearrangements`, `overallSynteny`, `conservedGenes`, `genomeSpecificGenes1/2` (core/dispensable, Tettelin 2005).
- **[C# API]** `ComparativeGenomics.CompareGenomes(...)`. Wraps (b)+(c); for >2 genomes / Heaps' law use `bio-metagenomics`.

### (e) Conserved gene clusters across ≥2 genomes
1. **[MCP]** `find_conserved_clusters`(genomes[][], orthologGroups, minClusterSize?=3, maxGap?=2) → `clusters` (each a sorted list of ortholog-group ids; common intervals, Heber & Stoye 2001). Needs ≥2 genomes; strict gap-free model (`maxGap` retained for compat only).
- **[C# API]** `ComparativeGenomics.FindConservedClusters(...)`.

## Envelope

No comparative-genomics unit is **guarded** in `docs/Validation/LIMITATIONS.md` (nothing throws
below a `MinimumMode`) — no per-call STOP rule; the standard `bio-rigor` envelope applies. Report the
**model floors**: ANI is **ANIb** (ungapped fragment placement, ≥70% fragment coverage); orthology is
**5-mer Jaccard**, not full alignment; `reversal_distance` is a **lower bound** (unsigned breakpoint
bound), not the exact sorting distance; `detect_rearrangements` here classifies only **Inversion /
Transposition** (del/dup/translocation → `bio-chromosome`).

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** Given two bacterial genomes (raw sequence + called genes), decide whether they are the
same species and describe how their gene order rearranged.
1. `calculate_ani`(g1, g2, fragmentSize=1000) → `ani`; is it ≥ ~0.95? (`ComparativeGenomics.CalculateANI`)
2. `find_reciprocal_best_hits`(g1Genes, g2Genes, minIdentity=0.3) → RBH pairs; build `orthologMap`. (`.FindReciprocalBestHits`)
3. `find_syntenic_blocks`(g1Genes, g2Genes, orthologMap, minBlockSize=3) → collinear blocks + `isInverted`. (`.FindSyntenicBlocks`)
4. `detect_rearrangements`(g1Genes, g2Genes, orthologMap) → Inversion/Transposition breakpoints; then `reversal_distance`(perm1, perm2) → lower-bound reversal count. (`.DetectRearrangements` → `.CalculateReversalDistance`)
```
Provenance
1) calculate_ani(g1,g2,fragmentSize=1000) → ani (ANIb; ~0.95 species boundary)
2) find_reciprocal_best_hits(g1Genes,g2Genes,minIdentity=0.3) → RBH → orthologMap
3) find_syntenic_blocks(...,orthologMap,minBlockSize=3) → blocks (isInverted flags)
4) detect_rearrangements(...,orthologMap) → Inversion/Transposition; reversal_distance(perm1,perm2) → ⌈bp/2⌉ (lower bound)
Model floors: ANIb; Jaccard orthology; reversal distance is a lower bound; no del/dup/transloc here (→ bio-chromosome).
Caveat: alpha software; not for clinical use — independently validate before relying on any call.
```

## Reference

- **This family's tool map (all 10 — curated index; NOT in domain-map.json, so there is NO
  generated slice):** [`reference/tool-map.md`](reference/tool-map.md)
- **Fuller recipes + parameter guidance + shared-name disambiguation:** [`reference/pipelines.md`](reference/pipelines.md)
- **Algorithm background (invariants/formulas/citations — link, don't copy):**
  [`Average_Nucleotide_Identity.md`](../../../docs/algorithms/Comparative_Genomics/Average_Nucleotide_Identity.md) ·
  [`Ortholog_Identification.md`](../../../docs/algorithms/Comparative_Genomics/Ortholog_Identification.md) ·
  [`Reciprocal_Best_Hits.md`](../../../docs/algorithms/Comparative_Genomics/Reciprocal_Best_Hits.md) ·
  [`Synteny_Block_Detection.md`](../../../docs/algorithms/Comparative_Genomics/Synteny_Block_Detection.md) ·
  [`Genome_Rearrangement_Detection.md`](../../../docs/algorithms/Comparative_Genomics/Genome_Rearrangement_Detection.md) ·
  [`Reversal_Distance.md`](../../../docs/algorithms/Comparative_Genomics/Reversal_Distance.md) ·
  [`Conserved_Gene_Clusters.md`](../../../docs/algorithms/Comparative_Genomics/Conserved_Gene_Clusters.md) ·
  [`Dot_Plot_Generation.md`](../../../docs/algorithms/Comparative_Genomics/Dot_Plot_Generation.md) ·
  [`Genome_Comparison.md`](../../../docs/algorithms/Comparative_Genomics/Genome_Comparison.md)
- **Operating envelope:** [`LIMITATIONS.md`](../../../docs/Validation/LIMITATIONS.md) (no guarded comparative-genomics unit) ·
  **cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup)
- **Neighbours:** [`bio-chromosome`](../bio-chromosome/SKILL.md) (chromosome/assembly-scale synteny + del/dup/transloc SV; owns the same-named `find_synteny_blocks` / `detect_rearrangements`) ·
  [`bio-annotation`](../bio-annotation/SKILL.md) (name-drops this family; gene calling upstream) ·
  [`bio-metagenomics`](../bio-metagenomics/SKILL.md) (multi-genome pan-genome / core-accessory / Heaps' law)
