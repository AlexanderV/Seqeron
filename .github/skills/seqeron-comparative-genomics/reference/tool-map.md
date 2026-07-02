# seqeron-comparative-genomics — tool map (all 10)

Server: **analysis**. One backing class: `ComparativeGenomics.*`.
This skill is **not** in `domain-map.json`, so it has **no** generated `_generated/tools.md` —
**this curated map is the index.** Verify schemas in `docs/mcp/tools/analysis/<tool>.md`.

> Coordinates & units: gene `start`/`end` are **0-based**; ANI and identities are **fractions in
> [0, 1]** (not percent); dot-plot points are `{x=offset in seq1, y=occurrence offset in seq2}`.
> A `Gene` is `{ id, genomeId, start, end, strand, sequence? }` — `sequence` is **required** for
> ortholog/RBH tools. Always confirm exact I/O in the tool doc.

## Genome identity

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `calculate_ani` | `ComparativeGenomics.CalculateANI` | ANIb (Goris 2007): genome-1 cut into `fragmentSize`-nt fragments, each placed by best ungapped match on genome 2; a fragment counts if identity > `minFragmentIdentity` **and** ≥70% of it aligns. Mean of qualifying fragments → fraction; ~0.95 = species boundary; **0** if empty / nothing qualifies. | `calculate_ani.md` |

## Orthology

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_orthologs` | `ComparativeGenomics.FindOrthologs` | One-directional **best-hit** ortholog pairs by **5-mer Jaccard** similarity of gene sequences, kept above `minIdentity` (0.3) & `minCoverage` (0.5). → `{gene1Id,gene2Id,identity,coverage,alignmentLength}`. | `find_orthologs.md` |
| `find_reciprocal_best_hits` | `ComparativeGenomics.FindReciprocalBestHits` | **RBH** — stricter: a pair is kept only if each gene is the other's best hit above `minIdentity` (0.3). RBH ⊆ best-hits. | `find_reciprocal_best_hits.md` |

## Synteny & rearrangements (gene-order level)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_syntenic_blocks` | `ComparativeGenomics.FindSyntenicBlocks` | Collinear runs of orthologous genes from an **ortholog map** (`g1 id → g2 id`): adjacent (within `maxGap`=5) consistently ordered anchors chained; block needs ≥ `minBlockSize` (3). → `{genome1Id,start1,end1,genome2Id,start2,end2,isInverted,geneCount,identity}`. **Owns**: gene-order level. (Chromosome-scale `find_synteny_blocks`/`ChromosomeAnalyzer.FindSyntenyBlocks` → bio-chromosome.) | `find_syntenic_blocks.md` |
| `detect_rearrangements` | `ComparativeGenomics.DetectRearrangements` | Breakpoints of the **signed gene-order permutation** (Bafna & Pevzner 1998) via `orthologMap`; every non-identity adjacency is one event, classified **Inversion** (sign flip) or **Transposition** (orientation-preserving discontinuity). **Owns**: gene-order level. (Chromosome-scale `detect_rearrangements`/`ChromosomeAnalyzer.DetectRearrangements`, which also calls del/dup/transloc from synteny blocks → bio-chromosome.) | `detect_rearrangements.md` |
| `reversal_distance` | `ComparativeGenomics.CalculateReversalDistance` | **Lower bound** on reversal (inversion) distance of two equal-length permutations: `d = ⌈breakpoints/2⌉` (unsigned breakpoint bound). Identical orders → 0. Not the exact sorting distance. | `reversal_distance.md` |

## Clusters, dot-plot, pipeline

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_conserved_clusters` | `ComparativeGenomics.FindConservedClusters` | Conserved gene clusters as **common intervals** of ortholog-group permutations (Uno & Yagiura 2000; Heber & Stoye 2001). Genes → groups via `orthologGroups`; ungrouped genes break windows. Needs ≥2 genomes (else empty); strict gap-free (`maxGap` retained for compat only); `minClusterSize` ≥ 2 (default 3). → `clusters` (sorted group-id lists). | `find_conserved_clusters.md` |
| `generate_dot_plot` | `ComparativeGenomics.GenerateDotPlot` | For every shared `wordSize`-mer (default 10, slid by `stepSize` 1 along seq1, located in seq2 via suffix tree) emit `{x,y}`. Diagonal runs reveal collinear similarity. | `generate_dot_plot.md` |
| `compare_genomes` | `ComparativeGenomics.CompareGenomes` | Full **two-genome** pipeline: RBH orthologs + syntenic blocks + rearrangements + summary → `syntenicBlocks`, `orthologs`, `rearrangements`, `overallSynteny`, `conservedGenes`, `genomeSpecificGenes1/2` (core/dispensable, Tettelin 2005). For >2 genomes / Heaps' law → bio-metagenomics `PanGenomeAnalyzer.*`. | `compare_genomes.md` |

## Shared tool names — who owns what

| Name | This skill (Analysis) | bio-chromosome (Chromosome) |
|---|---|---|
| synteny | `find_syntenic_blocks` / `ComparativeGenomics.FindSyntenicBlocks` (ortholog map over Gene records) | `find_synteny_blocks` / `ChromosomeAnalyzer.FindSyntenyBlocks` (ortholog pairs by chromosome, gaps in Mb); `find_syntenic_blocks_assemblies` / `GenomeAssemblyAnalyzer.FindSyntenicBlocks` (k-mer, assemblies) |
| rearrangements | `detect_rearrangements` / `ComparativeGenomics.DetectRearrangements` (signed-permutation breakpoints: Inversion/Transposition) | `detect_rearrangements` / `ChromosomeAnalyzer.DetectRearrangements` (inv/transloc/del/dup from synteny blocks) |

## Envelope

- No comparative-genomics unit is **guarded** in [`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md)
  (nothing throws below a `MinimumMode`) — no per-call STOP rule. Standard `bio-rigor` envelope
  applies. Report the model floors (ANIb; 5-mer-Jaccard orthology; reversal distance is a lower
  bound; only Inversion/Transposition here).

## Not a tool here (route elsewhere)

- **Chromosome / assembly-scale synteny + del/dup/translocation SV** → [`bio-chromosome`](../../bio-chromosome/SKILL.md).
- **Multi-genome pan-genome, core/accessory, Heaps' law** → [`bio-metagenomics`](../../bio-metagenomics/SKILL.md) (`PanGenomeAnalyzer.*`).
- **Gene calling / ORF prediction (to produce the `Gene` inputs)** → [`bio-annotation`](../../bio-annotation/SKILL.md).
- **Pairwise/MSA alignment identity between two sequences** → [`bio-alignment`](../../bio-alignment/SKILL.md) (ANI here is fragment-based, not a global alignment).
