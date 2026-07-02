# seqeron-comparative-genomics — fuller pipelines & parameter guidance

Dual-mode recipes for the genome-comparison family (analysis server, `ComparativeGenomics.*`).
Rigor (tool-only, provenance, cross-check, alpha caveat) is delegated to
[`bio-rigor`](../../bio-rigor/SKILL.md). Schemas: `docs/mcp/tools/analysis/<tool>.md`. Tool index:
[`tool-map.md`](tool-map.md).

## Parameter defaults (verify per tool doc)

- `calculate_ani`: `fragmentSize` = 1000 nt (> 0) · `minFragmentIdentity` = 0.7. ANI is a
  **fraction**; ~0.95 ≈ 70% DDH species boundary; returns **0** on empty / no qualifying fragment.
- `find_orthologs`: `minIdentity` = 0.3 · `minCoverage` = 0.5 (5-mer Jaccard).
- `find_reciprocal_best_hits`: `minIdentity` = 0.3 (no coverage arg).
- `find_syntenic_blocks`: `minBlockSize` = 3 · `maxGap` = 5 (anchor gap, in genes).
- `reversal_distance`: two **equal-length** integer permutations; result = `⌈breakpoints/2⌉`.
- `find_conserved_clusters`: `minClusterSize` = 3 (≥ 2) · `maxGap` = 2 (retained for compat; model
  is strict gap-free). Needs ≥ 2 genomes.
- `generate_dot_plot`: `wordSize` = 10 (≥ 1) · `stepSize` = 1 (≥ 1).
- `compare_genomes`: `minOrthologIdentity` = 0.3 · `minSyntenicBlockSize` = 3.
- A `Gene` = `{ id, genomeId, start, end, strand, sequence? }`; `start`/`end` **0-based**.
  `sequence` is **required** for `find_orthologs` / `find_reciprocal_best_hits` / `compare_genomes`.

## Shared tool names — disambiguate before calling (READ THIS)

Two tool names live on both Analysis and Chromosome. Pick by `Method ID` + input shape:

| You have… | Call (this skill) | Not (bio-chromosome) |
|---|---|---|
| Two gene sets + `orthologMap` (gene id → gene id), gene-order level | `find_syntenic_blocks` / `ComparativeGenomics.FindSyntenicBlocks` | — |
| Ortholog pairs keyed by chromosome, gaps in **Mb**, chromosome scale | — | `find_synteny_blocks` / `ChromosomeAnalyzer.FindSyntenyBlocks` |
| Two **assemblies**, no genes, k-mer anchoring | — | `find_syntenic_blocks_assemblies` / `GenomeAssemblyAnalyzer.FindSyntenicBlocks` |
| Signed gene-order permutation breakpoints (Inversion/Transposition) | `detect_rearrangements` / `ComparativeGenomics.DetectRearrangements` | — |
| SV calls (inv/transloc/**del/dup**) from **synteny blocks** | — | `detect_rearrangements` / `ChromosomeAnalyzer.DetectRearrangements` |

If unsure which scale the user means, ask, or route to [`bio-chromosome`](../../bio-chromosome/SKILL.md)
when chromosomes / assemblies / SV types (del/dup/transloc) are in play.

## 1. Genome identity — ANI + dot-plot corroboration

**Goal.** Quantify how similar two genome sequences are and sanity-check the number visually.

1. **[MCP]** `calculate_ani`(genome1Sequence, genome2Sequence, fragmentSize=1000, minFragmentIdentity=0.7)
   → `ani` (fraction). Interpret: ≥ ~0.95 → same species; lower → below the boundary.
2. **[MCP]** `generate_dot_plot`(sequence1, sequence2, wordSize=10, stepSize=1) → `points`. A strong
   unbroken main diagonal corroborates high ANI; broken/off-diagonal runs hint at rearrangement.
3. Cross-check: high `ani` with a sparse/diffuse dot-plot is contradictory — reconcile before reporting.

- **[C# API]** `ComparativeGenomics.CalculateANI(...)` · `.GenerateDotPlot(...)`.

```
Provenance
1) calculate_ani(g1,g2,fragmentSize=1000,minFragmentIdentity=0.7) → ani (ANIb, Goris 2007)
2) generate_dot_plot(g1,g2,wordSize=10,stepSize=1) → points (diagonal ↔ ani)
Model floor: ANIb = ungapped fragment placement, ≥70% fragment coverage; fraction not percent.
Caveat: alpha — not for clinical use.
```

## 2. Orthology — RBH first, best-hits as a looser corroborator

**Goal.** Get a trustworthy ortholog set between two gene sets.

1. **[MCP]** `find_reciprocal_best_hits`(genome1Genes, genome2Genes, minIdentity=0.3)
   → RBH pairs `{gene1Id,gene2Id,identity,coverage,alignmentLength}` (strict).
2. **[MCP]** `find_orthologs`(genome1Genes, genome2Genes, minIdentity=0.3, minCoverage=0.5)
   → best-hit pairs (looser, one-directional). RBH should be a **subset** — a RBH pair missing from
   best-hits is a red flag.

- **[C# API]** `ComparativeGenomics.FindReciprocalBestHits(...)` · `.FindOrthologs(...)`.

```
Provenance
1) find_reciprocal_best_hits(g1Genes,g2Genes,minIdentity=0.3) → RBH pairs
2) find_orthologs(g1Genes,g2Genes,minIdentity=0.3,minCoverage=0.5) → best-hits (RBH ⊆ this)
Model floor: similarity is 5-mer Jaccard, not full alignment.
Caveat: alpha.
```

## 3. Synteny → rearrangements → reversal distance (gene-order level)

**Goal.** From orthology, describe collinearity and the rearrangements that scramble it.

1. Build `orthologMap` (g1 id → g2 id) from pipeline 2's RBH pairs.
2. **[MCP]** `find_syntenic_blocks`(genome1Genes, genome2Genes, orthologMap, minBlockSize=3, maxGap=5)
   → blocks with coordinates in both genomes, `isInverted`, `geneCount`, `identity`.
3. **[MCP]** `detect_rearrangements`(genome1Genes, genome2Genes, orthologMap)
   → breakpoint events, each `Inversion` (sign flip) or `Transposition` (order-preserving discontinuity).
   Identical gene order → no events.
4. **[MCP]** `reversal_distance`(permutation1, permutation2) → lower-bound reversal count
   `⌈breakpoints/2⌉`. Derive the permutations by relabelling the shared markers to genome-2 rank.
5. Cross-check: an all-forward, single-block synteny (step 2) should coincide with zero events
   (step 3) and distance 0 (step 4).

- **[C# API]** `.FindSyntenicBlocks(...)` → `.DetectRearrangements(...)` → `.CalculateReversalDistance(...)`.
- **STOP-and-route:** need deletion / duplication / translocation calls, or chromosome-scale blocks?
  Those are **[`bio-chromosome`](../../bio-chromosome/SKILL.md)** (`ChromosomeAnalyzer.*`), not this
  skill — do not force del/dup/transloc out of the signed-permutation model here.

```
Provenance
1) orthologMap ← RBH (pipeline 2)
2) find_syntenic_blocks(g1Genes,g2Genes,orthologMap,minBlockSize=3,maxGap=5) → blocks (isInverted)
3) detect_rearrangements(g1Genes,g2Genes,orthologMap) → Inversion/Transposition (Bafna & Pevzner 1998)
4) reversal_distance(perm1,perm2) → ⌈breakpoints/2⌉ (lower bound, not exact sorting distance)
Model floor: only Inversion/Transposition here; del/dup/transloc → bio-chromosome.
Caveat: alpha.
```

## 4. One-shot two-genome comparison

**Goal.** Run RBH + synteny + rearrangements + core/dispensable in one call.

1. **[MCP]** `compare_genomes`(genome1Genes, genome2Genes, minOrthologIdentity=0.3, minSyntenicBlockSize=3)
   → `syntenicBlocks`, `orthologs` (RBH), `rearrangements`, `overallSynteny`, `conservedGenes`
   (core), `genomeSpecificGenes1`, `genomeSpecificGenes2` (dispensable, Tettelin 2005).

- **[C# API]** `ComparativeGenomics.CompareGenomes(...)`.
- This is a **two-genome** convenience wrapper of pipelines 2+3. For > 2 genomes, core/accessory
  spectra, or Heaps' law → [`bio-metagenomics`](../../bio-metagenomics/SKILL.md) (`PanGenomeAnalyzer.*`).

```
Provenance
1) compare_genomes(g1Genes,g2Genes,minOrthologIdentity=0.3,minSyntenicBlockSize=3)
   → syntenicBlocks,orthologs,rearrangements,overallSynteny,conservedGenes,genomeSpecificGenes1/2
Model floor: two genomes only; RBH orthology; Tettelin core/dispensable partition.
Caveat: alpha.
```

## 5. Conserved gene clusters (common intervals) across ≥2 genomes

**Goal.** Find gene sets that stay contiguous in every genome (operon/cluster conservation).

1. **[MCP]** `find_conserved_clusters`(genomes[][], orthologGroups, minClusterSize=3, maxGap=2)
   → `clusters`, each a sorted list of ortholog-group ids (common intervals; genes with no group
   break windows; needs ≥ 2 genomes; strict gap-free — `maxGap` retained for compatibility only).

- **[C# API]** `ComparativeGenomics.FindConservedClusters(...)`.

```
Provenance
1) find_conserved_clusters(genomes,orthologGroups,minClusterSize=3,maxGap=2) → clusters
Model floor: strict common intervals (gap-free); ≥2 genomes required (else empty).
Caveat: alpha.
```

## Scope reminders

- **Chromosome / assembly-scale synteny + del/dup/translocation SV** → [`bio-chromosome`](../../bio-chromosome/SKILL.md)
  (owns the same-named `find_synteny_blocks` / `detect_rearrangements`).
- **Multi-genome pan-genome, core/accessory, Heaps' law** → [`bio-metagenomics`](../../bio-metagenomics/SKILL.md).
- **Gene calling to produce `Gene` inputs** → [`bio-annotation`](../../bio-annotation/SKILL.md).
- No comparative-genomics unit is **guarded** in
  [`../../../../docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md) — no
  per-call STOP rule; just report the model floors above.
