---
name: bio-phylo-popgen
description: >-
  Build phylogenetic trees and compute population-genetics statistics with
  Seqeron (MCP tools OR the C# API). Phylogenetics: distance matrices
  (p-distance / Jukes-Cantor / Kimura-2P / Hamming), tree construction
  (Neighbor-Joining, UPGMA), Newick parse/serialize, tree stats (depth, length,
  leaves, MRCA, patristic + Robinson-Foulds distance), bootstrap support.
  Population genetics: allele/genotype frequencies, minor allele frequency (MAF)
  + MAF filtering, Hardy-Weinberg equilibrium test, Fst / F-statistics / pairwise
  Fst, nucleotide diversity (π), Watterson's θ, Tajima's D, heterozygosity,
  linkage disequilibrium (D'/r²) + haplotype blocks, runs of homozygosity /
  inbreeding, iHS, ancestry, selection scans. Triggers: "build a tree", "make a
  phylogeny / dendrogram from these sequences", "neighbor-joining / UPGMA tree",
  "distance matrix", "compute Fst", "test Hardy-Weinberg / HWE", "nucleotide
  diversity of this population", "Tajima's D", "linkage disequilibrium", "MAF
  filter". Servers: phylogenetics + population.
allowed-tools: Read, Bash, Grep, Glob
---

# bio-phylo-popgen — phylogeny + population genetics

Routing + orchestration skill for the **Phylogenetics** (13) and **Population** (18) servers = **31 tools**.
It picks the right tool for a tree-building or popgen question and gives a **dual-mode** recipe
(MCP tool calls **and** the equivalent `Seqeron.Genomics` C# `Method ID`s).

- **Rigor is delegated.** Parse-with-a-tool, envelope, provenance, cross-check, units/0-based
  coordinates, and the alpha / not-for-clinical-use caveat are owned by **`bio-rigor`** — it applies
  here by default; do not restate its rules.
- **Don't know the tool name?** Use **`seqeron-discovery`**
  (`python3 scripts/skills/find-tool.py <kw> --server phylogenetics|population`) — never guess.
- **Point, don't duplicate.** Full I/O schemas live in `docs/mcp/tools/{phylogenetics,population}/*.md`;
  algorithm invariants in `docs/algorithms/{Phylogenetics,Population_Genetics,PopGen}/*.md`. This skill
  links, it does not copy.

## Decision guide — which tool for which question

| Question | Tool ([MCP] / `Method ID`) |
|---|---|
| Build a tree **from aligned sequences** (distance + NJ/UPGMA in one step) | `build_phylogenetic_tree` / `PhylogeneticAnalyzer.BuildTree` |
| Build a tree **from a precomputed distance matrix** | `build_tree_from_matrix` / `PhylogeneticAnalyzer.BuildTreeFromMatrix` |
| Pairwise **distance matrix** of aligned seqs | `distance_matrix` / `PhylogeneticAnalyzer.CalculateDistanceMatrix` |
| Evolutionary **distance between two** aligned seqs | `pairwise_distance` / `PhylogeneticAnalyzer.CalculatePairwiseDistance` |
| **Serialize / parse** Newick | `to_newick` / `.ToNewick` · `parse_newick` / `.ParseNewick` |
| Tree **stats** (height / total length / leaves) | `tree_depth` `.GetTreeDepth` · `tree_length` `.CalculateTreeLength` · `tree_leaves` `.GetLeaves` |
| **MRCA** / path length between taxa | `mrca` `.FindMRCA` · `patristic_distance` `.PatristicDistance` |
| **Compare two trees** (topology distance) | `robinson_foulds_distance` / `.RobinsonFouldsDistance` |
| **Clade support** (resampling) | `bootstrap_support` / `.Bootstrap` |
| **Allele frequencies** from genotype counts | `allele_frequencies` / `PopulationGeneticsAnalyzer.CalculateAlleleFrequencies` |
| **MAF** from a genotype vector · **filter** variants by MAF | `minor_allele_frequency` `.CalculateMAF` · `filter_variants_by_maf` `.FilterByMAF` |
| **Hardy-Weinberg** equilibrium test | `hardy_weinberg_test` / `.TestHardyWeinberg` |
| **Fst** (two pops) · **F-statistics** (Fis/Fit/Fst) · **pairwise Fst matrix** | `fst` `.CalculateFst` · `f_statistics` `.CalculateFStatistics` · `pairwise_fst` `.CalculatePairwiseFst` |
| **Nucleotide diversity π** · **Watterson's θ** · all diversity stats in one pass | `nucleotide_diversity` `.CalculateNucleotideDiversity` · `wattersons_theta` `.CalculateWattersonTheta` · `diversity_statistics` `.CalculateDiversityStatistics` |
| **Tajima's D** (neutrality) | `tajimas_d` / `.CalculateTajimasD` |
| **Linkage disequilibrium** (D'/r²) · **haplotype blocks** | `linkage_disequilibrium` `.CalculateLD` · `haplotype_blocks` `.FindHaplotypeBlocks` |
| **Runs of homozygosity** · **inbreeding F_ROH** | `runs_of_homozygosity` `.FindROH` · `inbreeding_from_roh` `.CalculateInbreedingFromROH` |
| **iHS** · **ancestry** · **selection scan** | `integrated_haplotype_score` `.CalculateIHS` · `estimate_ancestry` `.EstimateAncestry` · `scan_selection_signals` `.ScanForSelection` |

Rule of thumb: **UPGMA** → rooted, ultrametric (assumes a molecular clock); **Neighbor-Joining** → unrooted,
no clock assumption (prefer for real divergence data). Distance model: `Hamming`/`PDistance` for closely
related seqs; `JukesCantor`/`Kimura2Parameter` correct for multiple hits but return `+Infinity` at
saturation (`p ≥ 0.75`). Default distance model is `JukesCantor`, default tree method is `UPGMA`.

## Canonical dual-mode pipelines

Method IDs on the `PopulationGeneticsAnalyzer` / `PhylogeneticAnalyzer` classes. Full param/return
schemas are in the linked per-tool docs — do not guess parameters.

### (a) Aligned sequences → distance matrix → NJ/UPGMA tree → Newick + stats
1. **[MCP]** `build_phylogenetic_tree`(sequences={name→alignedSeq}, distanceMethod=`JukesCantor`, treeMethod=`NeighborJoining`) → `newick`, `taxa`, `distanceMatrix`, `method`.
2. **[MCP]** `tree_length`(newick) + `tree_depth`(newick) + `tree_leaves`(newick) → tree summary.
- **[C# API]** `PhylogeneticAnalyzer.BuildTree(sequences, distanceMethod, treeMethod)` → `.CalculateTreeLength / .GetTreeDepth / .GetLeaves`.
- **Cross-check:** run `distance_matrix`(alignedSequences, method) independently and confirm it equals the `distanceMatrix` echoed by `build_phylogenetic_tree`; or feed that matrix to `build_tree_from_matrix` and compare the Newick.
```
Provenance
1) build_phylogenetic_tree(distanceMethod=JukesCantor, treeMethod=NeighborJoining) → newick, distanceMatrix
2) tree_length/tree_depth/tree_leaves(newick) → stats
Cross-check: distance_matrix(same method) == echoed distanceMatrix.
Envelope: none guarded. Caveat: alpha — validate before any decision use.
```

### (b) Genotype counts → allele/genotype frequencies → Hardy-Weinberg test
1. **[MCP]** `allele_frequencies`(homozygousMajor, heterozygous, homozygousMinor) → `majorFreq p`, `minorFreq q`.
2. **[MCP]** `hardy_weinberg_test`(variantId, observedAA, observedAa, observedaa, significanceLevel=0.05) → `expectedAA/Aa/aa`, `chiSquare`, `pValue`, `inEquilibrium`.
- **[C# API]** `.CalculateAlleleFrequencies(nAA,nAa,naa)` → `.TestHardyWeinberg(id,AA,Aa,aa,alpha)`.
- **Cross-check:** HWE-expected `2pqn` should use the same `p` from step 1; `inEquilibrium` ⇔ `pValue ≥ α` (1 df χ²).

### (c) Two populations' allele frequencies → Fst / F-statistics
1. **[MCP]** `fst`(population1=[{alleleFreq,sampleSize}…], population2=[…]) → `fst` (Wright, in [0,1]).
2. **[MCP]** (optional) `f_statistics`(…) → `Fis`, `Fit`, `Fst`; for >2 pops use `pairwise_fst` → symmetric matrix.
- **[C# API]** `.CalculateFst(pop1,pop2)` · `.CalculateFStatistics(...)` · `.CalculatePairwiseFst(...)`.
- Per-locus counts must match between the two populations (error 1001 otherwise). Fst=0 panmixia, 1 fixed difference.

### (d) Aligned population sample → nucleotide diversity / Tajima's D (neutrality)
1. **[MCP]** `diversity_statistics`(sequences=[aligned…]) → `nucleotideDiversity π`, `wattersonTheta`, `tajimasD`, `segregatingSites S`, `sampleSize n`, observed/expected heterozygosity — one pass.
- **[C# API]** `.CalculateDiversityStatistics(sequences)`.
- **Cross-check / split path:** `nucleotide_diversity` (π only) and `wattersons_theta`(S,n,L); then `tajimas_d`(averagePairwiseDifferences=k̂=π·L, segregatingSites=S, sampleSize=n) — Tajima's D takes **k̂ (NOT per-site)**. D<0 excess rare variants (expansion/purifying); D>0 deficit (balancing/contraction).

### (e) Genotype data → MAF filtering / linkage disequilibrium
1. **[MCP]** `minor_allele_frequency`(genotypes=[0/1/2…]) → `maf` ∈ [0,0.5]; `filter_variants_by_maf`(variants, minMAF=0.01, maxMAF=0.5) → variants passing (filters on stored `alleleFrequency`, does NOT recompute).
2. **[MCP]** `linkage_disequilibrium`(variant1Id, variant2Id, genotypes=[{geno1,geno2}…], distance) → `dPrime`, `rSquared`; chain into `haplotype_blocks` for block structure.
- **[C# API]** `.CalculateMAF(genotypes)` · `.FilterByMAF(variants,min,max)` · `.CalculateLD(...)` · `.FindHaplotypeBlocks(...)`.

## End-to-end grounded example (extends `docs/mcp/README.md`)

**Task.** Four aligned orthologs `A,B,C,D` (equal length). (1) Build an NJ tree, (2) report total tree
length + leaf branch lengths, (3) test whether the sample departs from neutrality, (4) corroborate the
distances independently.

Tool / `Method ID` chain (MCP names; C# path in parentheses):
1. `build_phylogenetic_tree`(sequences={A,B,C,D}, distanceMethod="JukesCantor", treeMethod="NeighborJoining")
   → `newick`, `distanceMatrix`, `taxa`. (`PhylogeneticAnalyzer.BuildTree`)
2. `tree_length`(newick) → total branch length; `tree_leaves`(newick) → per-taxon branch lengths.
   (`.CalculateTreeLength`, `.GetLeaves`)
3. `diversity_statistics`(sequences=[A,B,C,D]) → `π`, `wattersonTheta`, `tajimasD`, `S`, `n`.
   (`.CalculateDiversityStatistics`)
4. Cross-check: `distance_matrix`(alignedSequences, method="JukesCantor") must equal the `distanceMatrix`
   echoed in step 1. (`.CalculateDistanceMatrix`)

Expected-shape output (values illustrative; **compute them with the tools, do not eyeball**):
```
| tree_method | leaves | tree_length | pi     | wattersonTheta | tajimasD | S |
|-------------|-------:|------------:|-------:|---------------:|---------:|--:|
| NJ          |      4 |        …    |   …    |       …        |    …     | … |

Provenance
1) build_phylogenetic_tree(distanceMethod=JukesCantor, treeMethod=NeighborJoining) → newick, distanceMatrix
2) tree_length(newick), tree_leaves(newick) → branch lengths
3) diversity_statistics([A,B,C,D]) → pi, wattersonTheta, tajimasD, S, n
Cross-check: distance_matrix(JukesCantor) == echoed distanceMatrix.
Envelope: none guarded (Phylogenetics + Population tools within contract).
Caveat: alpha software; not for clinical use — independently validate before relying on any result.
```

## Reference

- **Full domain tool index (all 31, generated — do NOT hand-edit):** [`_generated/tools.md`](_generated/tools.md)
  (produced by `scripts/skills/gen-catalog.py`; if absent, run `seqeron-discovery`).
- **Fuller recipes + parameter/model guidance:** [`reference/pipelines.md`](reference/pipelines.md)
- **Tool map (all 31 by sub-task, one-liners + Method ID):** [`reference/tool-map.md`](reference/tool-map.md)
- **Algorithm background (invariants/formulas — link, don't copy):**
  [`docs/algorithms/Phylogenetics/Tree_Construction.md`](../../../docs/algorithms/Phylogenetics/Tree_Construction.md) ·
  [`Distance_Matrix.md`](../../../docs/algorithms/Phylogenetics/Distance_Matrix.md) ·
  [`Newick_Format.md`](../../../docs/algorithms/Phylogenetics/Newick_Format.md) ·
  [`Tree_Statistics.md`](../../../docs/algorithms/Phylogenetics/Tree_Statistics.md) ·
  [`Tree_Comparison.md`](../../../docs/algorithms/Phylogenetics/Tree_Comparison.md) ·
  [`Bootstrap_Analysis.md`](../../../docs/algorithms/Phylogenetics/Bootstrap_Analysis.md) ·
  [`Population_Genetics/Diversity_Statistics.md`](../../../docs/algorithms/Population_Genetics/Diversity_Statistics.md) ·
  [`F_Statistics.md`](../../../docs/algorithms/Population_Genetics/F_Statistics.md) ·
  [`Hardy_Weinberg_Test.md`](../../../docs/algorithms/Population_Genetics/Hardy_Weinberg_Test.md) ·
  [`Allele_Frequency.md`](../../../docs/algorithms/Population_Genetics/Allele_Frequency.md) ·
  [`Linkage_Disequilibrium.md`](../../../docs/algorithms/Population_Genetics/Linkage_Disequilibrium.md)
- **Cross-cutting:** [`bio-rigor`](../bio-rigor/SKILL.md) (rigor guardrail) · [`seqeron-discovery`](../seqeron-discovery/SKILL.md) (tool lookup).
