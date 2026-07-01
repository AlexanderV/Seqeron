# bio-phylo-popgen tool map (31 tools)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters.
Servers: **Phylogenetics** (13), **Population** (18).

All phylogenetics Method IDs are on `PhylogeneticAnalyzer`; all population Method IDs on
`PopulationGeneticsAnalyzer` (abbreviated `.Method` below).

## Phylogenetics — distances

| Tool | Method ID | Purpose |
|---|---|---|
| [`pairwise_distance`](../../../../docs/mcp/tools/phylogenetics/pairwise_distance.md) | `.CalculatePairwiseDistance` | Evolutionary distance between two aligned seqs (Hamming / PDistance / JukesCantor / Kimura2Parameter). |
| [`distance_matrix`](../../../../docs/mcp/tools/phylogenetics/distance_matrix.md) | `.CalculateDistanceMatrix` | Symmetric `n×n` distance matrix for aligned seqs; diagonal zero; JC69/K2P may return +Inf at saturation. |

## Phylogenetics — tree construction

| Tool | Method ID | Purpose |
|---|---|---|
| [`build_phylogenetic_tree`](../../../../docs/mcp/tools/phylogenetics/build_phylogenetic_tree.md) | `.BuildTree` | From aligned sequences: distance model + UPGMA/NJ → Newick, taxa, matrix, method. |
| [`build_tree_from_matrix`](../../../../docs/mcp/tools/phylogenetics/build_tree_from_matrix.md) | `.BuildTreeFromMatrix` | From a precomputed symmetric distance matrix + taxa → UPGMA/NJ Newick. |

## Phylogenetics — Newick I/O

| Tool | Method ID | Purpose |
|---|---|---|
| [`to_newick`](../../../../docs/mcp/tools/phylogenetics/to_newick.md) | `.ToNewick` | Serialize a tree to canonical Newick format. |
| [`parse_newick`](../../../../docs/mcp/tools/phylogenetics/parse_newick.md) | `.ParseNewick` | Parse a Newick string → structural summary. |

## Phylogenetics — tree statistics

| Tool | Method ID | Purpose |
|---|---|---|
| [`tree_depth`](../../../../docs/mcp/tools/phylogenetics/tree_depth.md) | `.GetTreeDepth` | Tree height: max edges from root to any leaf. |
| [`tree_length`](../../../../docs/mcp/tools/phylogenetics/tree_length.md) | `.CalculateTreeLength` | Sum of all branch lengths. |
| [`tree_leaves`](../../../../docs/mcp/tools/phylogenetics/tree_leaves.md) | `.GetLeaves` | Enumerate leaf (taxon) nodes with branch lengths. |
| [`mrca`](../../../../docs/mcp/tools/phylogenetics/mrca.md) | `.FindMRCA` | Most Recent Common Ancestor of two taxa (rooted tree). |
| [`patristic_distance`](../../../../docs/mcp/tools/phylogenetics/patristic_distance.md) | `.PatristicDistance` | Sum of branch lengths along the path between two taxa. |

## Phylogenetics — tree comparison & support

| Tool | Method ID | Purpose |
|---|---|---|
| [`robinson_foulds_distance`](../../../../docs/mcp/tools/phylogenetics/robinson_foulds_distance.md) | `.RobinsonFouldsDistance` | Symmetric clade-difference (topology) distance between two rooted trees. |
| [`bootstrap_support`](../../../../docs/mcp/tools/phylogenetics/bootstrap_support.md) | `.Bootstrap` | Non-parametric bootstrap support for clades. |

## Population — allele / genotype frequencies

| Tool | Method ID | Purpose |
|---|---|---|
| [`allele_frequencies`](../../../../docs/mcp/tools/population/allele_frequencies.md) | `.CalculateAlleleFrequencies` | Major/minor allele freqs (p,q) from AA/Aa/aa genotype counts. |
| [`minor_allele_frequency`](../../../../docs/mcp/tools/population/minor_allele_frequency.md) | `.CalculateMAF` | Folded MAF ∈ [0,0.5] from a 0/1/2 genotype vector. |
| [`filter_variants_by_maf`](../../../../docs/mcp/tools/population/filter_variants_by_maf.md) | `.FilterByMAF` | Keep variants whose stored `alleleFrequency`-derived MAF is in `[minMAF,maxMAF]`. |

## Population — Hardy-Weinberg

| Tool | Method ID | Purpose |
|---|---|---|
| [`hardy_weinberg_test`](../../../../docs/mcp/tools/population/hardy_weinberg_test.md) | `.TestHardyWeinberg` | χ² goodness-of-fit HWE test (1 df) from AA/Aa/aa; `inEquilibrium` iff pValue ≥ α. |

## Population — differentiation (F-statistics)

| Tool | Method ID | Purpose |
|---|---|---|
| [`fst`](../../../../docs/mcp/tools/population/fst.md) | `.CalculateFst` | Wright's Fst between two populations (per-locus `{alleleFreq,sampleSize}`). |
| [`f_statistics`](../../../../docs/mcp/tools/population/f_statistics.md) | `.CalculateFStatistics` | Wright's Fis / Fit / Fst between two populations. |
| [`pairwise_fst`](../../../../docs/mcp/tools/population/pairwise_fst.md) | `.CalculatePairwiseFst` | Symmetric pairwise Fst matrix for a set of populations. |

## Population — diversity & neutrality

| Tool | Method ID | Purpose |
|---|---|---|
| [`nucleotide_diversity`](../../../../docs/mcp/tools/population/nucleotide_diversity.md) | `.CalculateNucleotideDiversity` | π: average per-site pairwise differences over aligned seqs. |
| [`wattersons_theta`](../../../../docs/mcp/tools/population/wattersons_theta.md) | `.CalculateWattersonTheta` | Watterson's θ from segregating sites S, sample size n, length L. |
| [`tajimas_d`](../../../../docs/mcp/tools/population/tajimas_d.md) | `.CalculateTajimasD` | Tajima's D from k̂ (NOT per-site), S, n; neutrality test. |
| [`diversity_statistics`](../../../../docs/mcp/tools/population/diversity_statistics.md) | `.CalculateDiversityStatistics` | One pass: π, θ, Tajima's D, S, n, observed/expected heterozygosity. |

## Population — linkage & haplotypes

| Tool | Method ID | Purpose |
|---|---|---|
| [`linkage_disequilibrium`](../../../../docs/mcp/tools/population/linkage_disequilibrium.md) | `.CalculateLD` | D' and r² between two variants from 0/1/2 genotype pairs. |
| [`haplotype_blocks`](../../../../docs/mcp/tools/population/haplotype_blocks.md) | `.FindHaplotypeBlocks` | Detect haplotype blocks from adjacent-variant LD. |
| [`runs_of_homozygosity`](../../../../docs/mcp/tools/population/runs_of_homozygosity.md) | `.FindROH` | Identify runs of homozygosity from per-SNP genotype calls. |
| [`inbreeding_from_roh`](../../../../docs/mcp/tools/population/inbreeding_from_roh.md) | `.CalculateInbreedingFromROH` | Genomic inbreeding coefficient F_ROH from ROH segments. |

## Population — selection & ancestry (advanced)

| Tool | Method ID | Purpose |
|---|---|---|
| [`integrated_haplotype_score`](../../../../docs/mcp/tools/population/integrated_haplotype_score.md) | `.CalculateIHS(ehh0, ehh1, positions)` | iHS from precomputed EHH curves. |
| [`estimate_ancestry`](../../../../docs/mcp/tools/population/estimate_ancestry.md) | `.EstimateAncestry` | Per-individual ancestry proportions vs fixed reference populations. |
| [`scan_selection_signals`](../../../../docs/mcp/tools/population/scan_selection_signals.md) | `.ScanForSelection(regions, thresholds)` | Scan regions for selection using Tajima's D / Fst / iHS thresholds. |

> These three consume **precomputed** curves/references (EHH, reference-population panels, per-region
> statistics); they do not derive them from raw reads. Confirm you have those inputs before calling.

## Envelope

None of these 31 tools are among the 9 `LimitationPolicy`-guarded units
([`docs/Validation/LIMITATIONS.md`](../../../../docs/Validation/LIMITATIONS.md)) — no `MinimumMode` gate
applies to a normal phylo/popgen call. If a task chains into a guarded unit elsewhere, `bio-rigor`'s
envelope rule takes over: stop on a `SeqeronLimitationException` rather than forcing output.
