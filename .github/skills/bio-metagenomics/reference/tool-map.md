# bio-metagenomics tool map (19 tools)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters. Server: **Metagenomics** (19).
Method families: `MetagenomicsAnalyzer` (10) and `PanGenomeAnalyzer` (9).

## Classification (build a database, then classify)

| Tool | Method ID | Purpose |
|---|---|---|
| [`build_kmer_database`](../../../../docs/mcp/tools/metagenomics/build_kmer_database.md) | `MetagenomicsAnalyzer.BuildKmerDatabase` | Build Kraken canonical-k-mer → taxon-id database; shared k-mers collapse to taxonomy LCA (`k` default 31). |
| [`classify_reads`](../../../../docs/mcp/tools/metagenomics/classify_reads.md) | `MetagenomicsAnalyzer.ClassifyReads` | Kraken k-mer / LCA classification of reads → per-read taxon, RTL score, Kraken2 C/Q confidence. |

## Profiling & abundance

| Tool | Method ID | Purpose |
|---|---|---|
| [`taxonomic_profile`](../../../../docs/mcp/tools/metagenomics/taxonomic_profile.md) | `MetagenomicsAnalyzer.GenerateTaxonomicProfile` | Aggregate per-read classifications → relative abundance per rank + species-level Shannon/Simpson; excludes unclassified from denominator. |
| [`differential_abundance`](../../../../docs/mcp/tools/metagenomics/differential_abundance.md) | `MetagenomicsAnalyzer.DifferentialAbundance` | Per-taxon differential abundance between two condition groups. |

## Diversity

| Tool | Method ID | Purpose |
|---|---|---|
| [`alpha_diversity`](../../../../docs/mcp/tools/metagenomics/alpha_diversity.md) | `MetagenomicsAnalyzer.CalculateAlphaDiversity` | Within-sample: Shannon, Simpson, inverse Simpson, Chao1, observed species, Pielou evenness. Accepts counts or fractions. |
| [`beta_diversity`](../../../../docs/mcp/tools/metagenomics/beta_diversity.md) | `MetagenomicsAnalyzer.CalculateBetaDiversity` | Between two samples: Bray-Curtis, Jaccard (presence/absence), shared/unique species. UniFrac=0 (no tree). |
| [`functional_diversity`](../../../../docs/mcp/tools/metagenomics/functional_diversity.md) | `MetagenomicsAnalyzer.CalculateFunctionalDiversity` | Functional richness + Shannon functional diversity over predicted functions. |

## Function & resistance

| Tool | Method ID | Purpose |
|---|---|---|
| [`predict_functions`](../../../../docs/mcp/tools/metagenomics/predict_functions.md) | `MetagenomicsAnalyzer.PredictFunctions` | Predict functional annotations for proteins. |
| [`find_resistance_genes`](../../../../docs/mcp/tools/metagenomics/find_resistance_genes.md) | `MetagenomicsAnalyzer.FindResistanceGenes` | Search genes for antibiotic-resistance markers. |

## Binning (⚠ META-BIN-001 — guarded, min mode `Moderate`)

| Tool | Method ID | Purpose |
|---|---|---|
| [`bin_contigs`](../../../../docs/mcp/tools/metagenomics/bin_contigs.md) | `MetagenomicsAnalyzer.BinContigs` | k-means bin contigs into MAGs over GC/coverage/TETRA; reports completeness & contamination (domain-level CheckM). **Throws under `Strict`** — see STOP rule in `SKILL.md`. |

## Pan-genome (PanGenomeAnalyzer — multi-genome comparative)

| Tool | Method ID | Purpose |
|---|---|---|
| [`construct_pangenome`](../../../../docs/mcp/tools/metagenomics/construct_pangenome.md) | `PanGenomeAnalyzer.ConstructPanGenome` | Construct a pan-genome from a set of genomes. |
| [`cluster_genes`](../../../../docs/mcp/tools/metagenomics/cluster_genes.md) | `PanGenomeAnalyzer.ClusterGenes` | Cluster genes from multiple genomes into ortholog groups. |
| [`core_gene_clusters`](../../../../docs/mcp/tools/metagenomics/core_gene_clusters.md) | `PanGenomeAnalyzer.GetCoreGeneClusters` | Filter gene clusters down to the core set. |
| [`accessory_genes`](../../../../docs/mcp/tools/metagenomics/accessory_genes.md) | `PanGenomeAnalyzer.AnalyzeAccessoryGenes` | Summarise accessory gene clusters of a pan-genome. |
| [`core_genome_alignment`](../../../../docs/mcp/tools/metagenomics/core_genome_alignment.md) | `PanGenomeAnalyzer.CreateCoreGenomeAlignment` | Build a per-genome core-genome alignment block. |
| [`gene_presence_absence_matrix`](../../../../docs/mcp/tools/metagenomics/gene_presence_absence_matrix.md) | `PanGenomeAnalyzer.CreatePresenceAbsenceMatrix` | Build a per-genome gene presence/absence matrix. |
| [`find_genome_specific_genes`](../../../../docs/mcp/tools/metagenomics/find_genome_specific_genes.md) | `PanGenomeAnalyzer.FindGenomeSpecificGenes` | List cluster ids unique to each genome. |
| [`select_phylogenetic_markers`](../../../../docs/mcp/tools/metagenomics/select_phylogenetic_markers.md) | `PanGenomeAnalyzer.SelectPhylogeneticMarkers` | Select single-copy core clusters as phylogenetic markers. |
| [`fit_heaps_law`](../../../../docs/mcp/tools/metagenomics/fit_heaps_law.md) | `PanGenomeAnalyzer.FitHeapsLaw` | Fit Heaps' law to the pan-genome new-gene-discovery curve (open vs closed pan-genome). |

> For pairwise/MSA alignment of markers use **`bio-alignment`**; for phylogenetic tree building from
> selected markers use **`bio-phylo-popgen`**. This skill owns the metagenomics/community subset above.
