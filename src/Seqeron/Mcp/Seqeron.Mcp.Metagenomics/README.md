# Seqeron.Mcp.Metagenomics

MCP server — **Taxonomic classification, diversity, pan-genome and functional profiling.**

Exposes **19 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/metagenomics/`](../../../../docs/mcp/tools/metagenomics). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Metagenomics
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Metagenomics"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (19)

| Tool | Description |
|------|-------------|
| `accessory_genes` | Summarise accessory clusters (present in >1 but not all genomes) with their genome membership and frequency = genomeCount / totalGenomes. |
| `alpha_diversity` | Compute Shannon, Simpson, inverse Simpson, Chao1 (Chao 1984), observed species, and Pielou's evenness for a sample's abundance vector. |
| `beta_diversity` | Compute Bray–Curtis and Jaccard distances plus shared/unique species counts between two samples. |
| `bin_contigs` | Cluster contigs into metagenome-assembled genomes using k-means over GC content, normalized coverage and tetranucleotide frequency (TETRA… |
| `build_kmer_database` | Build a Kraken canonical-k-mer→taxon-id database from labeled reference sequences. |
| `classify_reads` | Classify metagenomic reads with the Kraken k-mer/LCA algorithm against a canonical-k-mer→taxon-id database and a taxonomy tree. |
| `cluster_genes` | Cluster genes from multiple genomes into ortholog groups using 7-mer Jaccard similarity. |
| `construct_pangenome` | Construct a pan-genome (core/accessory/unique partition, genome fluidity, open-vs-closed classification per Tettelin 2005) from a set of… |
| `core_gene_clusters` | Filter gene clusters down to the core set: those present in at least floor(threshold * totalGenomes) genomes. |
| `core_genome_alignment` | Concatenate a single genome's representative sequences for the supplied core clusters into a per-genome alignment block. |
| `differential_abundance` | Welch's t-test for differential taxon abundance between two condition groups. |
| `find_genome_specific_genes` | For each genome, list the cluster ids that occur only in that genome (singleton accessory clusters). |
| `find_resistance_genes` | Search nucleotide genes for antibiotic-resistance markers via motif containment against a resistance database; |
| `fit_heaps_law` | Fit Heaps' law n(N) = Intercept · N^(-Alpha) to the permuted new-gene-discovery curve (Tettelin 2008; |
| `functional_diversity` | Compute functional richness, Shannon functional diversity, and per-pathway hit counts from a set of functional annotations. |
| `gene_presence_absence_matrix` | Build a per-genome gene presence/absence matrix against a list of gene clusters. |
| `predict_functions` | Predict functional annotations for proteins by motif containment against a function database, returning function/pathway/KO and inferred… |
| `select_phylogenetic_markers` | Pick phylogenetic markers: single-copy core clusters (present in all genomes with exactly one gene each) with >= 1 parsimony-informative… |
| `taxonomic_profile` | Aggregate per-read classifications into rank-wise abundances plus Shannon/Simpson diversity (computed at species level). |
