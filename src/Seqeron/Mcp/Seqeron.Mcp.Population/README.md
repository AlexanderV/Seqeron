# Seqeron.Mcp.Population

MCP server — **Population genetics: allele frequencies, diversity, Fst, LD, selection scans, ROH.**

Exposes **18 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/population/`](../../../../docs/mcp/tools/population). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Population
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Population"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (18)

| Tool | Description |
|------|-------------|
| `allele_frequencies` | Calculate major and minor allele frequencies from diploid genotype counts (homozygous-major, heterozygous, homozygous-minor). |
| `diversity_statistics` | Compute combined diversity statistics (π, Watterson's θ, Tajima's D, segregating sites, observed/expected heterozygosity) from aligned se… |
| `estimate_ancestry` | Estimate per-individual ancestry proportions across reference populations via a simplified ADMIXTURE-like EM procedure. |
| `f_statistics` | Compute Wright's F-statistics (Fis, Fit, Fst) between two populations from per-variant heterozygosities and allele frequencies. |
| `filter_variants_by_maf` | Filter variants whose minor allele frequency lies within [minMAF, maxMAF]. |
| `fst` | Compute Wright's variance-based Fst between two populations from per-variant allele frequencies and sample sizes (Wright 1965). |
| `haplotype_blocks` | Detect haplotype blocks via adjacent-pair r² ≥ ldThreshold (simplified Gabriel et al. |
| `hardy_weinberg_test` | Test Hardy-Weinberg equilibrium for a biallelic variant via chi-square goodness-of-fit (1 df) on observed AA/Aa/aa counts. |
| `inbreeding_from_roh` | Estimate the genomic inbreeding coefficient F_ROH as Σ(ROH lengths) / genomeLength. |
| `integrated_haplotype_score` | Compute integrated haplotype score iHS = ln(iHH₁ / iHH₀) by trapezoidal integration of EHH curves over genomic positions. |
| `linkage_disequilibrium` | Compute pairwise linkage disequilibrium (D' and r²) between two variants from diploid genotype pairs (0/1/2 encoding). |
| `minor_allele_frequency` | Compute minor allele frequency (MAF) from a vector of diploid genotypes encoded as 0 (hom-ref), 1 (het), or 2 (hom-alt). |
| `nucleotide_diversity` | Compute nucleotide diversity π — average pairwise per-site differences across an aligned set of equal-length sequences. |
| `pairwise_fst` | Compute pairwise Fst for a set of populations; |
| `runs_of_homozygosity` | Identify runs of homozygosity (ROH) from per-SNP genotype calls (0/1/2) using minimum SNP count, minimum length, and maximum heterozygote… |
| `scan_selection_signals` | Scan genomic regions for selection signals using Tajima's D, Fst, and iHS thresholds. |
| `tajimas_d` | Compute Tajima's D (Tajima 1989) from average pairwise differences k̂, segregating sites S, and sample size n (n ≥ 3 required). |
| `wattersons_theta` | Compute Watterson's θ estimator from segregating site count, sample size, and sequence length. |
