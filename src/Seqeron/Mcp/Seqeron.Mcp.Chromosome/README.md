# Seqeron.Mcp.Chromosome

MCP server — **Chromosome-level analysis: karyotype, telomeres, centromere, assembly stats, synteny.**

Exposes **32 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/chromosome/`](../../../../docs/mcp/tools/chromosome). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Chromosome
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Chromosome"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (32)

| Tool | Description |
|------|-------------|
| `analyze_centromere` | Locate the centromere region by alpha-satellite-like repeat content and classify its type per Levan et al. |
| `analyze_karyotype` | Analyze karyotype from chromosome data; |
| `analyze_scaffolds` | Decompose scaffolds into contigs and gaps; |
| `analyze_telomeres` | Analyze 5'/3' telomere repeat tracts on a chromosome sequence; |
| `arm_ratio` | Compute chromosome arm ratio (p/q) from centromere position and total length. |
| `assembly_statistics` | Compute comprehensive assembly statistics: total length, N50/L50/N90/L90, largest/smallest, mean/median, GC, and gap metrics. |
| `assess_completeness` | Assess assembly completeness by aligning marker genes (BUSCO-like, k-mer based); |
| `au_n` | Compute auN (area under Nx curve) — a length-weighted contiguity metric robust to outliers. |
| `classify_chromosome_by_arm_ratio` | Classify a chromosome (Metacentric/Submetacentric/Acrocentric/Telocentric) from its p/q arm ratio per Levan et al. |
| `compare_assemblies` | Compare two assemblies by shared k-mer content; |
| `detect_aneuploidy` | Detect copy-number states from binned read-depth data; |
| `detect_ploidy` | Detect ploidy level from normalized read depth values; |
| `detect_rearrangements` | Detect chromosomal rearrangements (inversions, translocations, deletions, duplications) from synteny blocks. |
| `estimate_cell_divisions_from_telomere_length` | Estimate the number of cell divisions from current telomere length given birth length and per-division loss rate. |
| `estimate_completeness_from_kmers` | Estimate genome completeness, error rate and genome size from a k-mer count spectrum. |
| `estimate_telomere_length_from_ts_ratio` | Estimate telomere length in bp from a qPCR Telomere/Single-copy gene (T/S) ratio against a reference. |
| `extract_contigs` | Extract contigs (gap-free runs) from scaffolds with a minimum length filter. |
| `find_gaps` | Find gaps (N-runs) in assembled sequences; |
| `find_heterochromatin_regions` | Identify heterochromatin regions by k-mer repeat content; |
| `find_repetitive_regions` | Identify repetitive regions across assembled sequences using k-mer copy-number frequency. |
| `find_suspicious_regions` | Flag potentially misassembled regions by GC deviation, low complexity, and high N content. |
| `find_syntenic_blocks_assemblies` | Find syntenic blocks between two assemblies via k-mer anchor clustering; |
| `find_synteny_blocks` | Identify collinear synteny blocks between two genomes from a list of ortholog gene pairs. |
| `find_tandem_repeats` | Identify tandem repeats; |
| `gap_distribution` | Summarize a list of gaps: count, mean/median/max length, and counts by gap-length type. |
| `identify_whole_chromosome_aneuploidy` | Identify whole-chromosome aneuploidies (Monosomy/Trisomy/Tetrasomy/etc.) from per-bin copy-number states. |
| `length_distribution` | Bucket sequence lengths into bins (default: 100..1,000,000 powers). |
| `local_quality` | Compute per-window local quality metrics (GC content, N count, linguistic complexity). |
| `nx_curve` | Compute Nx/Lx for multiple thresholds (default 10..90 step 10). |
| `nx_statistics` | Compute Nx and Lx for a single threshold given pre-sorted (descending) sequence lengths and total length. |
| `predict_g_bands` | Predict cytogenetic G-band pattern from sequence GC content (gpos100/gpos50/gneg, simplified). |
| `repeat_content` | Compute total repeat length, repeat percentage and per-class lengths from repeat annotations. |
