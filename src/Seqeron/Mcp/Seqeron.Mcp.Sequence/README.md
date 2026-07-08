# Seqeron.Mcp.Sequence

MCP server — **DNA/RNA/protein validation, statistics, complexity, k-mers, translation.**

Exposes **35 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/sequence/`](../../../../docs/mcp/tools/sequence). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Sequence
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Sequence"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (35)

| Tool | Description |
|------|-------------|
| `amino_acid_composition` | Calculate amino acid composition, molecular weight, and other properties of a protein sequence. |
| `complement_base` | Get the Watson-Crick complement of a single nucleotide base (A↔T, C↔G for DNA; |
| `complexity_compression_ratio` | Estimate sequence complexity using compression ratio. |
| `complexity_dust_score` | Calculate DUST score for low-complexity filtering (as used in BLAST). |
| `complexity_kmer_entropy` | Calculate k-mer based Shannon entropy for DNA complexity analysis. |
| `complexity_linguistic` | Calculate DNA linguistic complexity as ratio of observed to possible subwords. |
| `complexity_mask_low` | Mask low-complexity regions in a DNA sequence using the DUST algorithm. |
| `complexity_shannon` | Calculate DNA Shannon entropy (bits per base). |
| `dna_reverse_complement` | Get the reverse complement of a DNA sequence. |
| `dna_validate` | Validate a DNA sequence. |
| `gc_content` | Calculate the GC content (percentage of G and C nucleotides) of a DNA/RNA sequence. |
| `hydrophobicity` | Calculate the grand average of hydropathy (GRAVY) index of a protein sequence. |
| `is_valid_dna` | Quick check if a sequence contains only valid DNA characters (A, T, G, C). |
| `is_valid_rna` | Quick check if a sequence contains only valid RNA characters (A, U, G, C). |
| `isoelectric_point` | Calculate the isoelectric point (pI) of a protein sequence. |
| `iupac_code` | Get the IUPAC ambiguity code that represents a set of nucleotide bases. |
| `iupac_match` | Check if two IUPAC codes can represent the same nucleotide base. |
| `iupac_matches` | Check if a specific nucleotide matches an IUPAC ambiguity code. |
| `kmer_analyze` | Comprehensive k-mer analysis including statistics about frequency distribution, entropy, and unique k-mers. |
| `kmer_count` | Count k-mer (substring of length k) frequencies in a sequence. |
| `kmer_distance` | Calculate k-mer based distance between two sequences using Euclidean distance of k-mer frequencies. |
| `kmer_entropy` | Calculate Shannon entropy based on k-mer frequencies. |
| `linguistic_complexity` | Calculate linguistic complexity of a sequence based on k-mer diversity. |
| `melting_temperature` | Calculate the melting temperature (Tm) of a DNA sequence using Wallace rule or GC formula. |
| `molecular_weight_nucleotide` | Calculate the molecular weight of a DNA or RNA sequence in Daltons (Da). |
| `molecular_weight_protein` | Calculate the molecular weight of a protein sequence in Daltons (Da). |
| `nucleotide_composition` | Calculate nucleotide composition (A, T, G, C, U counts) and GC content of a DNA/RNA sequence. |
| `protein_validate` | Validate a protein sequence. |
| `rna_from_dna` | Transcribe DNA to RNA by replacing T (thymine) with U (uracil). |
| `rna_validate` | Validate an RNA sequence. |
| `shannon_entropy` | Calculate Shannon entropy of a sequence. |
| `summarize_sequence` | Generate comprehensive summary statistics for a DNA/RNA sequence including composition, GC content, entropy, complexity, and Tm. |
| `thermodynamics` | Calculate thermodynamic properties (ΔH, ΔS, ΔG, Tm) of a DNA duplex using the nearest-neighbor method. |
| `translate_dna` | Translate a DNA sequence to protein using the standard genetic code. |
| `translate_rna` | Translate an RNA sequence to protein using the standard genetic code. |
