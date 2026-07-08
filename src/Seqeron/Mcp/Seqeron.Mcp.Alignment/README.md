# Seqeron.Mcp.Alignment

MCP server — **Pairwise and multiple sequence alignment, overlap/assembly and approximate matching.**

Exposes **22 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/alignment/`](../../../../docs/mcp/tools/alignment). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Alignment
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Alignment"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (22)

| Tool | Description |
|------|-------------|
| `alignment_statistics` | Computes match / mismatch / gap counts and percent identity, similarity, and gap percent for a previously produced pairwise alignment. |
| `assemble_de_bruijn` | Assembles input reads using a de Bruijn graph: shreds reads into k-mers, builds the (k-1)-mer node graph, traces non-branching paths into… |
| `assemble_olc` | Assembles input reads into contigs using the overlap–layout–consensus approach: detects pairwise suffix–prefix overlaps above thresholds,… |
| `assembly_stats` | Computes assembly statistics (N50, longest contig, total length, assembled-reads accounting) for a precomputed list of contigs. |
| `calculate_coverage` | Maps each read to its best-matching position on the reference (≥ minOverlap matching bases) and returns per-base coverage depth as an int… |
| `compute_consensus` | Builds a consensus sequence from a list of pre-aligned reads (same length, '-' / 'N' ignored) by majority vote per column. |
| `error_correct_reads` | Corrects single-base errors in reads using k-mer frequency: any k-mer occurring fewer than minKmerFrequency times is corrected by substit… |
| `find_all_overlaps` | Computes all suffix-of-i / prefix-of-j overlaps between read pairs above minOverlap length and minIdentity ratio. |
| `find_best_match` | Returns the single best (minimum-Hamming-distance) fixed-length window of pattern inside sequence. |
| `find_overlap` | Detects the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. |
| `find_with_edits` | Finds all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions) wit… |
| `find_with_mismatches` | Finds all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). |
| `format_alignment` | Renders a human-readable visual alignment (BLAST-style three-line block with '|' for matches, '.' for mismatches, ' ' for gaps), wrapped… |
| `frequent_kmers_with_mismatches` | Finds the most-frequent k-mers within sequence allowing up to d mismatches (counts each k-mer plus all its DNA neighbors at Hamming dista… |
| `global_align` | Performs global pairwise alignment using the Needleman–Wunsch dynamic programming algorithm (aligns full length of both sequences end-to-… |
| `local_align` | Performs local pairwise alignment using the Smith–Waterman algorithm (finds best-scoring substring alignment, with zero floor). |
| `merge_contigs` | Concatenates two contigs, collapsing the specified suffix/prefix overlap of overlapLength bases. |
| `multiple_align` | Anchor-based progressive multiple sequence alignment: picks a center sequence by 4-mer cosine similarity, builds a suffix tree on it, and… |
| `quality_trim_reads` | Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. |
| `scaffold_contigs` | Builds scaffolds by joining contigs using paired-end link records, inserting gapCharacter between linked contigs to span the indicated ga… |
| `semi_global_align` | Performs semi-global (fitting / glocal) pairwise alignment with free end gaps in sequence2; |
| `sequence_identity` | Computes percent-identity (fraction in [0,1]) between two equal-length sequences via case-insensitive position-by-position comparison. |
