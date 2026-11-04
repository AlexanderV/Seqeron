# Seqeron.Mcp.Alignment

MCP server — **Pairwise and multiple sequence alignment, overlap/assembly and approximate matching.**

Exposes **34 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/alignment/`](../../../../docs/mcp/tools/alignment). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Alignment
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Alignment"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (34)

| Tool | Description |
|------|-------------|
| `alignment_statistics` | Computes match / mismatch / gap counts and percent identity, similarity, and gap percent for a previously produced pairwise alignment. |
| `assemble_de_bruijn` | Assembles input reads using a de Bruijn graph: shreds reads into k-mers, builds the (k-1)-mer node graph, traces non-branching paths into… |
| `assemble_olc` | Assembles input reads into contigs using the overlap–layout–consensus approach: detects pairwise suffix–prefix overlaps above thresholds,… |
| `assembly_stats` | Computes assembly statistics (N50, longest contig, total length, assembled-reads accounting) for a precomputed list of contigs. |
| `calculate_coverage` | Maps each read to its best-matching position on the reference (≥ minOverlap matching bases) and returns per-base coverage depth as an int… |
| `compute_consensus` | Builds a consensus sequence from a list of pre-aligned reads (same length, '-' / 'N' ignored) by majority vote per column. |
| `damerau_alignment` | Optimal transposition-aware edit script turning sequence1 into sequence2 (Lowrance & Wagner 1975 trace). |
| `damerau_levenshtein_distance` | Edit distance with adjacent transpositions. |
| `edit_alignment` | Optimal global unit-cost (Levenshtein) alignment of query against target in edlib's convention: operations '=' match, 'X' mismatch, 'I' que… |
| `error_correct_reads` | Corrects single-base errors in reads using k-mer frequency: any k-mer occurring fewer than minKmerFrequency times is corrected by substit… |
| `find_all_overlaps` | Computes all suffix-of-i / prefix-of-j overlaps between read pairs above minOverlap length and minIdentity ratio. |
| `find_best_match` | Returns the single best (minimum-Hamming-distance) fixed-length window of pattern inside sequence. |
| `find_edit_end_positions` | Sellers (1980) k-differences search: reports every 0-based end position j in sequence at which some substring ending at j is within maxEdit… |
| `find_overlap` | Detects the longest suffix-of-sequence1 / prefix-of-sequence2 overlap satisfying minOverlap and minIdentity thresholds. |
| `find_with_edits` | Finds all approximate matches of pattern in sequence allowing up to maxEdits Levenshtein edits (insertions, deletions, substitutions) wit… |
| `find_with_mismatches` | Finds all occurrences of pattern inside sequence allowing up to maxMismatches substitutions (Hamming-style, fixed-length window). |
| `format_alignment` | Renders a human-readable visual alignment (BLAST-style three-line block with '|' for matches, '.' for mismatches, ' ' for gaps), wrapped… |
| `frequent_kmers_with_mismatches` | Finds the most-frequent k-mers within sequence allowing up to d mismatches (counts each k-mer plus all its DNA neighbors at Hamming dista… |
| `frequent_kmers_with_mismatches_and_revcomp` | Frequent words with mismatches and reverse complements (ROSALIND BA1J): all DNA k-mers P maximising Count_d(sequence, P) + Count_d(sequence… |
| `global_align` | Global pairwise alignment (Needleman-Wunsch): align the full length of both sequences end-to-end. Linear gap model: every gap position scor… |
| `global_align_affine` | Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (BLAST+ / Biopython convention… |
| `local_align` | Local pairwise alignment (Smith-Waterman, zero floor). Linear gap model: every gap position scores gapExtend; gapOpen is IGNORED (use local… |
| `local_align_affine` | Local pairwise alignment with affine gap costs (Smith-Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; … |
| `merge_contigs` | Concatenates two contigs, collapsing the specified suffix/prefix overlap of overlapLength bases. |
| `multiple_align` | Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total o… |
| `multiple_align_anchored` | Approximate, faster center-star MSA: 4-mer count cosine center, suffix-tree exact-match anchors (MUMmer-style seed, weighted colinear chain… |
| `multiple_align_consistency` | Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): primary library of global (Needleman-Wunsch) and local (Smith-Waterman… |
| `multiple_align_iterative` | Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge spli… |
| `multiple_align_progressive` | Guide-tree progressive MSA (Feng & Doolittle 1987 / CLUSTAL): ClustalW pairwise distances from all Needleman-Wunsch alignments (1 - identit… |
| `quality_trim_reads` | Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. |
| `scaffold_contigs` | Builds scaffolds by joining contigs using paired-end link records, inserting gapCharacter between linked contigs to span the indicated ga… |
| `semi_global_align` | Semi-global (fitting) pairwise alignment with free end gaps in sequence2. Linear gap model: every gap position scores gapExtend; gapOpen is… |
| `semi_global_align_affine` | Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing seq… |
| `sequence_identity` | Computes percent-identity (fraction in [0,1]) between two equal-length sequences via case-insensitive position-by-position comparison. |
