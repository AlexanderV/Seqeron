# Seqeron.Mcp.Alignment

MCP server — **Pairwise and multiple sequence alignment, overlap/assembly and approximate matching.**

Exposes **40 tools** — the same validated `Seqeron.Genomics` algorithms as the C# API, callable over
MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/alignment/`](../../../../docs/mcp/tools/alignment). Rollout status:
[`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Alignment
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Alignment"]`). New to MCP? The [hub guide](../../../../docs/mcp/README.md) lists all 11 servers and how to wire them up.

## Tools (35)

| Tool | Description |
|------|-------------|
| `alignment_statistics` | Identity / similarity / gap statistics of a pairwise alignment (EMBOSS needle convention: percentages over the full alignment length incl. gap columns); similarity = positive-score columns, optionally with `matrix` = BLOSUM62 / PAM250 / EDNAFULL; also `positives` and `gapOpens`. |
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
| `format_alignment` | Renders a pairwise alignment in the EMBOSS needle/water three-line style (markup `\|` identity, `:` similar, `.` mismatch, space gap), optionally with `matrix` = BLOSUM62 / PAM250 / EDNAFULL, wrapped to a configurable line width. |
| `frequent_kmers_with_mismatches` | Finds the most-frequent k-mers within sequence allowing up to d mismatches (counts each k-mer plus all its DNA neighbors at Hamming dista… |
| `frequent_kmers_with_mismatches_and_revcomp` | Frequent words with mismatches and reverse complements (ROSALIND BA1J): all DNA k-mers P maximising Count_d(sequence, P) + Count_d(sequence… |
| `general_gap_align` | Pairwise alignment with an arbitrary gap-cost function (Waterman–Smith–Beyer 1976): the score of a gap of length k is gapScores[k-1] (concave, convex, capped, step, ...); global / local / semi_global (endGaps), scalar or substitution-matrix scoring; Biopython callable gap_score oracle. |
| `global_align` | Global pairwise alignment (Needleman-Wunsch): align the full length of both sequences end-to-end. Linear gap model: every gap position scor… |
| `global_align_affine` | Global pairwise alignment with affine gap costs (Gotoh 1982): a gap of length k scores gapOpen + k·gapExtend (BLAST+ / Biopython convention… |
| `global_align_co_optimal` | ALL distinct co-optimal GLOBAL alignments (Biopython `PairwiseAligner(mode='global').align` set; linear or `affine` gaps, scalar or `matrix` scoring, any `endGaps` policy), capped at `maxAlignments`. |
| `local_align` | Local pairwise alignment (Smith-Waterman, zero floor). Linear gap model: every gap position scores gapExtend; gapOpen is IGNORED (use local… |
| `local_align_affine` | Local pairwise alignment with affine gap costs (Smith-Waterman 1981 with Gotoh 1982 gaps): a gap of length k scores gapOpen + k·gapExtend; … |
| `local_align_top_n` | Waterman-Eggert (1987) top-N NON-INTERSECTING local alignments (EMBOSS matcher -alternatives / lalign declumping): best alignment, then the best one sharing no aligned residue pair with it, ... up to `maxAlignments` of score >= `minScore`; linear or `affine` gaps, scalar or `matrix` scoring. |
| `local_align_co_optimal` | ALL distinct co-optimal local alignments (Biopython `PairwiseAligner(mode='local').align` set), capped at `maxAlignments`; linear or `affine` gaps, scalar or `matrix` scoring. |
| `merge_contigs` | Concatenates two contigs, collapsing the specified suffix/prefix overlap of overlapLength bases. |
| `multiple_align` | Center-star multiple sequence alignment (Gusfield 1993) with exact pairwise alignments: the center is the sequence with the highest total o… |
| `multiple_align_anchored` | Approximate, faster center-star MSA: 4-mer count cosine center, suffix-tree exact-match anchors (MUMmer-style seed, weighted colinear chain… |
| `multiple_align_clustalw` | ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN), ported step by step and equal row-for-row to the clustalw binary: pairwise stage (full DP or fast k-tuple), ClustalW neighbour-joining guide tree, tree-derived sequence weights, delayed divergent sequences, position-specific gap penalties, IUB matrix with `transitionWeight`, Myers-Miller profile alignment (`-GAPOPEN -GAPEXT -TRANSWEIGHT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE`). |
| `multiple_align_consistency` | Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): primary library of global (Needleman-Wunsch) and local (Smith-Waterman… |
| `multiple_align_exact` | EXACT sum-of-pairs optimal MSA (linear or natural affine gaps): k-dimensional dynamic program with the Carrillo-Lipman bound, at most 8 (linear) / 6 (affine) sequences, complexity guard `maxStateCells`; k = 2 is the pairwise optimum. |
| `multiple_align_iterative` | Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge spli… |
| `multiple_align_progressive` | Guide-tree progressive MSA (Feng & Doolittle 1987 / CLUSTAL): ClustalW pairwise distances from all Needleman-Wunsch alignments (1 - identit… |
| `quality_trim_reads` | Trims Phred+33 quality-encoded reads from both ends, dropping bases whose decoded quality is below minQuality. |
| `scaffold_contigs` | Builds scaffolds by joining contigs using paired-end link records, inserting gapCharacter between linked contigs to span the indicated ga… |
| `semi_global_align` | Semi-global (fitting) pairwise alignment with free end gaps in sequence2. Linear gap model: every gap position scores gapExtend; gapOpen is… |
| `semi_global_align_affine` | Semi-global (fitting) pairwise alignment with affine gap costs: sequence1 is aligned end to end inside sequence2 (free leading/trailing seq… |
| `sequence_identity` | Computes percent-identity (fraction in [0,1]) between two equal-length sequences via case-insensitive position-by-position comparison. |

**Protein / matrix multiple alignment (B09 F61, audit A1-14).** Every `multiple_align*` tool except `multiple_align_clustalw` takes `type = dna | protein` and `matrix = BLOSUM62 | PAM250 | EDNAFULL` (additive; `type=protein` defaults to BLOSUM62, EDNAFULL gives IUPAC DNA) and then scores residue pairs with the published substitution table and real-valued (fractional, EMBOSS-style 0.5) gap scores, linear or affine (`gapModel`); `gapOpen` / `gapExtend` are numbers (integers unless a matrix is used). The protein path is the same code as the DNA path (the string / `ProteinSequence` overloads of `SequenceAligner.MultipleAlign*` with `MatrixScoring`); results carry `exactTotalScore` (the real sum of pairs; `totalScore` is it rounded half away from zero) and, for `multiple_align_exact`, `scoreScale`. Verified against independent Python implementations of every method (238 protein / IUPAC-DNA sets, rows equal), Biopython 1.88 pairwise optima, exhaustive enumeration of tiny sets and the `t_coffee` program (160 protein sets for the progressive stage).

**ClustalW protein pairwise stages (B09 F62, audit A1-15b b).** With `type=protein`, `distance=clustalw` (full DP) and `distance=clustalw-quick` (`-QUICKTREE`) of `multiple_align_progressive`, `multiple_align_iterative` and `multiple_align_consistency` use ClustalW 2.1's PROTEIN pairwise stages (exact ports; the printed `Score:` equals the installed `clustalw -TYPE=PROTEIN` on 4 204 + 3 774 pairs): matrix `pwMatrix = gonnet (default) | blosum | pam | id`, `pwGapOpen` 10 / `pwGapExtend` 0.1 (length-dependent opening), fast stage k-tuple 1-2 with defaults 1 / 5 / 5 / 3 (`ktuple`, `window`, `topDiagonals`, `pairGap`). `multiple_align_clustalw` itself (the profile stage) remains DNA only.

**Memory and length limits (B09 F46).** The six pairwise tools (`global_align`, `local_align`, `semi_global_align` and their `_affine` variants, scalar and `matrix` scoring, every `endGaps` policy) switch automatically to exact linear-space programs (Hirschberg / Myers–Miller / Huang–Miller locate-then-solve) above 25 000 000 dynamic-programming cells; the output is identical to the full-matrix result, the memory is O(len1 + len2) (20 kb × 20 kb affine: 111 MB peak, 100 kb × 100 kb: 185 MB, 6.5 min on one core). `general_gap_align` with a non-affine gap table is limited to 4 000 000 cells (clear error above); linear/affine tables are delegated and unlimited. Progressive / iterative / consistency MSA tools throw above 50 000 000 profile-column cells; `multiple_align` (exact center-star) is linear-space.
