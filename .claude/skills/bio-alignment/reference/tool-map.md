# bio-alignment tool map (44 tools)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters. Servers: **Alignment** (26), **Core** (15).

## Pairwise alignment (Alignment server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`global_align`](../../../../docs/mcp/tools/alignment/global_align.md) | `SequenceAligner.GlobalAlign` | Needleman-Wunsch; align full length of both seqs end-to-end. |
| [`local_align`](../../../../docs/mcp/tools/alignment/local_align.md) | `SequenceAligner.LocalAlign` | Smith-Waterman; best-scoring shared substring, flanks ignored. |
| [`semi_global_align`](../../../../docs/mcp/tools/alignment/semi_global_align.md) | `SequenceAligner.SemiGlobalAlign` | Fitting/glocal; free end gaps in seq2 — fit a query into a reference / overlap. |

## Alignment reporting

| Tool | Method ID | Purpose |
|---|---|---|
| [`alignment_statistics`](../../../../docs/mcp/tools/alignment/alignment_statistics.md) | `SequenceAligner.CalculateStatistics` | matches/mismatches/gaps + identity/similarity/gap % (denom = full aln length incl. gaps). |
| [`format_alignment`](../../../../docs/mcp/tools/alignment/format_alignment.md) | `SequenceAligner.FormatAlignment` | Three-line human block (`\|`=identity, `:`=similar, ` `=gap/mismatch), wrapped. |
| [`sequence_identity`](../../../../docs/mcp/tools/alignment/sequence_identity.md) | `SequenceAssembler.CalculateIdentity` | % identity of two **equal-length** seqs (gapless denom); 0 if lengths differ. |

## Multiple sequence alignment & consensus

| Tool | Method ID | Purpose |
|---|---|---|
| [`multiple_align`](../../../../docs/mcp/tools/alignment/multiple_align.md) | `SequenceAligner.MultipleAlign` | Anchor-based progressive MSA → aligned rows, consensus, sum-of-pairs score. |
| [`compute_consensus`](../../../../docs/mcp/tools/alignment/compute_consensus.md) | `SequenceAssembler.ComputeConsensus` | Majority-vote consensus from pre-aligned equal-length reads (ties → `N`). |

## Approximate / fuzzy matching (Alignment server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`find_with_mismatches`](../../../../docs/mcp/tools/alignment/find_with_mismatches.md) | `ApproximateMatcher.FindWithMismatches` | All matches of a pattern up to `maxMismatches` (Hamming, fixed window). |
| [`find_with_edits`](../../../../docs/mcp/tools/alignment/find_with_edits.md) | `ApproximateMatcher.FindWithEdits` | All matches up to `maxEdits` Levenshtein edits (variable-length windows). |
| [`find_best_match`](../../../../docs/mcp/tools/alignment/find_best_match.md) | `ApproximateMatcher.FindBestMatch` | Single best (min Hamming) window of pattern in sequence; leftmost tie. |
| [`frequent_kmers_with_mismatches`](../../../../docs/mcp/tools/alignment/frequent_kmers_with_mismatches.md) | `ApproximateMatcher.FindFrequentKmersWithMismatches` | Most-frequent k-mers allowing up to d mismatches (neighborhood tally). |
| [`frequent_kmers_with_mismatches_and_revcomp`](../../../../docs/mcp/tools/alignment/frequent_kmers_with_mismatches_and_revcomp.md) | `ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements` | BA1J: most-frequent k-mers counting d-mismatch hits of the k-mer and its reverse complement. |
| [`find_edit_end_positions`](../../../../docs/mcp/tools/alignment/find_edit_end_positions.md) | `ApproximateMatcher.FindEditEndPositions` | Sellers k-differences: every end position within `maxEdits` + its minimum edit distance (optional weighted costs). |

## Similarity & distance (Core server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`edit_distance`](../../../../docs/mcp/tools/core/edit_distance.md) | `ApproximateMatcher.EditDistance` | Levenshtein distance (any lengths; Myers bit-parallel); optional insertion/deletion/substitution costs → weighted distance (rapidfuzz weights). |
| [`edit_alignment`](../../../../docs/mcp/tools/alignment/edit_alignment.md) | `ApproximateMatcher.GetEditAlignment` | Optimal Levenshtein alignment: `=`/`X`/`I`/`D` ops, edlib CIGAR, gapped strings (`linearSpace` → Hirschberg; optional weighted costs). Alignment server. |
| [`damerau_levenshtein_distance`](../../../../docs/mcp/tools/alignment/damerau_levenshtein_distance.md) | `ApproximateMatcher.DamerauLevenshteinDistance` | Edit distance with adjacent transpositions (`unrestricted` DL or `osa`); optional weighted ins/del/sub/transposition costs (Lowrance–Wagner). Alignment server. |
| [`damerau_alignment`](../../../../docs/mcp/tools/alignment/damerau_alignment.md) | `ApproximateMatcher.GetDamerauLevenshteinAlignment` / `GetOptimalStringAlignment` | Transposition-aware edit script (Lowrance–Wagner trace: `=`/`X`/`I`/`D`/`T` blocks with positions and costs), unit or weighted costs. Alignment server. |
| [`hamming_distance`](../../../../docs/mcp/tools/core/hamming_distance.md) | `ApproximateMatcher.HammingDistance` | Mismatch count for **equal-length** seqs. |
| [`calculate_similarity`](../../../../docs/mcp/tools/core/calculate_similarity.md) | `GenomicAnalyzer.CalculateSimilarity` | k-mer Jaccard similarity in [0,1] (optional `kmerSize`). |
| [`count_approximate_occurrences`](../../../../docs/mcp/tools/core/count_approximate_occurrences.md) | `ApproximateMatcher.CountApproximateOccurrences` | Count approximate occurrences of a pattern allowing mismatches. |

## Common regions & repeats (Core server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`find_longest_common_region`](../../../../docs/mcp/tools/core/find_longest_common_region.md) | `GenomicAnalyzer.FindLongestCommonRegion` | Longest common substring of two DNA seqs **+ positions in both**. |
| [`suffix_tree_lcs`](../../../../docs/mcp/tools/core/suffix_tree_lcs.md) | `SuffixTree.LongestCommonSubstring` | Longest common substring of two raw texts (`substring`, `length`). |
| [`find_longest_repeat`](../../../../docs/mcp/tools/core/find_longest_repeat.md) | `GenomicAnalyzer.FindLongestRepeat` | Longest repeated region in a DNA sequence. |
| [`suffix_tree_lrs`](../../../../docs/mcp/tools/core/suffix_tree_lrs.md) | `SuffixTree.LongestRepeatedSubstring` | Longest repeated substring in a text. |
| [`suffix_tree_all_lrs`](../../../../docs/mcp/tools/core/suffix_tree_all_lrs.md) | `SuffixTree.FindAllLongestRepeatedSubstrings` | Every tied longest repeated substring with all positions. |
| [`suffix_tree_find_mems`](../../../../docs/mcp/tools/core/suffix_tree_find_mems.md) | `SuffixTree.FindMaximalExactMatches` | All maximal exact matches ≥ minLength between reference and query (MUMmer `-maxmatch`). |
| [`suffix_tree_find_mums`](../../../../docs/mcp/tools/core/suffix_tree_find_mums.md) | `SuffixTree.FindMaximalUniqueMatches` | Maximal unique matches (MUMmer `-mum` / `-mumreference`). |
| [`suffix_tree_all_lcs`](../../../../docs/mcp/tools/core/suffix_tree_all_lcs.md) | `SuffixTree.FindAllDistinctLongestCommonSubstrings` | Every tied longest common substring of two texts with all positions in both. |
| [`suffix_tree_k_common_substrings`](../../../../docs/mcp/tools/core/suffix_tree_k_common_substrings.md) | `SuffixTree.FindLongestCommonSubstrings` | Longest substring(s) common to k strings or ≥ minSupport of them (Rosalind LCSM; l(q) for every q). |
| [`suffix_tree_maximal_repeats`](../../../../docs/mcp/tools/core/suffix_tree_maximal_repeats.md) | `SuffixTree.FindMaximalRepeatedPairs` | All maximal repeated pairs ≥ minLength in one text (MUMmer `repeat-match -f`; optional unique symbols such as N). |

## Exact substring queries (Core suffix tree)

| Tool | Method ID | Purpose |
|---|---|---|
| [`suffix_tree_contains`](../../../../docs/mcp/tools/core/suffix_tree_contains.md) | `SuffixTree.Contains` | Does pattern occur in text? |
| [`suffix_tree_count`](../../../../docs/mcp/tools/core/suffix_tree_count.md) | `SuffixTree.CountOccurrences` | How many times does pattern occur? |
| [`suffix_tree_find_all`](../../../../docs/mcp/tools/core/suffix_tree_find_all.md) | `SuffixTree.FindAllOccurrences` | All 0-based positions of pattern in text. |
| [`suffix_tree_stats`](../../../../docs/mcp/tools/core/suffix_tree_stats.md) | `SuffixTree.Properties` | Structural stats of a suffix tree built from text. |

## Assembly-adjacent overlap tools (Alignment server — see `bio-assembly` for full use)

These live on the Alignment server but belong to the assembly workflow; listed for completeness.

| Tool | Method ID | Purpose |
|---|---|---|
| [`find_overlap`](../../../../docs/mcp/tools/alignment/find_overlap.md) | `SequenceAssembler.FindOverlap` | Suffix-prefix overlap between two reads. |
| [`find_all_overlaps`](../../../../docs/mcp/tools/alignment/find_all_overlaps.md) | `SequenceAssembler.FindAllOverlaps` | All pairwise overlaps in a read set. |
| [`assemble_de_bruijn`](../../../../docs/mcp/tools/alignment/assemble_de_bruijn.md) | `SequenceAssembler.AssembleDeBruijn` | de Bruijn assembly. |
| [`assemble_olc`](../../../../docs/mcp/tools/alignment/assemble_olc.md) | `SequenceAssembler.AssembleOLC` | Overlap-layout-consensus assembly. |
| [`assembly_stats`](../../../../docs/mcp/tools/alignment/assembly_stats.md) | `SequenceAssembler.CalculateStats` | N50 and assembly metrics. |
| [`calculate_coverage`](../../../../docs/mcp/tools/alignment/calculate_coverage.md) | `SequenceAssembler.CalculateCoverage` | Per-position / mean coverage. |
| [`error_correct_reads`](../../../../docs/mcp/tools/alignment/error_correct_reads.md) | `SequenceAssembler.ErrorCorrectReads` | k-mer spectrum read correction. |
| [`quality_trim_reads`](../../../../docs/mcp/tools/alignment/quality_trim_reads.md) | `SequenceAssembler.QualityTrimReads` | Quality-based read trimming. |
| [`merge_contigs`](../../../../docs/mcp/tools/alignment/merge_contigs.md) | `SequenceAssembler.MergeContigs` | Merge overlapping contigs. |
| [`scaffold_contigs`](../../../../docs/mcp/tools/alignment/scaffold_contigs.md) | `SequenceAssembler.Scaffold` | Order/orient contigs into scaffolds. |

> For k-mer/assembly workflows use **`bio-assembly`**; this skill owns the alignment/similarity subset above.
