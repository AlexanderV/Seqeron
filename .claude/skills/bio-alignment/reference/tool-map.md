# bio-alignment tool map (86 tools)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters. Servers: **Alignment** (66), **Core** (18), **Analysis** (2).

## Pairwise alignment (Alignment server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`global_align`](../../../../docs/mcp/tools/alignment/global_align.md) | `SequenceAligner.GlobalAlign` | Needleman-Wunsch; align full length of both seqs end-to-end. |
| [`local_align`](../../../../docs/mcp/tools/alignment/local_align.md) | `SequenceAligner.LocalAlign` | Smith-Waterman; best-scoring shared substring, flanks ignored. |
| [`local_align_top_n`](../../../../docs/mcp/tools/alignment/local_align_top_n.md) | `SequenceAligner.LocalAlignTopN` | Waterman-Eggert top-N non-intersecting local alignments (repeats / multiple domains; EMBOSS matcher -alternatives, lalign); `maxAlignments`, `minScore`, `affine`, `matrix`. |
| [`local_align_co_optimal`](../../../../docs/mcp/tools/alignment/local_align_co_optimal.md) | `SequenceAligner.LocalAlignCoOptimal` | All co-optimal local alignments (Biopython `alignments` set), capped. |
| [`global_align_co_optimal`](../../../../docs/mcp/tools/alignment/global_align_co_optimal.md) | `SequenceAligner.GlobalAlignCoOptimal` | All co-optimal global / ends-free alignments (Biopython `align()` set), capped; `affine`, `matrix`, `endGaps`. |
| [`general_gap_align`](../../../../docs/mcp/tools/alignment/general_gap_align.md) | `SequenceAligner.GlobalAlignGeneralGap` | Arbitrary gap-cost function (Waterman–Smith–Beyer): concave / convex / capped / step gaps via a gapScores table; global, local or semi_global. |
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
| [`multiple_align`](../../../../docs/mcp/tools/alignment/multiple_align.md) | `SequenceAligner.MultipleAlign` | Center-star MSA (exact pairwise NW) → aligned rows, consensus, sum-of-pairs score; siblings `multiple_align_anchored`, `_progressive`, `_iterative`, `_consistency`, `_clustalw`. `type=protein` / `matrix=BLOSUM62|PAM250|EDNAFULL` (+ fractional gaps) on every `multiple_align*` tool except `_clustalw` scores residues with the substitution matrix (protein MSA, B09 F61). |
| [`multiple_align_exact`](../../../../docs/mcp/tools/alignment/multiple_align_exact.md) | `SequenceAligner.MultipleAlignExact` | EXACT sum-of-pairs optimal MSA (linear / natural affine gaps; k-dim DP + Carrillo-Lipman bound; <= 8 / 6 sequences, `maxStateCells` guard) - use when the heuristics' score must be proven optimal. |
| [`multiple_align_clustalw`](../../../../docs/mcp/tools/alignment/multiple_align_clustalw.md) | `SequenceAligner.MultipleAlignClustalW` | ClustalW 2.1 MSA of DNA or `type=protein` (user `-MATRIX` / `-PWMATRIX` files via `matrixFile` / `pwMatrixFile`), row-for-row equal to `clustalw -ALIGN` (NJ tree, tree weights, delayed divergent sequences, position-specific gaps; protein: matrix series `matrix=gonnet|blosum|pam|id`, residue-specific / hydrophilic / variable-region gap penalties, `gapDist`); `distance=clustalw-quick` = `-QUICKTREE`. |
| [`multiple_align_muscle`](../../../../docs/mcp/tools/alignment/multiple_align_muscle.md) | `SequenceAligner.MultipleAlignMuscle` | MUSCLE 3.8.1551 progressive MSA (k-mer distance, UPGMB tree, ClustalW weights, log-expectation / SPN profile DP; `maxIters` 1-2): port, row-for-row oracle (Docker hand-off) awaiting |
| [`multiple_align_tcoffee`](../../../../docs/mcp/tools/alignment/multiple_align_tcoffee.md) | `SequenceAligner.MultipleAlignTCoffeeProgram` | T-Coffee 13.41 program end to end (proba_pair / slow_pair / lalign_id_pair libraries, k-tuple NJ tree, relaxation + linked pair-wise DP): rows = `t_coffee` on 1 480 sets |
| [`multiple_align_kalign`](../../../../docs/mcp/tools/alignment/multiple_align_kalign.md) | `SequenceAligner.KalignRun` | Kalign 3.3.5 MSA (Myers bit-parallel distances, UPGMA / bisecting k-means tree, affine profile DP in linear space) of DNA / RNA or protein; rows equal the installed program on 526 / 526 runs. |
| [`kalign_command`](../../../../docs/mcp/tools/alignment/kalign_command.md) | `KalignProgram.Run` | Kalign 3.3.5 as a command line (`--type --gpo --gpe --tgpe -n -i -o -f`, getopt_long_only forms, FASTA files as text); equal to the program on 83 / 83 command lines. |
| [`multiple_align_mafft`](../../../../docs/mcp/tools/alignment/multiple_align_mafft.md) | `SequenceAligner.MultipleAlignMafft` | MAFFT 7.505 alignment of TWO sequences (stage 1 of the port: BLOSUM / Kimura scoring, `--op` / `--ep` / `--bl`, `--anysymbol`, free end gaps); rows equal the installed program on 884 / 884 pairs; T-Coffee `mafft_pair` uses it. |
| [`stretcher_align`](../../../../docs/mcp/tools/alignment/stretcher_align.md) | `SequenceAligner.StretcherAlign` | EMBOSS 6.6 stretcher: Myers-Miller linear-space global alignment, integer scores. |
| [`needle_align`](../../../../docs/mcp/tools/alignment/needle_align.md) | `SequenceAligner.NeedleAlign` | EMBOSS 6.6 needle: Needleman-Wunsch-Gotoh in the program's single-precision arithmetic, `-endweight` / `-endopen` / `-endextend`; rows and score equal the program on 887 / 887 pairs. |
| [`water_align`](../../../../docs/mcp/tools/alignment/water_align.md) | `SequenceAligner.WaterAlign` | EMBOSS 6.6 water: Smith-Waterman-Gotoh as the program reports it (ties, start coordinates, float score); equal to the program on 860 / 860 pairs. |
| [`write_alignment`](../../../../docs/mcp/tools/alignment/write_alignment.md) | `EmbossAlignmentWriter.Format` | Writes an alignment in an EMBOSS `-aformat` format (pair, markx, srspair, msf, clustal, mega, nexus, phylip, sam, ...). |
| [`dual_affine_align`](../../../../docs/mcp/tools/alignment/dual_affine_align.md) | `SequenceAligner.GlobalAlignDualAffine` / `LocalAlignDualAffine` | ksw2 / minimap2-style two-piece affine gaps, global or local, optional band. |
| [`supermatcher_align`](../../../../docs/mcp/tools/alignment/supermatcher_align.md) | `SequenceAligner.Supermatcher` | EMBOSS 6.6 supermatcher: Rabin-Karp word seeds + banded Smith-Waterman of large sequences. |
| [`needleall_align`](../../../../docs/mcp/tools/alignment/needleall_align.md) | `SequenceAligner.NeedleAll` | EMBOSS 6.6 needleall: every sequence of one set against every sequence of another (needle). |
| [`ksw2_extend`](../../../../docs/mcp/tools/alignment/ksw2_extend.md) | `Ksw2Aligner` | ksw2 / minimap2 extension alignment: affine or dual-affine gaps, band, Z-drop, end bonus, N penalty, CIGAR. |
| [`clustalw_pairwise_align`](../../../../docs/mcp/tools/alignment/clustalw_pairwise_align.md) | `SequenceAligner.ClustalWPairwiseAlign` | ClustalW 2.1 full pairwise stage (DNA / protein): local alignment, percent identity, guide-tree distance. |
| [`clustalw_quick_pairwise_align`](../../../../docs/mcp/tools/alignment/clustalw_quick_pairwise_align.md) | `SequenceAligner.ClustalWQuickPairwiseAlign` | ClustalW 2.1 fast k-tuple pairwise score / distance (-QUICKTREE). |
| [`clustalw_command`](../../../../docs/mcp/tools/alignment/clustalw_command.md) | `ClustalWProgram.Run` | ClustalW 2.1 as a command line: alignment, profiles, guide-tree files, -OUTORDER, seven output formats, structure masks. |
| [`clustalw_guide_tree`](../../../../docs/mcp/tools/alignment/clustalw_guide_tree.md) | `SequenceAligner.ClustalWGuideTreeText` | ClustalW 2.1 neighbour-joining guide-tree file (`.dnd`) from pairwise distances. |
| [`clustalw_profile_align`](../../../../docs/mcp/tools/alignment/clustalw_profile_align.md) | `SequenceAligner.AlignClustalWProfiles` | ClustalW 2.1 profile-profile alignment of two existing alignments (-PROFILE). |
| [`clustalw_sequences_to_profile`](../../../../docs/mcp/tools/alignment/clustalw_sequences_to_profile.md) | `SequenceAligner.AlignClustalWSequencesToProfile` | ClustalW 2.1 alignment of sequences to an existing profile (-SEQUENCES). |
| [`tcoffee_slow_pair_align`](../../../../docs/mcp/tools/alignment/tcoffee_slow_pair_align.md) | `SequenceAligner.TCoffeeSlowPairAlign` | T-Coffee slow_pair global pair alignment. |
| [`tcoffee_lalign_top_alignments`](../../../../docs/mcp/tools/alignment/tcoffee_lalign_top_alignments.md) | `SequenceAligner.TCoffeeLalignTopAlignments` | T-Coffee lalign_id_pair top-N non-intersecting local alignments. |
| [`tcoffee_primary_library`](../../../../docs/mcp/tools/alignment/tcoffee_primary_library.md) | `SequenceAligner.TCoffeePrimaryLibraryEntries` | T-Coffee primary library (weighted residue pairs). |
| [`tcoffee_blast_msa_library`](../../../../docs/mcp/tools/alignment/tcoffee_blast_msa_library.md) | `SequenceAligner.TCoffeeBlastMsaLibraryEntries` | T-Coffee `blastp_msa` library from the XML of an all-against-all BLAST search. |
| [`tcoffee_tblastx_library`](../../../../docs/mcp/tools/alignment/tcoffee_tblastx_library.md) | `SequenceAligner.TCoffeeTblastxLibraryEntries` | T-Coffee `tblastx_msa` library from the XML of an all-against-all tblastx search. |
| [`tcoffee_guide_tree`](../../../../docs/mcp/tools/alignment/tcoffee_guide_tree.md) | `SequenceAligner.TCoffeeGuideTree` | T-Coffee guide tree (Newick). |
| [`tcoffee_ktup_similarity`](../../../../docs/mcp/tools/alignment/tcoffee_ktup_similarity.md) | `SequenceAligner.TCoffeeKtupSimilarity` | T-Coffee k-tuple similarity matrix. |
| [`sum_of_pairs_score`](../../../../docs/mcp/tools/alignment/sum_of_pairs_score.md) | `SequenceAligner.SumOfPairsScore` | Sum-of-pairs score of an alignment (linear / affine, matrix). |
| [`multiple_alignment_statistics`](../../../../docs/mcp/tools/analysis/multiple_alignment_statistics.md) | `MultipleAlignmentStatistics.Calculate` | Analysis server: per-sequence infoalign stats vs reference/consensus, per-column counts/gaps/ClustalW conservation, pairwise identity matrix |
| [`showalign`](../../../../docs/mcp/tools/analysis/showalign.md) | `MultipleAlignmentDisplay.ShowAlign` | Analysis server: EMBOSS showalign pretty display of an alignment against the consensus or a reference row (-show / -order / -similarcase / numbers / ruler) |
| [`alignment_distance_matrix`](../../../../docs/mcp/tools/analysis/alignment_distance_matrix.md) | `MultipleAlignmentStatistics.CalculateDistanceMatrix` | Analysis server: EMBOSS distmat uncorrected / Jukes-Cantor distance matrix of an MSA |
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
| [`multiple_align_probcons`](../../../../docs/mcp/tools/alignment/multiple_align_probcons.md) | `SequenceAligner.ProbConsRun` | ProbCons 1.12 MSA (pair-HMM posteriors, consistency transformation, expected-accuracy tree, maximum-expected-accuracy profile alignment, iterative refinement with glibc `rand`); `consistency` / `iterativeRefinement` / `preTraining` / `alignmentOrder` / `cutoff` / `parameters` / `annotate`; equal to the installed program on 407 / 407 command lines. |
| [`probcons_command`](../../../../docs/mcp/tools/alignment/probcons_command.md) | `ProbConsProgram.Run` | ProbCons 1.12 as a command line (`-clustalw -c -ir -pre -pairs -viterbi -v -annot -t -e -p -a -go -ge -co`, MFA / MSF / Clustal files as text, stdout / stderr / written files byte for byte); equal to the program on 407 / 407 command lines. |
| [`multiple_align_poa`](../../../../docs/mcp/tools/alignment/multiple_align_poa.md) | `SequenceAligner.PoaRun` | POA 2.0 partial-order MSA (graph DP with the program's affine gap model, fusion into align rings, local / `global` / `progressive` / `pairScores`, `fuseAll`, heaviest bundle consensus); returns the rows and the program's PO / PIR / CLUSTAL texts; equal to the installed program on 359 / 359 command lines. |
| [`poa_command`](../../../../docs/mcp/tools/alignment/poa_command.md) | `PoaProgram.Run` | POA 2.0 as a command line (`-read_fasta -do_global -do_progressive -read_pairscores -fuse_all -hb -hbmin -best -preserve_seqorder -tolower -toupper -pir -clustal -po -printmatrix -v MATRIX`, files as text); equal to the program on 359 / 359 command lines. |
| [`align_gap_profile`](../../../../docs/mcp/tools/alignment/align_gap_profile.md) | `SequenceAligner.GlobalAlignAffine(GapScoreProfile)` | Biopython per-direction / per-end gap scores (left / right / internal x insertion / deletion, open / extend) |
| [`merger`](../../../../docs/mcp/tools/alignment/merger.md) | `SequenceAligner.Merger` | EMBOSS merger (merge two overlapping sequences) |
| [`minimap2_chains`](../../../../docs/mcp/tools/alignment/minimap2_chains.md) | `Minimap2Chaining.MapChains` | minimap2 2.26 seeding and chaining (minimizer anchors, mg_lchain_dp) |
| [`dotmatcher`](../../../../docs/mcp/tools/alignment/dotmatcher.md) | `DotMatcher.Compute` | EMBOSS dotmatcher windowed dot-plot segments of two sequences |

> For k-mer/assembly workflows use **`bio-assembly`**; this skill owns the alignment/similarity subset above.
