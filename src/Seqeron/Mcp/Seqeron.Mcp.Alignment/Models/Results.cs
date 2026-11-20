namespace Seqeron.Mcp.Alignment.Tools;

// ================================
// SequenceAligner DTOs
// ================================

/// <summary>Pairwise alignment result.</summary>
public record AlignmentResultDto(
    string AlignedSequence1,
    string AlignedSequence2,
    int Score,
    string AlignmentType,
    int StartPosition1,
    int StartPosition2,
    int EndPosition1,
    int EndPosition2,
    double ExactScore = 0,
    string? Formatted = null);

/// <summary>List of alignments (local Waterman–Eggert top-N, local / global co-optimal enumeration), in the order found.</summary>
public record LocalAlignmentListResult(AlignmentResultDto[] Alignments, int Count);

/// <summary>Alignment summary statistics.</summary>
public record AlignmentStatisticsDto(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent,
    int Positives = 0,
    int GapOpens = 0,
    int BlastPositives = 0,
    int BlastIdentityPercent = 0,
    int BlastPositivesPercent = 0,
    int BlastGapPercent = 0,
    string BlastPident = "",
    string BlastReport = "");

/// <summary>BLAST-style formatted alignment text.</summary>
public record FormatAlignmentResult(string Formatted);

/// <summary>An alignment file: its text, the format written and the dialect (tool whose layout it reproduces).</summary>
public record WriteAlignmentResult(string Text, string Format, string Dialect);

/// <summary>ClustalW 2.1 command: the files the program wrote (name to text) and what it printed to the error stream while reading the command line.</summary>
public record ClustalWCommandResultDto(
    Dictionary<string, string> Files,
    string[] Messages);

/// <summary>Kalign 3.3.5 alignment: the rows of the sequences that were aligned, which input sequences those are and the alphabet the program detected.</summary>
public record KalignResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    double ExactTotalScore,
    int[] KeptIndices,
    string Alphabet,
    int[][] GuideTree);

/// <summary>Kalign 3.3.5 command: exit code, standard output, the files the program wrote (name to text) and what it printed to the error stream.</summary>
public record KalignCommandResultDto(
    int ExitCode,
    string Stdout,
    Dictionary<string, string> Files,
    string[] Messages);

/// <summary>ProbCons 1.12 alignment: rows in output order, the input index of each row, the guide tree the program prints and the optional column annotation.</summary>
public record ProbConsResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    double ExactTotalScore,
    int[] Order,
    string GuideTree,
    int[]? Annotation);

/// <summary>ProbCons 1.12 command: exit code, standard output, standard error, and the files the program wrote (name to text).</summary>
public record ProbConsCommandResultDto(
    int ExitCode,
    string Stdout,
    string Stderr,
    Dictionary<string, string> Files);

/// <summary>Multiple sequence alignment result.</summary>
public record MultipleAlignmentResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    double ExactTotalScore = 0,
    string? Formatted = null);

/// <summary>MUSCLE 3.8 run: the alignment and the files the program can write (-tree1, -tree2, -distmx1, -distmx2, -computeweights).</summary>
public record MuscleRunResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    double ExactTotalScore,
    string[] Names,
    int[] OutputOrder,
    string? Tree1,
    string? Tree2,
    string? DistanceMatrix1,
    string? DistanceMatrix2,
    string? WeightsFile,
    bool ReferenceBehaviourUndefined = false);

/// <summary>Result of a MUSCLE 3.8 command line on alignments: rows (profile, refine, refinew, sw), a score (spscore, ppscore) or a tree (maketree).</summary>
public record MuscleCommandResultDto(
    string[] AlignedSequences,
    double? Score,
    string? Tree);

/// <summary>Exact sum-of-pairs multiple alignment with the size of the problem solved (heuristic seed bound, lattice and visited cells).</summary>
public record ExactMultipleAlignmentResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore,
    int LowerBound,
    long LatticeCells,
    long VisitedCells,
    long StateValues,
    double ExactTotalScore = 0,
    int ScoreScale = 1);

// ================================
// ApproximateMatcher DTOs
// ================================

/// <summary>Single approximate match record.</summary>
public record ApproximateMatchDto(
    int Position,
    string MatchedSequence,
    int Distance,
    int[] MismatchPositions,
    string MismatchType);

/// <summary>Container for approximate match collections.</summary>
public record ApproximateMatchListResult(ApproximateMatchDto[] Items);

/// <summary>Best (minimum distance) approximate match, or null.</summary>
public record FindBestMatchResult(ApproximateMatchDto? Match);

/// <summary>One Sellers end position with the minimum edit distance of a substring ending there.</summary>
public record EditEndPositionItem(int EndPosition, int Distance);

/// <summary>Container for Sellers end positions (increasing end position).</summary>
public record EditEndPositionsResult(EditEndPositionItem[] Items);

/// <summary>Unit-cost (Levenshtein) alignment path in edlib's convention.</summary>
public record EditAlignmentDto(
    int Distance,
    string Operations,
    string Cigar,
    string StandardCigar,
    string AlignedQuery,
    string AlignedTarget,
    int[] SubstitutionPositions,
    bool HasIndels);

/// <summary>Damerau–Levenshtein (unrestricted) or optimal string alignment distance.</summary>
public record DamerauDistanceResult(int Distance, string Variant);

/// <summary>One step of a Damerau edit script (kind: match, substitution, insertion, deletion, transposition).</summary>
public record DamerauOperationDto(string Kind, int SourcePosition, int SourceLength, int TargetPosition, int TargetLength, int Cost);

/// <summary>Transposition-aware (Lowrance–Wagner) edit script from sequence1 to sequence2.</summary>
public record DamerauAlignmentDto(int Distance, string Variant, string Script, int TranspositionCount, DamerauOperationDto[] Operations);

/// <summary>Frequent k-mer item.</summary>
public record FrequentKmerItem(string Kmer, int Count);

/// <summary>Container for frequent k-mer results.</summary>
public record FrequentKmersResult(FrequentKmerItem[] Items);

// ================================
// SequenceAssembler DTOs
// ================================

/// <summary>Genome assembly result with N50 and totals.</summary>
public record AssemblyResultDto(
    string[] Contigs,
    int TotalReads,
    int AssembledReads,
    double N50,
    int LongestContig,
    int TotalLength);

/// <summary>Pairwise read overlap.</summary>
public record OverlapDto(
    int ReadIndex1,
    int ReadIndex2,
    int OverlapLength,
    int Position1,
    int Position2);

/// <summary>Container for overlap collections.</summary>
public record OverlapListResult(OverlapDto[] Items);

/// <summary>Suffix-prefix overlap descriptor.</summary>
public record FindOverlapInfo(int Length, int Position1, int Position2);

/// <summary>find_overlap result wrapper (null if no qualifying overlap).</summary>
public record FindOverlapResult(FindOverlapInfo? Overlap);

/// <summary>Sequence identity ratio in [0,1].</summary>
public record IdentityResult(double Identity);

/// <summary>Merged contig result.</summary>
public record MergedContigResult(string Merged);

/// <summary>Scaffold strings result.</summary>
public record ScaffoldsResult(string[] Scaffolds);

/// <summary>Per-base coverage depth array.</summary>
public record CoverageResult(int[] Coverage);

/// <summary>Consensus sequence result.</summary>
public record ConsensusResult(string Consensus);

/// <summary>Quality-trimmed reads result.</summary>
public record TrimmedReadsResult(string[] Trimmed);

/// <summary>Error-corrected reads result.</summary>
public record CorrectedReadsResult(string[] Corrected);

// ================================
// Input DTOs
// ================================

/// <summary>Paired-end link between two contigs with estimated gap size.</summary>
public record ContigLinkInput(int Contig1, int Contig2, int GapSize);

/// <summary>What an external structure aligner wrote for one pair of structures (input of <c>multiple_align_tcoffee</c>): the method (<c>TMalign_pair</c> or <c>mustang_pair</c>), the input indices of the two sequences in the order the program calls the aligner, and the output text (TMalign: its standard output, only the last four lines matter; mustang: the content of results.afasta).</summary>
public record StructureAlignerOutputInput(string Method, int First, int Second, string Output);

/// <summary>Read with sequence and matching Phred+33 quality string.</summary>
public record QualityReadInput(string Sequence, string Quality);


// ================================
// Ported program stages (ClustalW 2.1 pairwise, T-Coffee 13.41 library / tree, sum of pairs)
// ================================

/// <summary>ClustalW 2.1 pairwise stage: the local alignment, ClustalW's percent identity and guide-tree distance.</summary>
public record ClustalWPairwiseAlignResult(
    int Score, int Identities, double PercentIdentity, double Distance, int ClustalWScore,
    string AlignedSequence1, string AlignedSequence2, int Start1, int End1, int Start2, int End2);

/// <summary>ClustalW 2.1 fast (k-tuple, -QUICKTREE) pairwise stage.</summary>
public record ClustalWQuickPairwiseAlignResult(int FragmentScore, double Score, double Distance, bool Aligned, bool PartialAlignment);

/// <summary>A pair of aligned rows.</summary>
public record AlignedPairResult(string RowA, string RowB);

/// <summary>One local alignment of the T-Coffee lalign stage (rows, 0-based starts, score).</summary>
public record TCoffeeLocalAlignmentItem(string RowA, string RowB, int StartA, int StartB, int Score);

/// <summary>Result of <c>tcoffee_lalign_top_alignments</c>.</summary>
public record TCoffeeLocalAlignmentsResult(TCoffeeLocalAlignmentItem[] Alignments, int Count);

/// <summary>One primary-library entry: residue PosA of sequence SeqA with residue PosB of sequence SeqB (0-based), weight 10 x percent identity, times added.</summary>
public record TCoffeeLibraryEntryItem(int SeqA, int PosA, int SeqB, int PosB, int Weight, int Count);

/// <summary>Result of <c>tcoffee_primary_library</c>.</summary>
public record TCoffeePrimaryLibraryResult(TCoffeeLibraryEntryItem[] Entries, int Count);

/// <summary>Guide tree of T-Coffee: Newick text (leaf labels = names or s0, s1, ...) and the leaf indices left to right (the merge order is the post-order of the tree).</summary>
public record TCoffeeGuideTreeResult(string Newick, int[] LeafOrder);

/// <summary>Pairwise k-tuple similarity matrix of T-Coffee's guide-tree stage (percent x 100 as the program stores it).</summary>
public record TCoffeeKtupSimilarityResult(int[][] Similarity);

/// <summary>Sum-of-pairs score of a multiple alignment.</summary>
public record SumOfPairsScoreResult(double Score, string GapModel, string Scoring);

/// <summary>One pairwise alignment reported by supermatcher / needleall: the rows, the residues of each sequence before (Offset) and after (Offend) the aligned part, and the single-precision score.</summary>
public record PairAlignmentItem(string Name1, string Name2, string AlignedSequence1, string AlignedSequence2, int Offset1, int Offset2, int Offend1, int Offend2, double Score);

/// <summary>Result of <c>supermatcher_align</c> / <c>needleall_align</c>: the alignments in the order the program writes them and the error-file text.</summary>
public record EmbossProgramAlignResult(PairAlignmentItem[] Alignments, int Count, string ErrorText, string? Formatted = null);

/// <summary>Result of <c>ksw2_extend</c>: the fields of ksw_extz_t and the CIGAR.</summary>
public record Ksw2ExtendResult(int Max, bool ZDropped, int MaxQ, int MaxT, int Mqe, int MqeT, int Mte, int MteQ, int Score, bool ReachEnd, string Cigar);
