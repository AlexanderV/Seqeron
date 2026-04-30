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
    int EndPosition2);

/// <summary>Alignment summary statistics.</summary>
public record AlignmentStatisticsDto(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent);

/// <summary>BLAST-style formatted alignment text.</summary>
public record FormatAlignmentResult(string Formatted);

/// <summary>Multiple sequence alignment result.</summary>
public record MultipleAlignmentResultDto(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore);

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

/// <summary>Read with sequence and matching Phred+33 quality string.</summary>
public record QualityReadInput(string Sequence, string Quality);
