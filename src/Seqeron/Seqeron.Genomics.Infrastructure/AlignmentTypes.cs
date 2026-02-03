using System;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Scoring parameters for sequence alignment.
/// </summary>
public sealed record ScoringMatrix(
    int Match,
    int Mismatch,
    int GapOpen,
    int GapExtend);

/// <summary>
/// Type of alignment performed.
/// </summary>
public enum AlignmentType
{
    Global,
    Local,
    SemiGlobal
}

/// <summary>
/// Result of pairwise sequence alignment.
/// </summary>
public sealed record AlignmentResult(
    string AlignedSequence1,
    string AlignedSequence2,
    int Score,
    AlignmentType AlignmentType,
    int StartPosition1,
    int StartPosition2,
    int EndPosition1,
    int EndPosition2)
{
    public static AlignmentResult Empty => new("", "", 0, AlignmentType.Global, 0, 0, 0, 0);
}

/// <summary>
/// Statistics calculated from an alignment.
/// </summary>
public readonly record struct AlignmentStatistics(
    int Matches,
    int Mismatches,
    int Gaps,
    int AlignmentLength,
    double Identity,
    double Similarity,
    double GapPercent)
{
    public static AlignmentStatistics Empty => new(0, 0, 0, 0, 0, 0, 0);
}

/// <summary>
/// Result of multiple sequence alignment.
/// </summary>
public sealed record MultipleAlignmentResult(
    string[] AlignedSequences,
    string Consensus,
    int TotalScore)
{
    public static MultipleAlignmentResult Empty => new(Array.Empty<string>(), "", 0);
}
