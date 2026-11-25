// SPDX-License-Identifier: AGPL-3.0-only
// This file is part of Seqeron.Mcp.OverrepresentedWords, licensed under the terms stated in LICENSE (this directory).
namespace Seqeron.Mcp.OverrepresentedWords.Tools;

/// <summary>A k-mer (or reverse-complement pair) scored by RSAT oligo-analysis occurrence statistics.</summary>
public record OligoMotifItem(
    string Sequence,
    string? ReverseComplement,
    int Count,
    int[] Positions,
    double ExpectedFrequency,
    double ExpectedOccurrences,
    double? Ratio,
    double OccurrenceProbability,
    double OccurrenceEValue,
    double OccurrenceSignificance)
{
    /// <summary>Sequence index of each position (0 = <c>sequence</c>); RSAT-option runs only.</summary>
    public int[]? SequenceIndices { get; init; }

    /// <summary>Discarded overlapping windows (<c>-noov</c>, RSAT <c>ovl_occ</c>); RSAT-option runs only.</summary>
    public int? Overlaps { get; init; }

    /// <summary>occ / total occurrences (RSAT <c>obs_freq</c>); RSAT-option runs only.</summary>
    public double? ObservedFrequency { get; init; }

    /// <summary>RSAT <c>exp_var</c> (with zscore or a calibration).</summary>
    public double? ExpectedVariance { get; init; }

    /// <summary>RSAT <c>ovlp</c> overlap coefficient (with zscore).</summary>
    public double? OverlapCoefficient { get; init; }

    /// <summary>RSAT z-score (occ − exp_occ)/√exp_var (with zscore; null when the variance is ≤ 0).</summary>
    public double? ZScore { get; init; }

    /// <summary>Distribution behind occ_P: Binomial, Poisson or NegativeBinomial; RSAT-option runs only.</summary>
    public string? FittedDistribution { get; init; }

    /// <summary>Best lexicon segmentation "prefix|suffix" (background 'lexicon').</summary>
    public string? LexiconSegmentation { get; init; }
}

/// <summary>Result of <c>oligo_analysis</c>.</summary>
public record OligoAnalysisResultDto(
    OligoMotifItem[] Motifs,
    int OligoLength,
    string Strands,
    bool CountOverlapping,
    long TotalOccurrences,
    int TestedPatterns,
    double? PossibleOligos)
{
    /// <summary>Number of analysed sequences; RSAT-option runs only.</summary>
    public int? SequenceCount { get; init; }

    /// <summary>Degenerate-word mode ('none', 'oneN', 'onedeg'); RSAT-option runs only.</summary>
    public string? Degenerate { get; init; }

    /// <summary>Sequence type ('protein' or 'other', RSAT -seqtype); null for DNA.</summary>
    public string? SequenceType { get; init; }

    /// <summary>RSAT alphabet_size (20 for protein, number of distinct residues for other); null for DNA.</summary>
    public int? AlphabetSize { get; init; }
}

/// <summary>A spaced dyad (or reverse-complement pair) scored by RSAT dyad-analysis.</summary>
public record DyadItem(
    string Pattern,
    string FirstMonad,
    int Spacing,
    string SecondMonad,
    string? ReverseComplement,
    int Occurrences,
    int Overlaps,
    int[] SequenceIndices,
    int[] Positions,
    double ObservedFrequency,
    double? ExpectedFrequency,
    double? ExpectedOccurrences,
    double? ExpectedVariance,
    double OverlapCoefficient,
    double? ZScore,
    double Ratio,
    double? OccurrenceProbability,
    double? OccurrenceEValue,
    double? OccurrenceSignificance,
    bool IsDirectRepeat,
    bool IsReversePalindrome,
    bool ExpectedFromMonads);

/// <summary>Per-spacing positions and counts of <c>dyad_analysis</c>.</summary>
public record DyadSpacingItem(int Spacing, long PossiblePositions, long Occurrences, long Overlaps);

/// <summary>Result of <c>dyad_analysis</c>.</summary>
public record DyadAnalysisResultDto(
    DyadItem[] Dyads,
    int MonadLength,
    int MinSpacing,
    int MaxSpacing,
    string DyadType,
    string Strands,
    bool CountOverlapping,
    int SequenceCount,
    long MonadOccurrences,
    DyadSpacingItem[] Spacings,
    int TestedPatterns,
    double? PossibleDyads);

/// <summary>A k-mer (or pair) scored by RSAT oligo-analysis matching-sequence statistics.</summary>
public record SharedMotifSignificanceItem(
    string Sequence,
    string? ReverseComplement,
    int[] SequenceIndices,
    double Prevalence,
    double ExpectedFrequency,
    double ExpectedMatchingSequences,
    double MatchingSequenceProbability,
    double? MatchingSequenceEValue,
    double? MatchingSequenceSignificance);

/// <summary>Result of <c>shared_motifs_significance</c>.</summary>
public record SharedMotifSignificanceResult(
    SharedMotifSignificanceItem[] Motifs,
    int OligoLength,
    string Strands,
    int SequenceCount,
    long PossiblePositions,
    double? PossibleOligos)
{
    /// <summary>Degenerate-word mode ('oneN' or 'onedeg', RSAT -oneN / -onedeg); null for plain words.</summary>
    public string? Degenerate { get; init; }
}
