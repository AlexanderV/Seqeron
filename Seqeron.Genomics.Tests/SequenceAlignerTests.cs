using System;
using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceAlignerTests
{
    #region Semi-Global Alignment Tests


    [Test]
    public void SemiGlobalAlign_ShortInLong_FindsMatch()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("AAAATGCAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference);

        Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
    }

    [Test]
    public void SemiGlobalAlign_FreeEndGaps_NoGapPenalty()
    {
        var query = new DnaSequence("ATGC");
        var reference = new DnaSequence("ATGCAAAA");

        var result = SequenceAligner.SemiGlobalAlign(query, reference);

        // Semi-global should not heavily penalize trailing gaps
        Assert.That(result.Score, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.AlignedSequence1.Replace("-", ""), Is.EqualTo("ATGC"));
    }

    #endregion

    #region Scoring Matrix Tests

    [Test]
    public void SimpleDna_HasExpectedValues()
    {
        Assert.That(SequenceAligner.SimpleDna.Match, Is.EqualTo(1));
        Assert.That(SequenceAligner.SimpleDna.Mismatch, Is.EqualTo(-1));
        Assert.That(SequenceAligner.SimpleDna.GapOpen, Is.EqualTo(-2));
        Assert.That(SequenceAligner.SimpleDna.GapExtend, Is.EqualTo(-1));
    }

    [Test]
    public void BlastDna_HasExpectedValues()
    {
        Assert.That(SequenceAligner.BlastDna.Match, Is.EqualTo(2));
        Assert.That(SequenceAligner.BlastDna.Mismatch, Is.EqualTo(-3));
    }

    [Test]
    public void HighIdentityDna_HasExpectedValues()
    {
        Assert.That(SequenceAligner.HighIdentityDna.Match, Is.EqualTo(5));
        Assert.That(SequenceAligner.HighIdentityDna.GapOpen, Is.EqualTo(-10));
    }

    #endregion

    #region Alignment Statistics Tests

    [Test]
    public void CalculateStatistics_PerfectMatch_100Identity()
    {
        var result = new AlignmentResult(
            AlignedSequence1: "ATGC",
            AlignedSequence2: "ATGC",
            Score: 4,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 3, EndPosition2: 3);

        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.That(stats.Identity, Is.EqualTo(100));
        Assert.That(stats.Matches, Is.EqualTo(4));
        Assert.That(stats.Mismatches, Is.EqualTo(0));
        Assert.That(stats.Gaps, Is.EqualTo(0));
    }

    [Test]
    public void CalculateStatistics_WithMismatches_CountsCorrectly()
    {
        var result = new AlignmentResult(
            AlignedSequence1: "ATGC",
            AlignedSequence2: "ATTC",
            Score: 2,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 3, EndPosition2: 3);

        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.That(stats.Matches, Is.EqualTo(3));
        Assert.That(stats.Mismatches, Is.EqualTo(1));
        Assert.That(stats.Identity, Is.EqualTo(75));
    }

    [Test]
    public void CalculateStatistics_WithGaps_CountsGaps()
    {
        var result = new AlignmentResult(
            AlignedSequence1: "AT-GC",
            AlignedSequence2: "ATXGC",
            Score: 2,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 3, EndPosition2: 4);

        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.That(stats.Gaps, Is.EqualTo(1));
        Assert.That(stats.AlignmentLength, Is.EqualTo(5));
    }

    [Test]
    public void CalculateStatistics_EmptyAlignment_ReturnsEmpty()
    {
        var stats = SequenceAligner.CalculateStatistics(AlignmentResult.Empty);

        Assert.That(stats.AlignmentLength, Is.EqualTo(0));
        Assert.That(stats.Identity, Is.EqualTo(0));
    }

    #endregion

    #region Format Alignment Tests

    [Test]
    public void FormatAlignment_CreatesVisualOutput()
    {
        var result = new AlignmentResult(
            AlignedSequence1: "ATGC",
            AlignedSequence2: "ATGC",
            Score: 4,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 3, EndPosition2: 3);

        string formatted = SequenceAligner.FormatAlignment(result);

        Assert.That(formatted, Does.Contain("ATGC"));
        Assert.That(formatted, Does.Contain("||||")); // Match indicators
    }

    [Test]
    public void FormatAlignment_ShowsMismatches()
    {
        var result = new AlignmentResult(
            AlignedSequence1: "ATGC",
            AlignedSequence2: "ATTC",
            Score: 2,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 3, EndPosition2: 3);

        string formatted = SequenceAligner.FormatAlignment(result);

        Assert.That(formatted, Does.Contain(".")); // Mismatch indicator
    }

    [Test]
    public void FormatAlignment_EmptyAlignment_ReturnsEmpty()
    {
        string formatted = SequenceAligner.FormatAlignment(AlignmentResult.Empty);

        Assert.That(formatted, Is.Empty);
    }

    [Test]
    public void FormatAlignment_LongSequence_WrapsLines()
    {
        var aligned = new string('A', 100);
        var result = new AlignmentResult(
            AlignedSequence1: aligned,
            AlignedSequence2: aligned,
            Score: 100,
            AlignmentType: AlignmentType.Global,
            StartPosition1: 0, StartPosition2: 0,
            EndPosition1: 99, EndPosition2: 99);

        string formatted = SequenceAligner.FormatAlignment(result, lineWidth: 50);

        Assert.That(formatted.Split('\n').Length, Is.GreaterThan(4)); // Multiple blocks
    }

    #endregion

    #region Multiple Alignment Tests

    [Test]
    public void MultipleAlign_TwoSequences_Aligns()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.AlignedSequences.Length, Is.EqualTo(2));
        Assert.That(result.Consensus, Is.EqualTo("ATGC"));
    }

    [Test]
    public void MultipleAlign_ThreeSequences_CreatesConsensus()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.AlignedSequences.Length, Is.EqualTo(3));
        Assert.That(result.Consensus, Is.Not.Empty);
    }

    [Test]
    public void MultipleAlign_DifferentLengths_PadsWithGaps()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.AlignedSequences[0].Length, Is.EqualTo(result.AlignedSequences[1].Length));
    }

    [Test]
    public void MultipleAlign_SingleSequence_ReturnsSame()
    {
        var sequences = new[] { new DnaSequence("ATGC") };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.AlignedSequences.Length, Is.EqualTo(1));
        Assert.That(result.Consensus, Is.EqualTo("ATGC"));
    }

    [Test]
    public void MultipleAlign_Empty_ReturnsEmpty()
    {
        var result = SequenceAligner.MultipleAlign(Array.Empty<DnaSequence>());

        Assert.That(result, Is.EqualTo(MultipleAlignmentResult.Empty));
    }

    [Test]
    public void MultipleAlign_ReturnsTotalScore()
    {
        var sequences = new[]
        {
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC"),
            new DnaSequence("ATGC")
        };

        var result = SequenceAligner.MultipleAlign(sequences);

        Assert.That(result.TotalScore, Is.GreaterThan(0));
    }

    #endregion

    #region Edge Cases

    [Test]
    public void SemiGlobalAlign_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.SemiGlobalAlign(null!, new DnaSequence("ATGC")));
    }

    [Test]
    public void CalculateStatistics_NullAlignment_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.CalculateStatistics(null!));
    }

    [Test]
    public void FormatAlignment_NullAlignment_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.FormatAlignment(null!));
    }

    [Test]
    public void MultipleAlign_NullSequences_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.MultipleAlign(null!));
    }

    #endregion
}
