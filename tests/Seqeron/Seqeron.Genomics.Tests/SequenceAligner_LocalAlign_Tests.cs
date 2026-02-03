using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for ALIGN-LOCAL-001: Local Alignment (Smith–Waterman)
/// Evidence: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
/// </summary>
[TestFixture]
public class SequenceAligner_LocalAlign_Tests
{
    // Wikipedia example scoring: match +3, mismatch -3, gap -2
    private static readonly ScoringMatrix WikipediaScoring = new(
        Match: 3,
        Mismatch: -3,
        GapOpen: -2,
        GapExtend: -2);

    #region M1: Wikipedia Example - Evidence-based invariant validation

    /// <summary>
    /// M1: Wikipedia example with sequences TGTTACGG and GGTTGACTA.
    /// Evidence: Smith–Waterman algorithm example shows optimal local alignment.
    /// Validates INV-1 through INV-5.
    /// </summary>
    [Test]
    public void LocalAlign_WikipediaExample_ValidatesAllInvariants()
    {
        var seq1 = new DnaSequence("TGTTACGG");
        var seq2 = new DnaSequence("GGTTGACTA");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            // INV-1: Score ≥ 0 (zero floor property)
            Assert.That(result.Score, Is.GreaterThanOrEqualTo(0),
                "INV-1: Local alignment score must be non-negative (zero floor)");

            // INV-2: Aligned region is contiguous in both sequences (verified via positions)
            Assert.That(result.EndPosition1, Is.GreaterThanOrEqualTo(result.StartPosition1),
                "INV-2: EndPosition must be >= StartPosition in seq1");
            Assert.That(result.EndPosition2, Is.GreaterThanOrEqualTo(result.StartPosition2),
                "INV-2: EndPosition must be >= StartPosition in seq2");

            // INV-3: AlignmentType is Local
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Local),
                "INV-3: AlignmentType must be Local");

            // INV-4: Removing gaps from aligned sequences yields substrings of originals
            var gapsRemoved1 = RemoveGaps(result.AlignedSequence1);
            var gapsRemoved2 = RemoveGaps(result.AlignedSequence2);
            Assert.That(seq1.Sequence, Does.Contain(gapsRemoved1),
                "INV-4: Aligned seq1 with gaps removed must be substring of original");
            Assert.That(seq2.Sequence, Does.Contain(gapsRemoved2),
                "INV-4: Aligned seq2 with gaps removed must be substring of original");

            // INV-5: Positions are within sequence bounds
            Assert.That(result.StartPosition1, Is.GreaterThanOrEqualTo(0),
                "INV-5: StartPosition1 must be >= 0");
            Assert.That(result.EndPosition1, Is.LessThanOrEqualTo(seq1.Length),
                "INV-5: EndPosition1 must be <= seq1.Length");
            Assert.That(result.StartPosition2, Is.GreaterThanOrEqualTo(0),
                "INV-5: StartPosition2 must be >= 0");
            Assert.That(result.EndPosition2, Is.LessThanOrEqualTo(seq2.Length),
                "INV-5: EndPosition2 must be <= seq2.Length");
        });
    }

    #endregion

    #region M2: Score Non-Negative (Zero Floor)

    /// <summary>
    /// M2: Score is always non-negative due to zero floor property.
    /// Evidence: Smith–Waterman algorithm sets negative scores to 0.
    /// </summary>
    [Test]
    public void LocalAlign_AnyInput_ScoreIsNonNegative()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("TGCATGCA");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.That(result.Score, Is.GreaterThanOrEqualTo(0),
            "Zero floor: Local alignment score must be non-negative");
    }

    #endregion

    #region M3: AlignmentType is Local

    /// <summary>
    /// M3: AlignmentType property is set to Local for local alignments.
    /// Evidence: Implementation contract.
    /// </summary>
    [Test]
    public void LocalAlign_ValidInput_AlignmentTypeIsLocal()
    {
        var seq1 = new DnaSequence("ACGT");
        var seq2 = new DnaSequence("ACGT");

        var result = SequenceAligner.LocalAlign(seq1, seq2);

        Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Local));
    }

    #endregion

    #region M4: Position Validity

    /// <summary>
    /// M4: Start and End positions are valid within sequence bounds.
    /// Evidence: Smith–Waterman traceback from matrix positions.
    /// </summary>
    [Test]
    public void LocalAlign_ValidInput_PositionsAreWithinBounds()
    {
        var seq1 = new DnaSequence("AAATGCAAA");
        var seq2 = new DnaSequence("TGCTGC");

        var result = SequenceAligner.LocalAlign(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(result.StartPosition1, Is.GreaterThanOrEqualTo(0));
            Assert.That(result.EndPosition1, Is.GreaterThanOrEqualTo(result.StartPosition1));
            Assert.That(result.EndPosition1, Is.LessThanOrEqualTo(seq1.Length));

            Assert.That(result.StartPosition2, Is.GreaterThanOrEqualTo(0));
            Assert.That(result.EndPosition2, Is.GreaterThanOrEqualTo(result.StartPosition2));
            Assert.That(result.EndPosition2, Is.LessThanOrEqualTo(seq2.Length));
        });
    }

    #endregion

    #region M5: Gaps Removed Yields Substrings

    /// <summary>
    /// M5: Removing gaps from aligned sequences yields substrings of originals.
    /// Evidence: Smith–Waterman traceback rules.
    /// </summary>
    [Test]
    public void LocalAlign_WikipediaExample_GapsRemovedYieldsSubstrings()
    {
        var seq1 = new DnaSequence("TGTTACGG");
        var seq2 = new DnaSequence("GGTTGACTA");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        var gapsRemoved1 = RemoveGaps(result.AlignedSequence1);
        var gapsRemoved2 = RemoveGaps(result.AlignedSequence2);

        Assert.Multiple(() =>
        {
            Assert.That(seq1.Sequence, Does.Contain(gapsRemoved1),
                "Aligned seq1 with gaps removed must be substring of original");
            Assert.That(seq2.Sequence, Does.Contain(gapsRemoved2),
                "Aligned seq2 with gaps removed must be substring of original");
        });
    }

    #endregion

    #region M6: String Overload Parity

    /// <summary>
    /// M6: String overload returns same result as DnaSequence overload.
    /// Evidence: ASSUMPTION (wrapper parity).
    /// </summary>
    [Test]
    public void LocalAlign_StringOverload_MatchesDnaSequenceResult()
    {
        const string seq1 = "TGTTACGG";
        const string seq2 = "GGTTGACTA";

        var dnaResult = SequenceAligner.LocalAlign(new DnaSequence(seq1), new DnaSequence(seq2), WikipediaScoring);
        var stringResult = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        Assert.Multiple(() =>
        {
            Assert.That(stringResult.AlignmentType, Is.EqualTo(AlignmentType.Local));
            Assert.That(stringResult.Score, Is.EqualTo(dnaResult.Score));
            Assert.That(stringResult.AlignedSequence1, Is.EqualTo(dnaResult.AlignedSequence1));
            Assert.That(stringResult.AlignedSequence2, Is.EqualTo(dnaResult.AlignedSequence2));
        });
    }

    #endregion

    #region M7: Empty Input Returns Empty

    /// <summary>
    /// M7: Empty string input returns AlignmentResult.Empty.
    /// Evidence: ASSUMPTION (implementation behavior).
    /// </summary>
    [Test]
    public void LocalAlign_StringOverload_EmptyInput_ReturnsEmpty()
    {
        var result = SequenceAligner.LocalAlign(string.Empty, "ACGT", WikipediaScoring);

        Assert.That(result, Is.EqualTo(AlignmentResult.Empty));
    }

    #endregion

    #region M8: Null Throws ArgumentNullException

    /// <summary>
    /// M8: Null DnaSequence throws ArgumentNullException.
    /// Evidence: ASSUMPTION (implementation behavior).
    /// </summary>
    [Test]
    public void LocalAlign_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.LocalAlign((DnaSequence)null!, new DnaSequence("ACGT")));
    }

    #endregion

    #region S1: Identical Sequences Full Match

    /// <summary>
    /// S1: Identical sequences produce full-length alignment with high score.
    /// Evidence: ASSUMPTION.
    /// </summary>
    [Test]
    public void LocalAlign_IdenticalSequences_FullLengthMatch()
    {
        var seq1 = new DnaSequence("ACGTACGT");
        var seq2 = new DnaSequence("ACGTACGT");

        var result = SequenceAligner.LocalAlign(seq1, seq2);

        Assert.Multiple(() =>
        {
            Assert.That(result.Score, Is.GreaterThan(0), "Identical sequences should have positive score");
            Assert.That(result.AlignedSequence1, Is.EqualTo("ACGTACGT"), "Should match full sequence");
            Assert.That(result.AlignedSequence2, Is.EqualTo("ACGTACGT"), "Should match full sequence");
        });
    }

    #endregion

    #region S2: Dissimilar Sequences - Zero Floor

    /// <summary>
    /// S2: Completely dissimilar sequences produce score >= 0 (zero floor).
    /// Evidence: Smith–Waterman zero floor property.
    /// </summary>
    [Test]
    public void LocalAlign_DissimilarSequences_ScoreNonNegative()
    {
        var seq1 = new DnaSequence("AAAA");
        var seq2 = new DnaSequence("TTTT");

        var result = SequenceAligner.LocalAlign(seq1, seq2, WikipediaScoring);

        // Zero floor ensures score is never negative
        Assert.That(result.Score, Is.GreaterThanOrEqualTo(0),
            "Zero floor: Score must be non-negative even for dissimilar sequences");
    }

    #endregion

    #region Helper Methods

    private static string RemoveGaps(string alignedSequence)
    {
        return new string(alignedSequence.Where(c => c != '-').ToArray());
    }

    #endregion
}
