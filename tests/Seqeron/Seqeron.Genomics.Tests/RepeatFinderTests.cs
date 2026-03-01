using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for RepeatFinder: TandemRepeatSummary + InvertedRepeat validation.
/// Microsatellite detection tests consolidated into RepeatFinder_Microsatellite_Tests.cs (REP-STR-001).
/// </summary>
[TestFixture]
public class RepeatFinderTests
{
    #region Tandem Repeat Summary Tests

    [Test]
    public void GetTandemRepeatSummary_MixedRepeats_CorrectSummary()
    {
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACACATGATGATG");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.TotalRepeats, Is.GreaterThan(0));
        Assert.That(summary.TotalRepeatBases, Is.GreaterThan(0));
        Assert.That(summary.PercentageOfSequence, Is.GreaterThan(0));
    }

    [Test]
    public void GetTandemRepeatSummary_NoRepeats_ZeroSummary()
    {
        var sequence = new DnaSequence("ACGT");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.TotalRepeats, Is.EqualTo(0));
        Assert.That(summary.TotalRepeatBases, Is.EqualTo(0));
    }

    [Test]
    public void GetTandemRepeatSummary_MononucleotideCount_Correct()
    {
        var sequence = new DnaSequence("AAAAAATTTTTGGGGGCCCCC");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.MononucleotideRepeats, Is.EqualTo(4)); // A, T, G, C runs
    }

    [Test]
    public void GetTandemRepeatSummary_LongestRepeat_Identified()
    {
        var sequence = new DnaSequence("AAACAGCAGCAGCAGCAGCAGAAA"); // 6x CAG = 18bp
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.LongestRepeat, Is.Not.Null);
        Assert.That(summary.LongestRepeat!.Value.TotalLength, Is.EqualTo(18));
    }

    #endregion

    #region Edge Cases and Validation

    [Test]
    public void FindInvertedRepeats_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindInvertedRepeats((DnaSequence)null!, 4, 10, 3).ToList());
    }

    #endregion
}
