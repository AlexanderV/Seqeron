namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for global, local, semi-global, and multiple alignment.
/// Verifies structural invariants of alignment results.
///
/// Test Units: ALIGN-GLOBAL-001, ALIGN-LOCAL-001, ALIGN-SEMI-001, ALIGN-MULTI-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Alignment")]
public class AlignmentProperties
{
    /// <summary>
    /// Global alignment: both aligned sequences must have equal length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GlobalAlign_AlignedSequences_HaveEqualLength()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT");
        Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
    }

    /// <summary>
    /// Local alignment: both aligned sequences have equal length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LocalAlign_AlignedSequences_HaveEqualLength()
    {
        var result = SequenceAligner.LocalAlign("ACGTACGTACGT", "TACGTAC");
        Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
    }

    /// <summary>
    /// Semi-global alignment: aligned sequences have equal length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void SemiGlobalAlign_AlignedSequences_HaveEqualLength()
    {
        var s1 = new DnaSequence("ACGTACGTACGTACGT");
        var s2 = new DnaSequence("TACGTACG");
        var result = SequenceAligner.SemiGlobalAlign(s1, s2);
        Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
    }

    /// <summary>
    /// Identity alignment (same sequence) should have score ≥ 0 and perfect match.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GlobalAlign_IdenticalSequences_MaxScore()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGTACGT");
        Assert.That(result.Score, Is.GreaterThan(0));
        Assert.That(result.AlignedSequence1, Is.EqualTo(result.AlignedSequence2));
    }

    /// <summary>
    /// Local alignment score is non-negative.
    /// </summary>
    [Test]
    [Category("Property")]
    public void LocalAlign_Score_IsNonNegative()
    {
        var result = SequenceAligner.LocalAlign("AAAA", "TTTT");
        Assert.That(result.Score, Is.GreaterThanOrEqualTo(0));
    }

    /// <summary>
    /// Multiple alignment: all aligned sequences have equal length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MultipleAlign_AllSequences_HaveEqualLength()
    {
        var seqs = new[]
        {
            new DnaSequence("ACGTACGT"),
            new DnaSequence("ACGACGT"),
            new DnaSequence("ACGTAGT")
        };
        var result = SequenceAligner.MultipleAlign(seqs);
        int len = result.AlignedSequences[0].Length;

        for (int i = 1; i < result.AlignedSequences.Length; i++)
            Assert.That(result.AlignedSequences[i].Length, Is.EqualTo(len),
                $"Aligned sequence {i} length differs");
    }

    /// <summary>
    /// Multiple alignment consensus length equals aligned sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void MultipleAlign_ConsensusLength_EqualsAlignedLength()
    {
        var seqs = new[]
        {
            new DnaSequence("ACGTACGT"),
            new DnaSequence("ACGACGT"),
            new DnaSequence("ACGTACG")
        };
        var result = SequenceAligner.MultipleAlign(seqs);
        Assert.That(result.Consensus.Length, Is.EqualTo(result.AlignedSequences[0].Length));
    }

    /// <summary>
    /// Alignment statistics: identity + mismatch + gap percentages should account for full alignment.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Statistics_Matches_Plus_Mismatches_Plus_Gaps_EqualsLength()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGTT");
        var stats = SequenceAligner.CalculateStatistics(result);
        Assert.That(stats.Matches + stats.Mismatches + stats.Gaps, Is.EqualTo(stats.AlignmentLength));
    }

    /// <summary>
    /// Alignment statistics: identity is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void Statistics_Identity_InRange()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGTT");
        var stats = SequenceAligner.CalculateStatistics(result);
        Assert.That(stats.Identity, Is.InRange(0.0, 100.0));
    }
}
