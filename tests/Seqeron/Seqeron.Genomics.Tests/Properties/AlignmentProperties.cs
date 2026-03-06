using System.Threading;

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
    /// Global alignment score is symmetric: GlobalAlign(a,b).Score == GlobalAlign(b,a).Score.
    /// This follows from the NW recurrence — S(a,b) is symmetric and gap penalties
    /// are the same for both sequences.
    /// Source: Needleman–Wunsch algorithm symmetry (Wikipedia)
    /// </summary>
    [Test]
    [Category("Property")]
    public void GlobalAlign_ScoreSymmetry_ReversedInputs()
    {
        var resultAB = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT");
        var resultBA = SequenceAligner.GlobalAlign("ACGACGT", "ACGTACGT");
        Assert.That(resultAB.Score, Is.EqualTo(resultBA.Score),
            "NW score must be symmetric: score(A,B) == score(B,A)");
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
    /// CancellationToken overload produces the same result as the standard overload.
    /// Verifies the separate code path (non-pooled 2D array) is functionally equivalent.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GlobalAlign_CancellationOverload_SameResultAsStandard()
    {
        using var cts = new CancellationTokenSource();
        var standard = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT");
        var withToken = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGT", null, cts.Token);

        Assert.Multiple(() =>
        {
            Assert.That(withToken.Score, Is.EqualTo(standard.Score),
                "CancellationToken overload must produce same score");
            Assert.That(withToken.AlignedSequence1, Is.EqualTo(standard.AlignedSequence1));
            Assert.That(withToken.AlignedSequence2, Is.EqualTo(standard.AlignedSequence2));
        });
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
    /// Alignment statistics: percentage fields satisfy structural invariants.
    ///   - 0 ≤ Identity ≤ 100
    ///   - 0 ≤ Similarity ≤ 100
    ///   - 0 ≤ GapPercent ≤ 100
    ///   - Identity ≤ Similarity  (matches ≤ matches + mismatches)
    ///   - Similarity + GapPercent ≈ 100  (non-gap + gap = total)
    /// </summary>
    [Test]
    [Category("Property")]
    public void Statistics_PercentageFields_SatisfyInvariants()
    {
        var result = SequenceAligner.GlobalAlign("ACGTACGT", "ACGACGTT");
        var stats = SequenceAligner.CalculateStatistics(result);

        Assert.Multiple(() =>
        {
            Assert.That(stats.Identity, Is.InRange(0.0, 100.0), "Identity in [0,100]");
            Assert.That(stats.Similarity, Is.InRange(0.0, 100.0), "Similarity in [0,100]");
            Assert.That(stats.GapPercent, Is.InRange(0.0, 100.0), "GapPercent in [0,100]");
            Assert.That(stats.Identity, Is.LessThanOrEqualTo(stats.Similarity),
                "Identity ≤ Similarity (matches ≤ matches + mismatches)");
            Assert.That(stats.Similarity + stats.GapPercent, Is.EqualTo(100.0).Within(0.001),
                "Similarity + GapPercent = 100 (all positions accounted for)");
        });
    }
}
