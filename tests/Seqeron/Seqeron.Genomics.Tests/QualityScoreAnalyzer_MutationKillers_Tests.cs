using System.Linq;
using NUnit.Framework;
using static Seqeron.Genomics.IO.QualityScoreAnalyzer;
using QualityEncoding = Seqeron.Genomics.IO.QualityScoreAnalyzer.QualityEncoding;
using FastqRecord = Seqeron.Genomics.IO.QualityScoreAnalyzer.FastqRecord;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// QUALITY-PHRED/STATS-001 mutation killers: exact-value / boundary tests for the FASTQ quality
/// helpers the canonical fixtures only smoke-tested — Phred encoding auto-detection (Cock et al. 2010),
/// end/sliding-window quality trimming, read filtering, low-quality masking, max-quality consensus and
/// low-quality region detection. Q-score thresholds are inclusive (a base at exactly Qn is counted).
/// </summary>
[TestFixture]
public class QualityScoreAnalyzer_MutationKillers_Tests
{
    private const double Tol = 1e-9;

    // Build a Phred+33 quality string from explicit Phred scores.
    private static string Q(params int[] scores) => new string(scores.Select(s => (char)(s + 33)).ToArray());

    #region DetectEncoding boundaries (Cock et al. 2010)

    [TestCase("!I", QualityEncoding.Phred33)]   // minChar 33 (< 59) ⇒ Phred+33
    [TestCase("@K", QualityEncoding.Phred64)]   // min 64 (== offset), max 75 (> 74) ⇒ Phred+64
    [TestCase("@J", QualityEncoding.Phred33)]   // max 74 (NOT > 74) ⇒ Phred+33 (kills > → >=)
    [TestCase("@F", QualityEncoding.Phred33)]   // max 70 (NOT > 74) ⇒ Phred+33 (kills && → ||)
    public void DetectEncoding_Boundaries(string qual, QualityEncoding expected)
        => Assert.That(DetectEncoding(qual), Is.EqualTo(expected));

    #endregion

    #region QualityTrim — inclusive end thresholds

    [Test]
    public void QualityTrim_KeepsTerminalBasesExactlyAtThreshold()
    {
        // Both ends at exactly Q20 (== minQuality) must be kept (kills < → <=).
        var r = QualityTrim("ACG", Q(20, 30, 20), minQuality: 20);
        Assert.That(r.Sequence, Is.EqualTo("ACG"));
        Assert.That(r.TrimmedFromStart, Is.EqualTo(0));
        Assert.That(r.TrimmedFromEnd, Is.EqualTo(0));
    }

    [Test]
    public void QualityTrim_SingleHighQualityBaseIsKept()
    {
        // start == end == 0 ⇒ NOT empty (kills start > end → start >= end).
        var r = QualityTrim("A", Q(30), minQuality: 20);
        Assert.That(r.FinalLength, Is.EqualTo(1));
        Assert.That(r.Sequence, Is.EqualTo("A"));
    }

    [Test]
    public void QualityTrim_TrimsLowQualityEnds()
    {
        var r = QualityTrim("AACGTA", Q(5, 30, 30, 30, 30, 5), minQuality: 20);
        Assert.That(r.Sequence, Is.EqualTo("ACGT"));
        Assert.That(r.TrimmedFromStart, Is.EqualTo(1));
        Assert.That(r.TrimmedFromEnd, Is.EqualTo(1));
        Assert.That(r.FinalLength, Is.EqualTo(4));
    }

    #endregion

    #region SlidingWindowTrim

    [Test]
    public void SlidingWindowTrim_TrimsTailWhereWindowAverageDrops()
    {
        // First 4 bases Q30, last 4 Q5; window 4, minAvg 20. The window [1..4] (avg 23.75) is the
        // last whose average clears the threshold ⇒ cutoff = 5.
        var r = SlidingWindowTrim("AAAAAAAA", Q(30, 30, 30, 30, 5, 5, 5, 5), windowSize: 4, minAverageQuality: 20);
        Assert.That(r.FinalLength, Is.EqualTo(5));
        Assert.That(r.Sequence, Is.EqualTo("AAAAA"));
        Assert.That(r.TrimmedFromEnd, Is.EqualTo(3));
    }

    [Test]
    public void SlidingWindowTrim_AllHighQualityKeepsFullRead()
    {
        var r = SlidingWindowTrim("AAAAAAAA", Q(30, 30, 30, 30, 30, 30, 30, 30), windowSize: 4, minAverageQuality: 20);
        Assert.That(r.FinalLength, Is.EqualTo(8));
        Assert.That(r.TrimmedFromEnd, Is.EqualTo(0));
    }

    #endregion

    #region FilterReads

    private static FastqRecord Read(string id, int length, int phred)
        => new(id, new string('A', length), Q(Enumerable.Repeat(phred, length).ToArray()));

    [Test]
    public void FilterReads_LengthBoundsAreInclusive()
    {
        var reads = new[]
        {
            Read("L4", 4, 30), Read("L5", 5, 30), Read("L10", 10, 30), Read("L11", 11, 30),
        };
        var kept = FilterReads(reads, minLength: 5, maxLength: 10).Select(r => r.Id).ToList();
        Assert.That(kept, Is.EqualTo(new[] { "L5", "L10" })); // 4 and 11 excluded
    }

    [Test]
    public void FilterReads_MeanQualityThresholdIsInclusiveLowerBound()
    {
        // mean == minMeanQuality must pass (kills mean < thr → mean <= thr would drop it).
        var reads = new[] { Read("Q20", 10, 20) };
        var kept = FilterReads(reads, minMeanQuality: 20).Select(r => r.Id).ToList();
        Assert.That(kept, Is.EqualTo(new[] { "Q20" }));
    }

    #endregion

    #region MaskLowQualityBases

    [Test]
    public void MaskLowQualityBases_KeepsBasesAtExactlyThreshold()
        // All bases at exactly Q20 are kept, not masked (kills >= → >).
        => Assert.That(MaskLowQualityBases("ACGT", Q(20, 20, 20, 20), minQuality: 20), Is.EqualTo("ACGT"));

    [Test]
    public void MaskLowQualityBases_MasksBasesBelowThreshold()
        => Assert.That(MaskLowQualityBases("ACGT", Q(30, 30, 5, 30), minQuality: 20), Is.EqualTo("ACNT"));

    #endregion

    #region CalculateConsensusQuality

    [Test]
    public void ConsensusQuality_TakesMaximumPerPosition()
    {
        // Per position the most confident (highest) Phred is chosen: max(Q40,Q20) = Q40 ⇒ 'I'.
        Assert.That(CalculateConsensusQuality(new[] { Q(40, 40), Q(20, 20) }), Is.EqualTo(Q(40, 40)));
    }

    [Test]
    public void ConsensusQuality_HandlesRaggedLengths()
    {
        // Shorter strings contribute only where present; longer positions keep their own score.
        Assert.That(CalculateConsensusQuality(new[] { Q(40, 40, 40), Q(20) }), Is.EqualTo(Q(40, 40, 40)));
    }

    #endregion

    #region FindLowQualityRegions

    [Test]
    public void FindLowQualityRegions_DetectsInternalRegion()
    {
        // Phred [5,5,5,30,30], window 2, maxQuality 15: windows at i=0,1 are low (avg 5) then i=2
        // (avg 17.5) ends the region ⇒ one region (start 0, end i+window−1 = 3, mean 5).
        var regions = FindLowQualityRegions(Q(5, 5, 5, 30, 30), windowSize: 2, maxQuality: 15).ToList();
        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].start, Is.EqualTo(0));
        Assert.That(regions[0].end, Is.EqualTo(3));
        Assert.That(regions[0].meanQuality, Is.EqualTo(5.0).Within(Tol));
    }

    [Test]
    public void FindLowQualityRegions_DetectsRegionExtendingToEnd()
    {
        // Phred [30,30,5,5,5], window 2, maxQuality 15: low windows at i=2,3 run to the end ⇒
        // final region (start 2, end length−1 = 4, mean 5).
        var regions = FindLowQualityRegions(Q(30, 30, 5, 5, 5), windowSize: 2, maxQuality: 15).ToList();
        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].start, Is.EqualTo(2));
        Assert.That(regions[0].end, Is.EqualTo(4));
        Assert.That(regions[0].meanQuality, Is.EqualTo(5.0).Within(Tol));
    }

    #endregion
}
