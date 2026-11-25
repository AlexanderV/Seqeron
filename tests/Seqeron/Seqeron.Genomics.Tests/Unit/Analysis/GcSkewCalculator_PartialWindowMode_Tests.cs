using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B10 R1-WP1 (A15-2): additive <c>GcSkewCalculator.CalculateWindowedGcContent(seq, w, step, fraction, mode, includePartialWindow)</c>
/// and the <c>SequenceStatistics.CalculateGcContentProfile(…, includePartialWindow)</c> projections.
/// Oracle: tools/review/oracles/B10/gc_profile_partial_oracle.py — Biopython 1.88
/// <c>[gc_fraction(seq[i:i+w], ambiguous=mode) * 100 for i in range(0, len(seq), step)]</c> (GC_skew windowing).
/// </summary>
[TestFixture]
public class GcSkewCalculator_PartialWindowMode_Tests
{
    private const string Seq = "ACGTNSWRYacgtnGGCC";
    private const int W = 5, Step = 4;

    private static readonly object[][] Oracle =
    {
        new object[] { SequenceExtensions.GcAmbiguityMode.Remove, new[] { 50.0, 50.0, 50.0, 75.0, 100.0 } },
        new object[] { SequenceExtensions.GcAmbiguityMode.Ignore, new[] { 40.0, 20.0, 40.0, 60.0, 100.0 } },
        new object[] { SequenceExtensions.GcAmbiguityMode.Weighted, new[] { 50.0, 50.0, 50.0, 70.0, 100.0 } },
    };

    [TestCaseSource(nameof(Oracle))]
    public void CalculateWindowedGcContent_PartialWindow_Mode_MatchesBiopython(
        SequenceExtensions.GcAmbiguityMode mode, double[] expected)
    {
        var points = GcSkewCalculator.CalculateWindowedGcContent(Seq, W, Step, false, mode, includePartialWindow: true).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points.Select(p => p.GcContent).ToArray(), Is.EqualTo(expected).Within(1e-12));
            Assert.That(points.Select(p => p.WindowStart).ToArray(), Is.EqualTo(new[] { 0, 4, 8, 12, 16 }));
            Assert.That(points[^1].WindowEnd, Is.EqualTo(Seq.Length - 1), "last window truncated at the sequence end");
        });
    }

    [TestCaseSource(nameof(Oracle))]
    public void SequenceStatistics_GcProfile_PartialWindow_Mode_ProjectsSameValues(
        SequenceExtensions.GcAmbiguityMode mode, double[] expected)
    {
        Assert.That(
            SequenceStatistics.CalculateGcContentProfile(Seq, W, Step, false, mode, true).ToArray(),
            Is.EqualTo(expected).Within(1e-12));
    }

    [TestCaseSource(nameof(Oracle))]
    public void CalculateWindowedGcContent_NoPartialWindow_Mode_EqualsExistingOverload(
        SequenceExtensions.GcAmbiguityMode mode, double[] _)
    {
        Assert.That(
            GcSkewCalculator.CalculateWindowedGcContent(Seq, W, Step, true, mode, includePartialWindow: false).Select(p => p.GcContent).ToArray(),
            Is.EqualTo(GcSkewCalculator.CalculateWindowedGcContent(Seq, W, Step, true, mode).Select(p => p.GcContent).ToArray()));
    }

    [Test]
    public void CalculateWindowedGcContent_PartialWindow_Mode_RejectsBadWindowAndStep()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedGcContent(Seq, 0, 1, false, SequenceExtensions.GcAmbiguityMode.Remove, true));
        Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.CalculateWindowedGcContent(Seq, 1, 0, false, SequenceExtensions.GcAmbiguityMode.Remove, true));
    }
}
