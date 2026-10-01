using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_approximate_tandem_repeats</c> MCP tool. Expected rows = compiled TRF 4.10.0
/// (`trf craft.fa 2 7 7 80 10 50 500 -h -d -ngs`, 1-based in TRF), locked in RepeatFinder_TrfParameters_Tests
/// (Evidence REP-APPROX-001 §WP6). NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindApproximateTandemRepeatsTests
{
    private const string U1 =
        "CTATCCTAACCCGACCCTAGGAGCGGTTGGCGTGTATGCCGTGAATTTTCTCATTTCCGCTCATTGGTCATTGGTCANTGGTCATTGGTCATTGGTCATTNGTCATTGGTCATTGGTCATTGGTAGACATAATCGTTCTGCCTATATCTGGACAACATCCCGGCGACTTAGGCGACCCACAGA";
    private const string U3 =
        "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC";

    [Test]
    public void FindApproximateTandemRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindApproximateTandemRepeats(U3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, matchProbability: 70));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, maxPeriod: 2001));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, minPeriod: 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindApproximateTandemRepeats(U3, mismatchPenalty: 3, examineUpToMaxPeriodOnly: true));
    }

    [Test]
    public void FindApproximateTandemRepeats_Binding_RecommendedParameters_ReproduceTrfRow()
    {
        // TRF U1 row: 61 124 7 9.1 7 92 0 110 14 14 26 42 1.83 TCATTGG
        var r = AnalysisTools.FindApproximateTandemRepeats(U1).Items.Single();
        Assert.Multiple(() =>
        {
            Assert.That((r.Start + 1, r.Start + r.SpanLength, r.Period, r.ConsensusSize), Is.EqualTo((61, 124, 7, 7)));
            Assert.That(r.CopyNumber, Is.EqualTo(9.1).Within(0.05 + 1e-9));
            Assert.That(((int)(r.PercentMatches + 1e-9), (int)(r.PercentIndels + 1e-9), r.AlignmentScore), Is.EqualTo((92, 0, 110)));
            Assert.That(((int)(r.PercentA + 1e-9), (int)(r.PercentC + 1e-9), (int)(r.PercentG + 1e-9), (int)(r.PercentT + 1e-9)),
                Is.EqualTo((14, 14, 26, 42)));
            Assert.That(r.EntropyTrf, Is.EqualTo(1.829258111162015).Within(1e-12));
            Assert.That(r.Consensus, Is.EqualTo("TCATTGG"));
            Assert.That(r.AlignedSequence, Is.Not.Null);
        });
    }

    [Test]
    public void FindApproximateTandemRepeats_CustomWeights_ReproduceTrfRow()
    {
        // `trf 2 3 5 80 10 40 200` → U3 1 60 2 30.0 2 89 0 105 CA
        var r = AnalysisTools.FindApproximateTandemRepeats(U3, maxPeriod: 200, minScore: 40, mismatchPenalty: 3, indelPenalty: 5).Items.Single();
        Assert.That((r.Start + 1, r.Start + r.SpanLength, r.Period, (int)(r.PercentMatches + 1e-9), r.AlignmentScore, r.Consensus),
            Is.EqualTo((1, 60, 2, 89, 105, "CA")));
    }

    [Test]
    public void FindApproximateTandemRepeats_NoRedundancyElimination_ReportsTrfMultiples()
    {
        // `trf ... -r`: U1 reports periods 7, 14, 21 over 61..124 (all score 110).
        var items = AnalysisTools.FindApproximateTandemRepeats(U1, eliminateRedundancy: false).Items;
        Assert.That(items.Select(r => (r.Start + 1, r.Start + r.SpanLength, r.Period, r.AlignmentScore)),
            Is.EqualTo(new[] { (61, 124, 7, 110), (61, 124, 14, 110), (61, 124, 21, 110) }));
    }
}
