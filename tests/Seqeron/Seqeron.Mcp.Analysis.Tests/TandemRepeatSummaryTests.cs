using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>tandem_repeat_summary</c> MCP tool.
/// Expected values derived from the microsatellite aggregation on a single (CAG)3 STR
/// and an STR-free sequence, NOT the wrapper output.
/// </summary>
[TestFixture]
public class TandemRepeatSummaryTests
{
    [Test]
    public void TandemRepeatSummary_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.TandemRepeatSummary("CAGCAGCAG", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary(null!, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TandemRepeatSummary("XYZ", 3));
    }

    [Test]
    public void TandemRepeatSummary_Binding_InvokesSuccessfully()
    {
        // (CAG)3 -> one trinucleotide STR spanning the whole 9-mer.
        var s = AnalysisTools.TandemRepeatSummary("CAGCAGCAG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(1));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(9));
            Assert.That(s.PercentageOfSequence, Is.EqualTo(100.0).Within(1e-9));
            Assert.That(s.TrinucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.DinucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.MostFrequentUnit, Is.EqualTo("CAG"));
        });

        // No STRs.
        var none = AnalysisTools.TandemRepeatSummary("ACGT", 3);
        Assert.Multiple(() =>
        {
            Assert.That(none.TotalRepeats, Is.EqualTo(0));
            Assert.That(none.TotalRepeatBases, Is.EqualTo(0));
            Assert.That(none.PercentageOfSequence, Is.EqualTo(0.0).Within(1e-9));
        });
    }
}
