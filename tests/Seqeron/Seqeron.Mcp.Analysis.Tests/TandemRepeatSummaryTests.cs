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
            Assert.That(s.LongestRepeat, Is.EqualTo(new MicrosatelliteItem(0, "CAG", 3, 9, "Trinucleotide")));
            Assert.That(s.PentanucleotideRepeats, Is.EqualTo(0));
            Assert.That(s.HexanucleotideRepeats, Is.EqualTo(0));
        });

        // No STRs.
        var none = AnalysisTools.TandemRepeatSummary("ACGT", 3);
        Assert.Multiple(() =>
        {
            Assert.That(none.TotalRepeats, Is.EqualTo(0));
            Assert.That(none.TotalRepeatBases, Is.EqualTo(0));
            Assert.That(none.PercentageOfSequence, Is.EqualTo(0.0).Within(1e-9));
            // No STR → no longest repeat (previously a default Position-0 item with a null unit).
            Assert.That(none.LongestRepeat, Is.Null);
            Assert.That(none.MostFrequentUnit, Is.Null);
        });
    }

    [Test]
    public void TandemRepeatSummary_PentaHexa_DocExample3()
    {
        // docs/mcp/tools/analysis/tandem_repeat_summary.md Example 3; values from the independent Python
        // brute-force maximal-run reference: (AAAGA)4@0, (TTAGGG)4@22 and 8 homopolymer runs; 71 repeat bases,
        // union coverage [0,20) U [22,46) = 44/46.
        var s = AnalysisTools.TandemRepeatSummary("AAAGAAAAGAAAAGAAAAGACCTTAGGGTTAGGGTTAGGGTTAGGG", 3);
        Assert.Multiple(() =>
        {
            Assert.That(s.TotalRepeats, Is.EqualTo(10));
            Assert.That(s.TotalRepeatBases, Is.EqualTo(71));
            Assert.That(s.PercentageOfSequence, Is.EqualTo(95.65217391304348).Within(1e-12));
            Assert.That(s.MononucleotideRepeats, Is.EqualTo(8));
            Assert.That(s.PentanucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.HexanucleotideRepeats, Is.EqualTo(1));
            Assert.That(s.LongestRepeat, Is.EqualTo(new MicrosatelliteItem(22, "TTAGGG", 4, 24, "Hexanucleotide")));
            Assert.That(s.MostFrequentUnit, Is.EqualTo("A"));
        });
    }
}
