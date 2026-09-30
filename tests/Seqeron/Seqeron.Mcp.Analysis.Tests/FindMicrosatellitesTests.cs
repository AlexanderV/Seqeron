using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_microsatellites</c> MCP tool.
/// Expected values from the STR definition (unit x count, type classification) on
/// canonical (CA)4 and (CAG)3 repeats, NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindMicrosatellitesTests
{
    [Test]
    public void FindMicrosatellites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindMicrosatellites("CACACACA", 2, 6, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMicrosatellites("", 2, 6, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMicrosatellites(null!, 2, 6, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindMicrosatellites("CACACACA", 0, 6, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindMicrosatellites("CACACACA", 2, 6, 1));
    }

    [Test]
    public void FindMicrosatellites_Binding_InvokesSuccessfully()
    {
        // (CA)4 -> one dinucleotide STR: unit CA, count 4, span 8.
        var ca = AnalysisTools.FindMicrosatellites("CACACACA", 2, 6, 3).Items;
        var dinuc = ca.Single(i => i.RepeatUnit == "CA");
        Assert.Multiple(() =>
        {
            Assert.That(dinuc.Position, Is.EqualTo(0));
            Assert.That(dinuc.RepeatCount, Is.EqualTo(4));
            Assert.That(dinuc.TotalLength, Is.EqualTo(8));
            Assert.That(dinuc.RepeatType, Is.EqualTo("Dinucleotide"));
        });

        // (CAG)3 -> trinucleotide STR: unit CAG, count 3, span 9.
        var cag = AnalysisTools.FindMicrosatellites("CAGCAGCAG", 3, 6, 3).Items;
        var trinuc = cag.Single(i => i.RepeatUnit == "CAG");
        Assert.Multiple(() =>
        {
            Assert.That(trinuc.RepeatCount, Is.EqualTo(3));
            Assert.That(trinuc.TotalLength, Is.EqualTo(9));
            Assert.That(trinuc.RepeatType, Is.EqualTo("Trinucleotide"));
        });
    }

    [Test]
    public void FindMicrosatellites_MisaThresholdsAndCompounds_MatchMisaPl()
    {
        // perl misa.pl (default misa.ini 1-10 2-6 3-5 4-5 5-5 6-5, interruptions 100):
        // "c (TA)6tccgt(GA)7ttttt(A)12 48 4 51".
        var r = AnalysisTools.FindMicrosatellites(
            "ACGTATATATATATATccgtGAGAGAGAGAGAGAtttttAAAAAAAAAAAAT", misaThresholds: true, maxCompoundInterruption: 100);
        Assert.Multiple(() =>
        {
            Assert.That(r.Items.Select(i => (i.Position, i.RepeatUnit, i.RepeatCount)),
                Is.EqualTo(new[] { (39, "A", 12), (3, "TA", 6), (20, "GA", 7) }));
            Assert.That(r.Compounds, Has.Length.EqualTo(1));
            Assert.That(r.Compounds![0].Notation, Is.EqualTo("(TA)6tccgt(GA)7ttttt(A)12"));
            Assert.That(r.Compounds[0].Type, Is.EqualTo("c"));
            Assert.That((r.Compounds[0].Start, r.Compounds[0].End, r.Compounds[0].Length), Is.EqualTo((3, 51, 48)));
            Assert.That(r.Compounds[0].Components, Has.Length.EqualTo(3));
        });

        // Defaults unchanged: no compounds requested → null; CCCCCCCCC (9) is an STR only without MISA thresholds.
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.FindMicrosatellites("CACACACA", 2, 6, 3).Compounds, Is.Null);
            Assert.That(AnalysisTools.FindMicrosatellites("CCCCCCCCC").Items, Has.Length.EqualTo(1));
            Assert.That(AnalysisTools.FindMicrosatellites("CCCCCCCCC", misaThresholds: true).Items, Is.Empty);
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindMicrosatellites("CACACACA", 7, 8, 3, misaThresholds: true));
        });
    }
}
