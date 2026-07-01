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
}
