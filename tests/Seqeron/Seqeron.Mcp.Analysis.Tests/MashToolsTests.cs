using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>B09 A5-10: the Mash text tools (<c>mash_bounds</c>, <c>mash_dist</c>, <c>mash_triangle</c>); the text is that of the installed mash 2.3 (locked byte for byte by the 660-row oracle of the library).</summary>
[TestFixture]
public class MashToolsTests
{
    private const string G1 = ">g1 first\nACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCGACCCTTTAGGGATCGATTTGACCAGGACTTAGCAATGCCGATTAGGCATCG\n";
    private const string G2 = ">g2 second\nACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCGACCCTTAAGGGATCGATTTGACCAGGACTTAGCAATGCCGATTAGGCATCG\n";

    [Test]
    public void MashBounds_FirstLines()
    {
        var r = AnalysisTools.MashBounds();
        Assert.That(r.Text, Does.Contain("\nSketch\t0.05\t0.1\t0.15\t0.2\t0.25\t0.3\t0.35\t0.4\n100\t0.0270708\t0.0867606\tinf"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashBounds(33));
    }

    [Test]
    public void MashDist_And_Triangle_AgreeOnTheDistance()
    {
        var d = AnalysisTools.MashDist([G1], [G2], ["g1.fa"], ["g2.fa"], k: 9, sketchSize: 20);
        var t = AnalysisTools.MashTriangle([G1, G2], ["g1.fa", "g2.fa"], k: 9, sketchSize: 20);
        string distance = d.Text.Split('\t')[2];
        Assert.That(d.Text, Does.StartWith("g1.fa\tg2.fa\t"));
        Assert.That(t.Text, Is.EqualTo($"\t2\ng1.fa\ng2.fa\t{distance}\n"));
        Assert.That(AnalysisTools.MashDist([G1], [G2], k: 9, sketchSize: 20, maxDistance: 0.0).Text, Is.Empty);        // -d 0: only identical pairs pass
        Assert.That(AnalysisTools.MashDist([G1], [G2], k: 9, sketchSize: 20, table: true).Text, Does.StartWith("#query\treference1\nquery1\t"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MashDist([], [G2]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MashTriangle([G1, G2], seed: -1));
    }
}
