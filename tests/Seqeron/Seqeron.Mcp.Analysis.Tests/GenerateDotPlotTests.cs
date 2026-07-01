using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>generate_dot_plot</c> MCP tool.
/// Expected coordinates from the exact word-match definition (each shared wordSize-mer
/// is a point (i, j)), NOT the wrapper output.
/// </summary>
[TestFixture]
public class GenerateDotPlotTests
{
    [Test]
    public void GenerateDotPlot_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateDotPlot("ATGC", "ATGC", 4, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDotPlot("", "ATGC", 4, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDotPlot("ATGC", "", 4, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDotPlot(null!, "ATGC", 4, 1));
    }

    [Test]
    public void GenerateDotPlot_Binding_InvokesSuccessfully()
    {
        // "ATGCATGC" vs itself, word 4 step 4: ATGC at i=0,4 matches j=0,4.
        var pts = AnalysisTools.GenerateDotPlot("ATGCATGC", "ATGCATGC", 4, 4).Points;
        Assert.That(pts.Select(p => (p.X, p.Y)),
            Is.EquivalentTo(new[] { (0, 0), (0, 4), (4, 0), (4, 4) }));

        // Disjoint sequences share no 4-mer.
        var none = AnalysisTools.GenerateDotPlot("ATGC", "TTTT", 4, 1).Points;
        Assert.That(none, Is.Empty);
    }
}
