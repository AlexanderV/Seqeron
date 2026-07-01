using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>predict_coiled_coils</c> MCP tool.
/// Expected values from ProteinMotifFinder's own unit test
/// (ProteinMotifFinder_PredictCoiledCoils_Tests, (LAALAAA)x5 -> region 0..34 score 1.0;
/// poly-G -> none), NOT the wrapper output.
/// </summary>
[TestFixture]
public class PredictCoiledCoilsTests
{
    private static string Repeat(string unit, int n) => string.Concat(Enumerable.Repeat(unit, n));

    [Test]
    public void PredictCoiledCoils_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PredictCoiledCoils(Repeat("LAALAAA", 5)));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictCoiledCoils(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PredictCoiledCoils(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.PredictCoiledCoils("LAAL", 0));
    }

    [Test]
    public void PredictCoiledCoils_Binding_InvokesSuccessfully()
    {
        // (LAALAAA)x5 = 35 aa perfect heptad -> one region 0..34, score 1.0.
        var regions = AnalysisTools.PredictCoiledCoils(Repeat("LAALAAA", 5)).Items;
        Assert.Multiple(() =>
        {
            Assert.That(regions, Has.Length.EqualTo(1));
            Assert.That(regions[0].Start, Is.EqualTo(0));
            Assert.That(regions[0].End, Is.EqualTo(34));
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10));
        });

        // Poly-G has no heptad periodicity.
        var none = AnalysisTools.PredictCoiledCoils(new string('G', 40)).Items;
        Assert.That(none, Is.Empty);
    }
}
