using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>bulge_loop_energy</c> MCP tool.
/// Expected values taken from the algorithm's own unit tests
/// (RnaSecondaryStructureTests, NNDB Turner 2004 bulge.html), NOT the wrapper's output.
/// </summary>
[TestFixture]
public class BulgeLoopEnergyTests
{
    [Test]
    public void BulgeLoopEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.BulgeLoopEnergy(1, "A", "G", "C", "G", "C"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BulgeLoopEnergy(1, "", "G", "C", "G", "C"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BulgeLoopEnergy(1, "AA", "G", "C", "G", "C"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BulgeLoopEnergy(1, "A", null!, "C", "G", "C"));
    }

    [Test]
    public void BulgeLoopEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // n=1 'A' between GC/GC: init(1)=3.8 + stacking GG/CC=-3.26 = 0.54.
            Assert.That(AnalysisTools.BulgeLoopEnergy(1, "A", "G", "C", "G", "C").Energy,
                Is.EqualTo(0.54).Within(0.01));
            // n=1 special-C bonus: init(1)=3.8 + stacking GC/CG=-3.42 + (-0.9) = -0.52.
            Assert.That(AnalysisTools.BulgeLoopEnergy(1, "C", "G", "C", "C", "G").Energy,
                Is.EqualTo(-0.52).Within(0.01));
            // n=3 between GC/AU: init(3)=3.2 + AU terminal penalty 0.45 = 3.65.
            Assert.That(AnalysisTools.BulgeLoopEnergy(3, "A", "G", "C", "A", "U").Energy,
                Is.EqualTo(3.65).Within(0.01));
            // n=2 both AU/GU ends: init(2)=2.8 + 0.45 + 0.45 = 3.7.
            Assert.That(AnalysisTools.BulgeLoopEnergy(2, "A", "A", "U", "G", "U").Energy,
                Is.EqualTo(3.7).Within(0.01));
        });
    }
}
