using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>dangling_end_energy</c> MCP tool.
/// Expected values taken directly from the NNDB Turner 2004 dangling-end tables embedded in
/// RnaSecondaryStructure (DanglingEnd3 / DanglingEnd5), keyed closing5' + dangling + closing3'.
/// NOT from the wrapper's output; swapping the 3'/5' flag would change the answer.
/// </summary>
[TestFixture]
public class DanglingEndEnergyTests
{
    [Test]
    public void DanglingEndEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DanglingEndEnergy("G", "C", "A", true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DanglingEndEnergy("", "C", "A", true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DanglingEndEnergy("GG", "C", "A", true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DanglingEndEnergy("G", null!, "A", true));
    }

    [Test]
    public void DanglingEndEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // 3' dangle, key "GAC" -> DanglingEnd3["GAC"] = -1.1
            Assert.That(AnalysisTools.DanglingEndEnergy("G", "C", "A", true).Energy, Is.EqualTo(-1.1).Within(1e-12));
            // 5' dangle, key "GAC" -> DanglingEnd5["GAC"] = -0.5
            Assert.That(AnalysisTools.DanglingEndEnergy("G", "C", "A", false).Energy, Is.EqualTo(-0.5).Within(1e-12));
            // 3' dangle, key "CAG" -> DanglingEnd3["CAG"] = -1.7
            Assert.That(AnalysisTools.DanglingEndEnergy("C", "G", "A", true).Energy, Is.EqualTo(-1.7).Within(1e-12));
            // Unknown key -> 0.0 (no table entry for a non-pair closing context).
            Assert.That(AnalysisTools.DanglingEndEnergy("A", "A", "A", true).Energy, Is.EqualTo(0.0).Within(1e-12));
        });
    }
}
