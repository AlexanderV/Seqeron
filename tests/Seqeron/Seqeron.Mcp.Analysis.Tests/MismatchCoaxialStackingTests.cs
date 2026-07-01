using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>mismatch_coaxial_stacking</c> MCP tool.
/// Expected values from RnaSecondaryStructure's NNDB coax unit test
/// (GC/AA -> -3.6 = tm(-1.1) + base(-2.1) + WC bonus(-0.4)), NOT the wrapper output.
/// </summary>
[TestFixture]
public class MismatchCoaxialStackingTests
{
    [Test]
    public void MismatchCoaxialStacking_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MismatchCoaxialStacking("G", "C", "A", "A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MismatchCoaxialStacking("", "C", "A", "A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MismatchCoaxialStacking("GG", "C", "A", "A"));
    }

    [Test]
    public void MismatchCoaxialStacking_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // GC closing, AA mismatch -> tm(-1.1) + base(-2.1) + WC bonus(-0.4) = -3.6.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("G", "C", "A", "A").Energy, Is.EqualTo(-3.6).Within(0.01));
            // GU closing, AA mismatch -> tm(-0.3) + base(-2.1) + GU bonus(-0.2) = -2.6.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("G", "U", "A", "A").Energy, Is.EqualTo(-2.6).Within(0.01));
        });
    }
}
