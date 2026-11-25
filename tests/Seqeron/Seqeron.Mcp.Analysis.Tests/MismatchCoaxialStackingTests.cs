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
            // NNDB coaxial stacking §8.3 / turner_2004_coax: tm − 2.1, −0.4 / −0.2 only if the MISMATCH could pair.
            // GC closing, AA mismatch (cannot pair) -> tm(-1.1) + (-2.1) = -3.2.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("G", "C", "A", "A").Energy, Is.EqualTo(-3.2).Within(1e-9));
            // GU closing, AA mismatch -> tm(-0.3) + (-2.1) = -2.4.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("G", "U", "A", "A").Energy, Is.EqualTo(-2.4).Within(1e-9));
            // NNDB §8.4 example: AU followed by a GG mismatch = -0.8 - 2.1 = -2.9.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("A", "U", "G", "G").Energy, Is.EqualTo(-2.9).Within(1e-9));
            // GC closing, AU mismatch (Watson-Crick capable) -> tm(-1.5) - 2.1 - 0.4 = -4.0.
            Assert.That(AnalysisTools.MismatchCoaxialStacking("G", "C", "A", "U").Energy, Is.EqualTo(-4.0).Within(1e-9));
        });
    }
}
