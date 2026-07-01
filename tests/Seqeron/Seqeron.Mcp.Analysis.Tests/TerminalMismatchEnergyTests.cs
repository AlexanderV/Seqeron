using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>terminal_mismatch_energy</c> MCP tool.
/// Expected values from RnaSecondaryStructure's NNDB Turner-2004 unit test
/// (GetTerminalMismatchEnergy_MatchesNNDB), NOT the wrapper output.
/// </summary>
[TestFixture]
public class TerminalMismatchEnergyTests
{
    [Test]
    public void TerminalMismatchEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.TerminalMismatchEnergy("G", "C", "A", "A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TerminalMismatchEnergy("", "C", "A", "A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TerminalMismatchEnergy("GG", "C", "A", "A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.TerminalMismatchEnergy("G", null!, "A", "A"));
    }

    [Test]
    public void TerminalMismatchEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // GC closing, AA mismatch -> -1.1.
            Assert.That(AnalysisTools.TerminalMismatchEnergy("G", "C", "A", "A").Energy, Is.EqualTo(-1.1).Within(0.01));
            // CG closing, AA mismatch -> -1.5.
            Assert.That(AnalysisTools.TerminalMismatchEnergy("C", "G", "A", "A").Energy, Is.EqualTo(-1.5).Within(0.01));
            // GC closing, GA mismatch -> -1.6.
            Assert.That(AnalysisTools.TerminalMismatchEnergy("G", "C", "G", "A").Energy, Is.EqualTo(-1.6).Within(0.01));
        });
    }
}
