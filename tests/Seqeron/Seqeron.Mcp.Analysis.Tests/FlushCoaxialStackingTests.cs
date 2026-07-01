using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>flush_coaxial_stacking</c> MCP tool.
/// Expected value from RnaSecondaryStructure's NNDB coax unit test
/// (GC onto CG -> -3.42), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FlushCoaxialStackingTests
{
    [Test]
    public void FlushCoaxialStacking_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FlushCoaxialStacking("G", "C", "C", "G"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FlushCoaxialStacking("", "C", "C", "G"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FlushCoaxialStacking("GG", "C", "C", "G"));
    }

    [Test]
    public void FlushCoaxialStacking_Binding_InvokesSuccessfully()
    {
        // Flush stack of GC onto CG -> -3.42 (WC stacking table).
        var e = AnalysisTools.FlushCoaxialStacking("G", "C", "C", "G").Energy;
        Assert.That(e, Is.EqualTo(-3.42).Within(0.01));
    }
}
