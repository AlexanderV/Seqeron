using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_exact_motif</c> MCP tool.
/// Expected positions from the exact-match/overlapping definition (order unspecified),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindExactMotifTests
{
    [Test]
    public void FindExactMotif_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindExactMotif("ATGATG", "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindExactMotif("", "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindExactMotif(null!, "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindExactMotif("ATGATG", ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindExactMotif("ATGATG", null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindExactMotif("XYZ", "ATG"));
    }

    [Test]
    public void FindExactMotif_Binding_InvokesSuccessfully()
    {
        var atg = AnalysisTools.FindExactMotif("ATGATG", "ATG").Positions;
        Assert.That(atg, Is.EquivalentTo(new[] { 0, 3 }));

        var aa = AnalysisTools.FindExactMotif("AAAA", "AA").Positions;
        Assert.That(aa, Is.EquivalentTo(new[] { 0, 1, 2 }));

        var none = AnalysisTools.FindExactMotif("ATGATG", "CCC").Positions;
        Assert.That(none, Is.Empty);
    }
}
