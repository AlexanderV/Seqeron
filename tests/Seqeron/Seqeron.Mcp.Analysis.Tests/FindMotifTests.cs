using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_motif</c> MCP tool.
/// Expected positions from the exact-match / overlapping-occurrence definition
/// (suffix-tree order is unspecified), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindMotifTests
{
    [Test]
    public void FindMotif_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindMotif("ATGATG", "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotif("", "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotif(null!, "ATG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotif("ATGATG", ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotif("ATGATG", null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotif("XYZ123", "ATG"));
    }

    [Test]
    public void FindMotif_Binding_InvokesSuccessfully()
    {
        // ATG occurs at 0 and 3 (order from suffix tree is unspecified).
        var atg = AnalysisTools.FindMotif("ATGATG", "ATG").Positions;
        Assert.That(atg, Is.EquivalentTo(new[] { 0, 3 }));

        // Overlapping AA in AAAA -> {0,1,2}.
        var aa = AnalysisTools.FindMotif("AAAA", "AA").Positions;
        Assert.That(aa, Is.EquivalentTo(new[] { 0, 1, 2 }));

        // Absent motif -> empty.
        var none = AnalysisTools.FindMotif("ATGATG", "CCC").Positions;
        Assert.That(none, Is.Empty);
    }
}
