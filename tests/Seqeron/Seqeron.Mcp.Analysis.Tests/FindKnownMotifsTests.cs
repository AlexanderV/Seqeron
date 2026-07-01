using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_known_motifs</c> MCP tool.
/// Expected values from the exact set-matching definition (ascending positions,
/// misses omitted, overlapping reported), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindKnownMotifsTests
{
    [Test]
    public void FindKnownMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindKnownMotifs("ATGATGCC", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs(null!, new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("XYZ", new[] { "ATG" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindKnownMotifs("ATGATGCC", null!));
    }

    [Test]
    public void FindKnownMotifs_Binding_InvokesSuccessfully()
    {
        // "ATGATGCC": ATG at 0,3; CC at 6; GG absent (omitted).
        var m = AnalysisTools.FindKnownMotifs("ATGATGCC", new[] { "ATG", "CC", "GG" }).Matches;
        Assert.Multiple(() =>
        {
            Assert.That(m.ContainsKey("ATG"), Is.True);
            Assert.That(m["ATG"], Is.EqualTo(new[] { 0, 3 }));
            Assert.That(m.ContainsKey("CC"), Is.True);
            Assert.That(m["CC"], Is.EqualTo(new[] { 6 }));
            Assert.That(m.ContainsKey("GG"), Is.False);
        });

        // Overlapping AA in AAAA -> 0,1,2.
        var overlap = AnalysisTools.FindKnownMotifs("AAAA", new[] { "AA" }).Matches;
        Assert.That(overlap["AA"], Is.EqualTo(new[] { 0, 1, 2 }));
    }
}
