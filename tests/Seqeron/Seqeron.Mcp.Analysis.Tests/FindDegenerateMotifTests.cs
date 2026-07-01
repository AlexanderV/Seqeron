using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_degenerate_motif</c> MCP tool.
/// Expected results derived from the IUPAC matching rules in MotifFinder.FindDegenerateMotif
/// (R = A/G, etc.), matched-sequence and score 1.0. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindDegenerateMotifTests
{
    [Test]
    public void FindDegenerateMotif_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindDegenerateMotif("AGGTAG", "RGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateMotif("", "RGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateMotif("AUGC", "RGG")); // not DNA
        // Invalid IUPAC code in the pattern throws.
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateMotif("ACGT", "AZ"));
    }

    [Test]
    public void FindDegenerateMotif_Binding_InvokesSuccessfully()
    {
        // "AGGTAG", motif "RGG" (R = A/G): only position 0 "AGG" matches.
        var one = AnalysisTools.FindDegenerateMotif("AGGTAG", "RGG").Items;
        Assert.That(one, Has.Length.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(one[0].Position, Is.EqualTo(0));
            Assert.That(one[0].MatchedSequence, Is.EqualTo("AGG"));
            Assert.That(one[0].Pattern, Is.EqualTo("RGG"));
            Assert.That(one[0].Score, Is.EqualTo(1.0).Within(1e-12));
        });

        // "ACG", motif "N" (any base): matches every position 0,1,2.
        var any = AnalysisTools.FindDegenerateMotif("ACG", "N").Items;
        Assert.That(any.Select(m => m.Position), Is.EqualTo(new[] { 0, 1, 2 }));
    }
}
