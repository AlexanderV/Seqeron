using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>detect_rearrangements</c> MCP tool.
/// Expected values traced by hand through the signed-permutation breakpoint model in
/// ComparativeGenomics.DetectRearrangements (Bafna &amp; Pevzner 1998). Identity gene order
/// yields no breakpoints; a single strand flip (local inversion) creates two Inversion
/// breakpoints. NOT taken from the wrapper's output.
/// </summary>
[TestFixture]
public class DetectRearrangementsTests
{
    private static GeneInput G(string id, char strand, int start) => new(id, "G1", start, start + 10, strand);

    [Test]
    public void DetectRearrangements_Schema_ValidatesCorrectly()
    {
        var g1 = new[] { G("A", '+', 0), G("B", '+', 20) };
        var g2 = new[] { new GeneInput("A", "G2", 0, 10, '+'), new GeneInput("B", "G2", 20, 30, '+') };
        var map = new Dictionary<string, string> { ["A"] = "A", ["B"] = "B" };

        Assert.DoesNotThrow(() => AnalysisTools.DetectRearrangements(g1, g2, map));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DetectRearrangements(null!, g2, map));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DetectRearrangements(g1, Array.Empty<GeneInput>(), map));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DetectRearrangements(g1, g2, null!));
    }

    [Test]
    public void DetectRearrangements_IdenticalOrder_HasNoBreakpoints()
    {
        var g1 = new[] { G("A", '+', 0), G("B", '+', 20), G("C", '+', 40) };
        var g2 = new[]
        {
            new GeneInput("A", "G2", 0, 10, '+'),
            new GeneInput("B", "G2", 20, 30, '+'),
            new GeneInput("C", "G2", 40, 50, '+')
        };
        var map = new Dictionary<string, string> { ["A"] = "A", ["B"] = "B", ["C"] = "C" };

        var items = AnalysisTools.DetectRearrangements(g1, g2, map).Items;
        Assert.That(items, Is.Empty);
    }

    [Test]
    public void DetectRearrangements_LocalInversion_ProducesInversionBreakpoints()
    {
        // Gene B is on the opposite strand in genome 2 -> signed permutation [1,-2,3]
        // -> two breakpoints at boundaries (1->-2) and (-2->3), both classified Inversion.
        var g1 = new[] { G("A", '+', 0), G("B", '+', 20), G("C", '+', 40) };
        var g2 = new[]
        {
            new GeneInput("A", "G2", 0, 10, '+'),
            new GeneInput("B", "G2", 20, 30, '-'),
            new GeneInput("C", "G2", 40, 50, '+')
        };
        var map = new Dictionary<string, string> { ["A"] = "A", ["B"] = "B", ["C"] = "C" };

        var items = AnalysisTools.DetectRearrangements(g1, g2, map).Items;
        Assert.That(items, Has.Length.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(items[0].Type, Is.EqualTo("Inversion"));
            Assert.That(items[1].Type, Is.EqualTo("Inversion"));
            Assert.That(items[0].GenomeId, Is.EqualTo("G1"));
            // First breakpoint anchors at gene B (genome-1 index 1, Start 20);
            // second at gene C (index 2, Start 40).
            Assert.That(items[0].Position, Is.EqualTo(20));
            Assert.That(items[1].Position, Is.EqualTo(40));
        });
    }
}
