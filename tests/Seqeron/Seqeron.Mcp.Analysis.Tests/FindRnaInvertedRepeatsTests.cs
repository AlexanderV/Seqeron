using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_rna_inverted_repeats</c> MCP tool.
/// Expected values from RnaSecondaryStructure's own unit test
/// (FindInvertedRepeats on "GCGCAAAAAAGCGC" finds the GCGC...GCGC hairpin stem),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindRnaInvertedRepeatsTests
{
    [Test]
    public void FindRnaInvertedRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindRnaInvertedRepeats("GCGCAAAAAAGCGC", 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRnaInvertedRepeats("", 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRnaInvertedRepeats(null!, 4, 3));
    }

    [Test]
    public void FindRnaInvertedRepeats_Binding_InvokesSuccessfully()
    {
        // GCGC (5') pairs antiparallel with GCGC (3') across a 6-nt loop.
        var items = AnalysisTools.FindRnaInvertedRepeats("GCGCAAAAAAGCGC", 4, 3).Items;
        Assert.Multiple(() =>
        {
            Assert.That(items, Is.Not.Empty);
            // The left arm starts at 0 with a stem of at least 4 bp; the right arm lies in
            // the 3' GCGC region (Start2 >= 10) and the arms don't overlap the 6-nt loop.
            var stem = items.First(i => i.Start1 == 0 && i.Length >= 4);
            Assert.That(stem.End1, Is.EqualTo(stem.Start1 + stem.Length - 1));
            Assert.That(stem.Start2, Is.GreaterThanOrEqualTo(10));
            Assert.That(stem.End2, Is.EqualTo(stem.Start2 + stem.Length - 1));
        });

        // Poly-A has no complementary arms.
        var none = AnalysisTools.FindRnaInvertedRepeats("AAAAAAAAAA", 4, 3).Items;
        Assert.That(none, Is.Empty);
    }
}
