using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>detect_pseudoknots</c> MCP tool.
/// Expected results derived from the crossing condition documented in
/// RnaSecondaryStructure.DetectPseudoknots: normalized/ordered pairs (i,j),(k,l) cross iff
/// i &lt; k &lt; j &lt; l. NOT from the wrapper's output. Nested/disjoint pairs must NOT be reported.
/// </summary>
[TestFixture]
public class DetectPseudoknotsTests
{
    private static BasePairItem Bp(int p1, int p2) => new(p1, p2, 'G', 'C', "WatsonCrick");

    [Test]
    public void DetectPseudoknots_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DetectPseudoknots(Array.Empty<BasePairItem>()));
        Assert.DoesNotThrow(() => AnalysisTools.DetectPseudoknots(new[] { Bp(0, 5), Bp(3, 8) }));
        // Invalid pair-type string must throw.
        Assert.Throws<ArgumentException>(() =>
            AnalysisTools.DetectPseudoknots(new[] { new BasePairItem(0, 5, 'G', 'C', "Bogus"), Bp(3, 8) }));
    }

    [Test]
    public void DetectPseudoknots_Binding_InvokesSuccessfully()
    {
        // Crossing: (0,5) and (3,8): 0 < 3 < 5 < 8 -> one pseudoknot.
        var crossing = AnalysisTools.DetectPseudoknots(new[] { Bp(0, 5), Bp(3, 8) }).Items;
        Assert.That(crossing, Has.Length.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(crossing[0].Start1, Is.EqualTo(0));
            Assert.That(crossing[0].End1, Is.EqualTo(5));
            Assert.That(crossing[0].Start2, Is.EqualTo(3));
            Assert.That(crossing[0].End2, Is.EqualTo(8));
        });

        // Nested (0,8) and (3,5) do NOT cross.
        var nested = AnalysisTools.DetectPseudoknots(new[] { Bp(0, 8), Bp(3, 5) }).Items;
        Assert.That(nested, Is.Empty);

        // Disjoint (0,3) and (5,8) do NOT cross.
        var disjoint = AnalysisTools.DetectPseudoknots(new[] { Bp(0, 3), Bp(5, 8) }).Items;
        Assert.That(disjoint, Is.Empty);
    }
}
