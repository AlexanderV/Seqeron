using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>pseudoknot_orders</c> MCP tool (RNA-PSEUDOKNOT-001). Expected values: biotite
/// <c>structure.pseudoknots</c> docstring example ("ABCbac": orders [0 0 1], "(([))]") and ViennaRNA 2.7.2
/// <c>RNA.db_from_ptable</c> on the same pair table.
/// </summary>
[TestFixture]
public class PseudoknotOrdersTests
{
    private static BasePairItem Bp(int p1, int p2) => new(p1, p2, 'G', 'C', "WatsonCrick");

    [Test]
    public void PseudoknotOrders_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PseudoknotOrders(Array.Empty<BasePairItem>(), 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PseudoknotOrders(new[] { Bp(0, 5), Bp(5, 8) }, 10));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PseudoknotOrders(new[] { Bp(0, 5) }, 5));
    }

    [Test]
    public void PseudoknotOrders_Binding_InvokesSuccessfully()
    {
        var r = AnalysisTools.PseudoknotOrders(new[] { Bp(0, 4), Bp(1, 3), Bp(2, 5) }, 6);
        Assert.Multiple(() =>
        {
            Assert.That(r.Orders, Is.EqualTo(new[] { 0, 0, 1 }));
            Assert.That(r.OrderDotBracket, Is.EqualTo("(([))]"));
            Assert.That(r.DotBracket, Is.EqualTo("(([))]"));
        });

        // Larger 3' stem kept at order 0 (maximum nested subset); ViennaRNA's 5'→3' passes open the 5' pair first.
        var u = AnalysisTools.PseudoknotOrders(new[] { Bp(0, 6), Bp(2, 10), Bp(3, 9), Bp(4, 8) }, 12);
        Assert.Multiple(() =>
        {
            Assert.That(u.Orders, Is.EqualTo(new[] { 1, 0, 0, 0 }));
            Assert.That(u.OrderDotBracket, Is.EqualTo("[.(((.].)))."));
            Assert.That(u.DotBracket, Is.EqualTo("(.[[[.).]]]."));
        });
    }
}
