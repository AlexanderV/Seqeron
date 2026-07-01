using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class PairwiseFstTests
{
    private static PopulationItem P(string id, double freq) =>
        new(id, new[] { new AlleleItem(freq, 100) });

    [Test]
    public void PairwiseFst_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.PairwiseFst(
            new[] { P("Pop1", 0.5), P("Pop2", 0.6) }));
    }

    [Test]
    public void PairwiseFst_Binding_InvokesSuccessfully()
    {
        // Three single-locus populations 0.5, 0.6, 0.9:
        //   [0,1] = 1/99, [0,2] = 4/21, [1,2] = 3/25; diagonal 0; symmetric.
        var result = PopulationTools.PairwiseFst(
            new[] { P("Pop1", 0.5), P("Pop2", 0.6), P("Pop3", 0.9) });

        Assert.Multiple(() =>
        {
            Assert.That(result.PopulationIds, Is.EqualTo(new[] { "Pop1", "Pop2", "Pop3" }));

            Assert.That(result.Matrix[0][0], Is.EqualTo(0.0));
            Assert.That(result.Matrix[1][1], Is.EqualTo(0.0));
            Assert.That(result.Matrix[2][2], Is.EqualTo(0.0));

            Assert.That(result.Matrix[0][1], Is.EqualTo(1.0 / 99.0).Within(1e-10));
            Assert.That(result.Matrix[0][2], Is.EqualTo(4.0 / 21.0).Within(1e-10));
            Assert.That(result.Matrix[1][2], Is.EqualTo(3.0 / 25.0).Within(1e-10));

            // Symmetry.
            Assert.That(result.Matrix[1][0], Is.EqualTo(result.Matrix[0][1]).Within(1e-10));
            Assert.That(result.Matrix[2][0], Is.EqualTo(result.Matrix[0][2]).Within(1e-10));
            Assert.That(result.Matrix[2][1], Is.EqualTo(result.Matrix[1][2]).Within(1e-10));
        });
    }

    [Test]
    public void PairwiseFst_Binding_IdenticalPopulations_ZeroMatrix()
    {
        var result = PopulationTools.PairwiseFst(new[] { P("A", 0.4), P("B", 0.4) });
        Assert.Multiple(() =>
        {
            Assert.That(result.Matrix[0][1], Is.EqualTo(0.0).Within(1e-10));
            Assert.That(result.Matrix[1][0], Is.EqualTo(0.0).Within(1e-10));
        });
    }
}
