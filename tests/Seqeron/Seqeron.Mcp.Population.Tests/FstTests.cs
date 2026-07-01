using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class FstTests
{
    private static AlleleItem A(double freq, int n) => new(freq, n);

    [Test]
    public void Fst_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.Fst(
            new[] { A(0.8, 100) }, new[] { A(0.2, 100) }));

        // Mismatched per-locus counts must throw (no silent truncation).
        Assert.Throws<ArgumentException>(() => PopulationTools.Fst(
            new[] { A(0.8, 100), A(0.2, 100) }, new[] { A(0.2, 100) }));
    }

    [Test]
    public void Fst_Binding_InvokesSuccessfully()
    {
        // p1=0.8, p2=0.2, equal sizes: p̄=0.5, var=0.09, het=0.25 → Fst = 0.36.
        var single = PopulationTools.Fst(new[] { A(0.8, 100) }, new[] { A(0.2, 100) });
        Assert.That(single.Fst, Is.EqualTo(0.36).Within(1e-10));

        // Fixed differences → complete differentiation Fst = 1.0.
        var fixedDiff = PopulationTools.Fst(new[] { A(1.0, 100) }, new[] { A(0.0, 100) });
        Assert.That(fixedDiff.Fst, Is.EqualTo(1.0).Within(1e-10));

        // Identical populations → panmixia Fst = 0.
        var identical = PopulationTools.Fst(
            new[] { A(0.5, 100), A(0.3, 100) },
            new[] { A(0.5, 100), A(0.3, 100) });
        Assert.That(identical.Fst, Is.EqualTo(0.0).Within(1e-10));
    }

    [Test]
    public void Fst_Binding_MultiLocus_ExactValue()
    {
        // pop1=[(0.9,100),(0.8,100)], pop2=[(0.1,100),(0.2,100)]:
        //   Locus1 var=0.16, Locus2 var=0.09, het=0.25 each → Fst = 0.25/0.5 = 0.50.
        var fst = PopulationTools.Fst(
            new[] { A(0.9, 100), A(0.8, 100) },
            new[] { A(0.1, 100), A(0.2, 100) });

        Assert.That(fst.Fst, Is.EqualTo(0.50).Within(1e-10));
    }
}
