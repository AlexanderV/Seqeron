using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class EstimateAncestryTests
{
    private static IndividualItem Ind(string id, params int[] g) => new(id, g);
    private static RefPopItem Ref(string id, params double[] f) => new(id, f);

    [Test]
    public void EstimateAncestry_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.EstimateAncestry(
            new[] { Ind("IND1", 2, 0) },
            new[] { Ref("A", 0.8, 0.2), Ref("B", 0.2, 0.8) },
            maxIterations: 1));

        // Empty individuals or empty reference panels → empty result (not an error).
        Assert.DoesNotThrow(() => PopulationTools.EstimateAncestry(
            Array.Empty<IndividualItem>(),
            new[] { Ref("A", 0.8), Ref("B", 0.2) }));
        Assert.DoesNotThrow(() => PopulationTools.EstimateAncestry(
            new[] { Ind("IND1", 2, 0) },
            Array.Empty<RefPopItem>()));
    }

    [Test]
    public void EstimateAncestry_Binding_InvokesSuccessfully()
    {
        // One EM iteration on symmetric panel f_A=(0.8,0.2), f_B=(0.2,0.8), g=(2,0).
        // Eq. 4 (Alexander et al. 2009): q_A = 0.8, q_B = 0.2.
        var oneIter = PopulationTools.EstimateAncestry(
            new[] { Ind("IND1", 2, 0) },
            new[] { Ref("A", 0.8, 0.2), Ref("B", 0.2, 0.8) },
            maxIterations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(oneIter.Items, Has.Count.EqualTo(1));
            var p = oneIter.Items[0];
            Assert.That(p.IndividualId, Is.EqualTo("IND1"));
            Assert.That(p.Proportions["A"], Is.EqualTo(0.8).Within(1e-10));
            Assert.That(p.Proportions["B"], Is.EqualTo(0.2).Within(1e-10));
        });

        // Converged: diagnostic individual driven to its source population (q_A → 1).
        var converged = PopulationTools.EstimateAncestry(
            new[] { Ind("IND1", 2, 0) },
            new[] { Ref("A", 0.8, 0.2), Ref("B", 0.2, 0.8) },
            maxIterations: 1000);

        Assert.Multiple(() =>
        {
            Assert.That(converged.Items[0].Proportions["A"], Is.EqualTo(1.0).Within(1e-3));
            Assert.That(converged.Items[0].Proportions["B"], Is.EqualTo(0.0).Within(1e-3));
        });
    }

    [Test]
    public void EstimateAncestry_Binding_EdgeCases()
    {
        Assert.Multiple(() =>
        {
            // Mismatched genotype length → that individual is skipped.
            var mixed = PopulationTools.EstimateAncestry(
                new[] { Ind("SHORT", 2), Ind("OK", 2, 0) },
                new[] { Ref("A", 0.8, 0.2), Ref("B", 0.2, 0.8) },
                maxIterations: 1);
            Assert.That(mixed.Items, Has.Count.EqualTo(1));
            Assert.That(mixed.Items[0].IndividualId, Is.EqualTo("OK"));

            // Zero iterations → uniform prior 1/K.
            var zero = PopulationTools.EstimateAncestry(
                new[] { Ind("IND1", 2, 0) },
                new[] { Ref("A", 0.8, 0.2), Ref("B", 0.2, 0.8) },
                maxIterations: 0);
            Assert.That(zero.Items[0].Proportions["A"], Is.EqualTo(0.5).Within(1e-10));
            Assert.That(zero.Items[0].Proportions["B"], Is.EqualTo(0.5).Within(1e-10));
        });
    }
}
