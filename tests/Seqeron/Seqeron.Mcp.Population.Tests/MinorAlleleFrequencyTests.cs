using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class MinorAlleleFrequencyTests
{
    [Test]
    public void MinorAlleleFrequency_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.MinorAlleleFrequency(new[] { 0, 1, 2 }));
        // Empty genotype vector is valid input; returns 0.
        Assert.DoesNotThrow(() => PopulationTools.MinorAlleleFrequency(Array.Empty<int>()));
    }

    [Test]
    public void MinorAlleleFrequency_Binding_InvokesSuccessfully()
    {
        // Alt alleles = 4·1 + 1·2 = 6, total = 20, altFreq = 0.3 → MAF = 0.3.
        var low = PopulationTools.MinorAlleleFrequency(new[] { 0, 0, 0, 0, 0, 1, 1, 1, 1, 2 });
        Assert.That(low.Maf, Is.EqualTo(0.3).Within(1e-10));

        // Alt freq 0.7 folds to MAF 0.3.
        var high = PopulationTools.MinorAlleleFrequency(new[] { 0, 1, 1, 1, 1, 2, 2, 2, 2, 2 });
        Assert.That(high.Maf, Is.EqualTo(1.0 - 14.0 / 20.0).Within(1e-10));

        // Balanced polymorphism → MAF = 0.5.
        var balanced = PopulationTools.MinorAlleleFrequency(new[] { 0, 1, 1, 2 });
        Assert.That(balanced.Maf, Is.EqualTo(0.5).Within(1e-10));
    }

    [Test]
    public void MinorAlleleFrequency_Binding_MonomorphicAndEmpty()
    {
        Assert.Multiple(() =>
        {
            // Monomorphic reference / alternate → MAF 0.
            Assert.That(PopulationTools.MinorAlleleFrequency(new[] { 0, 0, 0, 0, 0 }).Maf, Is.EqualTo(0.0));
            Assert.That(PopulationTools.MinorAlleleFrequency(new[] { 2, 2, 2, 2, 2 }).Maf, Is.EqualTo(0.0));
            // Empty → 0.
            Assert.That(PopulationTools.MinorAlleleFrequency(Array.Empty<int>()).Maf, Is.EqualTo(0.0));
        });
    }
}
