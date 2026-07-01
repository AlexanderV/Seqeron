using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class NucleotideDiversityTests
{
    [Test]
    public void NucleotideDiversity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.NucleotideDiversity(new[] { "AAAA", "TTTT" }));
        // Single sequence is valid input; returns 0.
        Assert.DoesNotThrow(() => PopulationTools.NucleotideDiversity(new[] { "ACGT" }));
    }

    [Test]
    public void NucleotideDiversity_Binding_InvokesSuccessfully()
    {
        // Two sequences differing at every site → π = 1.0.
        Assert.That(PopulationTools.NucleotideDiversity(new[] { "AAAA", "TTTT" }).Pi,
            Is.EqualTo(1.0).Within(1e-10));

        // Identical sequences → π = 0.
        Assert.That(PopulationTools.NucleotideDiversity(new[] { "ACGT", "ACGT", "ACGT" }).Pi,
            Is.EqualTo(0.0).Within(1e-10));

        // Wikipedia Tajima's D dataset → π = 0.1.
        Assert.That(PopulationTools.NucleotideDiversity(new[]
        {
            "00000000000000000000",
            "00100000000010000010",
            "00000000000010000010",
            "00000010000000000010",
            "00000010000010000010",
        }).Pi, Is.EqualTo(0.1).Within(1e-4));
    }

    [Test]
    public void NucleotideDiversity_Binding_PartialPolymorphism()
    {
        // ACGT / ACGA / ACGT: pairwise diffs (1,2)=1,(1,3)=0,(2,3)=1 → total 2;
        // π = 2 / (C(3,2)=3 × 4) = 2/12 = 1/6.
        var pi = PopulationTools.NucleotideDiversity(new[] { "ACGT", "ACGA", "ACGT" }).Pi;
        Assert.That(pi, Is.EqualTo(2.0 / 12.0).Within(1e-4));
    }
}
