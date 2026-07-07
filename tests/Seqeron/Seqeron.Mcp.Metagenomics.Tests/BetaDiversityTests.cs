using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.CalculateBetaDiversity.
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_BetaDiversity_Tests
// (Wikipedia Bray-Curtis aquarium example; Bray & Curtis 1957; Jaccard 1901).
[TestFixture]
public class BetaDiversityTests
{
    private static AbundanceItem[] Ab(params (string Name, double Fraction)[] items)
        => Array.ConvertAll(items, t => new AbundanceItem(t.Name, t.Fraction));

    [Test]
    public void BetaDiversity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.BetaDiversity(
            "s1", Ab(("A", 1.0)), "s2", Ab(("B", 1.0))));

        // Empty vectors are defined input (distances 0), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.BetaDiversity(
            "s1", Array.Empty<AbundanceItem>(), "s2", Array.Empty<AbundanceItem>()));
    }

    [Test]
    public void BetaDiversity_Binding_InvokesSuccessfully()
    {
        // Wikipedia Bray-Curtis aquarium example.
        // Tank1: Goldfish=6, Guppy=7, Rainbow=4; Tank2: Goldfish=10, Guppy=0, Rainbow=6.
        // BC = 1 - 2*10/33 = 13/33; Jaccard = 1 - 2/3 = 1/3.
        var result = MetagenomicsTools.BetaDiversity(
            "Tank 1", Ab(("Goldfish", 6), ("Guppy", 7), ("Rainbow", 4)),
            "Tank 2", Ab(("Goldfish", 10), ("Guppy", 0), ("Rainbow", 6)));

        Assert.Multiple(() =>
        {
            Assert.That(result.Sample1, Is.EqualTo("Tank 1"));
            Assert.That(result.Sample2, Is.EqualTo("Tank 2"));
            Assert.That(result.BrayCurtis, Is.EqualTo(13.0 / 33.0).Within(1e-10),
                "Bray-Curtis = 1 - 20/33 (Wikipedia worked example).");
            Assert.That(result.JaccardDistance, Is.EqualTo(1.0 / 3.0).Within(1e-10),
                "Jaccard distance = 1 - 2/3 = 1/3.");
            Assert.That(result.SharedSpecies, Is.EqualTo(2), "Goldfish and Rainbow shared.");
            Assert.That(result.UniqueToSample1, Is.EqualTo(1), "Guppy unique to Tank 1 (0 in Tank 2).");
            Assert.That(result.UniqueToSample2, Is.EqualTo(0));
            Assert.That(result.UniFracDistance, Is.EqualTo(0),
                "UniFrac is 0 (no phylogenetic tree input).");
        });
    }
}
