using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class FilterVariantsByMafTests
{
    private static VariantItem V(string id, double af) =>
        new(id, "chr1", 100, "A", "G", af, 100);

    [Test]
    public void FilterVariantsByMaf_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.FilterVariantsByMaf(new[] { V("V1", 0.1) }));
        Assert.DoesNotThrow(() => PopulationTools.FilterVariantsByMaf(Array.Empty<VariantItem>()));
    }

    [Test]
    public void FilterVariantsByMaf_Binding_InvokesSuccessfully()
    {
        // V1 MAF = 0.005 < 0.01 → excluded; V2 MAF = 0.05 ≥ 0.01 → kept.
        var result = PopulationTools.FilterVariantsByMaf(
            new[] { V("V1", 0.005), V("V2", 0.05) },
            minMAF: 0.01);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1));
            Assert.That(result.Items[0].Id, Is.EqualTo("V2"));
        });
    }

    [Test]
    public void FilterVariantsByMaf_Binding_FoldsHighFrequencyAndPreservesOrder()
    {
        Assert.Multiple(() =>
        {
            // AF = 0.95 → MAF = min(0.95, 0.05) = 0.05, passes [0.01, 0.1].
            var high = PopulationTools.FilterVariantsByMaf(
                new[] { V("V1", 0.95) }, minMAF: 0.01, maxMAF: 0.1);
            Assert.That(high.Items, Has.Count.EqualTo(1));
            Assert.That(high.Items[0].Id, Is.EqualTo("V1"));

            // maxMAF excludes MAF = 0.45; input order preserved for survivors.
            var ordered = PopulationTools.FilterVariantsByMaf(
                new[] { V("A", 0.3), V("B", 0.45), V("C", 0.2) },
                minMAF: 0.01, maxMAF: 0.4);
            Assert.That(ordered.Items.Select(i => i.Id), Is.EqualTo(new[] { "A", "C" }));
        });
    }
}
