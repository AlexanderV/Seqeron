using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class HaplotypeBlocksTests
{
    private static VariantGenotypesItem V(string id, int pos, params int[] g) => new(id, pos, g);

    [Test]
    public void HaplotypeBlocks_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.HaplotypeBlocks(
            new[] { V("V1", 100, 0, 0, 1, 1, 2, 2), V("V2", 200, 0, 0, 1, 1, 2, 2) }));

        // Fewer than 2 variants → no blocks (not an error).
        Assert.DoesNotThrow(() => PopulationTools.HaplotypeBlocks(Array.Empty<VariantGenotypesItem>()));
        Assert.DoesNotThrow(() => PopulationTools.HaplotypeBlocks(new[] { V("V1", 100, 0, 1, 2) }));
    }

    [Test]
    public void HaplotypeBlocks_Binding_InvokesSuccessfully()
    {
        // Three identical-genotype variants → r² = 1.0 → one block spanning 100..300.
        var geno = new[] { 0, 0, 1, 1, 2, 2 };
        var result = PopulationTools.HaplotypeBlocks(
            new[]
            {
                new VariantGenotypesItem("V1", 100, geno),
                new VariantGenotypesItem("V2", 200, geno),
                new VariantGenotypesItem("V3", 300, geno),
            },
            ldThreshold: 0.7);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1));
            Assert.That(result.Items[0].Start, Is.EqualTo(100));
            Assert.That(result.Items[0].End, Is.EqualTo(300));
            Assert.That(result.Items[0].Variants, Is.EqualTo(new[] { "V1", "V2", "V3" }));
        });
    }

    [Test]
    public void HaplotypeBlocks_Binding_TwoBlocksSeparatedByLdBreak()
    {
        var genoA = new[] { 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        var genoB = new[] { 0, 1, 2, 0, 1, 2, 0, 1, 2 }; // r²=0 vs genoA

        var result = PopulationTools.HaplotypeBlocks(
            new[]
            {
                new VariantGenotypesItem("V1", 100, genoA),
                new VariantGenotypesItem("V2", 200, genoA),
                new VariantGenotypesItem("V3", 300, genoB),
                new VariantGenotypesItem("V4", 400, genoA),
                new VariantGenotypesItem("V5", 500, genoA),
            },
            ldThreshold: 0.7);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(2));
            Assert.That(result.Items[0].Variants, Is.EqualTo(new[] { "V1", "V2" }));
            Assert.That(result.Items[1].Variants, Is.EqualTo(new[] { "V4", "V5" }));
        });
    }
}
