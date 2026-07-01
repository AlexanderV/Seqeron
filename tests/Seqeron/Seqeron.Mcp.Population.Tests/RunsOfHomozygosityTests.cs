using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class RunsOfHomozygosityTests
{
    // n consecutive SNPs at the given spacing, homozygous (0) except the listed het indices.
    private static GenotypePositionItem[] Snps(int n, int spacing, params int[] hetIndices)
    {
        var het = new HashSet<int>(hetIndices);
        var list = new GenotypePositionItem[n];
        for (int i = 0; i < n; i++)
            list[i] = new GenotypePositionItem(i * spacing, het.Contains(i) ? 1 : 0);
        return list;
    }

    [Test]
    public void RunsOfHomozygosity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.RunsOfHomozygosity(Snps(100, 20_000)));
        Assert.DoesNotThrow(() => PopulationTools.RunsOfHomozygosity(Array.Empty<GenotypePositionItem>()));

        // minSnps < 1 must throw.
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            PopulationTools.RunsOfHomozygosity(Snps(10, 20_000), minSnps: 0));
    }

    [Test]
    public void RunsOfHomozygosity_Binding_SingleRun_ExactBounds()
    {
        // 100 homozygous SNPs @ 20 kb → span [0, 1,980,000], count 100.
        var result = PopulationTools.RunsOfHomozygosity(
            Snps(100, 20_000), minSnps: 50, minLength: 1_000_000);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1));
            Assert.That(result.Items[0].Start, Is.EqualTo(0));
            Assert.That(result.Items[0].End, Is.EqualTo(1_980_000));
            Assert.That(result.Items[0].SnpCount, Is.EqualTo(100));
        });
    }

    [Test]
    public void RunsOfHomozygosity_Binding_ShortRunFilteredAndHetSplit()
    {
        Assert.Multiple(() =>
        {
            // Fewer than minSnps → no runs.
            var shortRun = PopulationTools.RunsOfHomozygosity(Snps(20, 20_000), minSnps: 100);
            Assert.That(shortRun.Items, Is.Empty);

            // Zero tolerance, het at midpoint splits 100 SNPs into [0..49] and [51..99].
            var split = PopulationTools.RunsOfHomozygosity(
                Snps(100, 20_000, 50), minSnps: 40, minLength: 500_000, maxHeterozygotes: 0);
            Assert.That(split.Items, Has.Count.EqualTo(2));
            Assert.That(split.Items[0].SnpCount, Is.EqualTo(50));
            Assert.That(split.Items[0].End, Is.EqualTo(49 * 20_000));
            Assert.That(split.Items[1].Start, Is.EqualTo(51 * 20_000));
            Assert.That(split.Items[1].SnpCount, Is.EqualTo(49));
        });
    }
}
