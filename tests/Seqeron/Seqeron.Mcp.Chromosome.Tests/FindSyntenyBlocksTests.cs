using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>find_synteny_blocks</c>. ChromosomeAnalyzer.FindSyntenyBlocks groups ortholog pairs by
/// chromosome pair and reports collinear runs of >= minGenes as one block per run.
/// </summary>
[TestFixture]
public class FindSyntenyBlocksTests
{
    private static List<OrthologPair> Collinear(int n)
    {
        var list = new List<OrthologPair>();
        for (int i = 0; i < n; i++)
            list.Add(new OrthologPair("chr1", i * 1_000_000, i * 1_000_000 + 1000, "g" + i,
                                       "chrA", i * 1_000_000, i * 1_000_000 + 1000, "h" + i));
        return list;
    }

    [Test]
    public void FindSyntenyBlocks_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.FindSyntenyBlocks(Collinear(5)));
        // Fewer pairs than minGenes -> no blocks (documented), not an error.
        Assert.DoesNotThrow(() => ChromosomeTools.FindSyntenyBlocks(Collinear(1)));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenyBlocks(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.FindSyntenyBlocks(Collinear(5), minGenes: 0));
    }

    [Test]
    public void FindSyntenyBlocks_Binding_InvokesSuccessfully()
    {
        // Five perfectly collinear same-strand pairs -> one forward block of 5 genes.
        var result = ChromosomeTools.FindSyntenyBlocks(Collinear(5), 3, 10);

        Assert.That(result.Items, Has.Count.EqualTo(1));
        var b = result.Items[0];
        Assert.Multiple(() =>
        {
            Assert.That(b.Species1Chromosome, Is.EqualTo("chr1"));
            Assert.That(b.Species2Chromosome, Is.EqualTo("chrA"));
            Assert.That(b.GeneCount, Is.EqualTo(5));
            Assert.That(b.Strand, Is.EqualTo('+'));
        });
    }
}
