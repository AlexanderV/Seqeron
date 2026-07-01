using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;
using SourceCa = Seqeron.Genomics.Chromosome.ChromosomeAnalyzer;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>detect_rearrangements</c>. ChromosomeAnalyzer.DetectRearrangements inspects synteny
/// blocks for inversions/translocations/etc. A single forward collinear block has no rearrangement.
/// </summary>
[TestFixture]
public class DetectRearrangementsTests
{
    private static List<SourceCa.SyntenyBlock> OneForwardBlock()
    {
        var pairs = new List<OrthologPair>();
        for (int i = 0; i < 5; i++)
            pairs.Add(new OrthologPair("chr1", i * 1_000_000, i * 1_000_000 + 1000, "g" + i,
                                        "chrA", i * 1_000_000, i * 1_000_000 + 1000, "h" + i));
        return ChromosomeTools.FindSyntenyBlocks(pairs, 3, 10).Items.ToList();
    }

    [Test]
    public void DetectRearrangements_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.DetectRearrangements(OneForwardBlock()));
        Assert.DoesNotThrow(() => ChromosomeTools.DetectRearrangements(new List<SourceCa.SyntenyBlock>()));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.DetectRearrangements(null!));
    }

    [Test]
    public void DetectRearrangements_SingleForwardBlock_NoRearrangements()
    {
        var result = ChromosomeTools.DetectRearrangements(OneForwardBlock());
        Assert.That(result.Items, Is.Empty);
    }
}
