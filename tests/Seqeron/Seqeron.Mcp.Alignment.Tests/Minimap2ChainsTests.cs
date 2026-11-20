using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>MCP <c>minimap2_chains</c> (B09 A3-9c, F151): the seeding and chaining stage of minimap2 2.26; expected values are those printed by <c>minimap2 -k 11 -w 5 --print-seeds --print-chains</c>.</summary>
[TestFixture]
public class Minimap2ChainsTests
{
    private const string Reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCT";
    private const string Query = "CCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCTTTCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAA";

    [Test]
    public void Minimap2Chains_OneChainOf53Anchors_ScoreAsPafS1()
    {
        // minimap2 -k 11 -w 5: 53 seeds, one chain of 53 anchors, PAF cm:i:53 s1:i:172; first anchor "31 + 11 11", last "199 + 181 11"
        var r = AlignmentTools.Minimap2Chains([Reference], Query, 11, 5);
        Assert.That(r.SeedCount, Is.EqualTo(53));
        var chain = r.Chains.Single();
        Assert.That((chain.Score, chain.Strand, chain.Target, chain.Anchors.Length), Is.EqualTo((172, "forward", 0, 53)));
        Assert.That(chain.Anchors[0], Is.EqualTo(new Minimap2AnchorDto(0, "forward", 31, 11, 11)));
        Assert.That(chain.Anchors[^1], Is.EqualTo(new Minimap2AnchorDto(0, "forward", 199, 181, 11)));
    }

    [Test]
    public void Minimap2Chains_Validation_AndThresholds()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Chains([], Query));
        Assert.That(AlignmentTools.Minimap2Chains([Reference], "").Chains, Is.Empty);   // the pairwise tools accept an empty sequence
        Assert.That(AlignmentTools.Minimap2Chains([Reference], Query, 11, 5, minChainScore: 1000).Chains, Is.Empty);
        Assert.That(AlignmentTools.Minimap2Chains([Reference], Query, 11, 5, minAnchors: 60).Chains, Is.Empty);
    }
}
