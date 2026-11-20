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

    [Test]
    public void Minimap2Chains_PafMode_IsThePafLineOfMinimap2()
    {
        // minimap2 -k 11 -w 5 ref query (PAF without -c): q 182 1 182 + r 300 21 200 172 181 60 tp:A:P cm:i:53 s1:i:172 s2:i:0 dv:f:0.0017 rl:i:0 (A4-7)
        var r = AlignmentTools.Minimap2Chains([Reference], Query, 11, 5, mode: "paf", queryName: "q");
        var hit = r.Hits!.Single();
        Assert.That((hit.QueryStart, hit.QueryEnd, hit.Strand, hit.TargetStart, hit.TargetEnd, hit.MatchLength, hit.BlockLength, hit.MapQ), Is.EqualTo((1, 182, "forward", 21, 200, 172, 181, 60)));
        Assert.That((hit.Primary, hit.AnchorCount, hit.Score, hit.SecondaryScore, hit.Divergence), Is.EqualTo((true, 53, 172, 0, 0.0017)));
        Assert.That((r.RepeatLength, r.MidOcc), Is.EqualTo((0, 10)));
        Assert.That(r.Chains.Single().Score, Is.EqualTo(172));
        // -H: homopolymer-compressed k-mers (q 182 5 173 + r 300 25 191 164 170 60 cm:i:34 s1:i:163 dv:f:0.0018)
        var hpc = AlignmentTools.Minimap2Chains([Reference], Query, 11, 5, mode: "paf", queryName: "q", homopolymerCompressed: true).Hits!.Single();
        Assert.That((hpc.QueryStart, hpc.QueryEnd, hpc.TargetStart, hpc.TargetEnd, hpc.AnchorCount, hpc.Score, hpc.Divergence), Is.EqualTo((5, 173, 25, 191, 34, 163, 0.0018)));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Chains([Reference], Query, mode: "sam"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.Minimap2Chains([Reference], Query, allVsAll: true));
    }

    [Test]
    public void Minimap2Chains_FrequentMinimizers_AreDropped_AndTheRepeatLengthReported()
    {
        // a 60 kb reference with 3 000 poly-A (the one frequent key) and a query of 24 A + the first 200 bases: the program drops the poly-A minimizer; rl = the dropped query bases
        var rng = new Random(5);
        string reference = new string(Enumerable.Range(0, 60_000).Select(_ => "ACGT"[rng.Next(4)]).ToArray()) + new string('A', 3000);
        var r = AlignmentTools.Minimap2Chains([reference], new string('A', 24) + reference[..200]);
        // A5-4: the exact values of `minimap2 --print-seeds ref.fa q.fa` (2.26) on these inputs: "mid_occ = 10", "RS 24"
        Assert.That((r.RepeatLength, r.MidOcc), Is.EqualTo((24, 10)));
    }
}
