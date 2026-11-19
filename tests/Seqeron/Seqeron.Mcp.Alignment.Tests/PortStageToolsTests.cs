using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// B09 A2-36: the MCP tools of the public stage ports (ClustalW pairwise stages, T-Coffee pair methods / primary library / guide tree / k-tuple similarity, sum of pairs).
/// Every expected value is copied from the reference-locked library tests (clustalw 2.1 / t_coffee 13.41 hand runs and oracles) or computed by hand.
/// </summary>
[TestFixture]
public class PortStageToolsTests
{
    private const string HbA = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF";
    private const string HbB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV";
    private const string HbC = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHG";

    [Test]
    public void ClustalWPairwiseAlign_Protein_ReproducesTheClustalWScores()
    {
        // SequenceAligner_ClustalWProtein_Tests: clustalw 2.1 prints 40 / 38 / 31 for these globin fragments (default GONNET pairwise stage)
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.ClustalWPairwiseAlign(HbA, HbB, "protein").ClustalWScore, Is.EqualTo(40));
            Assert.That(AlignmentTools.ClustalWPairwiseAlign(HbA, HbC, "protein").ClustalWScore, Is.EqualTo(38));
            Assert.That(AlignmentTools.ClustalWPairwiseAlign(HbB, HbC, "protein").ClustalWScore, Is.EqualTo(31));
        });
        var r = AlignmentTools.ClustalWPairwiseAlign(HbA, HbB, "protein");
        Assert.That(r.AlignedSequence1, Has.Length.EqualTo(r.AlignedSequence2.Length));
        Assert.That(r.Distance, Is.EqualTo(1.0 - (double)r.Identities / Math.Min(HbA.Length, HbB.Length)).Within(1e-12));
    }

    [Test]
    public void ClustalWPairwiseAlign_Dna_DelegatesToTheLibrary()
    {
        const string a = "ACGTACGTTTGACCAGTAGGCATCGAT", b = "ACGTACGATTGACCAGTAGCCATCGAT";
        var tool = AlignmentTools.ClustalWPairwiseAlign(a, b);
        var lib = global::Seqeron.Genomics.Alignment.SequenceAligner.ClustalWPairwiseAlign(a, b);
        Assert.That((tool.Score, tool.Identities, tool.ClustalWScore, tool.AlignedSequence1, tool.AlignedSequence2, tool.Start1, tool.End2),
            Is.EqualTo((lib.Score, lib.Identities, lib.ClustalWScore, lib.AlignedSequence1, lib.AlignedSequence2, lib.Start1, lib.End2)));
        Assert.That(AlignmentTools.ClustalWPairwiseAlign(a, b, dnaMatrix: "clustalw").Score,
            Is.EqualTo(global::Seqeron.Genomics.Alignment.SequenceAligner.ClustalWPairwiseAlign(a, b,
                new global::Seqeron.Genomics.Infrastructure.ClustalWDnaPairwiseOptions { Matrix = global::Seqeron.Genomics.Infrastructure.ClustalWDnaMatrix.Clustalw }).Score));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWPairwiseAlign(a, b, "rna"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWPairwiseAlign(a, b, "protein", pwMatrix: "x"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWPairwiseAlign(null!, b));
    }

    [Test]
    public void ClustalWQuickPairwiseAlign_ReproducesTheHandRun()
    {
        // `clustalw -QUICKTREE -ALIGN` prints 92.5926 / 27.7778 / 27.7778 (SequenceAligner_ClustalWQuick_Tests)
        const string a = "ACGTACGTTTGACCAGTAGGCATCGAT", b = "ACGTACGATTGACCAGTAGCCATCGAT", c = "GGGGGGAAAATTTTCCCC";
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.ClustalWQuickPairwiseAlign(a, b, maxSequenceLength: 27).Score, Is.EqualTo(92.5926).Within(5e-5));
            Assert.That(AlignmentTools.ClustalWQuickPairwiseAlign(a, c, maxSequenceLength: 27).Score, Is.EqualTo(27.7778).Within(5e-5));
            Assert.That(AlignmentTools.ClustalWQuickPairwiseAlign(b, c, maxSequenceLength: 27).Distance, Is.EqualTo((100.0 - 27.7778) / 100.0).Within(1e-6));
        });
        var absolute = AlignmentTools.ClustalWQuickPairwiseAlign(a, b, percent: false, maxSequenceLength: 27);
        Assert.That(absolute.Score, Is.EqualTo(absolute.FragmentScore));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ClustalWQuickPairwiseAlign(a, b, "rna"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ClustalWQuickPairwiseAlign(a, b, ktuple: 0));
    }

    [Test]
    public void TCoffeeSlowPairAlign_PlacesTheGapInside_WorkedExample()
    {
        var r = AlignmentTools.TCoffeeSlowPairAlign("TTTCAACTC", "TTTTCAACTC");
        Assert.That((r.RowA, r.RowB), Is.EqualTo(("T-TTCAACTC", "TTTTCAACTC")));
        Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeeSlowPairAlign(null!, "ACGT"));
    }

    [Test]
    public void TCoffeeLalignTopAlignments_ScoresInOrder()
    {
        var r = AlignmentTools.TCoffeeLalignTopAlignments("TTTTCGAACTTC", "TTTTCAACTC", 4);
        Assert.That(r.Alignments.Select(x => x.Score), Is.EqualTo(new[] { 790, 500, 400, 370 }));
        Assert.That((r.Alignments[0].RowA, r.Alignments[0].RowB), Is.EqualTo(("TTTTCGAACT", "TTTTC-AACT")));
        Assert.That(r.Count, Is.EqualTo(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.TCoffeeLalignTopAlignments("ACGT", "ACGT", 0));
    }

    [Test]
    public void TCoffeeKtupSimilarity_And_GuideTree()
    {
        var s = AlignmentTools.TCoffeeKtupSimilarity(["AAAAAAA", "AAAAAAA", "CCCCCCC"]).Similarity;
        Assert.That((s[0][1], s[0][2], s[2][0], s[1][1]), Is.EqualTo((100, 70, 70, 100)));

        var tree = AlignmentTools.TCoffeeGuideTree(["ACGTACGTAC", "ACGTACGTAA", "TTTTGGGGCC", "TTTTGGGGCA"]);
        Assert.That(tree.LeafOrder.OrderBy(x => x), Is.EqualTo(new[] { 0, 1, 2, 3 }));
        Assert.That(tree.Newick, Does.StartWith("(").And.EndWith(");"));
        var lib = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeeGuideTree(["ACGTACGTAC", "ACGTACGTAA", "TTTTGGGGCC", "TTTTGGGGCA"]);
        Assert.That(tree.Newick, Is.EqualTo(Render(lib) + ";"));
        Assert.That(AlignmentTools.TCoffeeGuideTree(["ACGT", "ACGA"], names: ["x", "y"]).Newick, Is.EqualTo("(x,y);"));
        Assert.That(AlignmentTools.TCoffeeGuideTree(["ACGT"]).Newick, Is.EqualTo("s0;"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeeGuideTree([]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeeGuideTree(["ACGT", "ACGA"], names: ["x"]));

        static string Render(global::Seqeron.Genomics.Alignment.TCoffeeTreeNode n)
            => n.Left is null ? $"s{n.LeafIndex}" : "(" + Render(n.Left) + "," + Render(n.Right!) + ")";
    }

    [Test]
    public void TCoffeePrimaryLibrary_DelegatesAndOrdersEntries()
    {
        string[] seqs = ["TTTTCAACTC", "TTTCAACTC", "TTTTCGAACTTC"];
        var r = AlignmentTools.TCoffeePrimaryLibrary(seqs, ["slow_pair", "lalign_id_pair"]);
        var lib = global::Seqeron.Genomics.Alignment.SequenceAligner.TCoffeePrimaryLibraryEntries(seqs,
            [global::Seqeron.Genomics.Alignment.TCoffeePairMethod.SlowPair, global::Seqeron.Genomics.Alignment.TCoffeePairMethod.LalignIdPair]);
        Assert.That(r.Count, Is.EqualTo(lib.Count).And.GreaterThan(0));
        Assert.That(r.Entries.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count)), Is.EqualTo(lib.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))));
        Assert.That(r.Entries.All(e => e.SeqA < e.SeqB), Is.True);
        Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeePrimaryLibrary(seqs, ["nope"]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.TCoffeePrimaryLibrary([]));
    }

    [Test]
    public void SumOfPairsScore_LinearAffineAndMatrix()
    {
        // 3 identical rows x 2 columns: 3 pairs x 2 columns x match 1 = 6 (linear); BLOSUM62 A/A = 4 -> 24 (SequenceAligner_ScoringContract_Tests)
        string[] rows = ["AA", "AA", "AA"];
        Assert.That(AlignmentTools.SumOfPairsScore(rows).Score, Is.EqualTo(6));
        Assert.That(AlignmentTools.SumOfPairsScore(rows, matrix: "BLOSUM62", gapOpen: -1, gapExtend: -1).Score, Is.EqualTo(24.0));
        // AC- / AC- / ACA with match 2, mismatch -1, gap -2 per position (linear): pairs (0,1) 2+2+0(gap-gap) = 4; (0,2) 2+2-2 = 2; (1,2) 2; total 8
        Assert.That(AlignmentTools.SumOfPairsScore(["AC-", "AC-", "ACA"], match: 2, mismatch: -1, gapOpen: -3, gapExtend: -2).Score, Is.EqualTo(8));
        // affine: pair (0,2): the single gap column costs gapOpen + 1 x gapExtend = -3 + -2 = -5 -> 4 + (4 - 5) + (4 - 5) = 2
        var affine = AlignmentTools.SumOfPairsScore(["AC-", "AC-", "ACA"], match: 2, mismatch: -1, gapOpen: -3, gapExtend: -2, gapModel: "affine");
        Assert.That(affine.Score, Is.EqualTo(2));
        Assert.That(affine.GapModel, Is.EqualTo("Affine"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(["AA", "A"]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(rows, gapOpen: -0.5));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore(rows, gapModel: "x"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SumOfPairsScore([]));
    }
}
