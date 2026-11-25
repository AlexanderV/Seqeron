// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>
/// B09 A2-36: the MCP tools of the public stage ports (ClustalW pairwise stages, T-Coffee pair methods / primary library / guide tree / k-tuple similarity, sum of pairs).
/// Every expected value is copied from the reference-locked library tests (clustalw 2.1 / t_coffee 13.41 hand runs and oracles) or computed by hand.
/// </summary>
[TestFixture]
public class PortStageToolsTests
{
    private const string HbA = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF";
    private const string HbC = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHG";
    private const string HbB = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKV";

    [Test]
    public void ClustalWPairwiseAlign_Protein_ReproducesTheClustalWScores()
    {
        // SequenceAligner_ClustalWProtein_Tests: clustalw 2.1 prints 40 / 38 / 31 for these globin fragments (default GONNET pairwise stage)
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentTools.ClustalWPairwiseAlign(HbA, HbB, "protein").ClustalWScore, Is.EqualTo(40));
            Assert.That(MultipleAlignmentTools.ClustalWPairwiseAlign(HbA, HbC, "protein").ClustalWScore, Is.EqualTo(38));
            Assert.That(MultipleAlignmentTools.ClustalWPairwiseAlign(HbB, HbC, "protein").ClustalWScore, Is.EqualTo(31));
        });
        var r = MultipleAlignmentTools.ClustalWPairwiseAlign(HbA, HbB, "protein");
        Assert.That(r.AlignedSequence1, Has.Length.EqualTo(r.AlignedSequence2.Length));
        Assert.That(r.Distance, Is.EqualTo(1.0 - (double)r.Identities / Math.Min(HbA.Length, HbB.Length)).Within(1e-12));
    }

    [Test]
    public void ClustalWPairwiseAlign_Dna_DelegatesToTheLibrary()
    {
        const string a = "ACGTACGTTTGACCAGTAGGCATCGAT", b = "ACGTACGATTGACCAGTAGCCATCGAT";
        var tool = MultipleAlignmentTools.ClustalWPairwiseAlign(a, b);
        var lib = global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWPairwiseAlign(a, b);
        Assert.That((tool.Score, tool.Identities, tool.ClustalWScore, tool.AlignedSequence1, tool.AlignedSequence2, tool.Start1, tool.End2),
            Is.EqualTo((lib.Score, lib.Identities, lib.ClustalWScore, lib.AlignedSequence1, lib.AlignedSequence2, lib.Start1, lib.End2)));
        Assert.That(MultipleAlignmentTools.ClustalWPairwiseAlign(a, b, dnaMatrix: "clustalw").Score,
            Is.EqualTo(global::Seqeron.Genomics.ClustalW.ClustalWAligner.ClustalWPairwiseAlign(a, b,
                new global::Seqeron.Genomics.ClustalW.ClustalWDnaPairwiseOptions { Matrix = global::Seqeron.Genomics.ClustalW.ClustalWDnaMatrix.Clustalw }).Score));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWPairwiseAlign(a, b, "rna"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWPairwiseAlign(a, b, "protein", pwMatrix: "x"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWPairwiseAlign(null!, b));
    }

    [Test]
    public void ClustalWQuickPairwiseAlign_ReproducesTheHandRun()
    {
        // `clustalw -QUICKTREE -ALIGN` prints 92.5926 / 27.7778 / 27.7778 (SequenceAligner_ClustalWQuick_Tests)
        const string a = "ACGTACGTTTGACCAGTAGGCATCGAT", b = "ACGTACGATTGACCAGTAGCCATCGAT", c = "GGGGGGAAAATTTTCCCC";
        Assert.Multiple(() =>
        {
            Assert.That(MultipleAlignmentTools.ClustalWQuickPairwiseAlign(a, b, maxSequenceLength: 27).Score, Is.EqualTo(92.5926).Within(5e-5));
            Assert.That(MultipleAlignmentTools.ClustalWQuickPairwiseAlign(a, c, maxSequenceLength: 27).Score, Is.EqualTo(27.7778).Within(5e-5));
            Assert.That(MultipleAlignmentTools.ClustalWQuickPairwiseAlign(b, c, maxSequenceLength: 27).Distance, Is.EqualTo((100.0 - 27.7778) / 100.0).Within(1e-6));
        });
        var absolute = MultipleAlignmentTools.ClustalWQuickPairwiseAlign(a, b, percent: false, maxSequenceLength: 27);
        Assert.That(absolute.Score, Is.EqualTo(absolute.FragmentScore));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ClustalWQuickPairwiseAlign(a, b, "rna"));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentTools.ClustalWQuickPairwiseAlign(a, b, ktuple: 0));
    }

    [Test]
    public void TCoffeeSlowPairAlign_PlacesTheGapInside_WorkedExample()
    {
        var r = MultipleAlignmentTools.TCoffeeSlowPairAlign("TTTCAACTC", "TTTTCAACTC");
        Assert.That((r.RowA, r.RowB), Is.EqualTo(("T-TTCAACTC", "TTTTCAACTC")));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeSlowPairAlign(null!, "ACGT"));
    }

    [Test]
    public void TCoffeeLalignTopAlignments_ScoresInOrder()
    {
        var r = MultipleAlignmentTools.TCoffeeLalignTopAlignments("TTTTCGAACTTC", "TTTTCAACTC", 4);
        Assert.That(r.Alignments.Select(x => x.Score), Is.EqualTo(new[] { 790, 500, 400, 370 }));
        Assert.That((r.Alignments[0].RowA, r.Alignments[0].RowB), Is.EqualTo(("TTTTCGAACT", "TTTTC-AACT")));
        Assert.That(r.Count, Is.EqualTo(4));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentTools.TCoffeeLalignTopAlignments("ACGT", "ACGT", 0));
    }

    [Test]
    public void TCoffeeKtupSimilarity_And_GuideTree()
    {
        var s = MultipleAlignmentTools.TCoffeeKtupSimilarity(["AAAAAAA", "AAAAAAA", "CCCCCCC"]).Similarity;
        Assert.That((s[0][1], s[0][2], s[2][0], s[1][1]), Is.EqualTo((100, 70, 70, 100)));

        var tree = MultipleAlignmentTools.TCoffeeGuideTree(["ACGTACGTAC", "ACGTACGTAA", "TTTTGGGGCC", "TTTTGGGGCA"]);
        Assert.That(tree.LeafOrder.OrderBy(x => x), Is.EqualTo(new[] { 0, 1, 2, 3 }));
        Assert.That(tree.Newick, Does.StartWith("(").And.EndWith(");"));
        var lib = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeGuideTree(["ACGTACGTAC", "ACGTACGTAA", "TTTTGGGGCC", "TTTTGGGGCA"]);
        Assert.That(tree.Newick, Is.EqualTo(Render(lib) + ";"));
        Assert.That(MultipleAlignmentTools.TCoffeeGuideTree(["ACGT", "ACGA"], names: ["x", "y"]).Newick, Is.EqualTo("(x,y);"));
        Assert.That(MultipleAlignmentTools.TCoffeeGuideTree(["ACGT"]).Newick, Is.EqualTo("s0;"));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeGuideTree([]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeGuideTree(["ACGT", "ACGA"], names: ["x"]));

        static string Render(global::Seqeron.Genomics.TCoffee.TCoffeeTreeNode n)
            => n.Left is null ? $"s{n.LeafIndex}" : "(" + Render(n.Left) + "," + Render(n.Right!) + ")";
    }

    [Test]
    public void TCoffeePrimaryLibrary_DelegatesAndOrdersEntries()
    {
        string[] seqs = ["TTTTCAACTC", "TTTCAACTC", "TTTTCGAACTTC"];
        var r = MultipleAlignmentTools.TCoffeePrimaryLibrary(seqs, ["slow_pair", "lalign_id_pair"]);
        var lib = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeePrimaryLibraryEntries(seqs,
            [global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.SlowPair, global::Seqeron.Genomics.TCoffee.TCoffeePairMethod.LalignIdPair]);
        Assert.That(r.Count, Is.EqualTo(lib.Count).And.GreaterThan(0));
        Assert.That(r.Entries.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count)), Is.EqualTo(lib.Select(e => (e.SeqA, e.PosA, e.SeqB, e.PosB, e.Weight, e.Count))));
        Assert.That(r.Entries.All(e => e.SeqA < e.SeqB), Is.True);
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeePrimaryLibrary(seqs, ["nope"]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeePrimaryLibrary([]));
    }}
