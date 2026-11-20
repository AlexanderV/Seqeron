using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP <c>local_align_top_n</c> / <c>local_align_co_optimal</c> (B09 F47, audit A1-7): Waterman–Eggert non-intersecting local alignments and
/// all co-optimal local alignments. Expected values: independent full-recompute Waterman–Eggert reference + exhaustive enumeration
/// (scripts/alignment/generate_topn_oracle.py), EMBOSS matcher 6.6.0 -alternatives, Biopython 1.88 co-optimal set (see the core tests).
/// ACGTAC x3 in sequence 1 against one ACGTAC: three copies of score 12 (+2/-3).
/// </summary>
[TestFixture]
public class TopNAlignTests
{
    private const string Repeats = "ACGTACTTTACGTACGGGACGTAC";

    [Test]
    public void TopN_ReportsEveryCopyOfTheRepeat_AndFirstIsLocalAlign()
    {
        var r = AlignmentTools.LocalAlignTopN(Repeats, "ACGTAC", maxAlignments: 5, match: 2, mismatch: -3, gapExtend: -2);
        Assert.That(r.Count, Is.EqualTo(r.Alignments.Length));
        Assert.That(r.Alignments.Take(3).Select(a => a.Score), Is.EqualTo(new[] { 12, 12, 12 }));
        Assert.That(r.Alignments.Take(3).Select(a => a.StartPosition1), Is.EqualTo(new[] { 0, 9, 18 }));
        Assert.That(r.Alignments[0], Is.EqualTo(AlignmentTools.LocalAlign(Repeats, "ACGTAC", 2, -3, -2, -2)));
        Assert.That(r.Alignments.Select(a => a.Score), Is.Ordered.Descending);
    }

    [Test]
    public void TopN_MinScoreAndAffine_AndMatrix()
    {
        var strong = AlignmentTools.LocalAlignTopN(Repeats, "ACGTAC", maxAlignments: 50, minScore: 12, match: 2, mismatch: -3, gapExtend: -2);
        Assert.That(strong.Count, Is.EqualTo(3));
        var affine = AlignmentTools.LocalAlignTopN(Repeats, "ACGTAC", maxAlignments: 3, affine: true, match: 2, mismatch: -3, gapOpen: -5, gapExtend: -2);
        Assert.That(affine.Alignments[0], Is.EqualTo(AlignmentTools.LocalAlignAffine(Repeats, "ACGTAC", 2, -3, -5, -2)));
        Assert.That(affine.Count, Is.EqualTo(3));

        // protein repeat, BLOSUM62 with fractional gaps: the second copy scores like the first
        var p = AlignmentTools.LocalAlignTopN("MKVLAAGIVGLLLAQPSMKVLAAGIVGLLLA", "MKVLAAGIVGLLL", maxAlignments: 2, affine: true, gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62");
        Assert.That(p.Count, Is.EqualTo(2));
        Assert.That(p.Alignments[1].ExactScore, Is.EqualTo(p.Alignments[0].ExactScore));
        Assert.That(p.Alignments[0].StartPosition1, Is.EqualTo(0));
        Assert.That(p.Alignments[1].StartPosition1, Is.EqualTo(17));
    }

    [Test]
    public void TopN_TieOrderSim_ReproducesEmbossMatcher_WhereWatermanEggertOrderDiffers()
    {
        // topn_oracle.tsv line 46 (B09 F145, audit A2-28d): matcher -alternatives 6 -gapopen 8 -gapextend 2 on EDNAFULL reports 55,38,21,16,13,11.
        const string a = "GCTCTAAATTTTGTAGCTCTAATT", b = "TTCGCTCTAAATTTC";
        var sim = AlignmentTools.LocalAlignTopN(a, b, maxAlignments: 6, affine: true, match: 5, mismatch: -4, gapOpen: -6, gapExtend: -2, tieOrder: "Sim");
        Assert.That(sim.Alignments.Select(x => x.Score), Is.EqualTo(new[] { 55, 38, 21, 16, 13, 11 }));
        var we = AlignmentTools.LocalAlignTopN(a, b, maxAlignments: 6, affine: true, match: 5, mismatch: -4, gapOpen: -6, gapExtend: -2);
        Assert.That(we.Alignments.Select(x => x.Score), Is.Not.EqualTo(sim.Alignments.Select(x => x.Score)));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignTopN(a, b, tieOrder: "nonsense"));
    }

    [Test]
    public void TopN_NoPositivePair_EmptyList_AndValidation()
    {
        Assert.That(AlignmentTools.LocalAlignTopN("AAAA", "CCCC").Count, Is.EqualTo(0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.LocalAlignTopN("ACGT", "ACGT", maxAlignments: 0));
    }

    [Test]
    public void TopN_AboveTheCellLimit_RunsInLinearSpace_AndEqualsTheLocalAlignOptimum()
    {
        // B09 F77 (A1-7c): 5 100 x 5 100 = 26.0 M cells > 25 M: the declumper runs in linear space instead of rejecting the input.
        var rng = new Random(11);
        var a = new char[5100];
        var b = new char[5100];
        for (int i = 0; i < a.Length; i++) { a[i] = "ACGT"[rng.Next(4)]; b[i] = "ACGT"[rng.Next(4)]; }
        int[] starts1 = [100, 1500, 3000], starts2 = [2000, 200, 4000];   // three disjoint exact 300-nt copies
        for (int s = 0; s < 3; s++)
            for (int i = 0; i < 300; i++) b[starts2[s] + i] = a[starts1[s] + i];
        string sa = new(a), sb = new(b);
        var r = AlignmentTools.LocalAlignTopN(sa, sb, maxAlignments: 3, match: 2, mismatch: -3, gapExtend: -4);
        Assert.That(r.Count, Is.EqualTo(3));
        Assert.That(r.Alignments[0], Is.EqualTo(AlignmentTools.LocalAlign(sa, sb, 2, -3, -2, -4)));
        Assert.That(r.Alignments.Select(x => x.Score), Is.Ordered.Descending);
        Assert.That(r.Alignments.Take(3).All(x => x.Score >= 600), Is.True, "the three planted 300-nt copies");
    }

    [Test]
    public void CoOptimal_GapPlacements_AndCap()
    {
        // GAAC vs GAAAC (+5/-3, gap -1): one gap, three placements, 4 x 5 - 1 = 19
        var all = AlignmentTools.LocalAlignCoOptimal("GAAC", "GAAAC", maxAlignments: 100, match: 5, mismatch: -3, gapExtend: -1);
        Assert.That(all.Alignments.Select(a => a.AlignedSequence1 + "/" + a.AlignedSequence2),
            Is.EqualTo(new[] { "G-AAC/GAAAC", "GA-AC/GAAAC", "GAA-C/GAAAC" }));   // order: end cell, then diagonal > up > left from the end
        Assert.That(all.Alignments.Select(a => a.Score), Is.All.EqualTo(19));
        Assert.That(AlignmentTools.LocalAlignCoOptimal("GAAC", "GAAAC", maxAlignments: 2, match: 5, mismatch: -3, gapExtend: -1).Count, Is.EqualTo(2));
        var affine = AlignmentTools.LocalAlignCoOptimal("GAAC", "GAAAC", maxAlignments: 100, affine: true, match: 5, mismatch: -3, gapOpen: -2, gapExtend: -1);
        Assert.That(affine.Count, Is.EqualTo(3));
        Assert.That(affine.Alignments[0].Score, Is.EqualTo(17));   // 4 x 5 - (gapOpen + 1 x gapExtend = -3)
    }
}
