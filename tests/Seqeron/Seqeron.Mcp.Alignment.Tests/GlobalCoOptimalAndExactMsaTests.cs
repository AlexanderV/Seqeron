using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP <c>global_align_co_optimal</c> (B09 F59, audit A1-7b) and <c>multiple_align_exact</c> (B09 F60, audit A1-11). Expected values: Biopython 1.88
/// <c>PairwiseAligner(mode='global').align</c> set and exhaustive enumeration (core tests, scripts/alignment/generate_global_cooptimal_oracle.py);
/// exact optimum of two independent Python programs and exhaustive enumeration (scripts/alignment/generate_msa_exact_oracle.py).
/// </summary>
[TestFixture]
public class GlobalCoOptimalAndExactMsaTests
{
    [Test]
    public void GlobalCoOptimal_GapPlacements_FirstIsGlobalAlign_AndCap()
    {
        // AAAAA vs A (+1/-4, gap -1): one match, four gaps, the A of sequence 2 against any of the five A's: 5 alignments, score 1 - 4 = -3
        var all = AlignmentTools.GlobalAlignCoOptimal("AAAAA", "A", maxAlignments: 100, match: 1, mismatch: -4, gapExtend: -1);
        Assert.That(all.Count, Is.EqualTo(5));
        Assert.That(all.Alignments.Select(a => a.AlignedSequence2).ToHashSet(), Is.EquivalentTo(new[] { "----A", "---A-", "--A--", "-A---", "A----" }));
        Assert.That(all.Alignments.Select(a => a.Score), Is.All.EqualTo(-3));
        Assert.That(all.Alignments[0], Is.EqualTo(AlignmentTools.GlobalAlign("AAAAA", "A", 1, -4, -2, -1)));
        Assert.That(AlignmentTools.GlobalAlignCoOptimal("AAAAA", "A", maxAlignments: 2, mismatch: -4).Count, Is.EqualTo(2));
    }

    [Test]
    public void GlobalCoOptimal_Affine_EndGaps_AndMatrix()
    {
        var affine = AlignmentTools.GlobalAlignCoOptimal("GAAC", "GAAAC", maxAlignments: 100, affine: true, match: 5, mismatch: -3, gapOpen: -2, gapExtend: -1);
        Assert.That(affine.Count, Is.EqualTo(3));
        Assert.That(affine.Alignments.Select(a => a.Score), Is.All.EqualTo(17));   // 4 x 5 - (gapOpen + gapExtend)
        Assert.That(affine.Alignments[0], Is.EqualTo(AlignmentTools.GlobalAlignAffine("GAAC", "GAAAC", 5, -3, -2, -1)));

        // free end gaps (fitting): the shorter query may sit anywhere inside the reference; default is None = global
        var fit = AlignmentTools.GlobalAlignCoOptimal("AC", "ACAC", maxAlignments: 100, match: 1, mismatch: -1, gapExtend: -1, endGaps: "Fitting");
        Assert.That(fit.Count, Is.EqualTo(2));
        Assert.That(fit.Alignments.Select(a => a.AlignedSequence1), Is.EquivalentTo(new[] { "AC--", "--AC" }));
        Assert.That(fit.Alignments, Has.All.Matches<AlignmentResultDto>(a => a.AlignmentType == "SemiGlobal" && a.Score == 2));
        Assert.That(AlignmentTools.GlobalAlignCoOptimal("AC", "ACAC", maxAlignments: 100, endGaps: "None").Alignments[0].AlignmentType, Is.EqualTo("Global"));

        // protein, BLOSUM62, fractional gaps
        var p = AlignmentTools.GlobalAlignCoOptimal("MKVLA", "MKVA", maxAlignments: 100, affine: true, gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62");
        Assert.That(p.Count, Is.GreaterThanOrEqualTo(1));
        Assert.That(p.Alignments[0].ExactScore, Is.EqualTo(AlignmentTools.GlobalAlignAffine("MKVLA", "MKVA", gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62").ExactScore));
    }

    [Test]
    public void GlobalCoOptimal_Validation()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.GlobalAlignCoOptimal("ACGT", "ACGT", maxAlignments: 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignCoOptimal("ACGT", "ACGT", endGaps: "Sideways"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlignCoOptimal("ACGT", "ACGT", gapExtend: -0.5));   // fractional gap without a matrix
    }

    [Test]
    public void MultipleAlignExact_TwoSequences_AreGlobalAlign_AndHandExample()
    {
        var two = AlignmentTools.MultipleAlignExact(["TTTTTTTTTTACGT", "ACGT"], 2, -3, -5, -2, "affine");
        var g = AlignmentTools.GlobalAlignAffine("TTTTTTTTTTACGT", "ACGT", 2, -3, -5, -2);
        Assert.That(two.AlignedSequences, Is.EqualTo(new[] { g.AlignedSequence1, g.AlignedSequence2 }));
        Assert.That(two.TotalScore, Is.EqualTo(-17));

        // three identical 8-mers: 3 pairs x 8 columns x match
        Assert.That(AlignmentTools.MultipleAlignExact(["ACGTACGT", "ACGTACGT", "ACGTACGT"]).TotalScore, Is.EqualTo(24));
        Assert.That(AlignmentTools.MultipleAlignExact([]).AlignedSequences, Is.Empty);
    }

    [Test]
    public void MultipleAlignExact_Guard_AndValidation()
    {
        string[] seqs = ["ACGTACGTACGTACGT", "TTGACCATTGACCATT", "GGCATTACGGCATTAC"];
        var ex = Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignExact(seqs, carrilloLipman: false, maxStateCells: 1000));
        Assert.That(ex!.Message, Does.Contain("complexity limit"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignExact(["ACGT", "ACXT", "ACGT"]));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignExact(["ACGT", "ACGT", "ACGT"], gapModel: "convex"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignExact(Enumerable.Repeat("ACGT", 9).ToArray()));
    }
}
