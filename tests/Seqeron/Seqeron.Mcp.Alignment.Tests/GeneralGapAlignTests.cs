using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// MCP <c>general_gap_align</c> (B09 F45, audit A1-5): Waterman–Smith–Beyer alignment with a gap-score table. Expected values:
/// Biopython 1.88 <c>PairwiseAligner(gap_score=callable)</c> (match 2, mismatch -3): AAAATTTT / AAAACCCCTTTT with
/// w(k) = -(2 + 3 ln k) -> 9.841116916640328 (one run of four); w(k) = -k^2 -> 5.
/// </summary>
[TestFixture]
public class GeneralGapAlignTests
{
    private static double[] Table(Func<int, double> w, int n = 12) => Enumerable.Range(1, n).Select(w).ToArray();

    [Test]
    public void Global_ConcaveLogGap_OneLongRun_BiopythonValue()
    {
        var r = AlignmentTools.GeneralGapAlign("AAAATTTT", "AAAACCCCTTTT", Table(k => -(2 + 3 * Math.Log(k))), "global", 2, -3);
        Assert.Multiple(() =>
        {
            Assert.That(r.ExactScore, Is.EqualTo(9.841116916640328).Within(1e-12));
            Assert.That(r.Score, Is.EqualTo(10));
            Assert.That(r.AlignedSequence1, Is.EqualTo("AAAA----TTTT"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("AAAACCCCTTTT"));
            Assert.That(r.AlignmentType, Is.EqualTo("Global"));
        });
    }

    [Test]
    public void Global_ConvexGap_BiopythonValue5()
    {
        var r = AlignmentTools.GeneralGapAlign("AAAATTTT", "AAAACCCCTTTT", Table(k => -(double)k * k), "global", 2, -3);
        Assert.That(r.ExactScore, Is.EqualTo(5));
    }

    [Test]
    public void AffineTable_EqualsAffineTool()
    {
        var general = AlignmentTools.GeneralGapAlign("ACGTACGTTTACGT", "ACGTACGTACGTAA", Table(k => -3.0 - k, 14), "global", 2, -3);
        var affine = AlignmentTools.GlobalAlignAffine("ACGTACGTTTACGT", "ACGTACGTACGTAA", 2, -3, -3, -1);
        Assert.That(general, Is.EqualTo(affine));
        var local = AlignmentTools.GeneralGapAlign("ACGTACGTTTACGT", "ACGTACGTACGTAA", Table(k => -3.0 - k, 14), "local", 2, -3);
        Assert.That(local, Is.EqualTo(AlignmentTools.LocalAlignAffine("ACGTACGTTTACGT", "ACGTACGTACGTAA", 2, -3, -3, -1)));
        var semi = AlignmentTools.GeneralGapAlign("ACGTACGT", "TTACGTACGTAA", Table(k => -3.0 - k), "semi_global", 2, -3, endGaps: "OverlapSuffix1Prefix2");
        Assert.That(semi, Is.EqualTo(AlignmentTools.SemiGlobalAlignAffine("ACGTACGT", "TTACGTACGTAA", 2, -3, -3, -1, endGaps: "OverlapSuffix1Prefix2")));
    }

    [Test]
    public void Matrix_Mode_ScoresProtein_FractionalAffineTable_EqualsEmbossNumbers()
    {
        // EMBOSS -gapopen 10 -gapextend 0.5 => w(k) = -(9.5 + 0.5 k); same value as global_align_affine with matrix
        var table = Table(k => -(9.5 + 0.5 * k), 20);
        var general = AlignmentTools.GeneralGapAlign("MKVLAAGIVGLLLAQ", "MKVLAGIVGLLLAPQ", table, "global", matrix: "BLOSUM62");
        var affine = AlignmentTools.GlobalAlignAffine("MKVLAAGIVGLLLAQ", "MKVLAGIVGLLLAPQ", gapOpen: -9.5, gapExtend: -0.5, matrix: "BLOSUM62");
        Assert.That(general.ExactScore, Is.EqualTo(affine.ExactScore).Within(1e-9));
    }

    [Test]
    public void Validation_ShortTable_UnknownMode_NaN_NullSequence()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign("ACGTACGT", "ACGT", new[] { -1.0, -2.0 }));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign("ACGT", "ACGT", Table(k => -k, 4), "sideways"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign("ACGT", "ACGT", new[] { -1.0, double.NaN, -3.0, -4.0 }));
        Assert.DoesNotThrow(() => AlignmentTools.GeneralGapAlign("ACXT", "ACGT", Table(k => -k, 4)));      // R14-2: any letters, like global_align
        Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign(null!, "ACGT", Table(k => -k, 4)));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GeneralGapAlign("ACGT", "ACGT", Table(k => -k, 4), "semi_global", endGaps: "Nope"));
    }
}
