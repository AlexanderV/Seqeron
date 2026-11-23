using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-11 (B09 F304): the <c>parasail_score</c> tool. Expected values: parasail 1.3.4 (<c>parasail.sw_stats("ACGTACGTAC", "TTACGTACGTAA", 5, 2, matrix_create("ACGT", 2, -3))</c> and the semi-global / global variants).
/// </summary>
[TestFixture]
public class ParasailScoreToolTests
{
    private const string Query = "ACGTACGTAC";
    private static readonly string[] References = ["TTACGTACGTAA", "GGGG", "ACGTACGTAC"];

    [Test]
    public void SwStats_ManyReferences_EqualParasail()
    {
        var r = AlignmentTools.ParasailScore(Query, References, "sw", 5, 2, statistics: true);
        Assert.That(r.Function, Is.EqualTo("sw"));
        Assert.That(r.Count, Is.EqualTo(3));
        // reference 1: ACGTACGTA aligned inside TTACGTACGTAA: nine matches (score 18), query 0..8, reference 2..10
        Assert.That((r.Results[0].Score, r.Results[0].EndQuery, r.Results[0].EndReference, r.Results[0].Matches, r.Results[0].Similar, r.Results[0].Length), Is.EqualTo((18, 8, 10, 9, 9, 9)));
        Assert.That(r.Results[2].Score, Is.EqualTo(20));
    }

    [Test]
    public void Global_SemiGlobal_TablesAndRowCol_AreConsistent()
    {
        var nw = AlignmentTools.ParasailScore("ACGT", ["ACGT"], "nw", 5, 2);
        Assert.That((nw.Results[0].Score, nw.Results[0].EndQuery, nw.Results[0].EndReference), Is.EqualTo((8, 3, 3)));
        var sg = AlignmentTools.ParasailScore("ACGT", ["TTACGTTT"], "sg_dx", 5, 2);
        Assert.That((sg.Results[0].Score, sg.Results[0].EndQuery, sg.Results[0].EndReference), Is.EqualTo((8, 3, 5)));
        var t = AlignmentTools.ParasailScore("ACGT", ["ACGT"], "nw", 5, 2, table: true, rowcol: true, statistics: true);
        Assert.That(t.Results[0].ScoreTable![3][3], Is.EqualTo(8));
        Assert.That(t.Results[0].ScoreRow![3], Is.EqualTo(8));
        Assert.That(t.Results[0].ScoreColumn![3], Is.EqualTo(8));
        Assert.That(t.Results[0].MatchesTable![3][3], Is.EqualTo(4));
        var parallel = AlignmentTools.ParasailScore(Query, References, "sg", 5, 2, parallel: true);
        var serial = AlignmentTools.ParasailScore(Query, References, "sg", 5, 2);
        Assert.That(parallel.Results.Select(x => x.Score), Is.EqualTo(serial.Results.Select(x => x.Score)));
    }

    [Test]
    public void Matrix_AndRefusals()
    {
        var protein = AlignmentTools.ParasailScore("MKTAYIAKQR", ["MKTAYIAKQR"], "nw", 11, 1, matrix: "BLOSUM62");
        Assert.That(protein.Results[0].Score, Is.EqualTo(49));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ParasailScore("ACGT", [], "nw"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ParasailScore("ACGT", ["ACGT"], "xx"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ParasailScore("ACGT", ["ACGT"], "nw", open: -1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.ParasailScore("ACGT", ["ACGT"], "nw", matrix: "GONNET1992"));
    }
}
